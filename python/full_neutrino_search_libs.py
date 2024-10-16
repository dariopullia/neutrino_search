import argparse
import sys
import numpy as np
import time
# import numba as nb
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# DAQ libraries
import daqdataformats
import detdataformats
import fddetdataformats
import trgdataformats
from hdf5libs import HDF5RawDataFile

# function to convert a TriggerPrimitive to a numpy array
def tp_to_numpy(tp):
    return np.array([tp.time_start/32, tp.time_over_threshold/32, tp.time_peak/32, tp.channel, tp.adc_integral, tp.adc_peak, tp.detid, tp.type, tp.algorithm, tp.version, tp.flag])


def search_neutrinos(filename, time_limit, adc_integral_cut):
    hdf5_file = HDF5RawDataFile(filename)
    # Get all records (TimeSlice)
    records = hdf5_file.get_all_record_ids()
    # Set the number of records to process
    print(f'Input file: {filename}')
    print(f'Total number of records: {len(records)}')

    # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
    tp_size = trgdataformats.TriggerPrimitive.sizeof()
    record_id=0
    all_sections = []

    for record in records:
        time_start = time.time()
        all_tps = []
        print(f'Processing record {record}') # TODO can add verbose option and put all these in it
        # Get all data sources
        datasets = hdf5_file.get_fragment_dataset_paths(record)
        # These are arrays since there can be multiple fragments per record
        n_tps = []
        n_tps_remaining = []
        next_tp_in_datasets = []
        fragments = []
        for dataset in datasets:
            # Get the fragment
            frag = hdf5_file.get_frag(dataset)
            fragments.append(frag)
            # Get the number of TPs
            n_tps_remaining.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            n_tps.append(int(frag.get_data_size()/trgdataformats.TriggerPrimitive.sizeof()))
            # Get the data
            next_tp_in_datasets.append(trgdataformats.TriggerPrimitive(frag.get_data(0)).time_start)
        # this is here to order the TPs by time_start and not just as they are written in the hdf5 file
        while len(n_tps_remaining) > 0:
            index = np.argmin(next_tp_in_datasets)
            tp = trgdataformats.TriggerPrimitive(fragments[index].get_data(tp_size*(n_tps[index]-n_tps_remaining[index])))
            # all_tps.append(tp_to_numpy(tp))
            tp_np = tp_to_numpy(tp)
            if not (tp_np[3]//2560 == 0 or tp_np[3]//2560 == 2):
                all_tps.append(tp_np)

            n_tps_remaining[index] -= 1
            if n_tps_remaining[index] == 0:
                del n_tps_remaining[index]
                del next_tp_in_datasets[index]
                del n_tps[index]
                del fragments[index]
            else:
                next_tp_in_datasets[index] = trgdataformats.TriggerPrimitive(fragments[index].get_data(tp_size*(n_tps[index]-n_tps_remaining[index]))).time_start
        record_id+=1
        all_tps = np.array(all_tps, dtype=int)
        current_record_interesting_sections = extract_interesting_sections(all_tps, time_limit, adc_integral_cut)
        print(f'Number of sections in record {record}: {len(current_record_interesting_sections)}')
        del all_tps
        all_sections.extend(current_record_interesting_sections)
        print(f'Time to process record {record}: {time.time()-time_start}')

    return all_sections



# @nb.njit
def extract_interesting_sections(all_tps, time_limit, adc_integral_cut):
    print(all_tps.shape)
    # remove second plane
    all_tps = all_tps[all_tps[:,3]//2560 != 2]
    all_tps = all_tps[all_tps[:,2].argsort()]

    section_start_idx = 0
    section_adc_integral = 0
    section_time_start = all_tps[0][2]

    all_sections = []
    for tp_index in range(len(all_tps)):
        tp = all_tps[tp_index]
        if tp[2] - section_time_start > time_limit:
            if section_adc_integral > adc_integral_cut:
                section = all_tps[section_start_idx:tp_index]
                all_sections.append(section)
            section_start_idx = tp_index
            section_adc_integral = 0
            section_time_start = tp[2]
        else:
            section_adc_integral += tp[4]
    return all_sections


def channel_to_z(channel):
    apa_length_in_cm = 230
    wire_pitch_in_cm_collection = 0.479
    offset_between_apa_in_cm = 2.4

    z_apa_offset = (channel // (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm)
    z_channel_offset = (((channel % 2560) - 1600) % 480) * wire_pitch_in_cm_collection
    z = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset
    return z


def is_groundshake(tps):
    total_charge = np.sum(tps[:,4])
    tot_charge_in_first_20cm = np.sum(tps[channel_to_z(tps[:,3]) < 20][:,4])
    if tot_charge_in_first_20cm > 0.02 * total_charge:
        return True
    return False



def is_vertical_cosmic(tps):
    x = tps[:,2]*0.08
    x_max = np.max(x)
    x_min = np.min(x)
    ntps_1_third = len(tps[x < x_min + (x_max - x_min)/3])
    ntps_2_third = len(tps[x < x_min + 2*(x_max - x_min)/3])
    ntps_2_third = ntps_2_third - ntps_1_third
    ntps_3_third = len(tps) - ntps_2_third - ntps_1_third

    # if ntps are similar in the 3 thirds of the image, it is a vertical cosmic
    if abs(ntps_1_third - ntps_2_third) < 0.15*len(tps) and abs(ntps_2_third - ntps_3_third) < 0.15*len(tps) and abs(ntps_1_third - ntps_3_third) < 0.15*len(tps):
        return True
    return False

def do_plots(tps, img_counter, filename, output_folder):
    time_offset = np.min(tps[:,0])
    tps = tps[tps[:,3]//2560 != 2]
    tps[:,0] = tps[:,0] - time_offset
    tps[:,2] = tps[:,2] - time_offset

    x = channel_to_z(tps[:,3])
    y = tps[:,2] * 0.08
    c = tps[:,4]

    # compute the total charge in a nxn subdivision of the image
    total_charge = np.sum(c)
    x_max = np.max(x)
    x_min = np.min(x)
    y_max = np.max(y)
    y_min = np.min(y)
    n = 3

    charge_matrix = np.zeros((n,n)) 
    for i in range(n):
        for j in range(n):
            x_low = x_min + i*(x_max-x_min)/n
            x_high = x_min + (i+1)*(x_max-x_min)/n
            y_low = y_min + j*(y_max-y_min)/n
            y_high = y_min + (j+1)*(y_max-y_min)/n
            mask = (x >= x_low) & (x < x_high) & (y >= y_low) & (y < y_high)
            charge_matrix[i,j] = np.sum(c[mask])


    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(4, 4)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    # Scatter plot
    scatter = ax_main.scatter(x, y, c=c, s=1)

    # ax_main.plot(z_avg, y_avg, color='red', linewidth=1)
    # ax_main.axvline(x=z_min, color='red', linestyle='--')
    # ax_main.axvline(x=z_max, color='red', linestyle='--')

    ax_main.set_xlabel('Z (cm)')
    ax_main.set_ylabel('X (cm)')
    ax_main.grid(True)
    # plot the charge matrix on top of the image as a heatmap
    for i in range(n):
        for j in range(n):
            ax_main.text(x_min + (i+0.5)*(x_max-x_min)/n, y_min + (j+0.5)*(y_max-y_min)/n, f'{charge_matrix[i,j]/total_charge:.2f}',
            ha='center', va='center', color='black', fontsize=8)


    # Histograms
    ax_xhist.hist(x, bins=20, color='gray', weights=c)
    ax_yhist.hist(y, bins=20, orientation='horizontal', color='gray', weights=c)

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    plt.setp(ax_xhist.get_xticklabels(), visible=False)
    plt.setp(ax_yhist.get_yticklabels(), visible=False)

    # Colorbar
    cbar = plt.colorbar(scatter, ax=ax_yhist, orientation='vertical')
    cbar.set_label('TP Value')
    cbar.ax.yaxis.set_label_position('right')
    cbar.ax.yaxis.set_ticks_position('right')

    # add title
    ax_main.set_title(f'Cluster {img_counter}. Number of TPs: {len(tps)}. Total charge: {total_charge/1E6:.2f} E6')
    plt.savefig(output_folder + '/cluster_' + filename + "_" + str(img_counter) + '.png')
    plt.close(fig)





