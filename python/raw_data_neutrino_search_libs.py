import os
import numpy as np
import argparse
import json
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# DAQ libraries
from hdf5libs import HDF5RawDataFile
import detdataformats
import fddetdataformats
from rawdatautils.unpack.wibeth import *
from rawdatautils.utilities.wibeth import *
import rawdatautils.utilities 
import detchannelmaps


def select_channels(wf, channel_limit_mins, channel_limit_maxs):
    '''
    Select the channels in the range channel_limit_mins and channel_limit_maxs
    '''
    if len(channel_limit_maxs) == 0:
        return wf
        
    survived_wf = []
    for i in range(len(channel_limit_maxs)):
        wf_over = wf[wf[:,0] >= channel_limit_mins[i]]
        wf_under = wf_over[wf_over[:,0] < channel_limit_maxs[i]]
        survived_wf.append(wf_under)

    survived_wf = np.concatenate(survived_wf, axis=0)
    return survived_wf


def from_hdf5_to_wf(filename, channel_map_name, det, channel_limit_mins=[], channel_limit_maxs=[], monitor_charge=[], monitor_time=[]):
    print("Processing file: ", filename)
    print(channel_limit_mins)
    print(channel_limit_maxs)
    h5_file = HDF5RawDataFile(filename)
    records = h5_file.get_all_record_ids()
    channel_map = detchannelmaps.make_map(channel_map_name)
    n_records_processed = 0
    records_adcs = []
    for r in records:
    # for r in records[:2]:
        print("Processing record: ", r)
        total_adcs = []
        tot_chan = []
        wib_geo_ids = h5_file.get_geo_ids_for_subdetector(r,detdataformats.DetID.string_to_subdetector(det))

        for gid in wib_geo_ids:
            frag_first = h5_file.get_frag(r,gid)
            break
        # find the arguments of the monitor matrix that are closest to the start of the fragment
        frag_min_time = frag_first.get_window_begin()
        monitor_index = np.argmin(np.abs(monitor_time - frag_min_time))
        print("Monitor index: ", monitor_index)
        print("Monitor time: ", monitor_time[monitor_index])
        print("Monitor charge: ", monitor_charge[monitor_index])
        if monitor_charge[monitor_index] < 2.1E6:
            print("Skipping event, charge too low. Charge: ", monitor_charge[monitor_index])
            continue
         


        for gid in wib_geo_ids:
            frag = h5_file.get_frag(r,gid)
            wf = fddetdataformats.WIBEthFrame(frag.get_data())
            adcs = np_array_adc(frag)
            # print(adcs)
            dh = wf.get_daqheader()
            N_CHANNELS_PER_FRAME = adcs.shape[-1]
            channels = [ channel_map.get_offline_channel_from_crate_slot_stream_chan(dh.crate_id, dh.slot_id, dh.stream_id, c) for c in range(N_CHANNELS_PER_FRAME) ]

            adcs = adcs.T
            # instert channel before adcs   
            adcs = np.insert(adcs, 0, channels, axis=1)
            adcs = select_channels(adcs, channel_limit_mins, channel_limit_maxs)
            total_adcs.append(adcs)
            tot_chan.append(channels)

        # check if all the lengths are the same, if not, skip the event
        list_of_lens = [x.shape[1] for x in total_adcs]
        if len(np.unique(list_of_lens)) != 1:
            print("Skipping event")
            print("Lengths: ", np.unique(list_of_lens))
            continue
        total_adcs = np.concatenate(total_adcs, axis=0)
        total_adcs = total_adcs.astype(np.int16)
        medians = np.median(total_adcs[:, 1:], axis=1)
        total_adcs[:, 1:] = total_adcs[:, 1:] - medians[:, None]
        print(total_adcs.shape)
        n_records_processed += 1
        # section = find_interesting_section(total_adcs, 800, 10E6)
        # if section is not None:
        #     records_adcs.append(section)
        #     print("Interesting section found")
        if contains_interesting_section(total_adcs, 800, 10E6):
            records_adcs.append(total_adcs)
            print("Interesting section found")

    return records_adcs


def contains_interesting_section(adcs, time_limit, adc_integral_cut):
    '''
    Check if the section contains a neutrino signal
    '''
    n_ticks = adcs.shape[1]
    for i in range(1, n_ticks-time_limit):
        if np.sum(adcs[:,i:i+time_limit]) > adc_integral_cut:
            return True
    return False

def find_interesting_section(adcs, time_limit, adc_integral_cut):
    '''
    Find the interesting section
    '''
    n_ticks = adcs.shape[1]
    i_start = 0
    max_integral = 0
    for i in range(1, n_ticks-time_limit):
        if np.sum(adcs[:,i:i+time_limit]) > adc_integral_cut:
            if np.sum(adcs[:,i:i+time_limit]) > max_integral:
                max_integral = np.sum(adcs[:,i:i+time_limit])
                i_start = i
    if max_integral == 0:
        return None
    return adcs[:,i_start:i_start+time_limit]




    
# def do_plot_wf(adc, img_counter, filename, output_folder):
#     '''
#     Plot the waveforms
#     '''
#     plt.figure(figsize=(10,10))
#     # adc = np.where(adc < 5, np.nan, adc)
#     adc = adc.astype(np.uint16)
#     image = np.zeros([adc.shape[1]-1, adc.shape[0]])
#     adc = adc[adc[:,0].argsort()]
#     for j in range(adc.shape[0]):
#         image[:,j] = adc[j,1:]
#     # plt.imshow(image, aspect=adc.shape[0]/adc.shape[1])
#     plt.imshow(image, aspect='auto')
#     plt.title(f'Total_charge: {np.sum(adc[:,1:])/1E6}E6') 
#     # plt.colorbar(fraction=adc.shape[0]/adc.shape[1], pad=0.04)
#     plt.savefig(f'{output_folder}/{filename}_{img_counter}.png')
#     plt.close()

def is_groundshake(adcs):
    '''
    Check if the event is a groundshake
    A groundshake is defined as an event where all the channels are active at the same time
    '''
    total_charge = np.sum(adcs[:,1:])
    for i in range(1, adcs.shape[1]-2):
        row = adcs[:,i]
        if np.where(np.abs(row) > 10)[0].shape[0] >= adcs.shape[0]*0.85:
            return True
    return False

def correlation_across_consecutive_channels(adcs):
    '''
    Calculate the correlation across consecutive channels
    '''
    correlations = []
    for i in range(adcs.shape[0]-1):
        column = adcs[i,1:]
        column_next = adcs[i+1,1:]
        corr = np.corrcoef(column, column_next)
        correlations.append(corr[0,1])

    return np.mean(correlations)



def do_plot_wf(adc, img_counter, filename, output_folder):
    if is_groundshake(adc):
        if not os.path.exists(f'{output_folder}/groundshake'):
            os.makedirs(f'{output_folder}/groundshake')
        output_folder = f'{output_folder}/groundshake'
    elif correlation_across_consecutive_channels(adc) > 0.7:
        if not os.path.exists(f'{output_folder}/high_corr'):
            os.makedirs(f'{output_folder}/high_corr')
        output_folder = f'{output_folder}/high_corr'
    else:   
        section = find_interesting_section(adc, 800, 10E6)
        chargeE7 = np.sum(adc[:,1:])/10E6
        chargeE7 = int(chargeE7)
        if not os.path.exists(f'{output_folder}/{chargeE7}E7'):
            os.makedirs(f'{output_folder}/{chargeE7}E7')
        output_folder = f'{output_folder}/{chargeE7}E7'

    print(correlation_across_consecutive_channels(adc))

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(4, 4)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)

    # Scatter plot
    adc_contrast = adc.astype(np.uint16)
    image = np.zeros([adc.shape[1]-1, adc.shape[0]])
    adc_contrast = adc_contrast[adc_contrast[:,0].argsort()]
    for j in range(adc_contrast.shape[0]):
        image[:,j] = adc_contrast[j,1:]
    ax_main.imshow(image, aspect='auto')
    ax_main.set_xlabel('Z (cm)')
    ax_main.set_ylabel('X (cm)')

    # Histograms
    ax_xhist.hist(np.arange(adc.shape[0]), bins=20, color='gray', weights=np.sum(adc[:,1:], axis=1))
    ax_yhist.hist(np.arange(adc.shape[1]-1), bins=20, orientation='horizontal', color='gray', weights=np.sum(adc[:,1:], axis=0))

    # Hide x labels and tick labels for top plots and y ticks for right plots.
    plt.setp(ax_xhist.get_xticklabels(), visible=False)
    plt.setp(ax_yhist.get_yticklabels(), visible=False)

    # add title
    plt.suptitle(f'Cluster {img_counter}. Total charge: {np.sum(adc[:,1:])/1E6:.2f} E6. Groundshake: {is_groundshake(adc)}. Correlation: {correlation_across_consecutive_channels(adc):.2f}', fontsize=16)  
    plt.savefig(f'{output_folder}/{filename}_{img_counter}.png')
    plt.close(fig)




