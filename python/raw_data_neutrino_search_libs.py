import os
import numpy as np
import argparse
import json
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import random 

# DAQ libraries
from hdf5libs import HDF5RawDataFile
import detdataformats
import fddetdataformats
from rawdatautils.unpack.wibeth import *
from rawdatautils.utilities.wibeth import *
import rawdatautils.utilities 
import detchannelmaps


from scipy.optimize import curve_fit

n_bins_x = 50
n_bins_y = 150
additional_angle = 7*np.pi/180


def gaussian(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def linear(x, a, b, x0):
    return a*(x-x0) + b



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


def from_hdf5_to_wf(filename, channel_map_name, det, channel_limit_mins=[], channel_limit_maxs=[], monitor_charge=[], monitor_time=[], beam_on=1, monitor_threshold=2.1E6, ask_for_close_beam=1):
    print("Processing file: ", filename)
    h5_file = HDF5RawDataFile(filename)
    records = h5_file.get_all_record_ids()
    channel_map = detchannelmaps.make_map(channel_map_name)
    n_records_processed = 0
    records_adcs = []
    records_ids = []
    for r in records:
    # for r in records[:6]:
        # if r[0] != 57570:
        # #     continue
        total_adcs = []
        tot_chan = []
        wib_geo_ids = h5_file.get_geo_ids_for_subdetector(r,detdataformats.DetID.string_to_subdetector(det))

        for gid in wib_geo_ids:
            frag_first = h5_file.get_frag(r,gid)
            break
        # find the arguments of the monitor matrix that are closest to the start of the fragment
        frag_min_time = frag_first.get_window_begin()
        monitor_index = np.argmin(np.abs(monitor_time - frag_min_time))
            
        if beam_on:
            if monitor_charge[monitor_index] < monitor_threshold:
                print("Skipping event, beam off. Charge: ", monitor_charge[monitor_index])
                continue
        else:       
            if monitor_charge[monitor_index] > monitor_threshold:
                print("Skipping event, beam on. Charge: ", monitor_charge[monitor_index])
                continue
            monitor_interval = monitor_charge[monitor_index-30:monitor_index+30]
            if np.max(monitor_interval) < monitor_threshold and ask_for_close_beam:
                print("Skipping event, maybe dump. Charge: ", monitor_charge[monitor_index])
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
        total_adcs = total_adcs[total_adcs[:,0].argsort()]
        medians = np.median(total_adcs[:, 1:], axis=1)
        total_adcs[:, 1:] = total_adcs[:, 1:] - medians[:, None]
        total_adcs = total_adcs[:,:3000]
        print(total_adcs.shape)
        n_records_processed += 1
        if contains_interesting_section(total_adcs, 800, 10E6):
            records_adcs.append(total_adcs)
            records_ids.append(r[0])
            print("Interesting section found")

    return records_adcs, records_ids


def contains_interesting_section(adcs, time_limit, adc_integral_cut):
    '''
    Check if the section contains a neutrino signal
    '''
    n_ticks = adcs.shape[1]
    for i in range(1, n_ticks-time_limit):
        if np.sum(adcs[:,i:i+time_limit]) > adc_integral_cut:
            return True
    return False

def find_most_interesting_section(adcs, time_limit):
    '''
    Find the interesting section
    '''
    n_ticks = adcs.shape[1]
    i_start = 0
    max_integral = 0
    for i in range(1, n_ticks-time_limit):
        if np.sum(adcs[:,i:i+time_limit]) > max_integral:
            max_integral = np.sum(adcs[:,i:i+time_limit])
            i_start = i
    # return adcs[:,i_start:i_start+time_limit]
    channels = adcs[:,0]
    i_end = i_start + time_limit
    i_start = np.max([0, i_start-100])
    i_end = np.min([n_ticks, i_end+100])

    final_adcs = adcs[:,i_start:i_end]
    final_adcs = np.insert(final_adcs, 0, channels, axis=1)
    final_adcs = final_adcs[final_adcs[:,0].argsort()]

    return final_adcs

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


def find_region_of_interest(adc, frac = 0.2):
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    frac_max_value_hist_y = np.max(hist_y[0]) * frac

    argmax_hist_y = np.argmax(hist_y[0])
    argmax_hist_y_min = np.argmax(hist_y[0])
    argmax_hist_y_max = np.argmax(hist_y[0])

    while hist_y[0][argmax_hist_y_max] > frac_max_value_hist_y and argmax_hist_y_max < len(hist_y[0])-1:
        argmax_hist_y_max += 1

    while hist_y[0][argmax_hist_y_min] > frac_max_value_hist_y and argmax_hist_y_min > 1:
        argmax_hist_y_min -= 1
    argmax_hist_y_min += 1

    adc_reduced = adc.copy()
    # take only yhe interval of interest
    channels = adc_reduced[:,0]
    adc_reduced = adc_reduced[:,int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])]
    adc_reduced = np.insert(adc_reduced, 0, channels, axis=1)

    hist_x = np.histogram(np.arange(adc_reduced.shape[0]), bins=n_bins_x, weights=np.sum(adc_reduced[:,1:], axis=1))
    frac_max_value_hist_x = np.max(hist_x[0]) * 0.2

    argmax_hist_x = np.argmax(hist_x[0])
    argmax_hist_x_min = np.argmax(hist_x[0])
    argmax_hist_x_max = np.argmax(hist_x[0])

    while hist_x[0][argmax_hist_x_max] > frac_max_value_hist_x and argmax_hist_x_max < len(hist_x[0])-1:
        argmax_hist_x_max += 1

    while hist_x[0][argmax_hist_x_min] > frac_max_value_hist_x and argmax_hist_x_min > 1:
        argmax_hist_x_min -= 1
    argmax_hist_x_min += 1

    return argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max


def fake_filters(adc):
    return True

def check_region_position(adc): 

    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc)
    # check if the region is at the left border
    if argmax_hist_x_min < 5:
        return False
    return True

def check_gaussian_width_y_coord(adc):
    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc)
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    initial_amplitude = np.max(hist_y[0][argmax_hist_y_min:argmax_hist_y_max])
    initial_mean = hist_y[1][argmax_hist_y_min] + (hist_y[1][argmax_hist_y_max] - hist_y[1][argmax_hist_y_min]) / 2
    initial_sigma = (hist_y[1][argmax_hist_y_max] - hist_y[1][argmax_hist_y_min]) / 4

    try:
        popt, pcov = curve_fit(gaussian, hist_y[1][argmax_hist_y_min:argmax_hist_y_max], hist_y[0][argmax_hist_y_min:argmax_hist_y_max], p0=[initial_amplitude, initial_mean, initial_sigma])

        if popt[2] > 100:
            return False
    except:
        if argmax_hist_y_max - argmax_hist_y_min > 550:
            return False

    return True

def check_if_there_is_a_muon_tail(adc):

    frac = 0.2
    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc, frac = frac)
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    hist_x = np.histogram(np.arange(adc.shape[0]), bins=n_bins_x, weights=np.sum(adc[:,1:], axis=1))
    adc_reduced = adc.copy()
    # take only yhe interval of interest
    channels = adc_reduced[:,0]
    adc_reduced = adc_reduced[:,int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])]

    # Find the vertex 
    weighted_average = []    
    corresponding_x = []
    for i in range(argmax_hist_x_min, argmax_hist_x_max):
        bin_values = np.sum(adc_reduced[int(hist_x[1][i]):int(hist_x[1][i+1]),1:], axis=0)
        weighted_average.append(np.average(np.arange(adc_reduced.shape[1]-1), weights=bin_values)+int(hist_y[1][argmax_hist_y_min]))
        corresponding_x.append(hist_x[1][i] + (hist_x[1][i+1] - hist_x[1][i])/2)

    corresponding_x = np.array(corresponding_x)
    weighted_average = np.array(weighted_average)
    
    # fit a line with the weighted average
    try:
        x_0 = corresponding_x[0]
        popt, pcov = curve_fit(lambda x, a, b: linear(x, a, b, x_0), corresponding_x, weighted_average)
    except:
        return False    
    
    image = np.zeros([adc.shape[1]-1, adc.shape[0]])
    adc = adc[adc[:,0].argsort()]
    for j in range(adc.shape[0]):
        image[:,j] = adc[j,1:]
    image = image[:, 0:int(hist_x[1][argmax_hist_x_min])]
    # everythin above the line should be zero
    n_muon_hits = 0
    for i in range(image.shape[1]):
        column = image[:,i]
        value1= linear(i, np.tan(np.arctan(popt[0]) - additional_angle), popt[1]+30, x_0)
        column[int(value1):] = 0
        value2= linear(i, np.tan(np.arctan(popt[0]) + additional_angle), popt[1]-30, x_0)
        column[:int(value2)] = 0

        # check if there is a muon hit, it means at least 5 consecutive channels with total value > 1000
        # print(np.sum(column))
        if find_muon_hit(column):
            n_muon_hits += 1
    # print(n_muon_hits)
    # print(image.shape[1])
    # print(n_muon_hits/(image.shape[1] -20))
    if n_muon_hits/(image.shape[1] - 20) > 0.85:
        return True
    return False


def find_muon_hit(column):
    '''
    Find a muon hit
    '''
    # for i in range(len(column)-10):
    #     if np.sum(column[i:i+10]) > 1500 and np.sum(column[i:i+5]) < 5000:
    #         return True
    if np.sum(column) > 1000:
        return True


    return False

def muon_tail_angle_filter(adc):
    frac = 0.2
    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc, frac = frac)
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    hist_x = np.histogram(np.arange(adc.shape[0]), bins=n_bins_x, weights=np.sum(adc[:,1:], axis=1))
    adc_reduced = adc.copy()
    # take only yhe interval of interest
    channels = adc_reduced[:,0]
    adc_reduced = adc_reduced[:,int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])]

    # Find the vertex 
    weighted_average = []    
    corresponding_x = []
    for i in range(argmax_hist_x_min, argmax_hist_x_max):
        bin_values = np.sum(adc_reduced[int(hist_x[1][i]):int(hist_x[1][i+1]),1:], axis=0)
        weighted_average.append(np.average(np.arange(adc_reduced.shape[1]-1), weights=bin_values)+int(hist_y[1][argmax_hist_y_min]))
        corresponding_x.append(hist_x[1][i] + (hist_x[1][i+1] - hist_x[1][i])/2)

    corresponding_x = np.array(corresponding_x)
    weighted_average = np.array(weighted_average)
    
    # fit a line with the weighted average
    try:
        x_0 = corresponding_x[0]
        popt, pcov = curve_fit(lambda x, a, b: linear(x, a, b, x_0), corresponding_x, weighted_average) 
    except:
        return True    
    
    if popt[0] < -1.3:
        return False
    if popt[0] > 1.3:
        return False

    return True

def pass_all_filters(adc):
    if not fake_filters(adc):
        return False
    if not check_region_position(adc):
        return False
    if not check_gaussian_width_y_coord(adc):
        return False
    if check_if_there_is_a_muon_tail(adc):
        return False
    if not muon_tail_angle_filter(adc):
        return False

    return True




def do_plot_wf(adc, img_counter, filename, output_folder):
    if is_groundshake(adc):
        return
        if not os.path.exists(f'{output_folder}/groundshake'):
            os.makedirs(f'{output_folder}/groundshake')
        output_folder = f'{output_folder}/groundshake'

    # elif correlation_across_consecutive_channels(adc) > 0.7:
    #     if not os.path.exists(f'{output_folder}/high_corr'):
    #         os.makedirs(f'{output_folder}/high_corr')
    #     output_folder = f'{output_folder}/high_corr'

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(4, 4)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    ax_text_space = fig.add_subplot(gs[0, 3])

    # Scatter plot
    image = np.zeros([adc.shape[1]-1, adc.shape[0]])
    for j in range(adc.shape[0]):
        image[:,j] = adc[j,1:]

    image = np.where(image > 1000, 1000, image)
    ax_main.imshow(image, aspect='auto')
    ax_main.set_xlabel('Z (cm)')
    ax_main.set_ylabel('X (cm)')
    ax_main.set_ylim(0, adc.shape[1]-1)
    # Histograms

    ax_yhist.hist(np.arange(adc.shape[1]-1), bins=n_bins_y, orientation='horizontal', color='gray', weights=np.sum(adc[:,1:], axis=0))

    frac = 0.2
    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc, frac = frac)
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    hist_x = np.histogram(np.arange(adc.shape[0]), bins=n_bins_x, weights=np.sum(adc[:,1:], axis=1))
    adc_reduced = adc.copy()
    # take only yhe interval of interest
    channels = adc_reduced[:,0]
    adc_reduced = adc_reduced[:,int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])]
    ax_xhist.hist(np.arange(adc_reduced.shape[0]), bins=n_bins_x, color='gray', weights=np.sum(adc_reduced[:,1:], axis=1))


    ax_xhist.axvline(hist_x[1][argmax_hist_x_min], color='b', linestyle='--')
    ax_xhist.axvline(hist_x[1][argmax_hist_x_max], color='b', linestyle='--')
    ax_yhist.axhline(hist_y[1][argmax_hist_y_min], color='b', linestyle='--')
    ax_yhist.axhline(hist_y[1][argmax_hist_y_max], color='b', linestyle='--')

    ax_main.axvline(hist_x[1][argmax_hist_x_min], color='w', linestyle='--')
    ax_main.axvline(hist_x[1][argmax_hist_x_max], color='w', linestyle='--')
    ax_main.axhline(hist_y[1][argmax_hist_y_min], color='w', linestyle='--')
    ax_main.axhline(hist_y[1][argmax_hist_y_max], color='w', linestyle='--')


    # Improved initial values for curve fitting
    initial_amplitude = np.max(hist_y[0][argmax_hist_y_min:argmax_hist_y_max])
    initial_mean = hist_y[1][argmax_hist_y_min] + (hist_y[1][argmax_hist_y_max] - hist_y[1][argmax_hist_y_min]) / 2
    initial_sigma = (hist_y[1][argmax_hist_y_max] - hist_y[1][argmax_hist_y_min]) / 4
    try:
        popt, pcov = curve_fit(gaussian, hist_y[1][argmax_hist_y_min:argmax_hist_y_max], hist_y[0][argmax_hist_y_min:argmax_hist_y_max], p0=[initial_amplitude, initial_mean, initial_sigma])
        ax_yhist.plot(gaussian(hist_y[1][argmax_hist_y_min:argmax_hist_y_max], *popt), hist_y[1][argmax_hist_y_min:argmax_hist_y_max], 'r-')
        # add text
        ax_yhist.text(0.5, 0.5, f'Mean: {popt[1]:.2f}', horizontalalignment='center', verticalalignment='center', transform=ax_yhist.transAxes, color='r')  
        ax_yhist.text(0.5, 0.4, f'Sigma: {popt[2]:.2f}', horizontalalignment='center', verticalalignment='center', transform=ax_yhist.transAxes, color='r')
    except:
        print("Error in Gaussian fitting")

    # set y limits
    ax_xhist.set_ylim(0, 3E6)
    ax_yhist.set_xlim(0, 3.5E6)

    # Find the vertex 
    weighted_average = []    
    corresponding_x = []
    for i in range(argmax_hist_x_min, argmax_hist_x_max):
        bin_values = np.sum(adc_reduced[int(hist_x[1][i]):int(hist_x[1][i+1]),1:], axis=0)
        weighted_average.append(np.average(np.arange(adc_reduced.shape[1]-1), weights=bin_values)+int(hist_y[1][argmax_hist_y_min]))
        corresponding_x.append(hist_x[1][i] + (hist_x[1][i+1] - hist_x[1][i])/2)

    corresponding_x = np.array(corresponding_x)
    weighted_average = np.array(weighted_average)

    ax_main.scatter(corresponding_x, weighted_average, color='r')
    # fit a line with the weighted average
    try:
        x_0 = corresponding_x[0]
        popt, pcov = curve_fit(lambda x, a, b: linear(x, a, b, x_0), corresponding_x, weighted_average)
        ax_main.plot(np.linspace(0, corresponding_x[-1], 100), linear(np.linspace(0, corresponding_x[-1], 100), popt[0], popt[1], x_0), 'r-', alpha=0.5)
        ax_main.plot(np.linspace(0, corresponding_x[0], 100), linear(np.linspace(0, corresponding_x[0], 100), np.tan(np.arctan(popt[0]) - additional_angle), popt[1]+30, x_0), 'r--', alpha=0.5)
        ax_main.plot(np.linspace(0, corresponding_x[0], 100), linear(np.linspace(0, corresponding_x[0], 100), np.tan(np.arctan(popt[0]) + additional_angle), popt[1]-30, x_0), 'r--', alpha=0.5)
        ax_main.text(0.5, 0.9, f'Coeff: {popt[0]:.2f}, Intercept: {popt[1]:.2f}', horizontalalignment='center', verticalalignment='center', transform=ax_main.transAxes, color='w')
    except:
        print("Error in linear fitting")


    # Hide x labels and tick labels for top plots and y ticks for right plots.
    plt.setp(ax_xhist.get_xticklabels(), visible=False)
    plt.setp(ax_yhist.get_yticklabels(), visible=False)

    # add text
    ax_text_space.text(
        0.05,
        0.9,
        f"Pos: {check_region_position(adc)}",
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax_text_space.transAxes,
        color="black",
    )
    ax_text_space.text(
        0.05,
        0.8,
        f"Gauss: {check_gaussian_width_y_coord(adc)}",
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax_text_space.transAxes,
        color="black",
    )
    ax_text_space.text(
        0.05,
        0.7,
        f"Muon tail: {check_if_there_is_a_muon_tail(adc)}",
        horizontalalignment="left",
        verticalalignment="center",
        transform=ax_text_space.transAxes,
        color="black",
    )
    ax_text_space.text(
        0.05,
        0.6,
        f"Muon tail angle: {muon_tail_angle_filter(adc)}",
        horizontalalignment="left",
        verticalalignment="center",
        transform = ax_text_space.transAxes,
        color = "black",
    )

    # use correct name
    chargeE7 = np.sum(adc[int(hist_x[1][argmax_hist_x_min]):int(hist_x[1][argmax_hist_x_max]), int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])])/1E7
    if not os.path.exists(f'{output_folder}/{int(chargeE7)}E7'):
        os.makedirs(f'{output_folder}/{int(chargeE7)}E7')
    output_folder = f'{output_folder}/{int(chargeE7)}E7'

    # add title
    plt.suptitle(f'Cluster {img_counter}. Total charge: {chargeE7*10:.2f} E6. Pass Filters: {pass_all_filters(adc)}.', fontsize=16)
    plt.savefig(f'{output_folder}/{filename}_{img_counter}.png')
    plt.close(fig)

def do_clean_plot_wf(adc, img_counter, filename, output_folder):
    if is_groundshake(adc):
        return

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(4, 4)

    ax_main = fig.add_subplot(gs[1:4, 0:3])
    ax_xhist = fig.add_subplot(gs[0, 0:3], sharex=ax_main)
    ax_yhist = fig.add_subplot(gs[1:4, 3], sharey=ax_main)
    ax_text_space = fig.add_subplot(gs[0, 3])

    # Scatter plot
    image = np.zeros([adc.shape[1]-1, adc.shape[0]])
    for j in range(adc.shape[0]):
        image[:,j] = adc[j,1:]

    image = np.where(image > 1000, 1000, image)
    ax_main.imshow(image, aspect='auto')
    ax_main.set_xlabel('Z (cm)')
    ax_main.set_ylabel('X (cm)')
    ax_main.set_ylim(0, adc.shape[1]-1)
    # Histograms
    ax_yhist.hist(np.arange(adc.shape[1]-1), bins=n_bins_y, orientation='horizontal', color='gray', weights=np.sum(adc[:,1:], axis=0))
    ax_xhist.hist(np.arange(adc.shape[0]), bins=n_bins_x, color='gray', weights=np.sum(adc[:,1:], axis=1))

    argmax_hist_x_min, argmax_hist_x_max, argmax_hist_y_min, argmax_hist_y_max = find_region_of_interest(adc, frac = 0.2)
    hist_y = np.histogram(np.arange(adc.shape[1]-1), bins=n_bins_y, weights=np.sum(adc[:,1:], axis=0))
    hist_x = np.histogram(np.arange(adc.shape[0]), bins=n_bins_x, weights=np.sum(adc[:,1:], axis=1))


    chargeE7 = np.sum(adc[int(hist_x[1][argmax_hist_x_min]):int(hist_x[1][argmax_hist_x_max]), int(hist_y[1][argmax_hist_y_min]):int(hist_y[1][argmax_hist_y_max])])/1E7
    if not os.path.exists(f'{output_folder}/{int(chargeE7)}E7'):
        os.makedirs(f'{output_folder}/{int(chargeE7)}E7')
    output_folder = f'{output_folder}/{int(chargeE7)}E7'
    # add title
    plt.suptitle(f'Cluster {img_counter}. Total charge: {chargeE7*10:.2f} E6. Pass Filters: {pass_all_filters(adc)}.', fontsize=16)
    plt.savefig(f'{output_folder}/{filename}_{img_counter}.png')
    plt.close(fig)
