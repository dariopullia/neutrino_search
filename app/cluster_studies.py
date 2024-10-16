import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import json
import argparse
from sklearn.linear_model import LinearRegression
    
sys.path.append('../python/') 
from cluster import *
from clustering_algorithm import *
from matplotlib import gridspec



####################################################################################################

def channel_to_z(channel):
    apa_length_in_cm = 230
    wire_pitch_in_cm_collection = 0.479
    offset_between_apa_in_cm = 2.4

    z_apa_offset = (channel // (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm)
    z_channel_offset = (((channel % 2560) - 1600) % 480) * wire_pitch_in_cm_collection
    z = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset
    return z

####################################################################################################


parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('--input_json', type=str, help='Input json file')
args = parser.parse_args()

input_json_file = args.input_json


# Read input json
with open(input_json_file) as f:
    input_json = json.load(f)

input_root_folder = input_json['input_root_folder']
output_folder = input_json['output_folder']
if not os.path.exists(output_folder):   
    os.makedirs(output_folder)
if not os.path.exists(output_folder + '/gs'):
    os.makedirs(output_folder + '/gs')
if not os.path.exists(output_folder + '/vc'):
    os.makedirs(output_folder + '/vc')
if not os.path.exists(output_folder + '/sc'):
    os.makedirs(output_folder + '/sc')
colours = [
    '#FF0000',  # Red
    '#00FF00',  # Green
    '#0000FF',  # Blue
    '#FFFF00',  # Yellow
    '#FF00FF',  # Magenta
    '#00FFFF',  # Cyan
    '#800000',  # Maroon
    '#808000',  # Olive
    '#008000',  # Dark Green
    '#800080',  # Purple
    '#008080',  # Teal
    '#000080',  # Navy
    '#FFA500',  # Orange
    '#A52A2A',  # Brown
    '#808080',  # Gray
    '#FFC0CB',  # Pink
    '#FFD700',  # Gold
    '#4B0082',  # Indigo
    '#D2691E',  # Chocolate
    '#00FF7F',  # Spring Green
]

img_counter = 0
# Read the root file

#extract list of folders in the input folder
folders = [f for f in os.listdir(input_root_folder) if os.path.isdir(os.path.join(input_root_folder, f))]

def unbalance_condition(cluster):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    z_min = np.min(z)
    z_max = np.max(z)
    half_length = (z_max - z_min) / 2
    z_center = z_min + half_length
    # count number of tps in the first half and second half
    n_first_half = len(tps[z < z_center])
    n_second_half = len(tps[z >= z_center])
    if n_first_half > n_second_half:
        return True

    return False


def second_half_flatness_condition(cluster):
    tps = cluster.get_tps()
    tps = tps[tps[:,4] > 2000]
    tps = tps[tps[:,4] < 4000]
    if tps.shape[0] < 10:
        return False

    z = channel_to_z(tps[:,3])
    z_min = np.min(z)
    z_max = np.max(z)
    half_length = (z_max - z_min) / 2
    z_center = z_min + half_length
    second_half = tps[z >= z_center]
    if len(second_half) == 0:
        return False
    second_half_z = channel_to_z(second_half[:,3])
    second_half_z_min = np.min(second_half_z)
    second_half_z_max = np.max(second_half_z)
    second_half_z_range = second_half_z_max - second_half_z_min
    
    second_half_x = second_half[:,2]*0.08

    second_half_x_min = np.min(second_half_x)
    second_half_x_max = np.max(second_half_x)
    second_half_x_range = second_half_x_max - second_half_x_min
    
    second_half_ratio = second_half_x_range / second_half_z_range
    if second_half_ratio < 0.2:
        return True
    
    return False

def is_groundshake(cluster):
    total_charge = np.sum(cluster.get_tps()[:,4])
    tot_charge_in_first_20cm = np.sum(cluster.get_tps()[channel_to_z(cluster.get_tps()[:,3]) < 20][:,4])
    if tot_charge_in_first_20cm > 0.03*total_charge:
        return True
    return False



def is_strong_candidate(cluster):
    tps = cluster.get_tps() 
    hist, _ = np.histogram(tps[:,2], bins=10)
    n_channels = len(np.unique(tps[:,3]))   
    if np.max(hist) > n_channels*0.95:
        return True
    return False   

def check_without_vertical_tracks(cluster, threshold_on_adc, threshold_total):
    tps = cluster.get_tps()
    tps = tps[tps[:,4] < threshold_on_adc]
    if np.sum(tps[:,4]) > threshold_total:
        return True
    return False

def large_charge_in_first_section(cluster):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    total_charge = np.sum(tps[:,4])
    tps = tps[z<50]
    charge_first_section = np.sum(tps[:,4])
    if charge_first_section > 0.05*total_charge:
        return True
    return False


def compute_avg_charge_pos(cluster):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    z_uni_list = []
    y_avg_list = []
    for z_unique in np.unique(z):
        mask = z == z_unique
        y_avg = np.average(tps[mask, 2], weights=tps[mask, 4])
        z_uni_list.append(z_unique)
        y_avg_list.append(y_avg)
    z_uni_list = np.array(z_uni_list)
    y_avg_list = np.array(y_avg_list)*0.08
    
    return y_avg_list, z_uni_list


def continuous_track_max_limits(y_avg, z_avg, max_displacement):
    z_avg = z_avg[z_avg.argsort()]
    y_avg = y_avg[z_avg.argsort()]
    start_z = z_avg[0]
    current_min_z = 0
    current_max_z = 0
    for i in range(1, len(z_avg)):
        # print(abs(y_avg[i] - y_avg[i-1]))
        if abs(y_avg[i] - y_avg[i-1]) < max_displacement:
            continue
        else:
            if z_avg[i-1] - start_z > current_max_z - current_min_z:
                current_max_z = z_avg[i-1]
                current_min_z = start_z
            start_z = z_avg[i]

    if z_avg[-1] - start_z > current_max_z - current_min_z:
        current_max_z = z_avg[-1]
        current_min_z = start_z
    return current_min_z, current_max_z

def flat_charge_mean(tps, max_displacement, min_z_length):
    # Initialize dictionaries for mean_time_per_z and total_charge_per_z
    mean_time_per_z = {}
    total_charge_per_z = {}

    # Iterate over the trigger primitives
    for tp in tps:
        z = channel_to_z(tp[3])
        charge = tp[4]
        time = tp[2] * 0.08

        # Update mean_time_per_z
        if z not in mean_time_per_z:
            mean_time_per_z[z] = 0
        mean_time_per_z[z] += charge * time

        # Update total_charge_per_z
        if z not in total_charge_per_z:
            total_charge_per_z[z] = 0
        total_charge_per_z[z] += charge

    # Normalize the mean time by dividing by the total charge per z
    for z, charge in total_charge_per_z.items():
        mean_time_per_z[z] /= charge

    # Get sorted z values
    z_values = sorted(total_charge_per_z.keys())

    # Check if the mean time is flat
    start_z = z_values[0]
    current_min_z = 0
    current_max_z = 0
    tollerance_counter = 0

    for i in range(1, len(z_values)):
        if abs(mean_time_per_z[z_values[i]] - mean_time_per_z[z_values[i-1]]) < max_displacement:
            tollerance_counter = 0
            continue
        else:
            if tollerance_counter < 3:
                tollerance_counter += 1
                continue
            if z_values[i-1] - start_z > current_max_z - current_min_z:
                current_max_z = z_values[i-1]
                current_min_z = start_z
            start_z = z_values[i]

    # Final check for the z length
    if z_values[-1] - start_z > current_max_z - current_min_z:
        current_max_z = z_values[-1]
        current_min_z = start_z

    time_in_interval = []
    for z in z_values:
        if z >= current_min_z and z <= current_max_z:
            time_in_interval.append(mean_time_per_z[z])

    return current_min_z, current_max_z, np.mean(time_in_interval)

def there_is_a_track(cluster, z_max, time_mean):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    tps = tps[z > z_max]
    tps = tps[tps[:,2]*0.08 < time_mean + 15]
    tps = tps[tps[:,2]*0.08 > time_mean - 15]

    space_left = 460 - z_max
    if len(tps) < 0.5*space_left:
        return False
    return True

def tps_in_first_50cm(cluster, time_mean):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    tps = tps[z < 50]
    tps = tps[tps[:,2]*0.08 < time_mean + 10]
    tps = tps[tps[:,2]*0.08 > time_mean - 10]
    if len(tps) > 20:
        return True
    return False

def is_vertical_cosmic(cluster):
    tps = cluster.get_tps()
    x = tps[:,2]*0.08
    x_max = np.max(x)
    x_min = np.min(x)
    ntps_1_third = len(tps[x < x_min + (x_max - x_min)/3])
    ntps_2_third = len(tps[x < x_min + 2*(x_max - x_min)/3])
    ntps_2_third = ntps_2_third - ntps_1_third
    ntps_3_third = len(tps) - ntps_2_third - ntps_1_third

    # if ntps are similar in the 3 thirds of the image, it is a vertical cosmic
    if abs(ntps_1_third - ntps_2_third) < 0.25*len(tps) and abs(ntps_2_third - ntps_3_third) < 0.25*len(tps):
        return True

    return False


# def do_plots(cluster, img_counter, output_folder, y_avg, z_avg):
# def do_plots(cluster, img_counter, output_folder, y_avg, z_avg, z_max, z_min):
# def do_plots(cluster, img_counter, output_folder, model):
def do_plots(cluster, img_counter, output_folder):
    tps = cluster.get_tps()
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
    plt.savefig(output_folder + '/cluster_' + str(img_counter) + '.png')
    plt.close(fig)




standard_deviations = []

# folders = ["/eos/user/d/dapullia/dune/protodune/neutrino_search/extract_neutrino_candidates/eliminami/"]

for folder_idx, folder in enumerate(folders[:]):
    print(f"Processing folder: {folder}.    {folder_idx+1}/{len(folders)}")
    # Extract list of root files in the folder
    root_files = [f for f in os.listdir(os.path.join(input_root_folder, folder)) if f.endswith('.root')]
    print(root_files)
    for root_file in root_files:
        print(f"Processing root file: {root_file}")
        input_file = os.path.join(input_root_folder, folder, root_file)

        clusters = read_root_file_to_clusters(input_file)
        print(f"Number of clusters: {len(clusters)}")


        for i, cluster in enumerate(clusters):

            # y_avg, z_avg = compute_avg_charge_pos(cluster)
            # z_min, z_max, time_mean = flat_charge_mean(cluster.get_tps(), max_displacement=10, min_z_length=50)

            # if np.sum(cluster.get_tps()[:,4]) < 7000000:
            #     continue
            
            # event_tps = cluster.get_tps()


            # hist, bin_edges = np.histogram(event_tps[:,2]*0.08, bins=20, weights=event_tps[:,4], range=(0, 450))
            # bin_centers = (bin_edges[:-1] + bin_edges[1:])/2
            # peak = False
            # for i in range(0, len(hist)-5):
            #     if np.sum(hist[i:i+5]) > 0.8*np.sum(hist):
            #         peak = True
            #         break
            # if not peak:
            #     continue
            
            # # select only the TPs that are in the peak
            # mask = (event_tps[:,2]*0.08 > bin_centers[i]) & (event_tps[:,2]*0.08 < bin_centers[i+5])
            # event_tps = event_tps[mask]
            # z = z[mask]
            
            # not_from_outside = True
            # hist, bin_edges = np.histogram(z, bins=20, weights=event_tps[:,4], range=(0, 460))
            # bin_centers = (bin_edges[:-1] + bin_edges[1:])/2

            # # get the maximum of the histogram
            # max_bin = np.argmax(hist)   
            # # check if the maximum is not in the first 2 or last 2 bins
            # if max_bin < 5:
            #     not_from_outside = False
            #     continue
            
            # sudden_rise = True
            # if not_from_outside:
            #     if np.sum(hist[:max_bin-4]) > 0.2*np.sum(hist):
            #         sudden_rise = False
            #         continue

            if is_groundshake(cluster):
                do_plots(cluster, img_counter, output_folder+'/gs/')
                img_counter += 1
                continue



            if if_vertical_cosmic(cluster):
                do_plots(cluster, img_counter, output_folder+'/vc/')
                img_counter += 1
                continue

            if np.sum(cluster.get_tps()[:,4]) > 17000000:
                do_plots(cluster, img_counter, output_folder+'/sc/')
                img_counter += 1
                continue
      
            do_plots(cluster, img_counter, output_folder, )
            # do_plots(cluster, img_counter, output_folder, y_avg, z_avg, z_max, z_min)

            img_counter += 1

print(f"Number of clusters: {img_counter}")#1201