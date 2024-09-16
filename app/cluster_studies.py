import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
import os
import json
import argparse

sys.path.append('../python/') 
from cluster import *
from clustering_algorithm import *

parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('--input_json', type=str, help='Input json file')
parser.add_argument('--output_folder', type=str, help='Output folder')
args = parser.parse_args()

input_json_file = args.input_json
output_folder = args.output_folder

if not os.path.exists(output_folder):   
    os.makedirs(output_folder)

# Read input json
with open(input_json_file) as f:
    input_json = json.load(f)

input_root_file = input_json['input_root_file']
input_tps = input_json['input_tps'] 

# Read the root file
clusters = read_root_file_to_clusters(input_root_file)
print(f"Number of clusters: {len(clusters)}")

tps = np.loadtxt(input_tps, max_rows=35000, skiprows=0)  
# print(tps[:10])
tps[:,0] = tps[:,0]/32
tps[:,1] = tps[:,1]/32
tps[:,2] = tps[:,2]/32

tps = tps[tps[:,3]%2560 >= 1600]
tps = tps[tps[:,1] > 10]

# clusters_py = cluster_maker(tps,ticks_limit=7,channel_limit=8,min_tps_to_cluster=1,adc_integral_cut=0)

# print(f"Number of clusters: {len(clusters_py)}")


tps = tps[tps[:,4] < 4000]
tps = tps[tps[:,4] > 2000]


# print(tps[:10])

# tps=tps[tps[:,1] > 10]

tps = tps[tps[:,3] > 2560]
tps = tps[tps[:,3]//2560 != 2]

def channel_to_z(channel):
    apa_length_in_cm = 230
    wire_pitch_in_cm_collection = 0.479
    offset_between_apa_in_cm = 2.4

    z_apa_offset = (channel // (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm)
    z_channel_offset = (((channel % 2560) - 1600) % 480) * wire_pitch_in_cm_collection
    z = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset
    return z

def filter_by_position(cluster, z_min, z_max, z_component_min, x_component_max):
    tps = cluster.get_tps()
    z = channel_to_z(tps[:,3])
    # not_from_outside_beam_side = False
    if np.min(z)<z_min:
        # not_from_outside_beam_side = True
        return False

    # goes_out_far_from_beam = False
    if np.max(z) < z_max:
        # goes_out_far_from_beam = True
        return False
    
    if (np.max(z) - np.min(z)) < z_component_min:
        return False

    x_max = np.max(tps[:,2])
    x_min = np.min(tps[:,2])
    x_length = (x_max - x_min)*0.08
    if x_length > x_component_max:
        return False



    return True

def extract_subclusters(cluster, adc_min, adc_max, ticks_limit=7, channel_limit=8, min_tps_to_cluster=1, adc_integral_cut=0):
    tps = cluster.get_tps()
    tps = tps[tps[:,4] < adc_max]
    tps = tps[tps[:,4] > adc_min]
    clusters = cluster_maker(tps, ticks_limit=ticks_limit, channel_limit=channel_limit, min_tps_to_cluster=min_tps_to_cluster, adc_integral_cut=adc_integral_cut)
    return clusters






print(f"Number of TPs: {len(tps)}")

plt.figure(figsize=(10,20))
# plt.scatter(tps[:,3], tps[:,2], c=tps[:,1], cmap='viridis', s=1)
plt.scatter(channel_to_z(tps[:,3]), tps[:,2], c=tps[:,4], cmap='viridis', s=1)
plt.colorbar()
plt.savefig(output_folder + '/tps.png')
plt.clf()


# colours = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
# create a list of colours for the clusters, at least 20
# colours = []
# for i in range(10):
#     colours.append((np.random.rand(), np.random.rand(), np.random.rand()))
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
tot_clusters = 0
# for i, cluster in enumerate(clusters):
#     if cluster.get_tps().shape[0] < 5:
#         continue
#     if cluster.get_tps()[:,3][0]//2560 == 2:
#         continue


#     plt.scatter(channel_to_z(cluster.get_tps()[:,3]), cluster.get_tps()[:,2], color=colours[i%len(colours)], s=1)
#     tot_clusters += 1

# plt.axvline(x=0, color='r', linestyle='--')
# plt.axvline(x=460, color='r', linestyle='--')
# # plt.xlim(9000, 10000)

# print(f"Number of clusters with more than 10 TPs: {tot_clusters}")
# plt.savefig(output_folder + '/cluster.png')
# plt.clf()

tot_clusters = 0
img_counter = 0
for i, cluster in enumerate(clusters):

    if cluster.get_tps().shape[0] < 10:
        continue
    if cluster.get_tps()[:,3][0]//2560 == 2:
        continue

    subclusters = extract_subclusters(cluster, adc_min=2000, adc_max=3500, ticks_limit=150, channel_limit=40, min_tps_to_cluster=1, adc_integral_cut=0)

    for j, subcluster in enumerate(subclusters):
        if subcluster.get_tps().shape[0] < 10:
            continue
        if not filter_by_position(subcluster, z_min=50, z_max=450, z_component_min=200, x_component_max=100):
            continue

        plt.scatter(channel_to_z(subcluster.get_tps()[:,3]), subcluster.get_tps()[:,2], c=subcluster.get_tps()[:,4], s=1)
        # plt.scatter(channel_to_z(subcluster.get_tps()[:,3]), subcluster.get_tps()[:,2], color=colours[(i+j)%len(colours)], s=1)
        tot_clusters += 1

        plt.axvline(x=0, color='r', linestyle='--')
        plt.axvline(x=460, color='r', linestyle='--')
        # print(f"Number of clusters with more than 10 TPs: {tot_clusters}")
        plt.colorbar()
        plt.savefig(output_folder + '/subcluster_' + str(img_counter) + '.png')
        print(f"Subcluster {img_counter}")
        img_counter += 1
        plt.clf()
# plt.axvline(x=0, color='r', linestyle='--')
# plt.axvline(x=460, color='r', linestyle='--')
# # plt.xlim(9000, 10000)
# # get the current y limits
# # ymin, ymax = plt.ylim()
# # print(ymin, ymax)
# # plt.ylim(ymax -20000, ymax-16000)
# print(f"Number of clusters with more than 10 TPs: {tot_clusters}")
# plt.savefig(output_folder + '/subcluster.png')
# plt.clf()
# tot_clusters = 0
# for i, cluster in enumerate(clusters_py):
#     if cluster.shape[0] < 5:
#         continue
#     if cluster[0,3]//2560 == 2:
#         continue


#     plt.scatter(channel_to_z(cluster[:,3]), cluster[:,2], color=colours[i%len(colours)], s=1)
#     tot_clusters += 1
# plt.axvline(x=0, color='r', linestyle='--')
# plt.axvline(x=460, color='r', linestyle='--')
# # plt.xlim(9000, 10000)
# print(f"Number of clusters with more than 10 TPs: {tot_clusters}")
# plt.savefig(output_folder + '/cluster_py.png')
# plt.clf()
