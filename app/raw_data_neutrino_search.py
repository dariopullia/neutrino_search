import os
import numpy as np
import argparse
import json
import sys
import matplotlib.pyplot as plt
# DAQ libraries
from hdf5libs import HDF5RawDataFile
import detdataformats
import fddetdataformats
from rawdatautils.unpack.wibeth import *
from rawdatautils.utilities.wibeth import *
import rawdatautils.utilities 
import detchannelmaps

sys.path.append('../python/') 
from raw_data_neutrino_search_libs import *



parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('--input_json', type=str, help='Input json file')
args = parser.parse_args()

input_json_file = args.input_json

filename = input_json_file.split('/')[-1].replace('.json', '')


# Read input json
with open(input_json_file) as f:
    input_json = json.load(f)

input_file = input_json['input_file']

output_folder = input_json['output_folder']
adc_integral_cut = input_json['adc_integral_cut']
time_limit = input_json['time_limit']
beam_status_file = input_json['beam_status_file']

monitor_matrix = np.loadtxt(beam_status_file, delimiter=",", dtype=str)
monitor_name = monitor_matrix[:,0]
monitor_charge = monitor_matrix[:,1].astype(float)
monitor_time = monitor_matrix[:,2].astype(float)

print(monitor_charge.shape)
print(monitor_time.shape)




if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# read the list of files
list_of_files = np.loadtxt(input_file, dtype=str)

records_adcs = []
for input_hdf5_file in list_of_files:
    records_adcs = from_hdf5_to_wf(input_hdf5_file, "PD2HDChannelMap", "HD_TPC", channel_limit_mins=[4160, 9280], channel_limit_maxs=[4640, 9760], monitor_charge=monitor_charge, monitor_time=monitor_time)
    for i in range(len(records_adcs)):
        print(records_adcs[i].shape)
        do_plot_wf(records_adcs[i], i, filename=input_hdf5_file.split('/')[-1].replace('.hdf5', ''), output_folder=output_folder)
