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

job_number = input_json_file.split('_')[-1].replace('.json', '')

filename = input_json_file.split('/')[-1].replace('.json', '')


# Read input json
with open(input_json_file) as f:
    input_json = json.load(f)

input_file = input_json['input_file']

output_folder = input_json['output_folder']
adc_integral_cut = input_json['adc_integral_cut']
time_limit = input_json['time_limit']
beam_status_file = input_json['beam_status_file']
beam_on = input_json['beam_on']
monitor_threshold = input_json['monitor_threshold']
ask_for_close_beam = input_json['ask_for_close_beam']
time_interval= input_json["time_interval"]
print("Input file: ", input_file)
print("Output folder: ", output_folder)
print("ADC integral cut: ", adc_integral_cut)
print("Time limit: ", time_limit)
print("Beam status file: ", beam_status_file)
print("Beam on: ", beam_on)
print("Monitor threshold: ", monitor_threshold)
print("Ask for close beam: ", ask_for_close_beam)
print("Time interval: ", time_interval)


if beam_on:
    output_folder = output_folder + "/beam_on/"
else:
    output_folder = output_folder + "/beam_off/"

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
passing_cosmics = 0
not_passing_cosmics = 0
passing_neutrinos = 0
not_passing_neutrinos = 0


for input_hdf5_file in list_of_files:
    records_adcs, records_ids = from_hdf5_to_wf(
        input_hdf5_file, 
        "PD2HDChannelMap", 
        "HD_TPC", 
        channel_limit_mins=[4160, 9280], 
        channel_limit_maxs=[4640, 9760], 
        time_interval=time_interval,
        monitor_charge=monitor_charge, 
        monitor_time=monitor_time, 
        beam_on=beam_on, 
        monitor_threshold=monitor_threshold, 
        ask_for_close_beam=ask_for_close_beam
    )
    # if the record id is in the list, remove from the list
    for i in range(len(records_adcs)):
        adcs = records_adcs[i]
        if pass_all_filters(adcs):
            passing_cosmics += 1
            do_clean_plot_wf(adcs, records_ids[i], filename=input_hdf5_file.split('/')[-1].replace('.hdf5', ''), output_folder=output_folder+'/passing/')
        else:
            not_passing_cosmics += 1
            do_plot_wf(adcs, records_ids[i], filename=input_hdf5_file.split('/')[-1].replace('.hdf5', ''), output_folder=output_folder+'/not_passing/')

    # break

print("Passing neutrinos: ", passing_neutrinos)
print("Not passing neutrinos: ", not_passing_neutrinos)
print("Passing cosmics: ", passing_cosmics)
print("Not passing cosmics: ", not_passing_cosmics)

if not os.path.exists(output_folder+'/summary/passing/'):
    os.makedirs(output_folder+'/summary/passing/')  
if not os.path.exists(output_folder+'/summary/not_passing/'):
    os.makedirs(output_folder+'/summary/not_passing/')

with open(output_folder+'/summary/passing/'+str(job_number)+'.txt', "w") as f:
    f.write(str(passing_cosmics))
with open(output_folder+'/summary/not_passing/'+str(job_number)+'.txt', "w") as f:
    f.write(str(not_passing_cosmics))
