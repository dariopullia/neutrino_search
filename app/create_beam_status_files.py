import numpy as np
import matplotlib.pyplot as plt
import time
import json
import argparse
import sys
import os

# DAQ libraries
import daqdataformats
import detdataformats
import fddetdataformats
import trgdataformats
from hdf5libs import HDF5RawDataFile

sys.path.append('../python/') 
from create_beam_status_files_libs import *

parser = argparse.ArgumentParser(description='Tranforms Trigger Primitives to images.')
parser.add_argument('--input_json', type=str, help='Input json file')
args = parser.parse_args()

input_json_file = args.input_json

# Read input json
with open(input_json_file) as f:
    input_json = json.load(f)

input_file = input_json['input_file']
output_folder = input_json['output_folder']
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


list_of_hdf5_paths = np.loadtxt(input_file, dtype=str)
name_out_file = output_folder + 'beam_status_' + input_file.split('/')[-1].split('.')[0] +'.csv'

hdf5_path = list_of_hdf5_paths[0]
n_tps_per_record, first_tp_time_start = extract_record_lengths(hdf5_path)
matrix = np.concatenate((np.array(n_tps_per_record).reshape(-1,1), np.array(first_tp_time_start).reshape(-1,1)), axis=1)
name_vector = np.array([hdf5_path]*len(n_tps_per_record)).reshape(-1,1)
full_matrix = np.concatenate((name_vector, matrix), axis=1)
print(full_matrix)

for hdf5_path in list_of_hdf5_paths[1:]:
    n_tps_per_record, first_tp_time_start = extract_record_lengths(hdf5_path)
    matrix = np.concatenate((np.array(n_tps_per_record).reshape(-1,1), np.array(first_tp_time_start).reshape(-1,1)), axis=1)
    name_vector = np.array([hdf5_path]*len(n_tps_per_record)).reshape(-1,1)
    matrix = np.concatenate((name_vector, matrix), axis=1)
    full_matrix = np.concatenate((full_matrix, matrix), axis=0)
    print(full_matrix)
np.savetxt(name_out_file, full_matrix, fmt='%s', delimiter=',')
