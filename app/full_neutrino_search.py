import os
import numpy as np
import argparse
import json
import sys
import matplotlib.pyplot as plt

sys.path.append('../python/') 
import full_neutrino_search_libs as fns

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

if not os.path.exists(output_folder):   
    os.makedirs(output_folder)
if not os.path.exists(output_folder + "/over_22/"):   
    os.makedirs(output_folder + "/over_22/")
if not os.path.exists(output_folder + "/20/"):   
    os.makedirs(output_folder + "/20/")
if not os.path.exists(output_folder + "/17/"):
    os.makedirs(output_folder + "/17/")
if not os.path.exists(output_folder + "/15/"):
    os.makedirs(output_folder + "/15/")
if not os.path.exists(output_folder + "/under_15/"):
    os.makedirs(output_folder + "/under_15/")
if not os.path.exists(output_folder + '/gs'):
    os.makedirs(output_folder + '/gs')
if not os.path.exists(output_folder + '/vc'):
    os.makedirs(output_folder + '/vc')


section = fns.search_neutrinos(input_file, time_limit, adc_integral_cut)

print(f'Number of sections: {len(section)}')

for index, s in enumerate(section):
    
    if fns.is_groundshake(s):
        fns.do_plots(s, index, filename, output_folder + '/gs')
        continue
    if fns.is_vertical_cosmic(s):
        fns.do_plots(s, index, filename, output_folder + '/vc')
        continue
    if np.sum(s[:,4]) > 22E6:
        fns.do_plots(s, index, filename, output_folder + "/over_22/")
    elif np.sum(s[:,4]) > 20E6:
        fns.do_plots(s, index, filename, output_folder + "/20/")
    elif np.sum(s[:,4]) > 17E6:
        fns.do_plots(s, index, filename, output_folder + "/17/")
    elif np.sum(s[:,4]) > 15E6:
        fns.do_plots(s, index, filename, output_folder + "/15/")
    else:
        fns.do_plots(s, index, filename, output_folder + "/under_15/")


