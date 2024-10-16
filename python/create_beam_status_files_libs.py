import numpy as np
import matplotlib.pyplot as plt
import time
# DAQ libraries
import daqdataformats
import detdataformats
import fddetdataformats
import trgdataformats
from hdf5libs import HDF5RawDataFile


def extract_record_lengths(filename):
    hdf5_file = HDF5RawDataFile(filename)
    # Get all records (TimeSlice)
    records = hdf5_file.get_all_record_ids()
    # Set the number of records to process
    print(f'Input file: {filename}')
    print(f'Total number of records: {len(records)}')

    # Get the size of a single TP in memory, useful to iterate over the TPs in the fragments
    tp_size = trgdataformats.TriggerPrimitive.sizeof()
    record_id=0

    n_tps_per_record = []
    min_time_frag = []
    size_of_TriggerPrimitive = trgdataformats.TriggerPrimitive.sizeof()
    for record in records[:]:
        time_start = time.time()
        all_tps = []
        print(f'Processing record {record}') # TODO can add verbose option and put all these in it
        # Get all data sources
        datasets = hdf5_file.get_fragment_dataset_paths(record)
        # These are arrays since there can be multiple fragments per record
        n_tps_per_record.append(0)
        min_time_frag.append(0)
        for dataset in datasets:
            # Get the fragment
            frag = hdf5_file.get_frag(dataset)
            # Get the number of TPs
            n_tps_remaining = (int(frag.get_data_size()/size_of_TriggerPrimitive))
            n_tps_per_record[-1] += n_tps_remaining
            time_frag = frag.get_window_begin()
            if time_frag < min_time_frag[-1] or min_time_frag[-1] == 0:
                min_time_frag[-1] = time_frag

    return (n_tps_per_record), (min_time_frag)
