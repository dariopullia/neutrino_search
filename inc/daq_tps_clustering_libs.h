// this is a .h file that contains the functions to read tps from files and create clusters
#ifndef DAQ_TPS_CLUSTERING_LIBS
#define DAQ_TPS_CLUSTERING_LIBS

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <climits>


#include "cluster.h"

// read the tps from the files and save them in a vector
std::vector<std::vector<double>> file_reader(std::vector<std::string> filenames, int plane=2, int max_tps_per_filename = INT_MAX, int min_tot_per_tp=0, int min_integral_per_tp=0, int n_skip_tps=0);


// create the clusters from the tps

bool channel_condition_with_pbc(double ch1, double ch2, int channel_limit);
std::vector<cluster> cluster_maker(std::vector<std::vector<double>>& all_tps, int ticks_limit=3, int channel_limit=1, int min_tps_to_cluster=1, int adc_integral_cut=0);

#endif // DAQ_TPS_CLUSTERING_LIBS
