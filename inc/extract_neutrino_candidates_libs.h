#ifndef EXTRACT_NEUTRINO_CANDIDATES_LIBS_H
#define EXTRACT_NEUTRINO_CANDIDATES_LIBS_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>

#include "daq_tps_clustering_libs.h"
#include "cluster.h"

std::vector<cluster> select_interesting_clusters(std::vector<cluster>& clusters, int sub_adc_min, int sub_adc_max, int sub_ticks_limit, int sub_channel_limit, int sub_min_tps_to_cluster, int sub_adc_integral_cut, int z_min, int z_max, int z_min_lenght, int x_max_length);

bool filter_by_position(cluster cl, int z_min, int z_max, int z_min_lenght, int x_max_length);

double channel_to_z(double channel);

std::vector<cluster> neutrino_explosion_finder(std::vector<std::vector<double>>& all_tps, int length_limit, int time_limit, int adc_integral_cut);

bool tps_pass_all_filters(std::vector<std::vector<double>>& tps);

bool flat_charge_mean(std::vector<std::vector<double>>& tps, double max_displacement, double min_z_length);

#endif // EXTRACT_NEUTRINO_CANDIDATES_LIBS_H