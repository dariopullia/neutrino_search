#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>

#include "daq_tps_clustering_libs.h"
#include "cluster.h"


double channel_to_z(double channel) {
    const double apa_length_in_cm = 230;
    const double wire_pitch_in_cm_collection = 0.479;
    const double offset_between_apa_in_cm = 2.4;

    // Calculate z_apa_offset
    double z_apa_offset = (static_cast<int>(channel) / (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm);

    // Calculate z_channel_offset
    double z_channel_offset = ((static_cast<int>(channel) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection;

    // Final z value
    double z = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset;
    
    return z;
}


bool filter_by_position(cluster cl, int z_min, int z_max, int z_min_lenght, int x_max_length){
    std::vector<std::vector<double>> tps = cl.get_tps();
    if (int(tps[0][3]/2560) == 2){
        return false;
    }
    // double this_z_min, this_z_max, this_x_min, this_x_max;
    double this_z_min = channel_to_z(tps[0][3]), this_z_max = channel_to_z(tps[0][3]);
    double this_x_min = tps[0][2]*0.08, this_x_max = tps[0][2]*0.08;
    
    for (auto& tp : tps){
        double z = channel_to_z(tp[3]);
        double x = tp[2]*0.08;
        if (z < this_z_min){
            this_z_min = z;
        }
        if (z > this_z_max){
            this_z_max = z;
        }
        if (x < this_x_min){
            this_x_min = x;
        }
        if (x > this_x_max){
            this_x_max = x;
        }
    }

    if (this_z_max - this_z_min < z_min_lenght or this_z_min < z_min or this_x_max - this_x_min > x_max_length or this_z_max < z_max){
        return false;
    }
    std::cout << "z_min: " << this_z_min << " z_max: " << this_z_max << " x_min: " << this_x_min << " x_max: " << this_x_max << std::endl;
    return true;
}

std::vector<cluster> select_interesting_clusters(std::vector<cluster>& clusters, int sub_adc_min, int sub_adc_max, int sub_ticks_limit, int sub_channel_limit, int sub_min_tps_to_cluster, int sub_adc_integral_cut, int z_min, int z_max, int z_min_lenght, int x_max_length){
    std::vector<cluster> clusters_selected;
    for (auto& cl : clusters){
        std::vector<std::vector<double>> tps = cl.get_tps();
        // select tps with the right adc
        std::vector<std::vector<double>> tps_selected;
        for (auto& tp : tps){
            if (tp[4] >= sub_adc_min && tp[4] <= sub_adc_max){
                tps_selected.push_back(tp);
            }
        }
        if (tps_selected.size() == 0){
            continue;
        }
        std::vector<cluster> clusters_sub = cluster_maker(tps_selected, sub_ticks_limit, sub_channel_limit, sub_min_tps_to_cluster, sub_adc_integral_cut);
        if (clusters_sub.size() == 0){
            continue;
        }
        for (auto& cl_sub : clusters_sub){
            if (filter_by_position(cl_sub, z_min, z_max, z_min_lenght, x_max_length)){
                clusters_selected.push_back(cl);
                break;
            }
        }
    }

    return clusters_selected;
}



bool flat_charge_mean(std::vector<std::vector<double>>& tps, double max_displacement, double min_z_length){
    // compute the mean time weighted by the charge per z
    std::map<double, double> mean_time_per_z;
    std::map<double, double> total_charge_per_z;

    for (auto& tp : tps){
        double z = channel_to_z(tp[3]);
        double charge = tp[4];
        double time = tp[2]*0.08;
        if (mean_time_per_z.find(z) == mean_time_per_z.end()){
            mean_time_per_z[z] = 0;
        }
        mean_time_per_z[z] += charge*time;
        if (total_charge_per_z.find(z) == total_charge_per_z.end()){
            total_charge_per_z[z] = 0;
        }
        total_charge_per_z[z] += charge;
    }

    for (auto& [z, charge] : total_charge_per_z){
        mean_time_per_z[z] /= charge;
    }

    // get a sorted vector of z values
    std::vector<double> z_values;
    for (auto& [z, charge] : total_charge_per_z){
        z_values.push_back(z);
    }
    std::sort(z_values.begin(), z_values.end());

    // check if the mean time is flat
    double start_z = z_values[0];
    double current_min_z = 0;
    double current_max_z = 0;
    int n_tps_in_flat_region = 0;

    for (int i = 1; i < z_values.size(); i++){
        if (z_values[i] < 50){
            continue;
        }

        if (std::abs(mean_time_per_z[z_values[i]] - mean_time_per_z[z_values[i-1]]) < max_displacement){
            n_tps_in_flat_region++;
            continue;
        }
        else{
            if (z_values[i-1] - start_z > current_max_z - current_min_z){
                if (n_tps_in_flat_region < 30){
                    n_tps_in_flat_region = 0;
                    start_z = z_values[i];
                    continue;
                }

                current_max_z = z_values[i-1];
                current_min_z = start_z;
            }
            n_tps_in_flat_region = 0;
            start_z = z_values[i];
        }
    }

    if (z_values.back() - start_z > current_max_z - current_min_z){
        if (n_tps_in_flat_region < 30){
            n_tps_in_flat_region = 0;
            start_z = z_values.back();
        }
        current_max_z = z_values.back();
        current_min_z = start_z;
    }

    if (current_max_z - current_min_z < min_z_length){
        return false;
    }


    // check if y max - y min is less than 1.5*max_displacement
    double time_max = 0;
    double time_min = 0;
    double time_mean = 0;
    int n_tps = 0;
    for (auto& z : z_values){
        if (z < current_min_z or z > current_max_z){
            continue;
        }
        if (time_max == 0){
            time_max = mean_time_per_z[z];
            time_min = mean_time_per_z[z];
        }

        if (mean_time_per_z[z] > time_max){
            time_max = mean_time_per_z[z];
        }
        if (mean_time_per_z[z] < time_min){
            time_min = mean_time_per_z[z];
        }
        n_tps++;
        time_mean += mean_time_per_z[z];

    }
    time_mean /= n_tps;

    if (time_max - time_min > 1.5*max_displacement){
        return false;
    }

    // extract tps after the flat region
    int n_tps_after_flat_region = 0;
    int n_tps_in_first_50_cm = 0;
    for (auto& tp : tps){
        double z = channel_to_z(tp[3]);
        if (z > current_max_z){
            // append the tp to the new vector if time within 1.5*max_displacement from the mean time
            if (std::abs(tp[2]*0.08 - time_mean) < 1.5*max_displacement){
                n_tps_after_flat_region++;
            }
        }
        if (z < 50){
            if (std::abs(tp[2]*0.08 - time_mean) < 1.5*max_displacement){
                n_tps_in_first_50_cm++;
            }
        }

    }

    double space_left_after_flat_region = 460 - current_max_z;

    if (n_tps_after_flat_region < space_left_after_flat_region*0.5){
        return false;
    }

    if (n_tps_in_first_50_cm > 20){
        return false;
    }
    
    return true;
}

bool tps_pass_all_filters(std::vector<std::vector<double>>& tps){
    // // -------------------------
    // double max_displacement = 10;
    // double min_z_length = 100;
    // // -------------------------
    // // check if the charge is flat
    // if (!flat_charge_mean(tps, max_displacement, min_z_length)){
    //     return false;
    // }
    return true;
}


std::vector<cluster> neutrino_explosion_finder(std::vector<std::vector<double>>& all_tps, int length_limit, int time_limit, int adc_integral_cut){
    // sort the tps by time
    std::sort(all_tps.begin(), all_tps.end(), [](std::vector<double>& a, std::vector<double>& b){
        return a[2] < b[2];
    });
    std::vector<cluster> clusters;
    int section_1_start_idx = 0;
    int section_2_start_idx = 0;
    double section_1_adc_integral = 0;
    double section_2_adc_integral = 0;
    double section_1_start_time = all_tps[0][2];  
    double section_2_start_time = all_tps[0][2] + 1500;

    for (int tp_index = 0; tp_index < all_tps.size(); tp_index++){
        std::vector<double> tp = all_tps[tp_index];
        // skip apa 2
        if (int(tp[3]/2560) == 2){
            continue;
        }

        if (tp[2] - section_1_start_time > time_limit){
            if (section_1_adc_integral > adc_integral_cut){
                std::vector<std::vector<double>> tps;
                for (int i = section_1_start_idx; i < tp_index; i++){
                    if (int(all_tps[i][3]/2560) == 2){
                        continue;
                    }
                    tps.push_back(all_tps[i]);
                    // std::cout << "Diff time: " << all_tps[i][2] - all_tps[section_1_start_idx][2] << std::endl;
                }
                if (tps_pass_all_filters(tps)){
                    clusters.push_back(cluster(tps));
                }
            }
            section_1_adc_integral = 0;
            section_1_start_idx = tp_index;
            section_1_start_time = tp[2];            
        }
        else{
            section_1_adc_integral += tp[4];
        }
    }

    return clusters;

}

