#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <climits>

// include root libraries
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMatrixD.h"

#include "daq_tps_clustering_libs.h"
#include "cluster.h"


std::vector<std::vector<double>> file_reader(std::vector<std::string> filenames, int plane, int max_tps_per_filename, int min_tot_per_tp, int min_integral_per_tp, int n_skip_tps) {
    std::vector<std::vector<double>> tps;
    std::string line;
    int skip_counter = 0;
    int n_events_offset = 0;
    int file_idx = 0;
    int excl_tot = 0;
    int this_plane=3;
    for (auto& filename : filenames) {
        // std::cout << filename << std::endl;
        std::ifstream infile(filename);
        // read and save the TPs
        while (std::getline(infile, line)) {

            std::istringstream iss(line);
            std::vector<double> tp;
            double val;
            while (iss >> val) {
                tp.push_back(val);
            }
            tp.push_back(file_idx);
            if (skip_counter < n_skip_tps) {
                skip_counter++;
                continue;
            }
            if (tps.size() >= max_tps_per_filename) {
                break;
            }
            if (tp[variables_to_index["time_over_threshold"]]/32 <= min_tot_per_tp) {
                excl_tot++;
                continue;
            }
            if (tp[variables_to_index["adc_integral"]] < min_integral_per_tp) {
                continue;
            }

            if (int(tp[variables_to_index["channel"]])%2560 < 800) {
                this_plane = 0;
            }
            else if (int(tp[variables_to_index["channel"]])%2560 < 1600) {
                this_plane = 1;
            }
            else {
                this_plane = 2;
            }

            if ( this_plane != plane) {
                continue;
            }

            tp[variables_to_index["time_start"]] = tp[variables_to_index["time_start"]]/32;
            tp[variables_to_index["time_over_threshold"]] = tp[variables_to_index["time_over_threshold"]]/32;
            tp[variables_to_index["time_peak"]] = tp[variables_to_index["time_peak"]]/32;
            tps.push_back(tp);
        }
    }
    // sort the TPs by time
    std::sort(tps.begin(), tps.end(), [](const std::vector<double>& a, const std::vector<double>& b) {
        return a[0] < b[0];
     });
    std::cout << "Number of TPs excluded due to TOT: " << excl_tot << std::endl;
    return tps;
}

bool channel_condition_with_pbc(double ch1, double ch2, int channel_limit) {
    if (int(ch1/2560) != int(ch2/2560)) {
        return false;
    }

    double diff = std::abs(ch1 - ch2);
    int n_chan;
    int mod_ch1 = int(ch1) % 2560;
    if (mod_ch1 >= 0 and mod_ch1 <800) {
        n_chan = 800;

    }
    else if (mod_ch1 >= 800 and mod_ch1 < 1600) {
        n_chan = 800;

    }
    else {
        n_chan = 960;
        if (diff <= channel_limit) {
            // std::cout << "chan" << std::endl;   
            return true;
        } else{
            return false;
        }

    }


    if (diff <= channel_limit) {
        // std::cout << "chan" << std::endl;   
        return true;
    }
    else if (diff >= n_chan - channel_limit) {
        std::cout << "chan with PCB!" << std::endl;
        std::cout << int(ch1) << " " << int(ch2) << std::endl;
        // std::cout << ch1 << " " << ch2 << std::endl;
        return true;
    }
    return false;
}

bool channel_condition_with_pos(double ch1, double ch2, int channel_limit) {
    
    if (int(ch1/2560)%2 != int(ch2/2560)%2) {
        return false;
    }

    const double apa_lenght_in_cm = 230;
    const double wire_pitch_in_cm_collection = 0.479;
    const double offset_between_apa_in_cm = 2.4;
    float z_apa_offset, z_channel_offset;

    z_apa_offset = int(ch1 / (2560*2)) * (apa_lenght_in_cm + offset_between_apa_in_cm);
    z_channel_offset = ((int(ch1) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection;
    float z1 = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset;

    z_apa_offset = int(ch2 / (2560*2)) * (apa_lenght_in_cm + offset_between_apa_in_cm);
    z_channel_offset = ((int(ch2) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection;
    float z2 = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset;

    if (std::abs(z1 - z2) <= channel_limit * wire_pitch_in_cm_collection) {
        return true;
    }
    return false;
}

std::vector<cluster> cluster_maker(std::vector<std::vector<double>>& all_tps, int ticks_limit, int channel_limit, int min_tps_to_cluster, int adc_integral_cut) {
    std::vector<std::vector<std::vector<double>>> buffer;
    std::vector<cluster> clusters;
    for (auto& tp : all_tps) {
        if (buffer.size() == 0) {
            std::vector<std::vector<double>> temp;
            temp.push_back(tp);
            buffer.push_back(temp);
        }
        else {
            std::vector<std::vector<std::vector<double>>> buffer_copy = buffer;
            buffer.clear();
            bool appended = false;
            int idx = 0;
            int idx_appended;
            for (auto& candidate : buffer_copy) {
                // get a the max containing the times of the TPs in the candidate
                double max_time = 0;
                for (auto& tp2 : candidate) {
                    max_time = std::max(max_time, tp2[0] + tp2[1]);
                }
                bool time_cond = (tp[0] - max_time) <= ticks_limit;
                if (time_cond) {
                    bool chan_cond = false;
                    for (auto& tp2 : candidate) {
                        if (channel_condition_with_pos(tp[3], tp2[3], channel_limit)) {
                        // if (channel_condition_with_pbc(tp[3], tp2[3], channel_limit)) {
                        // if (std::abs(tp[3] - tp2[3]) <= channel_limit) {
                            if (tp[0]-(tp2[0]+tp2[1]) <=ticks_limit) {
                                chan_cond = true;
                                break;
                            }
                        }
                    }
                    if (chan_cond) {
                        // std::cout << "chan" << std::endl;
                        if (!appended) {
                            candidate.push_back(tp);
                            buffer.push_back(candidate);
                            appended = true;
                            idx_appended = idx;
                            // std::cout << "appended" << std::endl;
                        }
                        else {
                            for (auto& tp2 : candidate) {
                                buffer[idx_appended].push_back(tp2);
                            }
                        }
                    }
                    else {
                        buffer.push_back(candidate);
                        ++idx;
                    }
                }
                else {
                    if (candidate.size() >= min_tps_to_cluster) {
                        int adc_integral = 0;
                        for (auto& tp2 : candidate) {
                            adc_integral += tp2[4];
                        }
                        if (adc_integral > adc_integral_cut) {
                            cluster g(candidate);
                            clusters.push_back(g);
                        }
                    }
                }
            }
            if (!appended) {
                std::vector<std::vector<double>> temp;
                temp.push_back(tp);
                buffer.push_back(temp);
            }
        }
    }
    if (buffer.size() > 0) {
        for (auto& candidate : buffer) {
            if (candidate.size() >= min_tps_to_cluster) {
                int adc_integral = 0;
                for (auto& tp : candidate) {
                    adc_integral += tp[4];
                }
                if (adc_integral > adc_integral_cut) {
                    cluster g(candidate);
                    clusters.push_back(g);
                }
            }
        }
    }

    return clusters;
}

// std::vector<cluster> cluster_maker(std::vector<std::vector<double>>& all_tps, int ticks_limit, int channel_limit, int min_tps_to_cluster, int adc_integral_cut) {
//     std::vector<std::vector<std::vector<double>>> buffer;
//     std::vector<cluster> clusters;
//     for (auto& tp : all_tps) {
//         if (buffer.size() == 0) {
//             std::vector<std::vector<double>> temp;
//             temp.push_back(tp);
//             buffer.push_back(temp);
//         }
//         else {
//             std::vector<std::vector<std::vector<double>>> buffer_copy = buffer;
//             buffer.clear();
//             bool appended = false;
//             int idx = 0;
//             int idx_appended;
//             for (auto& candidate : buffer_copy) {
//                 // get a the max containing the times of the TPs in the candidate
//                 double max_time = 0;
//                 for (auto& tp2 : candidate) {
//                     max_time = std::max(max_time, tp2[0] + tp2[1]);
//                 }
//                 bool time_cond = (tp[0] - max_time) <= ticks_limit;
//                 if (time_cond) {
//                     bool chan_cond = false;
//                     for (auto& tp2 : candidate) {
//                         if (channel_condition_with_pbc(tp[3], tp2[3], channel_limit)) {
//                         // if (std::abs(tp[3] - tp2[3]) <= channel_limit) {
//                             chan_cond = true;
//                             break;
//                         }
//                     }
//                     if (chan_cond) {
//                         // std::cout << "chan" << std::endl;
//                         if (!appended) {
//                             candidate.push_back(tp);
//                             buffer.push_back(candidate);
//                             appended = true;
//                             idx_appended = idx;
//                             // std::cout << "appended" << std::endl;
//                         }
//                         else {
//                             for (auto& tp2 : candidate) {
//                                 buffer[idx_appended].push_back(tp2);
//                             }
//                         }
//                     }
//                     else {
//                         buffer.push_back(candidate);
//                         ++idx;
//                     }
//                 }
//                 else {
//                     if (candidate.size() >= min_tps_to_cluster) {
//                         int adc_integral = 0;
//                         for (auto& tp2 : candidate) {
//                             adc_integral += tp2[4];
//                         }
//                         if (adc_integral > adc_integral_cut) {
//                             cluster g(candidate);
//                             clusters.push_back(g);
//                         }
//                     }
//                 }
//             }
//             if (!appended) {
//                 std::vector<std::vector<double>> temp;
//                 temp.push_back(tp);
//                 buffer.push_back(temp);
//             }
//         }
//     }
//     if (buffer.size() > 0) {
//         for (auto& candidate : buffer) {
//             if (candidate.size() >= min_tps_to_cluster) {
//                 int adc_integral = 0;
//                 for (auto& tp : candidate) {
//                     adc_integral += tp[4];
//                 }
//                 if (adc_integral > adc_integral_cut) {
//                     cluster g(candidate);
//                     clusters.push_back(g);
//                 }
//             }
//         }
//     }

//     return clusters;
// }
