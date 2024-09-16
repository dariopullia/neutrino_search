#include <vector>
#include <iostream>
#include "position_calculator.h"
#include "cluster.h"

std::vector<float> calculate_position(std::vector<double> tp) { // only works for plane 2
    float z;
    if (tp[variables_to_index["view"]] == 2) {
        float z_apa_offset = int(tp[variables_to_index["channel"]]) / (2560*2) * (apa_lenght_in_cm + offset_between_apa_in_cm);
        float z_channel_offset = ((int(tp[variables_to_index["channel"]]) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection;
        z = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset;
    } else {
        z = 0;
    }

    float y = 0;
    float x_signs = ((int(tp[variables_to_index["channel"]]) % 2560-2080<0) ? -1 : 1);
    float x = ((int(tp[variables_to_index["time_peak"]]) % EVENTS_OFFSET )* time_tick_in_cm + apa_width_in_cm/2) * x_signs;

    return {x, y, z};
}

std::vector<std::vector<float>> validate_position_calculation(std::vector<std::vector<double>> tps) {
    std::vector<std::vector<float>> positions;
    for (auto& tp : tps) {
        if (tp[0] == 2 and tp[variables_to_index["true_x"]] and tp[variables_to_index["true_y"]] and tp[variables_to_index["true_z"]]) {
            positions.push_back({int(tp[variables_to_index["true_x"]]), int(tp[variables_to_index["true_y"]]), int(tp[variables_to_index["true_z"]]), calculate_position(tp)[0], calculate_position(tp)[1], calculate_position(tp)[2]});
        }
    }
    return positions;
}

float distance(cluster cluster1, cluster cluster2) {
    float x1 = cluster1.get_reco_pos()[0];
    float y1 = cluster1.get_reco_pos()[1];
    float z1 = cluster1.get_reco_pos()[2];
    float x2 = cluster2.get_reco_pos()[0];
    float y2 = cluster2.get_reco_pos()[1];
    float z2 = cluster2.get_reco_pos()[2];
    return sqrt(pow(x1 - x2, 2) + pow(y1 - y2, 2) + pow(z1 - z2, 2));
}

float eval_y_knowing_z_U_plane(std::vector<std::vector<double>> tps, float z, float x_sign) {
    z = z - int(tps[0][3]) / (2560*2) * (apa_lenght_in_cm + offset_between_apa_in_cm);
    float ordinate;
    std::vector<float> Y_pred;
    for (auto& tp : tps) {
        if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 == 0) {
            if (x_sign < 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 < 400) {
                    if (z > (int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction + backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = apa_lenght_in_cm - z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = ((int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction - z) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 > 399) {
                    ordinate = (apa_lenght_in_cm - z + (int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            } else if (x_sign > 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 > 399) {
                    if (z < (799 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction - backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = (z - (799 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 < 400) {
                    ordinate = (z + (int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            }
        } else if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 == 1) {
            if (x_sign < 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 < 400) {
                    if (z < (399 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction - backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = (z - (399 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 > 399) {
                    ordinate = (z + (int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            } else if (x_sign > 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 > 399) {
                    if (z > (int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction + backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = apa_lenght_in_cm - z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = ((int(tp[variables_to_index["channel"]]) % 2560 - 400) * wire_pitch_in_cm_induction - z) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 < 400) {
                    ordinate = (apa_lenght_in_cm - z + (int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            }
        }
        // ordinate = (ordinate) - apa_height_in_cm if (tp[idx['channel']] / 2560) % 2 < 1 else apa_height_in_cm - (ordinate);
        if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 < 1) {
            ordinate = (ordinate) - apa_height_in_cm;
        } else {
            ordinate = apa_height_in_cm - (ordinate);
        }
        Y_pred.push_back(ordinate);
    }
    // mean
    float Y_pred_mean = 0;
    for (auto& y : Y_pred) {
        Y_pred_mean += y;
    }
    Y_pred_mean = Y_pred_mean / Y_pred.size();

    return Y_pred_mean;
}


float eval_y_knowing_z_V_plane(std::vector<std::vector<double>> tps, float z, float x_sign) {
    z = z - int(tps[0][3]) / (2560*2) * (apa_lenght_in_cm + offset_between_apa_in_cm);
    float ordinate;
    std::vector<float> Y_pred;
    for (auto& tp : tps) {
        if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 == 0) {
            if (x_sign < 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 < 1200) {
                    if (z < (1199 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction - backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = (z - (1199 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 > 1199) {
                    ordinate = (z + (int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            } else if (x_sign > 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 > 1199) {
                    if (z > (int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction + backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = apa_lenght_in_cm - z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = ((int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction - z) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 < 1200) {
                    ordinate = (apa_lenght_in_cm - z + (int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            } 
        } else if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 == 1) {
            if (x_sign < 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 < 1200) {
                    if (z > (int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction + backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = apa_lenght_in_cm - z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = ((int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction - z) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 > 1199) {
                    ordinate = (apa_lenght_in_cm - z + (int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            } else if (x_sign > 0) {
                if (int(tp[variables_to_index["channel"]]) % 2560 > 1199) {
                    if (z < (1599 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction - backtracker_error_margin) {
                        float until_turn = (int(tp[variables_to_index["channel"]]) % 2560 - 1200) * wire_pitch_in_cm_induction;
                        float all_the_way_behind = apa_lenght_in_cm;
                        float the_last_piece = z;
                        ordinate = (until_turn + all_the_way_behind + the_last_piece) * apa_angular_coeff;
                    } else {
                        ordinate = (z - (1599 - int(tp[variables_to_index["channel"]]) % 2560) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                    }
                } else if (int(tp[variables_to_index["channel"]]) % 2560 < 1200) {
                    ordinate = (z + (int(tp[variables_to_index["channel"]]) % 2560 - 800) * wire_pitch_in_cm_induction) * apa_angular_coeff;
                }
            }
        }
        if ((int(tp[variables_to_index["channel"]]) / 2560) % 2 < 1) {
            ordinate = (ordinate) - apa_height_in_cm;
        } else {
            ordinate = apa_height_in_cm - (ordinate);
        }

        Y_pred.push_back(ordinate);
    }
    // mean
    float Y_pred_mean = 0;
    for (auto& y : Y_pred) {
        Y_pred_mean += y;
    }
    Y_pred_mean = Y_pred_mean / Y_pred.size();

    return Y_pred_mean;
}




