#ifndef POSITION_CALCULATOR_H
#define POSITION_CALCULATOR_H

#include <vector>
#include <cmath>
#include "cluster.h"

const int EVENTS_OFFSET = 5000;
const double apa_lenght_in_cm = 230;
const double wire_pitch_in_cm_collection = 0.479;
const double wire_pitch_in_cm_induction_diagonal = 0.4669;
const double apa_angle = (90 - 35.7);
const double wire_pitch_in_cm_induction = wire_pitch_in_cm_induction_diagonal / sin(apa_angle * M_PI / 180);
const double offset_between_apa_in_cm = 2.4;
const double apa_height_in_cm = 598.4;
const double time_tick_in_cm = 0.0805;
const double apa_width_in_cm = 4.7;
const int backtracker_error_margin = 4;
const double apa_angular_coeff = tan(apa_angle * M_PI / 180);

std::vector<float> calculate_position(std::vector<double> tp);
std::vector<std::vector<float>> validate_position_calculation(std::vector<std::vector<double>> tps);
float distance(cluster cluster1, cluster cluster2);
float eval_y_knowing_z_U_plane(std::vector<std::vector<double>> tps, float z, float x_sign);
float eval_y_knowing_z_V_plane(std::vector<std::vector<double>> tps, float z, float x_sign);

#endif



