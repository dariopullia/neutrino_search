import numpy as np
import numba as nb
import sys

sys.path.append('../python/') 
from cluster import *


def channel_condition_with_pos(ch1, ch2, channel_limit):
    if int(ch1 / 2560) % 2 != int(ch2 / 2560) % 2:
        return False

    apa_length_in_cm = 230
    wire_pitch_in_cm_collection = 0.479
    offset_between_apa_in_cm = 2.4

    z_apa_offset = int(ch1 / (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm)
    z_channel_offset = ((int(ch1) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection
    z1 = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset

    z_apa_offset = int(ch2 / (2560 * 2)) * (apa_length_in_cm + offset_between_apa_in_cm)
    z_channel_offset = ((int(ch2) % 2560 - 1600) % 480) * wire_pitch_in_cm_collection
    z2 = wire_pitch_in_cm_collection + z_apa_offset + z_channel_offset

    if abs(z1 - z2) <= channel_limit * wire_pitch_in_cm_collection:
        return True
    return False



# @nb.njit
def cluster_maker(all_tps, ticks_limit, channel_limit, min_tps_to_cluster, adc_integral_cut):
    buffer = []
    clusters = []

    for tp in all_tps:
        if len(buffer) == 0:
            temp = [tp]
            buffer.append(temp)
        else:
            buffer_copy = buffer.copy()
            buffer.clear()
            appended = False
            idx = 0

            for candidate in buffer_copy:
                max_time = 0
                for tp2 in candidate:
                    max_time = max(max_time, tp2[0] + tp2[1])

                time_cond = (tp[0] - max_time) <= ticks_limit
                if time_cond:
                    chan_cond = False
                    for tp2 in candidate:
                        if channel_condition_with_pos(tp[3], tp2[3], channel_limit):
                            if abs(tp[2] - tp2[2]) <= ticks_limit:
                            # if tp[0] - (tp2[0] + tp2[1]) <= ticks_limit:
                                chan_cond = True
                                break

                    if chan_cond:
                        if not appended:
                            candidate.append(tp)
                            buffer.append(candidate)
                            appended = True
                            idx_appended = idx
                        else:
                            for tp2 in candidate:
                                buffer[idx_appended].append(tp2)
                    else:
                        buffer.append(candidate)
                        idx += 1
                else:
                    if len(candidate) >= min_tps_to_cluster:
                        adc_integral = 0
                        for tp2 in candidate:
                            adc_integral += tp2[4]
                        if adc_integral > adc_integral_cut:
                            clusters.append(cluster(np.array(candidate)))

            if not appended:
                temp = [tp]
                buffer.append(temp)

    if len(buffer) > 0:
        for candidate in buffer:
            if len(candidate) >= min_tps_to_cluster:
                adc_integral = 0
                for tp in candidate:
                    adc_integral += tp[4]
                if adc_integral > adc_integral_cut:
                    clusters.append(cluster(np.array(candidate)))
 
    return clusters
