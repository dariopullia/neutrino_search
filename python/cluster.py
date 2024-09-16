import numpy as np
import sys
import ROOT

class cluster: # Only good for HD geometry 1x2x6
    def __init__(self, tps):
        self.tps_ = tps
        self.reco_pos_ = None
        self.n_tps_ = tps.shape[0]
    def __len__(self):
        return self.n_tps_
    def __str__(self):
        return str(self.tps_)
    def __repr__(self):
        return str(self.tps_)
    def get_tps(self):
        return self.tps_
    def get_reco_pos(self):
        return self.reco_pos_
    def set_variables(self, reco_pos = np.array([0,0,0])):
        self.reco_pos_ = reco_pos


def read_root_file_to_clusters(filename):
        # Create a structured array with column names
    dt = np.dtype([('time_start', float), 
                    ('time_over_threshold', float),
                    ('time_peak', float),
                    ('channel', int),
                    ('adc_integral', int),
                    ('adc_peak', int),
                    ('detid', int),
                    ('type', int),
                    ('algorithm', int),
                    ('version', int),
                    ('flag', int)])

    file = ROOT.TFile.Open(filename, "READ")
    clusters = []
    for i in file.GetListOfKeys():
        # get the TTree
        tree = file.Get(i.GetName())
        print(f"Tree name: {tree.GetName()}")
        tree.Print()
        counter = 0
        for entry in tree:
            counter += 1
            if counter % 1000 == 0:
                print(f"cluster {counter}")
            nrows = entry.nrows
            # m = np.empty((nrows), dtype=dt)
            # for j in range(nrows):
            #     m[j]["time_start"] = entry.matrix[j][0]
            #     m[j]["time_over_threshold"] = entry.matrix[j][1]
            #     m[j]["time_peak"] = entry.matrix[j][2]
            #     m[j]["channel"] = entry.matrix[j][3]
            #     m[j]["adc_integral"] = entry.matrix[j][4]
            #     m[j]["adc_peak"] = entry.matrix[j][5]
            #     m[j]["detid"] = entry.matrix[j][6]
            #     m[j]["type"] = entry.matrix[j][7]
            #     m[j]["algorithm"] = entry.matrix[j][8]
            #     m[j]["version"] = entry.matrix[j][9]
            #     m[j]["flag"] = entry.matrix[j][10]

            m = np.empty((nrows, 11))
            for j in range(nrows):
                m[j, 0] = entry.matrix[j][0]
                m[j, 1] = entry.matrix[j][1]
                m[j, 2] = entry.matrix[j][2]
                m[j, 3] = entry.matrix[j][3]
                m[j, 4] = entry.matrix[j][4]
                m[j, 5] = entry.matrix[j][5]
                m[j, 6] = entry.matrix[j][6]
                m[j, 7] = entry.matrix[j][7]
                m[j, 8] = entry.matrix[j][8]
                m[j, 9] = entry.matrix[j][9]
                m[j, 10] = entry.matrix[j][10]
            this_cluster = cluster(m)
            reco_pos = np.array([entry.reco_pos_x, entry.reco_pos_y, entry.reco_pos_z])
            this_cluster.set_variables(reco_pos = reco_pos)
            clusters.append(this_cluster)
        break
    return clusters

