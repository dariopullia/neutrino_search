
#ifndef CLUSTER_H
#define CLUSTER_H

#include <vector>
#include <cmath>
#include <map>
#include <string>
// include root libraries
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TLeaf.h"
#include "TMatrixD.h"

extern std::map<std::string, int> variables_to_index;


class cluster {
public:
    cluster() {reco_pos_ = {0, 0, 0};}
    cluster(std::vector<std::vector<double>> tps) { tps_ = tps; reco_pos_ = {0, 0, 0}; update_cluster_info(); }
    ~cluster() { }

    void update_cluster_info();

    std::vector<std::vector<double>> get_tps() const { return tps_; }
    void set_tps(std::vector<std::vector<double>> tps) { tps_ = tps; update_cluster_info();}
    std::vector<double> get_tp(int i) { return tps_[i]; }
    int get_size() { return tps_.size(); }
    std::vector<float> get_reco_pos() { return reco_pos_; }
    void set_reco_pos(std::vector<float> pos) { reco_pos_ = pos; }
    void set_total_charge(float charge) { total_charge_ = charge; }
    float get_total_charge() { return total_charge_; }
private:
    std::vector<std::vector<double>> tps_;
    std::vector<float> reco_pos_;    
    float total_charge_;
};

#endif // CLUSTER_H

// write the clusters to a root file
void write_clusters_to_root(std::vector<cluster>& clusters, std::string root_filename);
