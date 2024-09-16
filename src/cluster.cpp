#include <iostream>
#include <vector>

#include "cluster.h"
#include "position_calculator.h" 

std::map<std::string, int> variables_to_index = {
    {"time_start", 0},
    {"time_over_threshold", 1},
    {"time_peak", 2},
    {"channel", 3},
    {"adc_integral", 4},
    {"adc_peak", 5},
    {"detid", 6},
    {"type", 7},
    {"algorithm", 8},
    {"version", 9},
    {"flag", 10},
};


void cluster::update_cluster_info() {
    float total_charge = 0;
    std::vector<float> reco_pos = {0, 0, 0};

    for (int i = 0; i < tps_.size(); i++) {
        total_charge += tps_[i][variables_to_index["adc_integral"]];
        std::vector<float> pos = calculate_position(tps_[i]);
        reco_pos[0] += pos[0];
        reco_pos[1] += pos[1];
        reco_pos[2] += pos[2];
    }

    total_charge_ = total_charge;
    int ntps = tps_.size();
    reco_pos[0] /= ntps;
    reco_pos[1] /= ntps;
    reco_pos[2] /= ntps;
    reco_pos_ = reco_pos;
}




void write_clusters_to_root(std::vector<cluster>& clusters, std::string root_filename) {
    // create folder if it does not exist
    std::string folder = root_filename.substr(0, root_filename.find_last_of("/"));
    std::string command = "mkdir -p " + folder;
    system(command.c_str());
    // create the root file
    TFile *f = new TFile(root_filename.c_str(), "recreate");
    TTree *tree = new TTree("tree", "tree");
    // prepare objects to save the data
    // std::vector<std::vector<int>> matrix;
    std::vector<std::vector<double>> matrix;
    int nrows;
    int reco_pos_x; 
    int reco_pos_y;
    int reco_pos_z;


    // create the branches
    tree->Branch("matrix", &matrix);
    tree->Branch("nrows", &nrows);
    tree->Branch("reco_pos_x", &reco_pos_x);
    tree->Branch("reco_pos_y", &reco_pos_y);
    tree->Branch("reco_pos_z", &reco_pos_z);

    // fill the tree
    for (auto& g : clusters) {
        matrix = g.get_tps();
        nrows = g.get_size();
        reco_pos_x = g.get_reco_pos()[0];
        reco_pos_y = g.get_reco_pos()[1];
        reco_pos_z = g.get_reco_pos()[2];
        tree->Fill();
    }
    // write the tree
    tree->Write();
    f->Close();

    return;   
}
