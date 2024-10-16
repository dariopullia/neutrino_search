#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <ctime>
#include <nlohmann/json.hpp>

#include "CmdLineParser.h"
#include "Logger.h"

#include "daq_tps_clustering_libs.h"
#include "cluster.h"
#include "position_calculator.h"
#include "extract_neutrino_candidates_libs.h"


LoggerInit([]{
  Logger::getUserHeader() << "[" << FILENAME << "]";
});

int main(int argc, char* argv[]) {
    CmdLineParser clp;

    clp.getDescription() << "> extract_neutrino_candidates app."<< std::endl;

    clp.addDummyOption("Main options");
    clp.addOption("json",    {"-j", "--json"}, "JSON file containing the configuration");

    clp.addDummyOption("Triggers");
    clp.addTriggerOption("verboseMode", {"-v"}, "RunVerboseMode, bool");

    clp.addDummyOption();
    // usage always displayed
    LogInfo << clp.getDescription().str() << std::endl;

    LogInfo << "Usage: " << std::endl;
    LogInfo << clp.getConfigSummary() << std::endl << std::endl;

    clp.parseCmdLine(argc, argv);

    LogThrowIf( clp.isNoOptionTriggered(), "No option was provided." );

    LogInfo << "Provided arguments: " << std::endl;
    LogInfo << clp.getValueSummary() << std::endl << std::endl;

    std::string json = clp.getOptionVal<std::string>("json");
    // read the configuration file
    std::ifstream i(json);
    nlohmann::json j;
    i >> j;
    std::string filename = j["filename"];
    std::cout << "Filename: " << filename << std::endl;
    std::string outfolder = j["output_folder"];
    std::cout << "Output folder: " << outfolder << std::endl;
    int plane = j["plane"];
    std::cout << "Plane: " << plane << std::endl;
    int max_tps_per_filename = j["max_tps_per_filename"];
    std::cout << "Max events per filename: " << max_tps_per_filename << std::endl;
    int adc_integral_cut = j["adc_integral_cut"];
    std::cout << "ADC integral cut: " << adc_integral_cut << std::endl;
    int time_limit = j["time_limit"];
    std::cout << "Time limit: " << time_limit << std::endl;
    int min_tot_per_tp = j["min_tot_per_tp"];
    std::cout << "Min TOT per TP: " << min_tot_per_tp << std::endl;
    int min_integral_per_tp = j["min_integral_per_tp"];
    std::cout << "Min integral per TP: " << min_integral_per_tp << std::endl;
    int n_skip_tps = j["n_skip_tps"];
    std::cout << "Number of TPs to skip: " << n_skip_tps << std::endl;


    std::vector<std::string> plane_names = {"U", "V", "X"};
    // start the clock
    std::clock_t start;
    // filename is the name of the file containing the filenames to read
    std::vector<std::string> filenames;
    // read the file containing the filenames and save them in a vector
    std::ifstream infile(filename);
    std::string line;
    std::cout<<"Opening file: "<< filename << std::endl;
    // read and save the TPs
    while (std::getline(infile, line)) {
        filenames.push_back(line);
    }
    std::cout << "Number of files: " << filenames.size() << std::endl;
    std::vector<std::vector<double>> tps = file_reader(filenames, plane, max_tps_per_filename, min_tot_per_tp, min_integral_per_tp, n_skip_tps);
    std::cout << "Number of tps: " << tps.size() << std::endl;
    // cluster the tps
    std::vector<cluster> clusters = neutrino_explosion_finder(tps, 0, time_limit, adc_integral_cut);
    std::cout << "Number of clusters: " << clusters.size() << std::endl;
    // check if I have at least one cluster
    if (clusters.size() == 0){
        std::cout << "No clusters found" << std::endl;
        return 0;
    }

    std::string root_filename = outfolder + "/clusters_selected.root";
    // std::string root_filename = outfolder + "/" + plane_names[plane] + "/clusters_tick_limits_" + std::to_string(ticks_limit) + "_channel_limits_" + std::to_string(channel_limit) + "_min_tps_to_cluster_" + std::to_string(min_tps_to_cluster) + ".root";
    write_clusters_to_root(clusters, root_filename);
    std::cout << "clusters written to " << root_filename << std::endl;
    // stop the clock
    std::clock_t end;
    end = std::clock();
    return 0;
}
