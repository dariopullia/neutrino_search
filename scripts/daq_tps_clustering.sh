#!bin/bash
INPUT_JSON=/afs/cern.ch/work/d/dapullia/public/dune/protodune/neutrino_search/json/daq_tps_clustering/standard.json
REPO_HOME=$(git rev-parse --show-toplevel)

# parse the input
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input_file)
            INPUT_JSON="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: ./daq_tps_clustering.sh --input_file <input_json>"
            exit 0
            ;;  
        *)
            shift
            ;;
    esac
done

echo "REPO_HOME: ${REPO_HOME}"
# compile
echo "Compiling..."
cd ${REPO_HOME}/build/
cmake ..
make -j 5
# if successful, run the app
if [ $? -ne 0 ]; then
    echo "Compilation failed"
fi

# Run the app
./app/daq_tps_clustering -j $INPUT_JSON

cd ${REPO_HOME}/scripts/
