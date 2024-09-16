#!bin/bash

REPO_HOME=$(git rev-parse --show-toplevel)

# move to the folder, run and come back to scripts
cd ${REPO_HOME}/app/

python cluster_studies.py --input_json ${REPO_HOME}/json/cluster_studies/standard.json --output_folder /eos/user/d/dapullia/dune/protodune/neutrino_search/cluster_studies/

cd ${REPO_HOME}/scripts/


