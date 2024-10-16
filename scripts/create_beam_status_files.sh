#!bin/bash

REPO_HOME=$(git rev-parse --show-toplevel)

# add the path /afs/cern.ch/work/d/dapullia/public/dune/protodune/neutrino_search/external_libs/lib64/python3.9/site-packages
export PYTHONPATH=${REPO_HOME}/external_libs/lib/python3.10/site-packages:$PYTHONPATH



# move to the folder, run and come back to scripts
cd ${REPO_HOME}/app/

python create_beam_status_files.py --input_json ${REPO_HOME}/json/create_beam_status_files/standard.json

cd ${REPO_HOME}/scripts/


