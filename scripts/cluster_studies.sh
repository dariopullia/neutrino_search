#!bin/bash

REPO_HOME=$(git rev-parse --show-toplevel)

# move to the folder, run and come back to scripts
cd ${REPO_HOME}/app/

python cluster_studies.py --input_json ${REPO_HOME}/json/cluster_studies/protodune-29002.json
python cluster_studies.py --input_json ${REPO_HOME}/json/cluster_studies/protodune-29006_1.json
python cluster_studies.py --input_json ${REPO_HOME}/json/cluster_studies/protodune-29006_2.json
python cluster_studies.py --input_json ${REPO_HOME}/json/cluster_studies/protodune-29006_3.json

cd ${REPO_HOME}/scripts/


