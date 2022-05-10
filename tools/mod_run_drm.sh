#!/bin/sh
# =======================================================================================
#
# This script will create, setup and build and run a multi-site ELM-FATES simulation  
# based on config.yaml.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../tools/yaml.sh"`
create_variables "$SCRIPT_DIR/../config.yaml"

# 2. To prepare shell script to submit on the back-end node:
sed -i "s/^#SBATCH -N.*/#SBATCH -N ${N_NODE} # number of nodes/g" run_drm.sh
sed -i "s/^#SBATCH -t.*/#SBATCH -t ${WALL_TIME}/g" run_drm.sh
sed -i "s/^#SBATCH -A.*/#SBATCH -A ${ACCOUNT}/g" run_drm.sh
