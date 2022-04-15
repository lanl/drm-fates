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

# 7. To prepare the parallel simulation command as per elm.py, run the following command:
python src/elm.py
MPICOMMAND=`cat $SCRIPT_DIR/../mpi_command.txt`
CASEDIR=`realpath $SCRIPT_DIR/../${CASE_DIR}`

sed -i "s/^#SBATCH -N.*/#SBATCH -N ${N_NODE} # number of nodes/g" src/run_elm.sh
sed -i "s/^#SBATCH -t.*/#SBATCH -t ${WALL_TIME}/g" src/run_elm.sh
sed -i "s|^casedir.*|casedir\=\'$CASEDIR\'|g" src/run_elm.sh
sed -i '/^mpi/d' src/run_elm.sh
sed -i "/jobid*/a $MPICOMMAND" src/run_elm.sh

# 8. To run parallel simulations on the back node, run sbatch:
rm slurm*
sbatch src/run_elm.sh
