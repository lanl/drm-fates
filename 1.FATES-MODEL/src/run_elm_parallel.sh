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

# 1.1 Unzip sample climate data:

unzip -u data/bci_0.1x0.1_v4.0i/bci_0.1x0.1_met.v5.1.zip -d data/bci_0.1x0.1_v4.0i/

# 2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv, run:

python src/generate.inputs.py
 
# 3. To generate a base case of ELM, run:

python src/create.basecase.py
 
# 4. To generate ELM clone cases each associated with an ensemble member, run:

python src/create.ELM.ensemble.py
 
# 5. To prepare the parallel simulation command as per elm.py, run the following command:
python src/elm.py
MPICOMMAND=`cat $SCRIPT_DIR/../mpi_command.txt`
CASEDIR=`realpath $SCRIPT_DIR/../${CASE_DIR}`

sed -i "s/^#SBATCH -N.*/#SBATCH -N ${N_NODE} # number of nodes/g" src/run_elm.sh
sed -i "s/^#SBATCH -t.*/#SBATCH -t ${WALL_TIME}/g" src/run_elm.sh
sed -i "s|^casedir.*|casedir\=\'$CASEDIR\'|g" src/run_elm.sh
sed -i '/^mpi/d' src/run_elm.sh
sed -i "/jobid*/a $MPICOMMAND" src/run_elm.sh

# 6. To run parallel simulations on the back node, run sbatch:
rm slurm*
sbatch src/run_elm.sh

# 7. Make sure cases have successfully finished, else re-run sbatch:
if grep -w "Finished on"  slurm*.out; then
   echo "sbatch has finished running"
else
echo "sbatch has not yet finished running"
fi
# The following test should ideally be run if sbatch above has finished running,
# but typically no harm if sbatch is submitted again. It will skip those cases that are completed:
# This will test if at least one case has finished successfully.
if grep -w "CASE.RUN HAS FINISHED" slurm*.out; then
   echo "Simulation ran successfully"
else
   echo "Rerunning batch";
   sbatch src/run_elm.sh
fi
