#!/bin/sh
# =======================================================================================
#
# This script will create, setup and build and run a multi-site ELM-FATES simulation  
# based on config.yaml.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../../tools/yaml.sh"`
create_variables "$SCRIPT_DIR/../../config.yaml"

# 2. Read python subprocess variables
export RESTART=$1

# 3. To prepare the parallel simulation command as per elm.py, run the following command:
python src/elm.py
MPICOMMAND=`cat $SCRIPT_DIR/../mpi_command.txt`
CASEDIR=`realpath $SCRIPT_DIR/../${CASE_DIR}`

sed -i "s|^setenv casedir.*|setenv casedir $CASEDIR|g" src/run_elm.sh # the script is now in tcsh after adding to be able to source ~/.tcshrc for the minimal conda env 
#sed -i "s|^casedir.*|casedir\=\'$CASEDIR\'|g" src/run_elm.sh # when the script is in bash
sed -i '/^srun/d' src/run_elm.sh
sed -i "/location*/a $MPICOMMAND" src/run_elm.sh

# 4. If this is a restart simulation run, modify cases to set CONTINUE_RUN variable to TRUE & update STOP_N:
if [ "$RESTART" == TRUE ] || [ "$RESTART" == true ]; then
	echo "Preparing FATES for a restart run"
	python src/restart.ELM.ensemble.py
else
        echo "Preparing FATES for an initial run"
fi

# 5. To run parallel simulations on the back node, run sbatch:
# CSXM (BGN)
#   1) at early time
#     a) sh     src/run_elm.sh cannot specify SBATCH properties (e.g., the project account and JOB_WALLCLOCK_TIME) 
#     b) sbatch src/run_elm.sh can    specify SBATCH properties (e.g., the project account and JOB_WALLCLOCK_TIME), but the next step will run without waiting for the completion of this step
#     c) we therefore keep using 1), but specify project account and JOB_WALLCLOCK_TIME in create.basecase.sh for FATES runs
#   2) at later time (https://stackoverflow.com/questions/71302330/run-multiple-files-consecutively-via-slurm-with-individual-timeout)
#     a) we found that using "sbatch --wait" can specify SBATCH properties (e.g., the project account and JOB_WALLCLOCK_TIME) while letting the next step to run until the completetion of this step
#     b) we therefore using "sbatch --wait src/run_elm.sh"
# CSXM (END)
sbatch --wait src/run_elm.sh # MSXM

#6. To find out if all simulations (cases) finished successfully, run this. It will also save which cases are complete (output/Filter.txt) and which are not (output/Missing.txt):
python src/create.filter.py

# 7. To extract FATES ensemble restart outputs (output/treelist.txt, output/litter.txt, output/litter.moisture.txt)
if [ "$HYDRO" == TRUE ] || [ "$HYDRO" == true ]; then
        python src/extract.moisture.py
fi
python src/extract.treelist.py
python src/extract.litter.py
