#!/bin/sh
# =======================================================================================
#
# This script will create, setup and build and run a multi-site ELM-FATES simulation  
# based on config.yaml and extracts outputs.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================
# 0. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
echo $SCRIPT_DIR
source "$SCRIPT_DIR/../tools/yaml.sh"
create_variables "$SCRIPT_DIR/../config.yaml"

# 1. To create and activate a conda environment with all the python and R modules needed, run this in a shell:

# ./src/create.conda.env.sh

# Alternatively, you could also create the environment from yml file:
# conda config --add envs_dirs $PROJECT_ROOT
## Name of the new conda env is conda_env in the environment.yml file. If such env elready xists, change the name in the file
# sed -i 's/conda_env/new_env/g' environment.yml
#conda env create -f environment.yml

# Activate the environment

#conda activate conda_env
 
# 2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv, run:

#python src/generate.inputs.py

 
# 3. To generate a base case of ELM, run:

#python src/create.basecase.py

 
# 4. To generate ELM clone cases each associated with an ensemble member, run:

#python src/create.ELM.ensemble.py

 
# 5. To run the parallel emsemble simulations as per elm.py, run the following command:
sed -i "s/^#SBATCH -N.*/#SBATCH -N ${N_NODE} # number of nodes/g" src/run_elm.sh
sed -i "s/^#SBATCH -t.*/#SBATCH -t ${WALL_TIME}/g" src/run_elm.sh
sbatch src/run_elm.sh
 
# 6. To find which cases are complete (output/Filter.txt) and which are not (output/Missing.txt), run:

#python src/create.filter.py

 
# 7. To extract outputs from ELM ensembles (output/exrtact/elm_daily_outputs.txt), run:

#python src/extract.output.py

