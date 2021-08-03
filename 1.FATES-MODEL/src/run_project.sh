#!/bin/sh
# =======================================================================================
#
# This script will create, setup and build and run a multi-site simulation  
# based on config.yaml.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================

# 1. To create and activate a conda environment with all the python and R modules needed, run this in a shell:

# ./create.conda.env.sh

# Alternatively, you could also create the environment by running

conda env create -f environment.yml

# Activate the environment

conda activate conda_env
 
# 2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv, run:

python src/generate.inputs.py

 
# 3. To generate a base case of ELM, run:

python src/create.basecase.py

 
# 4. To generate ELM clone cases each associated with an ensemble member, run:

python src/create.ELM.ensemble.py

 
# 5. To run the parallel emsemble simulations as per elm.py, run the following command:

sbatch src/run_elm.sh

 
# 6. To find which cases are complete (output/Filter.txt) and which are not (output/Missing.txt), run:

python src/create.filter.py

 
# 7. To extract outputs from ELM ensembles (output/exrtact/elm_daily_outputs.txt), run:

python src/extract.output.py

