#!/bin/sh
# =======================================================================================
# This will create a conda environment with all the python and R modules needed
# =======================================================================================

# To create the environment from a yml file, run this in a shell:
conda config --add envs_dirs $PROJECT_ROOT
## Name of the new conda env is conda_env in the environment.yml file. If such env elready xists, change the name in the file
sed -i 's/conda_env/toronto_env/g' environment.yml
sed -i '/prefix/d' environment.yml
sed -i -e '$a\' -e 'prefix:/turquoise/usr/projects/veg/rutuja/ACME_cases/TORONTO/toronto_env' environment.yml

# Alternatively, you could also create the environment from scratch:
# ./src/create.conda.env.sh



