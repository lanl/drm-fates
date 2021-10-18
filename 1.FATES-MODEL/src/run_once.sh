#!/bin/sh
# =======================================================================================
# This will create a conda environment with all the python and R modules needed
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
echo $SCRIPT_DIR
source "$SCRIPT_DIR/../tools/yaml.sh"
create_variables "$SCRIPT_DIR/../config.yaml"

PROJECT_ROOT="$SCRIPT_DIR/.."
echo $PROJECT_ROOT
# 2. To create the environment from a yml file, run this in a shell:
conda config --add envs_dirs $PROJECT_ROOT
## Name of the new conda env is conda_env in the environment.yml file. If such env elready xists, change the name in the file
sed -i 's/conda_env/elm_env/g' environment.yml
sed -i '/prefix/d' environment.yml
sed -i -e '$a\' -e 'prefix: $PROJECT_ROOT/elm_env' environment.yml
cd $PROJECT_ROOT
conda env create -f environment.yml

# Alternatively, you could also create the environment from scratch:
# ./src/create.conda.env.sh

## 3. Make sure you are on the desired ELM branch. 

#cd $PROJECT_ROOT
cd $ACME_ROOT
git clone git@github.com:rutujact/E3SM.git --branch rutuja/ELM_FATES_HYDRO_DFLT_pedotrf_HKSAT_ADJ
#git checkout rutuja/ELM_FATES_HYDRO_DFLT_pedotrf_HKSAT_ADJ

## If a submodule is not found, update them

# git submodule update --init --recursive

## 4. Make sure you are on the desired fates branch. 

#cd $PROJECT_ROOT/E3SM/components/clm/src/external_models
cd $FATES_ROOT
git clone git@github.com:rutujact/fates.git --branch rutuja/xuc/coastal_veg_ready_to_merge_July_2019_hybrid_static --single-branch

