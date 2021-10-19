#!/bin/sh
# =======================================================================================
# This will create a conda environment with all the python and R modules needed
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source "$SCRIPT_DIR/../tools/yaml.sh"
create_variables "$SCRIPT_DIR/../config.yaml"

PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
cd $PROJECT_ROOT
mkdir -p elm_env # create if not already present

# 2. To create the environment from a yml file, run this in a shell:
conda config --add envs_dirs $PROJECT_ROOT/elm_env
## Name of the new conda env is conda_env in the environment.yml file. If such env elready xists, change the name in the file
sed -i 's/conda_env/elm_env/g' environment.yml
sed -i '/prefix/d' environment.yml
sed -i -e '$a\' -e "prefix: $PROJECT_ROOT/elm_env" environment.yml
cd $PROJECT_ROOT
conda env create -f environment.yml

# Alternatively, you could also create the environment from scratch:
# ./src/create.conda.env.sh

## 3. Make sure you are on the desired ELM branch. 
cd $ACME_ROOT
git clone $ELM_REMOTE --branch $ELM_BRANCH --single-branch
git checkout $ELM_BRANCH

## If a submodule is not found, update them
# git submodule update --init --recursive

## 4. Make sure you are on the desired fates branch. 
cd $FATES_ROOT
git clone $FATES_REMOTE --branch $FATES_BRANCH --single-branch
git checkout $FATES_BRANCH

