#!/bin/sh
# =======================================================================================
# This will create a conda environment with all the python and R modules needed
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../tools/yaml.sh"`
create_variables "$SCRIPT_DIR/../config.yaml"

PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
cd $PROJECT_ROOT
mkdir -p elm_env # create if not already present

# 2. To create the environment from a yml file, run this in a shell:
# Let conda know the path to the new conda env and update it in environment.yml file
conda config --add envs_dirs $PROJECT_ROOT
sed -i '/prefix/d' environment.yml
sed -i -e '$a\' -e "prefix: $PROJECT_ROOT/elm_env" environment.yml
cd $PROJECT_ROOT
conda env create -f environment.yml

# Alternatively, you could also create the environment from scratch:
# ./src/create.conda.env.sh

## 3. Make sure you are on the desired ELM branch. 
mkdir -p $ACME_ROOT
cd $ACME_ROOT

if ! [ -d .git ]; then
  git clone $ELM_REMOTE .
fi
git remote -v | grep -w elm_repo || git remote add elm_repo $ELM_REMOTE
git remote set-url elm_repo $ELM_REMOTE
git checkout $ELM_BRANCH
# git checkout --track elm_repo/$ELM_BRANCH

## If a submodule is not found, update them
git submodule update --init --recursive

## Defaults paths are currently set to scratch3 by HPC folks, which is run-only. Change to scratch4
sed -i "s|scratch3|scratch4|g" /turquoise/usr/projects/veg/${USER}/ACME/cime/config/e3sm/machines/config_machines.xml

## 4. Make sure you are on the desired fates branch. 
cd $FATES_ROOT
cd ..
rm -rf fates
git clone $FATES_REMOTE
#git remote -v | grep -w fates_repo || git remote add fates_repo $FATES_REMOTE
cd $FATES_ROOT
git checkout $FATES_BRANCH

