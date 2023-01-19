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
conda config --set env_prompt envname
sed -i '/prefix/d' environment.yml
sed -i -e '$a\' -e "prefix: $PROJECT_ROOT/elm_env" environment.yml
cd $PROJECT_ROOT
conda env create -f environment.yml

# Alternatively, you could also create the environment from scratch:
# ./src/create.conda.env.sh

## 3. Make sure you are on the desired ELM branch. 
mkdir -p $E3SM_ROOT
cd $E3SM_ROOT

if ! [ -d .git ]; then
  git clone $ELM_REMOTE .
fi
git remote -v | grep -w elm_repo || git remote add elm_repo $ELM_REMOTE
git remote set-url elm_repo $ELM_REMOTE
git checkout $ELM_BRANCH
# git checkout --track elm_repo/$ELM_BRANCH

## If a submodule is not found, update them
git submodule update --init --recursive

## 4. Make sure you are on the desired fates branch. 
#cd $FATES_ROOT
#cd ..
#rm -rf fates
#git clone $FATES_REMOTE fates
#git remote -v | grep -w fates_repo || git remote add fates_repo $FATES_REMOTE
#cd $FATES_ROOT
#git checkout $FATES_BRANCH

