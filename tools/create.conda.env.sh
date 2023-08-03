#"""
# This file creates a conda env within the working directory according to LANL HPC instructions:
# https://hpc.lanl.gov/software/hpc-provided-software/anaconda-python.html#AnacondaPython-CreatingAnacondaPythonEnvironments
# Rutuja Chitra-Tarak
# July 15, 2021
#"""
##**************************************************************
# Uncomment this during the first instance of running this script
#module load python/3.8-anaconda-2020.07

#module avail python

#conda info

## To check installed conda environments
# conda env list
## OR
## conda info --envs

##**************************************************************

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from

PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
cd $PROJECT_ROOT

# add the project directory to the conda env dir
conda config --add envs_dirs $PROJECT_ROOT

# To remove:
# conda config --remove envs_dirs $PROJECT_ROOT 

# To remove older conda env
conda remove --prefix $PROJECT_ROOT/elm_env --all

# Create a conda environment with a specified name, elm_env, and a desired list of modules and their versions
conda create --prefix=elm_env pip #for conda to add the pathname of the project directory automatically, important that this command is run from $PROJECT_ROOT (this is set above)

# Activate the environment. Activation prepends to PATH.
conda activate elm_env # While ran from the prefix folder mentioned in conda create ($PROJECT ROOT), this invokes the elm_env referred to by the above prefix (even if there are other environments with the same name elm_env exist, with a different prefix)
# One can check whether the desired conda elm_env is activated, by checking with, conda info --envs, if the * marking the active env is marked at the correct path

# Add any other modules that were not found by the default conda channels
conda install python=3.9
python -m pip cache remove mpi4py
env CC="cc" MPICC="cc" python -m pip install --no-binary=mpi4py mpi4py
conda install matplotlib pandas pyyaml xarray
conda install -c conda-forge tzlocal scipy netcdf4 hdf5
conda install -c r r-base r-essentials r-ncdf4 r-slhd r-stringi
python -m pip cache remove rpy2
pip install rpy2 (install after other r-packages above, Jun-29-2023-modifications-for-cxu)
conda install -c conda-forge imageio
conda install -c r r-tidyverse

# Export your active environment to a new file:
conda env export > tools/environment.yml # this is the backup made from scratch from teh above script
conda env export > environment.yml # We will use this file to build the new environment from; and this is modifiable 
# Open and edit the names environment.yml to not have path but just the name elm_env

# ---------------------------------
# Useful commands
# -----------------------------------
# It is possible to further expand a conda environment by continuing to install python packages via the conda command. To do this, switch to the conda environment and issue the following: 
# conda install <space-separated list of packages>

# You will need to give execute command to run this file if running for the first time
# chmod +x create.conda.env.sh
# Run with
#./create.conda.env.sh

# With conda versions that end in yyyy-mo, path to the  elm_env would be automatically added. If not,
# Make sure the conda env is on the PATH
# echo $PATH
# set PATH=$currentdir'/elm_env/bin':$PATH
# echo $PATH

# Activate the environment
# conda activate elm_env
# Export your active environment to a new file:
# conda env export > tools/environment.yml

# Viewing a list of the packages in an environment
# If the environment is not activated
# conda list -n myenv
# If the environment is activated
# conda list
# -----End------
	
