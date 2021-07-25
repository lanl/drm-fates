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

## To initialise conda
# source "${LPYTHON_HOME}/etc/profile.d/conda.csh"

## To initialise conda at login
# module load python/3.8-anaconda-2020.07
# conda init tcsh
# You will need to make sure that ~/.tcshrc runs at the login
##**************************************************************

# add the curret directory to the conda env dir
set currentdir=`pwd`
conda config --add envs_dirs $currentdir
# To remove:
# conda config --remove envs_dirs $currentdir 

# To remove older conda env
conda remove --prefix $currentdir/conda_env --all

# Create a conda environment with a specified name, conda_env, and a desired list of modules and their versions
conda create --name conda_env python=3.6 r-base=3.6 r-essentials=3.6 rpy2 pandas r-ncdf4 mpi4py pyyaml

# Activate the environment. Activation prepends to PATH.
conda activate conda_env

# Add any other modules that were not found by the default conda channels
conda install -c conda-forge tzlocal 

# ---------------------------------
# Useful commands
# -----------------------------------
# It is possible to further expand a conda environment by continuing to install python packages via the conda command. To do this, switch to the conda environment and issue the following: 
# conda install <space-separated list of packages>

# You will need to give execute command to run this file if running for the first time
# chmod +x create.conda.env.sh
# Run with
#./create.conda.env.sh

# With conda versions that end in yyyy-mo, path to the  conda_env would be automatically added. If not,
# Make sure the conda env is on the PATH
# echo $PATH
# set PATH=$currentdir'/conda_env/bin':$PATH
# echo $PATH

# Activate the environment
# conda activate conda_env
# Export your active environment to a new file:
# conda env export > environment.yml

# Viewing a list of the packages in an environment
# If the environment is not activated
# conda list -n myenv
# If the environment is activated
# conda list
# -----End------
	
