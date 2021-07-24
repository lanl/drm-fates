#"""
# This file creates a conda env within the working directory according to LANL HPC instructions:
# https://hpc.lanl.gov/software/hpc-provided-software/anaconda-python.html#AnacondaPython-CreatingAnacondaPythonEnvironments
# Rutuja Chitra-Tarak
# July 15, 2021
#"""
module avail python
module load python/3.8-anaconda-2020.07

module avail python

conda info

# To check installed conda environments
conda env list
# OR
# conda info --envs

# To initialise conda
source "${LPYTHON_HOME}/etc/profile.d/conda.csh"

# To initialise conda at login
# module load python/3.8-anaconda-2020.07
# conda init
# move the contents of the .bashrc file to .tcshrc
# This will also allow for conda activate to work

# add the curret directory to the conda env dir
set currentdir=`pwd`
conda config --add envs_dirs $currentdir
# To remove:
# conda config --remove envs_dirs $currentdir 

# Create a conda environment with a specified name, conda_env, and a desired list of modules and their versions
conda create --name conda_env python=3.6 r-base=3.6 r-essentials=3.6 rpy2 pandas r-ncdf4

# Activate the environment
conda activate conda_env

# Add any other modules that were not found by the default conda channels
conda install -c conda-forge tzlocal 

# You will need to give execute command to run this file if running for the first time
# chmod +x create.conda.env.sh
# Run with
#./create.conda.env.sh
# you will be prompted if an exisiting command by the same name in teh current dir needs to be replaced
# And whether to proceed to load all the modules.

# Make sure the conda env is on the PATH
echo $PATH
set PATH=$currentdir'/conda_env/bin':$PATH
# Activate the environment
conda activate conda_env
# -----End------
	
