"""
# Creates clones of a base case via calling a shell script
# Clones end with IDs in HPU.Table
# Rutuja Chitra-Tarak
# July 25, 2021
"""
import os
import subprocess
import yaml
import argparse

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../config.yaml')
args = parser.parse_args()

# generate parameter data input
with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

SIM_ID_START = config_dict['SIM_ID_START']
SIM_ID_END = config_dict['SIM_ID_END']

# Set Case ID array
case_vec = list(range(SIM_ID_START, SIM_ID_END + 1))
#https://stackoverflow.com/questions/11392033/passing-python-array-to-bash-script-and-passing-bash-variable-to-python-functio
case_arr = ' '.join(str(v) for v in case_vec)

PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')

# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT + "/BASE_CASE_NAME.txt", "r")
BASE_CASE=ff.read()

# Set the CLONE ROOT Directory
CLONE_ROOT = PROJECT_ROOT + '/' + config_dict['CASE_DIR']

# Set the RUN ROOT Directory
runroot = os.environ["RUN_ROOT"]

# Set whether parameter files need to be varied across clones
MUTATE = str(config_dict['MUTATE'])

# Set the parameter basefile to clone
CLONE_BASE = config_dict['CLONE_BASE']
CLONE_TYPE = config_dict['CLONE_TYPE']
CLONE_FILE = str(config_dict['PARAM_FILE'][CLONE_TYPE][0])
# This will clone a base case ending with elements of case_arr
# For python > shell examples: https://stackoverflow.com/questions/32085956/pass-a-variable-from-python-to-shell-script
subprocess.call(['tcsh', './src/create.ELM.ensemble.csh', case_arr, BASE_CASE, CLONE_ROOT, runroot, MUTATE, CLONE_BASE, CLONE_FILE, CLONE_TYPE])
