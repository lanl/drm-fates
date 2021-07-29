"""
# Creates clones of a base case via calling a shell script
# Clones end with IDs in HPU.Table
# Rutuja Chitra-Tarak
# July 25, 2021
"""
import os
import subprocess
import shlex
import yaml
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config.yaml')
args = parser.parse_args()

# generate surface data input
with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

HPU_ID_START = config_dict['HPU_ID_START']
HPU_ID_END = config_dict['HPU_ID_END']

# Set Case ID array
case_arr = list(range(HPU_ID_START, HPU_ID_END + 1))
#https://stackoverflow.com/questions/11392033/passing-python-array-to-bash-script-and-passing-bash-variable-to-python-functio
os.putenv('case_arr', ' '.join(str(v) for v in case_arr))

# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff = open("BASE_CASE_NAME.txt", "r")
print(ff.read())
os.putenv('BASE_CASE', ff.read())

# Set the CLONE ROOT Directory
CLONE_ROOT = config_dict['PROJECT_ROOT']+ '/' + config_dict['CASE_DIR']
print(CLONE_ROOT)
os.putenv('CLONE_ROOT', CLONE_ROOT)

# This will clone a base case ending with elements of case_arr
subprocess.call('./src/create.ELM.ensemble.csh')

