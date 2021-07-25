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
                    default='3.config.yaml')
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

# Set the CLONE ROOT Directory
#dir_path = os.path.dirname(os.path.realpath(__file__))
#print(dir_path)
#CLONE_ROOT="$CASE_ROOT/ELM_Disease"
#os.putenv('CLONE_ROOT', str(dir_path))

# Set the BASE CASE Directory
BASE_CASE = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.met.v5.2016-2018"
os.putenv('BASE_CASE', BASE_CASE)

# This will clone a base case ending with elements of case_arr
subprocess.call('./Generate_Sen_Cases.csh')

