"""
# Runs a shell script to run ELM in parallel
# Rutuja Chitra-Tarak
# July 25, 2021
"""
import os
import subprocess
import yaml
import argparse

# To read CASE_DIR from yaml file

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

CLONE_ROOT = config_dict['PROJECT_ROOT']+ '/' + config_dict['CASE_DIR']
PY_SRC_PATH = config_dict['PROJECT_ROOT']+ '/src/elm.py'

subprocess.call(['sh', './src/run_elm.sh', CLONE_ROOT, PY_SRC_PATH])
