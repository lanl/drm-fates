import pandas as pd
#import random
import numpy as np
import yaml
#import pyarrow as pa
#import pyarrow.parquet as pq
import os
import argparse

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../../config.yaml')
args = parser.parse_args()

# generate surface data input
with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

#n = config_dict['DURATION']
SIM_ID_START = config_dict['SIM_ID_START']
SIM_ID_END = config_dict['SIM_ID_END']
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
LOG_PATH = PROJECT_ROOT+'/'+config_dict['LOG_DIR']
N_CPU_PER_NODE = config_dict['N_CPU_PER_NODE']
N_NODE = config_dict['N_NODE']
WALL_TIME = config_dict['WALL_TIME']
N_PROC = N_NODE * N_CPU_PER_NODE
TOTAL= SIM_ID_END - SIM_ID_START + 1
CASE_DIR=config_dict['CASE_DIR']

SCRIPT = PROJECT_ROOT+'/'+config_dict['ELM_RUN_PY']
ff = open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
BASE_CASE=base_case.strip()
RUN_ROOT = os.environ["RUN_ROOT"]
finalyear = int(config_dict['FINAL_TAG_YEAR'])
# First FATES simulation is set to run for 1 month
if int(config_dict['CYCLE_INDEX'])==0:
    FINALTAG = "elm.h0."+ str(finalyear) +"-01.nc"
else:
    FINALTAG = "elm.h0."+ str(finalyear) +"-12.nc"
print('finaltag check is ' + FINALTAG)
#command = "python " + str(SCRIPT) + " -c " + str(BASE_CASE) + "." + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(SIM_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH)

# command = "mpiexec -n " + str(TOTAL) + " python " + str(SCRIPT) + " -c " + "'"+ str(BASE_CASE) + "." + "'" + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(SIM_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH) # commented out by SXM
command = "srun -n " + str(TOTAL) + " python " + str(SCRIPT) + " -c " + "'"+ str(BASE_CASE) + "." + "'" + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(SIM_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH) # ASXM

file = open("mpi_command.txt", "w")
file.write(command + "\n")
file.close()

cc = open(PROJECT_ROOT+"/mpi_command.txt", "r")
stored_command=cc.read()

print('Command in "mpi_command.txt" is '+ stored_command)
exit(0)

# The following option could be used if this file needs to be passed onto src/run_elm.sh via run_elm.py
# And if front node parallel simulation needs to be run, one would then run python run_elm.py; after letting src/run_elm.sh accept variables from run_elm.py
# But since sbatch src/run_elm.sh works to parallel run on back nodes, this won't be required
#os.system(command) 
