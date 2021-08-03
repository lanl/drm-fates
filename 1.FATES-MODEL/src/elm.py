import pandas as pd
#import random
import numpy as np
import yaml
#import pyarrow as pa
#import pyarrow.parquet as pq
import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config.yaml')
args = parser.parse_args()

# generate surface data input
with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

#n = config_dict['DURATION']
HPU_ID_START = config_dict['HPU_ID_START']
HPU_ID_END = config_dict['HPU_ID_END']
LOG_PATH = config_dict['PROJECT_ROOT']+'/'+config_dict['LOG_DIR']
N_CPU_PER_NODE = config_dict['N_CPU_PER_NODE']
N_NODE = config_dict['N_NODE']
WALL_TIME = config_dict['WALL_TIME']
N_PROC = N_NODE * N_CPU_PER_NODE
TOTAL= HPU_ID_END - HPU_ID_START + 1

SCRIPT = config_dict['ELM_RUN_PY']
ff = open(config_dict['PROJECT_ROOT']+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
BASE_CASE=base_case.strip()
#BASE_CASE='BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.full.met.v6'
RUN_ROOT = config_dict['RUN_ROOT']
FINALTAG = "clm2.h0."+ str(config_dict['DATM_CLMNCEP_YR_END']) +"-12.nc"

command = "python ./" + str(SCRIPT) + " -c " + str(BASE_CASE) + "." + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(HPU_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH)

#command = "mpiexec -n " + str(N_PROC) + " python ./" + str(SCRIPT) + " -c " + str(BASE_CASE) + "." + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(HPU_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH)

print(command)
os.system(command)
