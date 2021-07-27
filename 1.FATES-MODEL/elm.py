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
HPU_PATH = config_dict['HPU_PATH']
LOG_PATH = config_dict['LOG_PATH']
TOTAL= HPU_ID_END - HPU_ID_START + 1

SCRIPT = "parallel.run.py"
BASE_CASE = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.met.v5.2016-2018"
RUN_ROOT = "/lustre/scratch3/turquoise/rutuja/ACME/cases"
FINALTAG = "clm2.h0.2020-12.nc"

command = "python ./" + str(SCRIPT) + " -c " + str(BASE_CASE) + "." + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(HPU_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH)

#command = "mpiexec -n " + str(TOTAL) + " python ./" + str(SCRIPT) + " -c " + str(BASE_CASE) + "." + " -r " + str(RUN_ROOT) + " -f " + str(FINALTAG) + " -s " +str(HPU_ID_START) + " -t " + str(TOTAL) + " -g " + str(LOG_PATH)

print(command)
os.system(command)
#os.system('./run_elm2.sh', 'command')
