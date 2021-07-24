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

n = config_dict['DURATION']
HPU_ID_START = config_dict['HPU_ID_START']
HPU_ID_END = config_dict['HPU_ID_END']
HPU_PATH = config_dict['HPU_PATH']
LOG_PATH = config_dict['LOG_PATH']
#HPU.table = read.csv(HPU_PATH)

t= HPU_ID_END - HPU_ID_START + 1
command = "python parallel.run.random.py -c 'BCI.ICLM45ED.wolf.intel.Cb14cb81-F812a621.' -r /lustre/scratch3/turquoise/cxu/ACME/cases -f clm2.h0.2003-12.nc"+ " -s " +str(HPU_ID_START) + " -t " + str(t)+ " -g " + str(LOG_PATH)
os.system(command)
