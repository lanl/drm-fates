"""
# Find ELM cases that did and did not run completely
# Clones end with IDs in HPU.Table
# Rutuja Chitra-Tarak
# July 25, 2021
"""
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os
import yaml
import argparse

# Defining the R script and loading the instance in Python
R_file = 'OutputExtract/create.filter.R'
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
r_function_filter = robjects.globalenv['function_filter']

# Preparing input parameters
outdir = 'OutputExtract'

parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
  config_dict = yaml.safe_load(in_file)

start_n = config_dict['HPU_ID_START']
stop_n = config_dict['HPU_ID_END']

filebase="BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.met.v5.2016-2018"
finaltag="clm2.h0.2020-12.nc"

# Invoking the R function and getting the result
df_result_r = r_function_filter(outdir, filebase, finaltag, start_n, stop_n)
# Converting it back to a pandas dataframe.
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.ri2py(df_result_r) # in later rpy2 versions use rpy2py
df_result
print(df_result, ' Cases Finished Successfully!')
exit(0)

