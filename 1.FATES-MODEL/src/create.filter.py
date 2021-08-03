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
import re

# Defining the R script and loading the instance in Python
dir_path = os.path.dirname(os.path.realpath(__file__))
R_file = dir_path+'/create.filter.R'
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
r_function_filter = robjects.globalenv['function_filter']

# Preparing input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
  config_dict = yaml.safe_load(in_file)

outdir = config_dict['PROJECT_ROOT']+'/'+config_dict['OUTPUT_DIR']

sam_start = config_dict['HPU_ID_START']
sam_stop = config_dict['HPU_ID_END']

start_year = config_dict['DATM_CLMNCEP_YR_START']
end_year = config_dict['DATM_CLMNCEP_YR_END']

# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(config_dict['PROJECT_ROOT']+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
#BASE_CASE=base_case.strip()
BASE_CASE='BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.full.met.v6.2017-2018'
timetag=str(start_year)+'-'+str(end_year)
filebase = re.sub(timetag, '', BASE_CASE)

finaltag = "clm2.h0."+ str(config_dict['DATM_CLMNCEP_YR_END']) +"-12.nc"

# Invoking the R function and getting the result
df_result_r = r_function_filter(outdir, filebase, finaltag, sam_start, sam_stop)
# Converting it back to a pandas dataframe.
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.ri2py(df_result_r) # in later rpy2 versions use rpy2py
df_result
print(df_result, ' Cases Finished Successfully!')
exit(0)

