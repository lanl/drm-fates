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
from pandas.core.index import Index as PandasIndex

# Defining the R script and loading the instance in Python
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
R_file = SCRIPT_DIR + '/create.filter.R'
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
r_function_filter = robjects.globalenv['function_filter']

# Preparing input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
  config_dict = yaml.safe_load(in_file)

PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
outdir = PROJECT_ROOT+'/'+config_dict['OUTPUT_DIR']
runroot = os.environ["RUN_ROOT"]

sam_start = config_dict['SIM_ID_START']
sam_stop = config_dict['SIM_ID_END']

nsam = sam_stop - sam_start + 1
print(nsam)
start_year = config_dict['DATM_CLMNCEP_YR_START']
end_year = config_dict['DATM_CLMNCEP_YR_END']

# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
filebase = base_case.strip()

finaltag = "clm2.h0."+ str(config_dict['DATM_CLMNCEP_YR_END']) +"-12.nc"

# Invoking the R function and getting the result
df_result_r = r_function_filter(outdir, runroot, filebase, finaltag, sam_start, sam_stop)
# Converting it back to a pandas dataframe.
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.ri2py(df_result_r) # in later rpy2 versions use rpy2py
df_result
print(int(df_result),' out of ',nsam,' cases Finished Successfully!')
exit(0)

