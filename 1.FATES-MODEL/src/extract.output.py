"""
# Extract ELM emsemble outputs
# Rutuja Chitra-Tarak
# July 26, 2021
"""
import pandas as pd
import os
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import argparse
import yaml
from pandas.core.index import Index as PandasIndex

# Determine input parameters
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

R_file = SCRIPT_DIR+'/extract.output.R'
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
print(PROJECT_ROOT)
outdir = PROJECT_ROOT+'/'+config_dict['OUTPUT_DIR']
runroot = os.environ["RUN_ROOT"]

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
#extractres_h0_r = robjects.globalenv['extractres_h0']
extractres_h1_r = robjects.globalenv['extractres_h1']
#extractres_h2_r = robjects.globalenv['extractres_h2']

sam_start = config_dict['HPU_ID_START']
sam_end = config_dict['HPU_ID_END']

start_year = config_dict['DATM_CLMNCEP_YR_START']
end_year = config_dict['DATM_CLMNCEP_YR_END']

# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
filebase=base_case.strip()

filterFile = "Filter.txt"
var_vec_h1 = ["RAIN", "QRUNOFF", "TWS", "BTRAN", "ZWT", "H2OSOI"] 
#['H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'GPP', 'TWS', 'ZWT', 'BTRAN', 'SOILPSI']
scale_vec_h1 = [1.0]*len(var_vec_h1) 

# Converting python objects into r objects for passing into r function

with localconverter(robjects.default_converter + pandas2ri.converter):
  var_vec_h1_r = robjects.vectors.StrVector(var_vec_h1) 
var_vec_h1_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  scale_vec_h1_r = robjects.vectors.FloatVector(scale_vec_h1) 
scale_vec_h1_r

#Invoking the R function and getting the result. Note that the sequence of arguments is critical
h1_result = extractres_h1_r(sam_start, sam_end, outdir, runroot, filebase, var_vec_h1_r, scale_vec_h1_r, filterFile, start_year, end_year)
if (h1_result):
    print('Daily outputs extracted successfully for ', var_vec_h1, 'at ', outdir+"/elm_daily_outputs.txt")
    exit(0)
else:
   print("Daily outputs not generated")
