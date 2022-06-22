"""
# Extract ELM ensemble treelist outputs
# Rutuja Chitra-Tarak
# MMay 14, 2022
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
                    default=SCRIPT_DIR+'/../../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

R_file = SCRIPT_DIR+'/extract.moisture.R'
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
outdir = PROJECT_ROOT+'/'+config_dict['OUTPUT_DIR']
runroot = os.environ["RUN_ROOT"]

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
extract_moisture_r = robjects.globalenv['extract_moisture']

sam_start = int(config_dict['SIM_ID_START'])
sam_end = int(config_dict['SIM_ID_END'])
finalyear = int(config_dict['FINAL_TAG_YEAR'])
fates_res = int(config_dict['FATES_RES'])
fates_c2b = 2 # Carbon to biomass 
fates_levscls = 13 # No. of FATES size classes in h0.cdl history output
# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
filebase=base_case.strip()

filterFile = "Filter.txt"
var_vec_re = ["FATES_LTH_SCPF", "FATES_STH_SCPF"] #["fates_litter_moisture_pa_nfsc"]

# Converting python objects into r objects for passing into r function

with localconverter(robjects.default_converter + pandas2ri.converter):
  var_vec_re_r = robjects.vectors.StrVector(var_vec_re)
var_vec_re_r

#Invoking the R function and getting the result. Note that the sequence of arguments is critical
moisture_result = extract_moisture_r(sam_start, sam_end, outdir, runroot, filebase, var_vec_re_r, filterFile, finalyear, fates_res, fates_c2b, fates_levscls)

if (moisture_result):
    print('Live fuel moisture extracted successfully at', outdir + "/livefuel.moisture.txt")
    exit(0)
else:
   print("Live fuel moisture not extracted")
