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

R_file = SCRIPT_DIR+'/update.restart.treelist.R'
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
outdir = PROJECT_ROOT+'/'+config_dict['OUTPUT_DIR']
VDM2FM = PROJECT_ROOT+'/'+config_dict['VDM2FM_DIR']
runroot = os.environ["RUN_ROOT"]

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
update_restart_treelist_r = robjects.globalenv['update_restart_treelist']

sam_start = int(config_dict['SIM_ID_START'])
sam_end = int(config_dict['SIM_ID_END'])
finalyear = int(config_dict['FINAL_TAG_YEAR'])
fire_res = int(config_dict['FIRE_RES'])
fates_res = int(config_dict['FATES_RES'])

grass_pft_index = 3
# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
filebase=base_case.strip()

filterFile = "Filter.txt"
# Converting python objects into r objects for passing into r function

#Invoking the R function and getting the result. Note that the sequence of arguments is critical
restart_result = update_restart_treelist_r(sam_start, sam_end, outdir, VDM2FM, runroot, filebase, filterFile, finalyear, fire_res, fates_res, grass_pft_index)
if (restart_result):
    print("Restart file successfully updated by removing trees that died in fire.")
    exit(0)
else:
   print("Restart file not updated.")
