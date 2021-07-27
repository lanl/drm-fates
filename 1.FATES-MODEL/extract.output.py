"""
# Extract ELM emsemble outputs
# Rutuja Chitra-Tarak
# July 26, 2021
"""
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os
import argparse
import yaml

dir_path = os.path.dirname(os.path.realpath(__file__))
R_file = 'extract.output.R'
outdir=dir_path + '/OutputExtract'
cloneroot = "/lustre/scratch3/turquoise/rutuja/ACME/cases"

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
extractres_h0_r = robjects.globalenv['extractres_h0']
extractres_h1_r = robjects.globalenv['extractres_h1']
extractres_h2_r = robjects.globalenv['extractres_h2']

# Determine input parameters
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default='config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)

sam_start = config_dict['HPU_ID_START']
sam_end = config_dict['HPU_ID_END']

#filebase = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.met.v5.2016-2018"
filebase = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.full.met.v6."

start_year = 2016
end_year = 2018

filterFile = "Filter.txt"
var_vec_h1 = ["H2OSOI", "QRUNOFF"] 
#['H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'GPP', 'TWS', 'ZWT', 'BTRAN', 'SOILPSI']
scale_vec_h1 = [1.0, 1.0] 
#[1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0]

# Converting python objects into r objects for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
  var_vec_h1_r = robjects.conversion.py2ri(var_vec_h1) # in later rpy2 versions use py2rpy
var_vec_h1_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  scale_vec_h1_r = robjects.conversion.py2ri(scale_vec_h1) # in later rpy2 versions use py2rpy
scale_vec_h1_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  sam_start_r = robjects.conversion.py2ri(sam_start) # in later rpy2 versions use py2rpy
sam_start_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  sam_end_r = robjects.conversion.py2ri(sam_end) # in later rpy2 versions use py2rpy
sam_end_r

#Invoking the R function and getting the result. Note that the sequence of arguments is critical
df_result = extractres_h1_r(sam_start_r, sam_end_r, outdir, cloneroot, filebase, var_vec_h1_r, scale_vec_h1_r, filterFile, start_year, end_year)

print(df_result)
