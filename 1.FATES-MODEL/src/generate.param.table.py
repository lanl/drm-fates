"""
# Generate parameter table for sensitivity analysis
# Rutuja Chitra-Tarak
# Dec 9, 2021
"""
import pandas as pd
import rpy2
from rpy2.robjects.packages import importr # import R's "base" package
base = importr('base')
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri # install any dependency package if you get error like "module not found"
pandas2ri.activate()
from rpy2.robjects.conversion import localconverter
import os 
import argparse
import yaml

print('rpy2 version is', rpy2.__version__)

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
PARAM_Table = PROJECT_ROOT+'/'+config_dict['BUILT_PARAM_TABLE_PATH']
outdir=config_dict['OUTPUT_DIR']
SLICES = config_dict['SLICES']

pd_df = pd.DataFrame({'par.names': config_dict['RANGE']['PARAMETERS'],
                      'min.param': config_dict['RANGE']['MINPARAM'],
                      'max.param': config_dict['RANGE']['MAXPARAM']})
R_file = PROJECT_ROOT+'/src/generate.param.table.R'

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
generate_param_table_function_r = robjects.globalenv['generate.param.table']

# Converting pandas object into r object for passing into r function
#df_r = pandas2ri.py2ri(df) # This works for rpy2 2.9.4
with localconverter(robjects.default_converter + pandas2ri.converter):
 df_r = robjects.conversion.py2ri(pd_df) #For rpy2 versions >2.9 use .py2rpy
#print(type(df_r))
#calling function under base package
#print(base.summary(r_df))

#Invoking the R function and getting the result
df_result_r = generate_param_table_function_r(df_r, SLICES, outdir, PARAM_Table)
# Converting it back to a pandas dataframe.
# df_result = pandas2ri.ri2py(df_result_r)

with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.ri2py(df_result_r) #For rpy2 versions >=2.9 use .rpy2py
#df_result

if(df_result):
    print('Successfully generated parameter table for sensitivity analysis.')
else:
    print("ERROR: Parameter table generation for sensitivity analysis  was unsuccessful.") 
    exit(0)

