"""
# Generate clones of a parameter file and substitutes desired parameter values
# Rutuja Chitra-Tarak
# July 7, 2021
"""
import pandas as pd
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import os 
import argparse
import yaml

SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
parser = argparse.ArgumentParser()
parser.add_argument('-c', '--config_file', action='store',
                    default=SCRIPT_DIR+'/../../config.yaml')
args = parser.parse_args()

with open(args.config_file, 'r') as in_file:
    config_dict = yaml.safe_load(in_file)
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
PARAM_DIR = config_dict['PARAM_DIR']
PARAM_PATH = PROJECT_ROOT +'/' + PARAM_DIR
CLONE_TYPE = config_dict['FILE_TYPE_TO_CLONE']['LIST']
PARAM_KEY = config_dict['PARAM_KEY']

clone_type = []
file_to_clone = []
clone_base = []
for i in range(0,len(CLONE_TYPE)):
   clone_type.append(CLONE_TYPE[i])
   file_to_clone.append(config_dict['PARAM_FILE'][CLONE_TYPE[i]])
   clone_base.append(config_dict['CLONE_BASE'][CLONE_TYPE[i]])

clone_pd_df = pd.DataFrame({'clone_type': clone_type,
		            'file_to_clone': file_to_clone,
			    'clone_base': clone_base})

if(config_dict['SENSITIVITY']):
    R_file = PROJECT_ROOT+'/src/generate.inputs_sen.R'
    if(config_dict['BUILD_PARAM_TABLE']):
        PARAM_Table = PROJECT_ROOT+'/'+config_dict['BUILT_PARAM_TABLE_PATH']
    else:
        PARAM_Table = PROJECT_ROOT+'/'+config_dict['GIVEN_PARAM_TABLE_PATH']
else:
        R_file = PROJECT_ROOT+'/src/generate.inputs.R'
        PARAM_Table = PROJECT_ROOT+'/'+config_dict['SITES_PARAM_TABLE_PATH']

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
generate_parameter_files_function_r = robjects.globalenv['generate.parameter.files']

# Reading and processing data
df = pd.read_csv(PARAM_Table)

# Converting it into r object for passing into r function
# df_r = pandas2ri.py2ri(df)
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_r = robjects.conversion.py2rpy(df) # For rpy2 versions >2.9 use py2rpy
df_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  clone_pd_df_r = robjects.conversion.py2rpy(clone_pd_df) # For rpy2 versions >2.9 use py2rpy
clone_pd_df_r

PARAM_KEY_r = rpy2.robjects.ListVector(PARAM_KEY)

#Invoking the R function and getting the result
df_result_r = generate_parameter_files_function_r(df_r, PARAM_PATH, clone_pd_df_r, PARAM_KEY_r)
# Converting it back to a pandas dataframe.
# df_result = pandas2ri.ri2py(df_result_r)
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.rpy2py(df_result_r) # For rpy2 versions >2.9 use rpy2py
df_result

print ("Parameter Table used for mutating parameter files is  "+PARAM_Table)

if(df_result):
    print('Parameter file clones generated. New Parameters Substituted Successfully!')
else:
    print("ERROR: Parameter substitution was unsuccessful.") 
    exit(0)

