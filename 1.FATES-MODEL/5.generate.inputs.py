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
dir_path = os.path.dirname(os.path.realpath(__file__))
R_file = '4.generate.inputs.R'
surf_basefile_rel_path = 'data-raw/surfdata_bci_panama_v1_c171113.nc'

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
generate_surface_files_function_r = robjects.globalenv['generate.surface.files']

# Reading and processing data
df = pd.read_csv("HPU.Table.csv")

# Converting it into r object for passing into r function
# df_r = pandas2ri.ri2py(df)
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_r = robjects.conversion.py2ri(df) # in later rpy2 versions use py2rpy
df_r

#Invoking the R function and getting the result
df_result_r = generate_surface_files_function_r(df_r, dir_path, surf_basefile_rel_path)
# Converting it back to a pandas dataframe.
# df_result = pandas2ri.py2ri(df_result_r)
with localconverter(robjects.default_converter + pandas2ri.converter):
  df_result = robjects.conversion.ri2py(df_result_r) # in later rpy2 versions use rpy2py
df_result

if(df_result):
    print('Surface Data file clones generated. New Parameters Substituted Successfully!')
    exit(0)

