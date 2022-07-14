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

R_file = SCRIPT_DIR+'/extract.treelist.R'
PROJECT_ROOT = os.path.abspath(SCRIPT_DIR+'/..')
outdir = PROJECT_ROOT+'/'+config_dict['OUTPUT_DIR']
VDM2FM = PROJECT_ROOT+'/'+config_dict['VDM2FM_DIR']
runroot = os.environ["RUN_ROOT"]

# Defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file)

# Loading the function we have defined in R.
extract_treelist_r = robjects.globalenv['extract_treelist']

sam_start = int(config_dict['SIM_ID_START'])
sam_end = int(config_dict['SIM_ID_END'])
finalyear = int(config_dict['FINAL_TAG_YEAR'])
fire_res = int(config_dict['FIRE_RES'])
fates_res = int(config_dict['FATES_RES'])
cycle_index = int(config_dict['CYCLE_INDEX'])

# HYDRO ON/OFF? If OFF (e.g. currently HYDRO doesn't work iwth grasses), moisture is set to 1.
HYDRO = int(config_dict['HYDRO'])
# FATES parameters
fates_CWD_frac_twig = 0.045
fates_c2b = 2
# fates_leaf_slatop = 0.00662, 0.0189200006425381 # For long-leaf pine & Turkey oak
# denleaf  = -2.3231_r8*sla/prt_params%c2b(ft) + 781.899_r8 L863 kg/m3
leafdensity = -2.3231*(0.00662+0.0189200006425381)/2/2 + 781.899 # kg/m3
# fates_wood_density: 0.58, 0.65, #g/cm3
# 0.615 g/cm3 = 615 kg/m3
wooddensity = (0.58 + 0.65)/2 # kg/m3

# Fire parameter; a version of fuel surface area to volume ratio
sizescale_pd_df = pd.DataFrame({'fates_pft': [1,2,3], # 1 = pine, 2 = braodleaf
                            'sizescale': [0.2,0.6,1]})
grass_pft_index = 3
# Set the BASE CASE name. This is generated from yaml and src/create.basecase.sh
ff=open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case=ff.read()
filebase=base_case.strip()

filterFile = "Filter.txt"
var_vec_re = ["fates_pft", "fates_dbh", "fates_height",  "fates_crown_depth",
               "fates_nplant", "fates_cohort_area", 
               "leaf_c_val_001", "sapw_c_val_001", "store_c_val_001","repro_c_val_001", "struct_c_val_001"]

#var_vec_re = ["fates_PatchesPerSite", "fates_CohortsPerPatch", "fates_patchdistturbcat", 
#               "fates_leaflittin_vec_001", "fates_livegrass", "fates_area", 
#               "fates_pft", "fates_nplant", "fates_size_class_lasttimestep", "fates_dbh", "fates_height",
#               "fates_crown_depth", "fates_cohort_area", "fates_scorch_ht_pa_pft", "fates_litter_moisture_pa_nfsc",
#               "fates_ag_cwd_vec_001", "fates_bg_cwd_vec_001", "fates_leaf_fines_vec_001", "fates_fnrt_fines_vec_001"]

# Converting python objects into r objects for passing into r function

with localconverter(robjects.default_converter + pandas2ri.converter):
  var_vec_re_r = robjects.vectors.StrVector(var_vec_re)
var_vec_re_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  sizescale_pd_df_r = robjects.conversion.py2ri(sizescale_pd_df) # For rpy2 versions >2.9 use py2rpy
sizescale_pd_df_r

#Invoking the R function and getting the result. Note that the sequence of arguments is critical
treelist_result = extract_treelist_r(sam_start, sam_end, outdir, VDM2FM, runroot, filebase, var_vec_re_r, filterFile, finalyear, fire_res, fates_res, fates_CWD_frac_twig, fates_c2b, leafdensity, wooddensity, sizescale_pd_df_r, grass_pft_index, HYDRO, cycle_index)
if (treelist_result):
    print('Treelist extracted successfully at', VDM2FM + "/treelist_VDM.dat")
    exit(0)
else:
   print("Treelist not extracted")
