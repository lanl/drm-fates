"""
# Extract ELM ensemble treelist outputs
# Rutuja Chitra-Tarak
# MMay 14, 2022
"""

# SXM Notes (BGN)
# 1. Deal with 3 PFTs as in the parameter file:
# fates_pftname =
#  "needleleaf_evergreen_extratrop_tree", --> long-leaf pine
#  "broadleaf_colddecid_extratrop_tree ", --> turkey oak
#  "c4_grass                           "  --> grass
# SXM Notes (END)

import pandas as pd
import os
import rpy2
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
import argparse
import yaml
from pandas.core.index import Index as PandasIndex

# determine input parameters
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

# defining the R script and loading the instance in Python
r = robjects.r
r['source'](R_file) #CSXM: to (re)use the functns defined in the R_file (https://www.statology.org/source-function-in-r/)

# loading the function we have defined in R
extract_treelist_r = robjects.globalenv['extract_treelist'] #CSXM: extract_treelist_r can be deemed as the name/handle of the extract_treelist function in extract.treelist.R

sam_start = int(config_dict['SIM_ID_START'])
sam_end = int(config_dict['SIM_ID_END'])
finalyear = int(config_dict['FINAL_TAG_YEAR'])
fire_res = int(config_dict['FIRE_RES'])
fates_res = int(config_dict['FATES_RES'])
cycle_index = int(config_dict['CYCLE_INDEX'])

# HYDRO ON/OFF? if OFF (e.g. currently HYDRO doesn't work iwth grasses), moisture is set to 1.
HYDRO = int(config_dict['HYDRO'])

# FATES parameters
fates_CWD_frac_twig = 0.045 #CSXM: CWD may indicate Coarse Woody Debris
fates_c2b = 2 #CSXM: carbon to biomass multiplier of bulk structural tissues (unit: ratio)
# fates_leaf_slatop = 0.00662, 0.0189200006425381 # for long-leaf pine & Turkey oak #CSXM: specific Leaf Area (SLA) at top of canopy, projected area basis (unit: m^2/gC)
# denleaf  = -2.3231_r8*sla/prt_params%c2b(ft) + 781.899_r8 L863 kg/m3 #CSXM: leaf dry mass per unit fresh leaf volume (unit: kg/m3) and this formula is from FatesPlantHydraulicsMod.F90
leafdensity = -2.3231*(0.00662+0.0189200006425381)/2/2 + 781.899 # kg/m3
# fates_wood_density: 0.58, 0.65, # g/cm3 #CSXM: mean density of woody tissue in plant (unit: g/cm3)
# 0.615 g/cm3 = 615 kg/m3
wooddensity = (0.58 + 0.65)/2 # kg/m3 #CSxM: the value seems for g/cm3

# fire parameter; a version of fuel surface area to volume ratio
sizescale_pd_df = pd.DataFrame({'fates_pft': [1,2,3], #CSXM: 1=long-leaf pine (needleleaf_evergreen_extratrop_tree), 2=turkey oak (broadleaf_colddecid_extratrop_tree), 3=grass (c4_grass)
                                'sizescale': [0.2,0.6,1]})
grass_pft_index = 3

# set the BASE CASE name. this is generated from yaml and src/create.basecase.sh
ff = open(PROJECT_ROOT+"/BASE_CASE_NAME.txt", "r")
base_case = ff.read()
filebase = base_case.strip()
filterFile = "Filter.txt" #CSXM: documents whether the FATES simus are completed (TRUE) or not (FALSE)

#CSXM (BGN)
# the following are from the FATES restarting file
# fates_pft:         plant functional type (index)
# fates_dbh:         diameter at breast height (cm)
# fates_height:      plant height (m)
# fates_crown_depth: plant crown depth fraction (fraction) 
# fates_nplant:      number of plants in the cohort (/patch)
# fates_cohort_area: area of the fates cohort (m2)
# leaf_c_val_001:    leaf         carbon, state var, position:001 (kg)
# sapw_c_val_001:    sapwood      carbon, state var, position:001 (kg)
# store_c_val_001:   storage      carbon, state var, position:001 (kg)
# repro_c_val_001:   reproductive carbon, state var, position:001 (kg)
# struct_c_val_001:  structural   carbon, state var, position:001 (kg)
#CSXM (END)
var_vec_re = ["fates_pft", "fates_dbh", "fates_height", "fates_crown_depth",
              "fates_nplant", "fates_cohort_area", 
              "leaf_c_val_001", "sapw_c_val_001", "store_c_val_001","repro_c_val_001", "struct_c_val_001"]

# converting python objects into r objects for passing into r function
with localconverter(robjects.default_converter + pandas2ri.converter):
  var_vec_re_r = robjects.vectors.StrVector(var_vec_re)
var_vec_re_r

with localconverter(robjects.default_converter + pandas2ri.converter):
  sizescale_pd_df_r = robjects.conversion.py2rpy(sizescale_pd_df) # For rpy2 versions >2.9 use py2rpy
sizescale_pd_df_r

# invoking the R function and getting the result. Note that the sequence of arguments is critical
treelist_result = extract_treelist_r(sam_start, sam_end, outdir, VDM2FM, runroot, filebase, var_vec_re_r, filterFile, finalyear, fire_res, fates_res, fates_CWD_frac_twig, fates_c2b, leafdensity, wooddensity, sizescale_pd_df_r, grass_pft_index, HYDRO, cycle_index)

if(treelist_result):
  print('**********----Treelist extracted successfully at', VDM2FM + "/treelist_VDM.dat----**********")
  exit(0)
else:
  print("**********----Treelist not extracted----**********")
