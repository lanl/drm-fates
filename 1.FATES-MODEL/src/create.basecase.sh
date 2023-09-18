#!/bin/sh
# =======================================================================================
# This script will create, setup and build a single-site simulation at 
# Barro Colorado Island, Panama
#
# A specialized domain and surface file was generated using Gautam Bisht's
# "matlab-script-for-clm-sparse-grid".
#
# Meteorological driving data was prepared using Ryan Knox's "ConvertMetCSVtoCLM".
# PLEASE SEE THE DRIVER DATA'S METADATA FOR ACKNOWLEDGMENTS AND ATTRIBUTION
#
# In this implementation, meteorological data will be cycled based on the data
# offsets. We make use of "CLM1PT" datm mode, and make a minor modification to the
# default stream--file to signal to CLM/ELM that no downwelling long-wave is available
# in our dataset.
##
# The base-version of this script, works off the assumption that driver data
# and surface/domain data have been unpaciked in the cime/scripts directory
# and can all be found in a parent directory with alias $SITE_DIR
#
# Ryan Knox (Mon Nov 13 13:53:03 PST 2017)
# Updated by Rutuja Chitra-Tarak

# modified by SXM at LANL (Apr 2023)

#++++++++++++++++++++
# READ YAML VARIABLES
#++++++++++++++++++++
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" # Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../../tools/yaml.sh"`
# parse_yaml "$SCRIPT_DIR/../config.yaml"
create_variables "$SCRIPT_DIR/../../config.yaml"
PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
RUNROOT=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`       # CSXM: set in .tcshrc as RUN_ROOT      and passed here through create.basecase.py as the 1st argument
ARCHIVEROOT=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`   # CSXM: set in .tcshrc as ARCHIVE_ROOT  and passed here through create.basecase.py as the 2nd argument
DIN_LOC_ROOTF=`echo $3 | sed 's/^[ \t]*//;s/[ \t]*$//'` # CSXM: set in .tcshrc as DIN_LOC_ROOTF and passed here through create.basecase.py as the 3rd argument

#++++++++++++++
# USER SETTINGS
#++++++++++++++
# USER MAY ALSO WANT TO ADJUST XML CHANGES, AND NAMELIST ARGUMENTS
export TAG=$TAG                                         # TAG             is set in config.yaml (to differentiate runs)
export COMPSET=$COMPSET                                 # COMPSET         is set in config.yaml
export MAC=$MACHINE                                     # MACHINE         is set in config.yaml
export COMPILER=$COMPILER                               # COMPILER        is set in config.yaml
export DOMAIN_ROOT=$PROJECT_ROOT/$DOMAIN_DIR            # DOMAIN_DIR      is set in config.yaml
export SURF_ROOT=$PROJECT_ROOT/$SURF_DIR                # SURF_DIR        is set in config.yaml
export PARAM_ELM_ROOT=$PROJECT_ROOT/$PARAM_ELM_DIR      # PARAM_ELM_DIR   is set in config.yaml
export PARAM_FATES_ROOT=$PROJECT_ROOT/$PARAM_FATES_DIR  # PARAM_FATES_DIR is set in config.yaml
export IC_DATA_ROOT=$PROJECT_ROOT/$IC_DATA_DIR          # IC_DATA_DIR     is set in config.yaml
export CLIM_DATA_ROOT=$PROJECT_ROOT/$CLIM_DATA_DIR      # CLIM_DATA_DIR   is set in config.yaml
export CASEROOT=$PROJECT_ROOT/$CASE_DIR                 # CASE_DIR        is set in config.yaml (where to the built is generated, probably on scratch)
export MAXPFT=$MAXPFT                                   # MAXPFT          is set in config.yaml
export CLIM_DATA_LINE=$CLIM_DATA_ROOT/$IC_DATA_DIR/CLM1PT_data

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# DEPENDENT PATHS AND VARIABLES (USER MIGHT CHANGE THESE..)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
export CLM_HASH=`cd ${E3SM_ROOT}/components/elm/;git log -n 1 --pretty=%h`                                    # E3SM_ROOT is set in .tcshrc
export FATES_HASH=`(cd ${E3SM_ROOT}/components/elm/src/external_models/fates;git log -n 1 --pretty=%h)`       # E3SM_ROOT is set in .tcshrc
export GIT_HASH=C${CLM_HASH}-F${FATES_HASH}                                                                   # set above
export RES=ELM_USRDAT                                                                                         # refers to resolution and part of rof_grid_supported in the $SrcCODE
export CASE_NAME=${TAG}.${COMPSET}.${MAC}.${COMPILER}.${GIT_HASH}.$DATM_CLMNCEP_YR_START-$DATM_CLMNCEP_YR_END # set either in config.yaml, above or .tcshrc
echo $CASE_NAME > BASE_CASE_NAME.txt                                                                          # record $CASE_NAME in a .txt file
rm -r ${CASEROOT}/${CASE_NAME}                                                                                # REMOVE EXISTING CASE IF PRESENT

#++++++++++++++++
# CREATE THE CASE
#++++++++++++++++
# CSXM (BGN)
# 1. --user-mods-dir
#   Full pathname to a directory containing any combination of user_nl_* files and a shell_commands script (typically containing xmlchange commands). 
#   The directory can also contain an SourceMods/ directory with the same structure as would be found in a case directory.
#   It can also contain a file named 'include_user_mods' which gives the path to one or more other directories that should be included.
#   Multiple directories can be given to the --user-mods-dirs argument, in which case changes from all of them are applied.
#   (If there are conflicts, later directories take precedence.)
#   (Care is needed if multiple directories include the same directory via 'include_user_mods': in this case, the included directory will be applied multiple times.)
# 2. see https://esmci.github.io/cime/versions/master/html/Tools_user/create_newcase.html for details. 
# CSXM (END)
${E3SM_ROOT}/cime/scripts/create_newcase --case ${CASEROOT}/${CASE_NAME} --res ${RES} --compset ${COMPSET} --mach ${MAC} --compiler ${COMPILER} --mpilib="mpi-serial" --user-mods-dir ${CASEROOT}/${CASE_NAME} --project ${ACCOUNT} --q standard --walltime ${WALL_TIME}
cd ${CASEROOT}/${CASE_NAME} 

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SET PATHS TO SCRATCH ROOT, DOMAIN AND MET DATA (USERS WILL PROB NOT CHANGE THESE)
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
echo `pwd`
./xmlchange --file env_run.xml --id ATM_DOMAIN_FILE --val ${DOMAIN_FILE}         # DOMAIN_FILE is set in config.yaml
./xmlchange --file env_run.xml --id ATM_DOMAIN_PATH --val ${DOMAIN_ROOT}         # set above
./xmlchange --file env_run.xml --id LND_DOMAIN_FILE --val ${DOMAIN_FILE}         # DOMAIN_FILE is set in config.yaml
./xmlchange --file env_run.xml --id LND_DOMAIN_PATH --val ${DOMAIN_ROOT}         # set above
./xmlchange --file env_run.xml --id DATM_MODE --val CLM1PT                       # data-atmosphere mode
./xmlchange --file env_run.xml --id ELM_USRDAT_NAME --val ${IC_DATA_DIR}         # IC_DATA_DIR is set in config.yaml 
./xmlchange --file env_run.xml --id DIN_LOC_ROOT_CLMFORC --val ${CLIM_DATA_ROOT} # set above 

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SPECIFY PE LAYOUT FOR SINGLE SITE RUN (USERS WILL PROB NOT CHANGE THESE)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
./xmlchange NTASKS_ATM=1
./xmlchange NTASKS_CPL=1
./xmlchange NTASKS_GLC=1
./xmlchange NTASKS_OCN=1
./xmlchange NTASKS_WAV=1
./xmlchange NTASKS_ICE=1
./xmlchange NTASKS_LND=1
./xmlchange NTASKS_ROF=1
./xmlchange NTASKS_ESP=1
./xmlchange ROOTPE_ATM=0
./xmlchange ROOTPE_CPL=0
./xmlchange ROOTPE_GLC=0
./xmlchange ROOTPE_OCN=0
./xmlchange ROOTPE_WAV=0
./xmlchange ROOTPE_ICE=0
./xmlchange ROOTPE_LND=0
./xmlchange ROOTPE_ROF=0
./xmlchange ROOTPE_ESP=0
./xmlchange NTHRDS_ATM=1
./xmlchange NTHRDS_CPL=1
./xmlchange NTHRDS_GLC=1
./xmlchange NTHRDS_OCN=1
./xmlchange NTHRDS_WAV=1
./xmlchange NTHRDS_ICE=1
./xmlchange NTHRDS_LND=1
./xmlchange NTHRDS_ROF=1
./xmlchange NTHRDS_ESP=1

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# SPECIFY RUN TYPE PREFERENCES (USERS WILL CHANGE THESE)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CSXM (BGN)
# 1. DATM_CLMNCEP_YR_ALIGN: I compsets only - simulation year corresponding to data starting year
# 2. DATM_CLMNCEP_YR_START: I compsets only - data model starting year to loop data over
# 3. DATM_CLMNCEP_YR_END:   I compsets only - data model ending year to loop data over
# Note: see http://esmci.github.io/cime/versions/ufs_release_v1.1/html/data_models/data-atm.html for details
# CSXM (END)
./xmlchange --file env_build.xml --id DEBUG --val FALSE
./xmlchange --file env_run.xml --id STOP_N --val $STOP_N                               # STOP_N                is set in config.yaml
./xmlchange --file env_run.xml --id RUN_STARTDATE --val $RUN_STARTDATE                 # RUN_STARTDATE         is set in config.yaml
./xmlchange --file env_run.xml --id STOP_OPTION --val nyears                           # nyears                is set in DRM_framework_coupling.py
./xmlchange --file env_run.xml --id REST_N --val $REST_N                               # REST_N                is set in config.yaml
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_START --val $DATM_CLMNCEP_YR_START # DATM_CLMNCEP_YR_START is set in config.yaml
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val $DATM_CLMNCEP_YR_END     # DATM_CLMNCEP_YR_END   is set in config.yaml

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# MACHINE SPECIFIC, AND/OR USER PREFERENCE CHANGES (USERS WILL CHANGE THESE)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#./xmlchange --file env_build.xml --id GMAKE --val make
#./xmlchange --file env_run.xml --id BATCHQUERY --val ''
#./xmlchange --file env_run.xml --id BATCHSUBMIT --val ''
#./xmlchange --file env_run.xml --id DOUT_S_SAVE_INTERIM_RESTART_FILES --val TRUE
#./xmlchange --file env_run.xml --id DOUT_S --val TRUE
./xmlchange --file env_run.xml --id DOUT_S_ROOT --val ${ARCHIVEROOT}              # set above
#./xmlchange --file env_run.xml --id RUNDIR --val ${RUN_ROOT}/${CASE_NAME}/run    # removed to use the default
#./xmlchange --file env_build.xml --id EXEROOT --val ${RUN_ROOT}/${CASE_NAME}/bld # removed to use the default
./xmlchange --file env_build.xml --id CIME_OUTPUT_ROOT --val ${RUN_ROOT}          # set above

#+++++++++++++++++++++++++++++++++++++++++++++
# SPECIFY INPUT DATA (USERS WILL CHANGE THESE)
#+++++++++++++++++++++++++++++++++++++++++++++
./xmlchange DIN_LOC_ROOT=${DIN_LOC_ROOTF} # set above

# ASXM (BGN)
# specify WALLCLOCK_TIME
./xmlchange --file env_workflow.xml --id JOB_WALLCLOCK_TIME --val ${WALL_TIME}
# ASXM (END)

#+++++++++++++++++++++++++++++++++++++++++++++++++
# MODIFY THE CLM NAMELIST (USERS MODIFY AS NEEDED)
#+++++++++++++++++++++++++++++++++++++++++++++++++
# CSXM: see iCloud "FATES: Namelists" for details
# commented out by SXM (BGN)
# cat >> user_nl_elm <<EOF
# fsurdat = '${ELM_SURFDAT_ROOT}/${SURF_FILE}'
# fates_paramfile = '${PARAM_FATES_ROOT}/${PARAM_FATES_FILE}'
# ! paramfile = '${ELM_PARAM_ROOT}/${PARAM_FILE_ELM}'
# maxpatch_pft = 17
# fates_spitfire_mode = 1
# use_fates_logging = .false.
# ! use_fates_planthydro = .true.
# ! use_fates_sp = .true.
# ! use_fates_ed_st3 = .false.
# ! use_var_soil_thick = .true.
# ! hist_empty_htapes = .true.
# use_fates_inventory_init = .false.
# ! fates_inventory_ctrl_filename = '${IC_DATA_ROOT}/eglin_inv_file_list.txt'
# hist_fincl2 = 'SOILWATER_10CM','H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'TWS', 'ZWT', 'BTRAN'
# hist_nhtfrq = 0, -24
# hist_mfilt = 1, 365
# EOF
# commented out by SXM (END)
# ASXM (BGN)
if [ ${MAXPFT} == 2 ]
then 
  export inventoryFile="eglin_inv_file_list_wo_grass.txt"
fi
if [ ${MAXPFT} == 3 ]
then
  export inventoryFile="eglin_inv_file_list.txt"
fi

#Update inventory file paths
export psscssFilepaths="${IC_DATA_ROOT}/Eglin_default.lat30.547lon-86.639.pss ${IC_DATA_ROOT}/Eglin_default.lat30.547lon-86.639.css"
# Replace third space with a unique string, then replace everythign after that string with the file paths
sed -i "2s|\s\+| XXX|3" ${IC_DATA_ROOT}/${inventoryFile}
sed -i "s|XXX.*|$psscssFilepaths|1" ${IC_DATA_ROOT}/${inventoryFile}

histDRM="'SOILWATER_10CM','H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'TWS', 'ZWT', 'BTRAN'"
histCX="'TLAI'"
histZJB="'FATES_NPLANT_SZ', 'FATES_CANOPYAREA_HT', 'FATES_BASALAREA_SZ', 'FATES_LEAFAREA_HT'"
histZJBS01="'FATES_LAI_CANOPY_SZ', 'FATES_LAI_USTORY_SZ', 'FATES_NPLANT_CANOPY_SZ', 'FATES_NPLANT_USTORY_SZ', 'FATES_VEGC_ABOVEGROUND_SZ'"
histZJBS02="'FATES_GPP_PF', 'FATES_NPLANT_PF', 'FATES_NPP_PF', 'FATES_STOREC_PF', 'FATES_VEGC_PF'"
hist_fincl2SXM="${histDRM}, ${histCX}, ${histZJB}, ${histZJBS01}, ${histZJBS02}"
cat >> user_nl_elm <<EOF
! domain file
! set in env_run.xml using ATM_DOMAIN_FILE, LND_DOMAIN_FILE and etc.

! surface data
fsurdat = '${SURF_ROOT}/${SURF_FILE}'

! initial condition
maxpatch_pft = ${MAXPFT}
use_fates_inventory_init = .true.
fates_inventory_ctrl_filename = '${IC_DATA_ROOT}/${inventoryFile}'

! parameter file
fates_paramfile = '${PARAM_FATES_ROOT}/${PARAM_FATES_FILE}'

! physics setups
fates_spitfire_mode = 1
use_fates_logging = .false.

! hist fields
hist_fincl2 = ${hist_fincl2SXM} 
hist_nhtfrq = 0, -24
hist_mfilt  = 1, 365

use_fates_planthydro = .true.
! options considered but not used by Rutuja
! use_fates_planthydro = .true.
! use_fates_sp = .true.
! use_fates_ed_st3 = .false.
! use_var_soil_thick = .true.
! hist_empty_htapes = .true.
EOF
# ASXM (END)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Useful user_nl_clm arguments: this couplet will enable hourly output
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# hist_mfilt   = 480      
# hist_nhtfrq  = -1  
# #hist_fincl1 ='NEP','NPP','GPP','TLAI','TSOI_10CM','QVEGT','EFLX_LH_TOT','AR','HR','ED_biomass','ED_bleaf','ED_balive','DDBH_S#CPF','BA_SCPF','NPLANT_SCPF','M1_SCPF','M2_SCPF','M3_SCPF','M4_SCPF','M5_SCPF','M6_SCPF','WIND','ZBOT','FSDS','RH','TBOT','P#BOT','QBOT','RAIN','FLDS'
#hist_fincl3 = 'H2OSOI', 'SOILPSI', 'ED_biomass', 'NEP', 'NPP', 'GPP', 'TV', 'C_LBLAYER', 'C_STOMATA','EFLX_LH_TOT', 'WIND', 'ZBOT', 'FSA', 'PARVEGLN', 'FSDS', 'FLDS', 'RH', 'TBOT', 'QBOT', 'PBOT', 'RAIN', 'QRUNOFF', 'QVEGE', 'QVEGT', 'QSOIL', 'TWS'
#hist_nhtfrq = 0, -24, -1
#hist_mfilt = 1, 365, 8760

# MODIFY THE DATM NAMELIST (DANGER ZONE - USERS BEWARE CHANGING)
# CSXM (BGN)
# 1. taxmode is the time axis mode. For CLM we usually have it set to cycle which means that once the end of the data is reached it will start over at the beginning (https://bb.cgd.ucar.edu/cesm/threads/error-about-taxmode.7246/)
# 2. the setup of taxmode depends on the number of datastreams. For instance, for four datastreams (https://www2.cesm.ucar.edu/models/cesm1.2/clm/models/lnd/clm/doc/UsersGuide/x12979.html)
#   taxmode = 'cycle','cycle','cycle','cycle'
# 3. for official description of taxmode, see $SrcCODE/components/data_comps/datm/cime_config/namelist_definition_datm.xml
# CSXM (END)
echo `pwd`
cat >> user_nl_datm <<EOF
taxmode = "cycle", "cycle", "cycle"
EOF

./case.setup

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# HERE WE NEED TO MODIFY THE STREAM FILE (DANGER ZONE - USERS BEWARE CHANGING)
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
./preview_namelists
cp ${RUN_ROOT}/${CASE_NAME}/run/datm.streams.txt.CLM1PT.ELM_USRDAT user_datm.streams.txt.CLM1PT.ELM_USRDAT
`sed -i '/FLDS/d' user_datm.streams.txt.CLM1PT.ELM_USRDAT`
# `sed -i "s/$CLIM_DATA_LINE/$CLIM_DATA_ROOT" user_datm.streams.txt.CLM1PT.ELM_USRDAT` # commented out by SXM
`sed -i "s|$CLIM_DATA_LINE|$CLIM_DATA_ROOT|" user_datm.streams.txt.CLM1PT.ELM_USRDAT`

./case.build
