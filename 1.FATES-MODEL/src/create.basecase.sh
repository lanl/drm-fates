#!/bin/sh
# =======================================================================================
#
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
# offsets.  We make use of "CLM1PT" datm mode, and make a minor modification to the
# default stream--file to signal to CLM/ELM that no downwelling long-wave is available
# in our dataset.
#
# The base-version of this script, works off the assumption that driver data
# and surface/domain data have been unpaciked in the cime/scripts directory
# and can all be found in a parent directory with alias $SITE_DIR
#
#
# Ryan Knox (Mon Nov 13 13:53:03 PST 2017)

# Updated by Rutuja Chitra-Tarak
# READ YAML VARIABLES
# =======================================================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../../tools/yaml.sh"`
#parse_yaml "$SCRIPT_DIR/../config.yaml"
create_variables "$SCRIPT_DIR/../../config.yaml"
PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
RUNROOT=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`
ARCHIVEROOT=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`

# USER SETTINGS
# USER MAY ALSO WANT TO ADJUST XML CHANGES, AND NAMELIST ARGUMENTS
# =======================================================================================

export SITE_NAME=$CLIM_DATA_DIR	                        	# Name of folder with site data
export SITE_BASE_DIR=$PROJECT_ROOT/$DATA_DIR                	# Where is the site folder located? (SITE_NAME)
export TAG=$TAG                                            	# User defined tag to differentiate runs
export COMPSET=$COMPSET                                    	# Compset (probably ICLM45ED or ICLM50ED)
export MAC=$MACHINE                                           	# Name your machine
export COMPILER=gnu	                                  	# Name your compiler
export CASEROOT=$PROJECT_ROOT/$CASE_DIR                 	# Where the build is generated (probably on scratch partition)
export ELM_USRDAT_DOMAIN="$PARAM_FILE_DOMAIN"   		# Name of domain file in scripts/${SITE_DIR}/
export ELM_USRDAT_SURDAT="$PARAM_FILE_SURF"			# Name of surface file in scripts/${SITE_DIR}/
export PARAM_BASE_DIR=$PARAM_DIR

# DEPENDENT PATHS AND VARIABLES (USER MIGHT CHANGE THESE..)
# =======================================================================================
export CLM_SURFDAT_DIR=$PROJECT_ROOT/$PARAM_BASE_DIR
export CLM_DOMAIN_DIR=$PROJECT_ROOT/$PARAM_BASE_DIR
export DIN_LOC_ROOT_FORCE=${SITE_BASE_DIR}
export FATES_PARAM_DIR=$PROJECT_ROOT/$PARAM_BASE_DIR  	#location of FATES parameter file
export ELM_PARAM_DIR=$PROJECT_ROOT/$PARAM_BASE_DIR   	#location of ELM parameter file

export CLM_HASH=`cd ${E3SM_ROOT}/components/elm/;git log -n 1 --pretty=%h`
export FATES_HASH=`(cd ${E3SM_ROOT}/components/elm/src/external_models/fates;git log -n 1 --pretty=%h)`
export GIT_HASH=C${CLM_HASH}-F${FATES_HASH}
export RES=ELM_USRDAT
export CASE_NAME=${TAG}.${COMPSET}.${MAC}.${COMPILER}.${GIT_HASH}.$DATM_CLMNCEP_YR_START-$DATM_CLMNCEP_YR_END
export FATES_PARAM="$PARAM_FILE_FATES" # Name of FATES parameter file in FATES_PARAM_DIR
export ELM_PARAM="$PARAM_FILE_ELM" # Name of ELM parameter file in ELM_PARAM_DIR

echo $CASE_NAME > BASE_CASE_NAME.txt

#REMOVE EXISTING CASE IF PRESENT
rm -r ${CASEROOT}/${CASE_NAME}

# CREATE THE CASE
${E3SM_ROOT}/cime/scripts/create_newcase -case ${CASEROOT}/${CASE_NAME} -res ${RES} -compset ${COMPSET} -mach ${MAC} -compiler ${COMPILER} -mpilib="mpi-serial" --user-mods-dir ${CASEROOT}/${CASE_NAME}

cd ${CASEROOT}/${CASE_NAME} 

echo `pwd`
# SET PATHS TO SCRATCH ROOT, DOMAIN AND MET DATA (USERS WILL PROB NOT CHANGE THESE)
# =================================================================================

./xmlchange --file env_run.xml --id ATM_DOMAIN_FILE --val ${ELM_USRDAT_DOMAIN}
./xmlchange --file env_run.xml --id ATM_DOMAIN_PATH --val ${CLM_DOMAIN_DIR}
./xmlchange --file env_run.xml --id LND_DOMAIN_FILE --val ${ELM_USRDAT_DOMAIN}
./xmlchange --file env_run.xml --id LND_DOMAIN_PATH --val ${CLM_DOMAIN_DIR}
./xmlchange --file env_run.xml --id DATM_MODE --val CLM1PT
./xmlchange --file env_run.xml --id ELM_USRDAT_NAME --val ${SITE_NAME}
./xmlchange --file env_run.xml --id DIN_LOC_ROOT_CLMFORC --val ${DIN_LOC_ROOT_FORCE}
#./xmlchange --file env_build.xml --id CESMSCRATCHROOT --val ${CASE_NAME}

# SPECIFY PE LAYOUT FOR SINGLE SITE RUN (USERS WILL PROB NOT CHANGE THESE)
# =================================================================================

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

# SPECIFY RUN TYPE PREFERENCES (USERS WILL CHANGE THESE)
# =================================================================================

./xmlchange --file env_build.xml --id DEBUG --val FALSE
./xmlchange --file env_run.xml --id STOP_N --val $STOP_N
./xmlchange --file env_run.xml --id RUN_STARTDATE --val $RUN_STARTDATE
./xmlchange --file env_run.xml --id STOP_OPTION --val nyears
./xmlchange --file env_run.xml --id REST_N --val $REST_N
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_START --val $DATM_CLMNCEP_YR_START
./xmlchange --file env_run.xml --id DATM_CLMNCEP_YR_END --val $DATM_CLMNCEP_YR_END

# ========
# MACHINE SPECIFIC, AND/OR USER PREFERENCE CHANGES (USERS WILL CHANGE THESE)
# =================================================================================

#./xmlchange --file env_build.xml --id GMAKE --val make
#./xmlchange --file env_run.xml --id BATCHQUERY --val ''
#./xmlchange --file env_run.xml --id BATCHSUBMIT --val ''
#./xmlchange --file env_run.xml --id DOUT_S_SAVE_INTERIM_RESTART_FILES --val TRUE
#./xmlchange --file env_run.xml --id DOUT_S --val TRUE

./xmlchange --file env_run.xml --id DOUT_S_ROOT --val ${ARCHIVEROOT}
#./xmlchange --file env_run.xml --id RUNDIR --val ${RUN_ROOT}/${CASE_NAME}/run # removed to use the default
#./xmlchange --file env_build.xml --id EXEROOT --val ${RUN_ROOT}/${CASE_NAME}/bld # removed to use the default
./xmlchange --file env_build.xml --id CIME_OUTPUT_ROOT --val ${RUN_ROOT}

# SPECIFY INPUT DATA (USERS WILL CHANGE THESE)
# =================================================================================
./xmlchange  DIN_LOC_ROOT=/usr/projects/higrad/e3sm_input_data/input_data

# ========
# MODIFY THE CLM NAMELIST (USERS MODIFY AS NEEDED)
# =================================================================================

cat >> user_nl_elm <<EOF
fsurdat = '${CLM_SURFDAT_DIR}/${ELM_USRDAT_SURDAT}'
fates_paramfile = '${FATES_PARAM_DIR}/${FATES_PARAM}'
!paramfile = '${ELM_PARAM_DIR}/${ELM_PARAM}'
maxpatch_pft = 17
fates_spitfire_mode = 1
use_fates_logging = .false.
!use_fates_planthydro = .true.
!use_fates_sp = .true.
!use_fates_ed_st3 = .false.
!use_var_soil_thick = .true.
!hist_empty_htapes = .true.
use_fates_inventory_init = .false.
!fates_inventory_ctrl_filename = '${SITE_BASE_DIR}/${SITE_NAME}/eglin_inv_file_list.txt'
hist_fincl2 = 'SOILWATER_10CM','H2OSOI', 'QRUNOFF', 'QOVER', 'QCHARGE', 'QDRAI', 'RAIN', 'QINTR', 'QDRIP', 'QVEGE', 'QVEGT', 'QSOIL', 'TWS', 'ZWT', 'BTRAN'
hist_nhtfrq = 0, -24
hist_mfilt = 1, 365
EOF

# Useful user_nl_clm arguments: 
# This couplet will enable hourly output
# hist_mfilt             = 480      
# hist_nhtfrq            = -1  
# #hist_fincl1='NEP','NPP','GPP','TLAI','TSOI_10CM','QVEGT','EFLX_LH_TOT','AR','HR','ED_biomass','ED_bleaf','ED_balive','DDBH_S#CPF','BA_SCPF','NPLANT_SCPF','M1_SCPF','M2_SCPF','M3_SCPF','M4_SCPF','M5_SCPF','M6_SCPF','WIND','ZBOT','FSDS','RH','TBOT','P#BOT','QBOT','RAIN','FLDS'

#hist_fincl3 = 'H2OSOI', 'SOILPSI', 'ED_biomass', 'NEP', 'NPP', 'GPP', 'TV', 'C_LBLAYER', 'C_STOMATA','EFLX_LH_TOT', 'WIND', 'ZBOT', 'FSA', 'PARVEGLN', 'FSDS', 'FLDS', 'RH', 'TBOT', 'QBOT', 'PBOT', 'RAIN', 'QRUNOFF', 'QVEGE', 'QVEGT', 'QSOIL', 'TWS'
#hist_nhtfrq = 0, -24, -1
#hist_mfilt = 1, 365, 8760

# MODIFY THE DATM NAMELIST (DANGER ZONE - USERS BEWARE CHANGING)
echo `pwd`

cat >> user_nl_datm <<EOF
taxmode = "cycle", "cycle", "cycle"
EOF

./case.setup
# HERE WE NEED TO MODIFY THE STREAM FILE (DANGER ZONE - USERS BEWARE CHANGING)
./preview_namelists
cp ${RUN_ROOT}/${CASE_NAME}/run/datm.streams.txt.CLM1PT.ELM_USRDAT user_datm.streams.txt.CLM1PT.ELM_USRDAT
`sed -i '/FLDS/d' user_datm.streams.txt.CLM1PT.ELM_USRDAT`
`sed -i "s/CLM1PT_data/$CLIM_DATA/" user_datm.streams.txt.CLM1PT.ELM_USRDAT`

./case.build
