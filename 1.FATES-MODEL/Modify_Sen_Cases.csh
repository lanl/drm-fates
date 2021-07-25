#!/bin/tcsh

#-----------------------------------------------------------------------------------------------

# Set the CLONE ROOT Directory
set CLONE_ROOT = $VEGHOME/ACME_cases/BCI;

# Set the BASE CASE Directory
set BASE_CASE = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.met.v6.1971.2020" 
set bestfit_arr = (3563, 3535, 2731, 3653, 3444, 1120, 2173, 1611, 3271, 2451, 2495, 2501, 2708, 1881, 3858, 2298, 1017, 1117, 1217, 1317, 1417, 1517, 4185, 2741, 3089, 3157, 2274, 2541, 3124, 4703, 3639, 1867, 2636, 2527, 1754, 1592, 3241, 2031, 2164, 2371, 2029, 2099, 2578, 2693, 3849, 1937, 3005, 3506, 4983, 1929, 2643, 1789, 4142, 2880, 3133, 2929, 3861, 1758, 4047, 2372, 3743, 1793, 2165, 1703, 2688, 1817, 2973, 3145, 4517, 2049, 3762, 3697, 2190, 2714, 1756, 4893, 2196, 1874, 3077, 3303, 2996, 3347, 2833, 1672, 3588, 3898, 2809, 2766, 2353, 2179, 2666, 3589, 4514, 3361, 2689, 1979, 3450, 2082, 2773, 2805)


# Now loop through and modify each case
foreach case_i (`seq 1 100`) 

  # Go to CASE scripts directory
  cd $ACME_ROOT/cime/scripts

  # Name the case
  set ACME_CASE_CLONE = "BCI.ICLM45ED.badger.intel.C700b46fec-F8c9cd1b0.full.met.v6.$case_i" 
  #rm -rf $CLONE_ROOT/$ACME_CASE_CLONE
  # Create the Clone Case
  #./create_clone -case $CLONE_ROOT/$ACME_CASE_CLONE -clone  $CLONE_ROOT/$BASE_CASE
  # Go to CLONE case directory
  cd $CLONE_ROOT/$ACME_CASE_CLONE
  #./case.setup --reset

  #set new_surfpar_i=$new_surfpar_arr[$case_i]
  #set bestfit_i = $bestfit_arr[$case_i]
  
  #sed -i 's/aitor/rutuja/g' user_nl_clm
  #sed -i 's/'c171113.1.nc'/'c171113.${bestfit_i}.nc'/g' user_nl_clm
  sed -i 's/',.nc'/'.nc'/g' user_nl_clm
  #sed -i /REST_N/s/10/5/ env_run.xml
  #set old_param_i=$old_param_arr[$case_i]
  #set new_param_i=$new_param_arr[$case_i]
  # changes files for both fates param and elm params
  #sed -i 's:'parameter_file_name1.nc':'parameter_file_name${bestfit_i}.nc':g' user_nl_clm
  #sed -i 's+elm.param.sam/elm.param/elm.param.sam.sam/+elm.param.sam/+g' user_nl_clm

  #./xmlchange -file env_run.xml -id RUN_STARTDATE -val '2008-01-01' 
  # Set up the Clone case
  #./case.setup
  #./xmlchange BUILD_COMPLETE="TRUE"
  #./xmlchange LND_DOMAIN_FILE="domain.lnd.1x1_Colorado${case_i}_FATES_1x1_Colorado${case_i}_FATES.161222.nc"
  #./xmlchange ATM_DOMAIN_FILE="domain.lnd.1x1_Colorado${case_i}_FATES_1x1_Colorado${case_i}_FATES.161222.nc"

  # Replace the parameter file for this run
  #set user_nl_clmi = "$CLONE_ROOT/user_nl_clm_files/user_nl_clm00$case_i" 
  #cp $user_nl_clmi $CLONE_ROOT/$ACME_CASE_CLONE/user_nl_clm

  #set depth_i=$soil_depth_arr[$case_i]

  #sed -i 's/'c171113_aveDTB10.nc'/'c171113_aveDTB${depth_i}.nc'/g' user_nl_clm
  #sed -i '/hist_fincl3/d' user_nl_clm
  
  #sed -i '/hist_nhtfrq/d' user_nl_clm
  #sed -i '/hist_mfilt/d' user_nl_clm
  #append with following lines\
  #sed -i '$ a hist_nhtfrq = 0, -24' user_nl_clm
  #sed -i '$ a hist_mfilt = 12, 365' user_nl_clm


  #sed -i 's/hydro${case_i}.nc/'hydro${case_i}.nc'/g' user_nl_clm
  #sed -i '$a hist_fincl4 = 'H2OSOI'' user_nl_clm
  #sed -i '$a hist_fincl3 =  'H2OSOI', 'QRUNOFF', 'QVEGE', 'QVEGT', 'QSOIL', 'GPP'' user_nl_clm
  
  # copy cesm.exe in base case run folder
  #cp $RUN_ROOT/$BASE_CASE/bld/e3sm.exe $RUN_ROOT/$ACME_CASE_CLONE/bld/
  #copy the env.specific.xml
  #cp $CLONE_ROOT/$BASE_CASE/env_mach_specific.xml $CLONE_ROOT/$ACME_CASE_CLONE/  
end
