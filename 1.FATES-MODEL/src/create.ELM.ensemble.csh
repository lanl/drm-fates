#!/bin/tcsh

#-----------------------------------------------------------------------------------------------

# Python script will pass the CLONE ROOT Directory and BASE_CASE
set case_arr=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set BASE_CASE=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set CLONE_ROOT=`echo $3 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set RUN_ROOT=`echo $4 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set MUTATE=`echo $5 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set clone_base=`echo $6 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set clone_file=`echo $7 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set clone_type=`echo $8 | sed 's/^[ \t]*//;s/[ \t]*$//'`

set clone_type_len=`echo $clone_type:q | sed 's/,/ /g'`
set arr_case=`echo $case_arr:q | sed 's/,/ /g'`
##======================================

# Now loop through and create each case 
foreach case_i (`seq 1 2)
#foreach case_i (`seq 1 $#arr_case`)
  # Go to CASE scripts directory
  cd $ACME_ROOT/cime/scripts
  # Name the case
  set ACME_CASE_CLONE = $BASE_CASE.$case_i 
  # Delete a case by the same name if previously present
  rm -rf $CLONE_ROOT/$ACME_CASE_CLONE
  # Create the Clone Case
 ./create_clone -case $CLONE_ROOT/$ACME_CASE_CLONE -clone  $CLONE_ROOT/$BASE_CASE
  # Go to CLONE case directory
  cd $CLONE_ROOT/$ACME_CASE_CLONE

  # Set up the Clone case
  if ($MUTATE == TRUE) then
        # Replace the parameter file for this run
        set case_num = $arr_case[$case_i]
           	foreach clone_type_i (`seq 1 $#clone_type_len`)
echo ${clone_file}[$clone_type_i]
echo clone_type_i
echo ${clone_type}[$clone_type_i]
echo ${clone_base}[$clone_type_i]${case_num}

        	#sed -i 's/'${clone_base}nc'/'${clone_base}${case_num}.nc'/g' user_nl_clm
#        	sed -i 's:'${clone_file}[$clone_type_i]':'${clone_type}[$clone_type_i].params/${clone_base}[$clone_type_i]${case_num}.nc':g' user_nl_clm
		#sed -i "s|'params/${clone_file}'|'params/${clone_type}.params/${clone_base}${case_num}.nc'|g" user_nl_clm
        	#sed -i 's:'parameter_file_name1.nc':'parameter_file_name${bestfit_i}.nc':g' user_nl_clm
  		endif
  endif
  ./case.setup
  ./xmlchange BUILD_COMPLETE="TRUE"
  #./xmlchange BATCH_SYSTEM="slurm" --file env_batch.xml
  #./xmlchange BATCH_SYSTEM="slurm" --file LockedFiles/env_batch.xml
  
  # copy cesm.exe in base case run folder
  # this is useful if FATES files are changed in the interim (they would have to be built for the base case though and then are copied over here)
  #rm -rf $RUN_ROOT/$ACME_CASE_CLONE
  cp $RUN_ROOT/$BASE_CASE/bld/e3sm.exe $RUN_ROOT/$ACME_CASE_CLONE/bld/ 
  echo "============From a base case successfully created a clone with a case name prefix = $case_i============"
end
echo "============In total created $#arr_case clones from a base case============"
