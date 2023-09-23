#!/bin/tcsh

#-----------------------------------------------------------------------------------------------

# Python script will pass the CLONE ROOT Directory and BASE_CASE
set case_arr=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set BASE_CASE=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set CLONE_ROOT=`echo $3 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set STOP_N=`echo $4 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set REST_N=`echo $4 | sed 's/^[ \t]*//;s/[ \t]*$//'`
set arr_case=`echo $case_arr:q | sed 's/,/ /g'`
##======================================

## Now loop through and create each case
echo "============In total, restarting $#arr_case cases============"

foreach case_i (`seq 1 $#arr_case`)
  # Go to CASE scripts directory
  cd $E3SM_ROOT/cime/scripts
  # Name the case
  set E3SM_CASE_CLONE = $BASE_CASE.$case_i
  # Go to CLONE case directory
  cd $CLONE_ROOT/$E3SM_CASE_CLONE
  ./xmlchange --file env_run.xml --id CONTINUE_RUN --val TRUE
  ./xmlchange --file env_run.xml --id STOP_N --val $STOP_N 
  ./xmlchange --file env_run.xml --id REST_N --val $REST_N
  ./xmlchange --file env_run.xml --id STOP_OPTION --val nyears
  ./xmlchange --file env_run.xml --id REST_OPTION --val nyears
  ./xmlchange --file env_run.xml --id REST_DATE --val -999
end
