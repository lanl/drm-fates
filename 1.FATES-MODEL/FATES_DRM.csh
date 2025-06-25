#!/bin/sh


eval "$(conda shell.bash hook)"
conda activate /users/zjrobbins/.conda/envs/mpi4pyEnv_GNU

export CLONE_ROOT="/lustre/scratch5/.mdt1/zjrobbins/E3SM_cases/Eglin"
export SRC="/lustre/scratch5/.mdt1/zjrobbins/drm-fates"
export RUN_ROOT="/lustre/scratch5/.mdt1/zjrobbins/E3SM/scratch"
export E3SM_CASE_CLONE="EGLIN_31824_new5.IELMFATES"
export STOP_N=5

rm $RUN_ROOT/$E3SM_CASE_CLONE/run/$E3SM_CASE_CLONE.elm.h0*
rm $RUN_ROOT/$E3SM_CASE_CLONE/run/$E3SM_CASE_CLONE.elm.r*

rm $RUN_5/drm-fates/1.FATES-MODEL/VDM2FM/Pre_fire*
rm $RUN_5/drm-fates/1.FATES-MODEL/FM2VDM/Afterfire*
rm $RUN_5/drm-fates/1.FATES-MODEL/VDM2FM/fuels*

#rm /lustre/scratch5/.mdt1/zjrobbins/E3SM/scratch/EGLIN_31824.IELMFATES/run/EGLIN_31824.IELMFATES.elm.h0.*
#rm /lustre/scratch5/.mdt1/zjrobbins/E3SM/scratch/EGLIN_31824.IELMFATES/run/EGLIN_31824.IELMFATES.elm.r.*


cd $CLONE_ROOT/$E3SM_CASE_CLONE
./xmlchange --file env_run.xml --id CONTINUE_RUN --val FALSE

./case.build --clean-all
./case.build
	
for case_i in $(seq 1 20)
do        
        echo "Running iteration $case_i"
	cd $CLONE_ROOT/$E3SM_CASE_CLONE
#	echo "tooo"
        ./case.submit --no-batch 
        cd $SRC
	python FATES_to_DRM.py -c ${RUN_ROOT}/${E3SM_CASE_CLONE}/run/ -s 1990 -IT $STOP_N -n $case_i -p ${E3SM_CASE_CLONE}
        python DRM_short.py
        cd $CLONE_ROOT/$E3SM_CASE_CLONE
        ./xmlchange --file env_run.xml --id CONTINUE_RUN --val TRUE
        ./xmlchange --file env_run.xml --id STOP_N --val $STOP_N
        ./xmlchange --file env_run.xml --id REST_N --val $STOP_N
        ./xmlchange --file env_run.xml --id STOP_OPTION --val nyears
        ./xmlchange --file env_run.xml --id REST_OPTION --val nyears
        cd $SRC
        python DRM_to_FATES.py -c ${RUN_ROOT}/${E3SM_CASE_CLONE}/run/ -s 1990 -IT $STOP_N -n $case_i -p ${E3SM_CASE_CLONE}

done


rm -r $SRC/outputs/$E3SM_CASE_CLONE/; mkdir $SRC/outputs/$E3SM_CASE_CLONE/

cp -r $SRC/7.QUICFIRE-MODEL/projects/Tester/Plots3/  $SRC/outputs/$E3SM_CASE_CLONE/
cp $RUN_5/drm-fates/1.FATES-MODEL/VDM2FM/fuels*  $SRC/outputs/$E3SM_CASE_CLONE/
cp $RUN_5/drm-fates/1.FATES-MODEL/VDM2FM/Pre_fire*  $SRC/outputs/$E3SM_CASE_CLONE/
cp $RUN_5/drm-fates/1.FATES-MODEL/FM2VDM/Afterfire*  $SRC/outputs/$E3SM_CASE_CLONE/

