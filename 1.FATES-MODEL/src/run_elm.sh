#!/usr/bin/env bash

casedir='/usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/CASE_DIR'
echo 'Started on '`date`' for '$casedir' starting in '$casedir

cd $casedir

location=`date "+%F-%T"`
mpiexec -n 4 python /usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/src/parallel.run.py -c 'drm.newqf.ch.test.IELMFATES.chicoma.gnu.Ce4e912868b-Fe663a6e6.1990-2021.' -r /lustre/scratch4/turquoise/.mdt5/rutuja/E3SM.DRM/runs -f clm2.h0.1990-12.nc -s 1 -t 4 -g /usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/log


echo 'Finished on '`date`
