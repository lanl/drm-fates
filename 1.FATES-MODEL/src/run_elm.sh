#!/usr/bin/env bash

casedir='/turquoise/usr/projects/higrad/rutuja/DRM/1.FATES-MODEL/CASE_DIR'
echo 'Started on '`date`' for '$casedir' starting in '$casedir

cd $casedir

location=`date "+%F-%T"`
mpiexec -n 2 python /turquoise/usr/projects/higrad/rutuja/DRM/1.FATES-MODEL/src/parallel.run.py -c 'drm.test.1.E3SMv2.IELMFATES.badger.intel.C60051649f4-F38708294.2010-2015.' -r /lustre/scratch4/turquoise/rutuja/E3SM/scratch -f clm2.h0.2012-12.nc -s 1 -t 2 -g /turquoise/usr/projects/higrad/rutuja/DRM/1.FATES-MODEL/log


echo 'Finished on '`date`
