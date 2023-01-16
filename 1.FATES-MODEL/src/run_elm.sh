#!/usr/bin/env bash

casedir='/usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/CASE_DIR'
echo 'Started on '`date`' for '$casedir' starting in '$casedir

cd $casedir

location=`date "+%F-%T"`
mpiexec -n 4 python /usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/src/parallel.run.py -c 'drm.newtrees.test.ch.IELMFATES.badger.intel.C60051649f4-F3b085a322.1990-2021.' -r /lustre/scratch4/turquoise/rutuja/E3SM/scratch -f clm2.h0.1999-12.nc -s 1 -t 4 -g /usr/projects/higrad/rutuja/newqf/1.FATES-MODEL/log


echo 'Finished on '`date`
