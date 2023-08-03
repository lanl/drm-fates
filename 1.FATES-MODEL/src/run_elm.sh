#!/usr/bin/env bash
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w22_fire
#SBATCH -N 1 # number of nodes
#SBATCH -t 00:05:00
#SBATCH -L scratch4@slurmdb

casedir='/usr/projects/higrad/rutuja/E3SM_cases/proj3/1.FATES-MODEL/CASE_DIR'
echo 'Started on '`date`' for '$casedir' starting in '$casedir

cd $casedir

location=`date "+%F-%T"`
srun -n 4 python /usr/projects/higrad/rutuja/E3SM_cases/proj3/1.FATES-MODEL/src/parallel.run.py -c 'drm.newqf.FATES.api.24.1.0.from.envyml.IELMFATES.chicoma.gnu.Ce4e912868b-Fa2f9f60b3.1990-1991.' -r /usr/projects/higrad/rutuja/E3SM.DRM/runs -f clm2.h0.1990-12.nc -s 1 -t 4 -g /usr/projects/higrad/rutuja/E3SM_cases/proj3/1.FATES-MODEL/log

echo 'Finished on '`date`
