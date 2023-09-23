#!/bin/tcsh
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w22_fire
#SBATCH -N 1 # number of nodes
#SBATCH -t 00:05:00

. ~/.tcshrc

#conda activate elm_env
conda activate mpi4pyEnv_GNU

setenv casedir /usr/projects/higrad/rutuja/E3SM_cases/proj1/1.FATES-MODEL/CASE_DIR

echo 'Started on '`date`' for '$casedir' starting in '$casedir

cd $casedir

setenv location `date "+%F-%T"`
srun -n 4 python /usr/projects/higrad/rutuja/E3SM_cases/proj1/1.FATES-MODEL/src/parallel.run.py -c 'drm.FATES.api24.12.1.1.0.5.IELMFATES.chicoma.gnu.Cc6963c159f-Fbeceb996d.1990-1991.' -r /usr/projects/higrad/rutuja/E3SM.DRM/runs -f elm.h0.1990-12.nc -s 1 -t 4 -g /usr/projects/higrad/rutuja/E3SM_cases/proj1/1.FATES-MODEL/log


echo 'Finished on '`date`
