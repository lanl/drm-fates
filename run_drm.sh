#!/usr/bin/env bash
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w22_fire
#SBATCH -N 1 # number of nodes
#SBATCH -t 00:20:00

echo 'Started on '`date`' for '$casedir' starting in '$casedir

jobid=$SLURM_JOB_ID

python DRM_framework_coupling.py

echo 'Finished on '`date`
