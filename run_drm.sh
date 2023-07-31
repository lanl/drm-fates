#!/usr/bin/env bash
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w22_fire
#SBATCH -N 1 # number of nodes
#SBATCH -t 01:00:00
#SBATCH -L scratch4@slurmdb

echo 'Started on '`date`' for '$casedir' starting in '$casedir

jobid=$SLURM_JOB_ID

python DRM_framework_coupling.py

echo 'Finished on '`date`
