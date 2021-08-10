#!/usr/bin/env bash
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w18_califorest
#SBATCH -N 1 # number of nodes
#SBATCH -t 1:00:00

casedir=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`
PY_SRC_PATH=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`

echo 'Started on '`date`' for '$casedir' starting in '$casedir

module purge
module use /usr/projects/cesm/software/local/modulefiles/all
module load  acmeenv/1.0.0
#module load  clmedenv/1.0.0
#module load  intel/15.0.5
module load  intel/17.0.1
#module load  openmpi/1.6.5
#module load  openmpi
#module load  mkl/10.3.11
module load mkl
#module load  cmake/3.5.2
module load  ncview
module load  nco/4.4.4
module load  ncl
#module load  python/2.7-anaconda-2.1.0
module load python/anaconda-2.7-climate
module load  mpi4py/3.0.0
module load netcdf/4.4.1

cd $casedir

date=`date "+%F-%T"`

jobid=$SLURM_JOB_ID

cd $casedir

python $PY_SRC_PATH

echo 'Finished on '`date`
