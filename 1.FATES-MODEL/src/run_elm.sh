#!/usr/bin/env bash
#SBATCH -o slurm%j.out
#SBATCH -e slurm%j.err
#SBATCH -A w18_califorest
#SBATCH -N 1 # number of nodes
#SBATCH -t 1:00:00

#casedir=`echo $1 | sed 's/^[ \t]*//;s/[ \t]*$//'`
#PY_SRC_PATH=`echo $2 | sed 's/^[ \t]*//;s/[ \t]*$//'`

casedir='/turquoise/usr/projects/veg/rutuja/ACME_cases/ELM_Disease/CASE_DIR'
echo 'Started on '`date`' for '$casedir' starting in '$casedir

module purge
module use /usr/projects/cesm/software/local/modulefiles/all
module load  acmeenv/1.0.0
module load  intel/17.0.1
module load mkl
module load  ncview
module load  nco/4.4.4
module load  ncl
module load python/anaconda-2.7-climate
module load  mpi4py/3.0.0
module load netcdf/4.4.1

cd $casedir

date=`date "+%F-%T"`

jobid=$SLURM_JOB_ID

cd $casedir

mpiexec -n 2 python /turquoise/usr/projects/veg/rutuja/ACME_cases/ELM_Disease/src/parallel.run.py -c 'HPU.TEST.ICLM45ED.badger.intel.C4a252d21f-F8c9cd1b0d.2017-2018'. -r /lustre/scratch3/turquoise/rutuja/ACME/cases -f clm2.h0.2018-12.nc -s 1 -t 2 -g /turquoise/usr/projects/veg/rutuja/ACME_cases/ELM_Disease/log

#python $PY_SRC_PATH

echo 'Finished on '`date`
