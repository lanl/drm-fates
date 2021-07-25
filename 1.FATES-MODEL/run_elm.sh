
currentdir=`pwd`'/'
echo 'Started on '`date`' for '$currentdir' starting in '$currentdir

cd $currentdir

date=`date "+%F-%T"`

jobid=$SLURM_JOB_ID

cd $currentdir

python 8.elm.py

echo 'Finished on '`date`
