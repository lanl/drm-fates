This project will create and run point (e.g. each HPU) simulations of ELM-FATES in parallel and extract outputs.

0. To set-up this project on LANL HPC, do this and follow instructions in README.md:
ssh user@wtrw.lanl.gov

ssh ba-fe

cd /usr/projects/veg/$user

mkdir -p ELM_cases

cd ELM_cases

! Before cloning repo, you may need to ssh auntheticate, if you haven't done that already:
! Follow instructions at: https://docs.github.com/en/authentication/connecting-to-github-with-ssh

rm -rf ELM_Workflow

git clone git@github.com:rutujact/ELM_Workflow.git

cd ELM_Workflow

!vim README.md

1. First set variables in config.yaml (Pre-configured for a toy model run, so no changes required.)

2. Set environmental variables:

source src/.tcshrc

3. cd into ELM_Disease, then run once to create and activate a conda environment (takes ~10-15 min) and load ELM and FATES branches:

sh src/run_once.sh

4. Activate conda environment

! sh src/activate_env.sh 
! Somehow the conda activate command throws an error (needs fixing), so run it outside, provided that the environment name matches the one created (xx_env):

conda activate elm_env

5. Then to set-up and run ELM simulations in parallel on the back node, run: 

sh src/run_elm_parallel.sh

6. Often the first attempt of running a parallel simulation is unsuccessful. Run this script until you get "sbatch has finished running" & "a simulation ran successfully":

sh src/run_success.sh

7. To find which simulations (cases) are complete (output/Filter.txt) and which are not (output/Missing.txt), run:

python src/create.filter.py

8.  If all are not complete, try submitting the job again (cases that were already complete, will be skipped):

sh src/run_elm.sh

9. Then to extract ELM outputs with:

sh src/extract_outputs.sh

