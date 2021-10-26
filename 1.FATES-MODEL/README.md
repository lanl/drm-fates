
This project will create and run a parallel ELM-FATES simulation and extract outputs.

1. First set variables in config.yaml (Well-configured for a toy model run, so no changes needed)

2. Set environmental variables:

source src/.tcshrc

3. cd into ELM_Disease, then run once to load ELM and FATES branches:

sh src/run_once.sh

4. Activate conda environment

! sh src/activate_env.sh 
! Somehow the conda activate command throws an error (needs fixing), so run it outside, provided that the environment name matches the one created (*_env):

conda activate elm_env

5. Then to set-up and run ELM simulations in parallel on the back node, run: 

sh src/run_elm_parallel.sh

6. You need to wait till simulations have run--depnding on years in config.yaml. Note: Currently sbatch run twice because first instance raises false alarm". You can confirm success if the second round of file of type slurmJOB#.out records "CASE.RUN HAS FINISHED"

7. Then extract ELM outputs with:

sh src/extract_outputs.sh

