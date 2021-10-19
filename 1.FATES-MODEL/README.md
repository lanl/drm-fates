
This project will create and run a parallel ELM-FATES simulation and extract outputs.

1. First set variables in config.yaml (Well-configured for a toy model run, so no changes needed)
2. Run once to load ELM and FATES branches:
sh src/run_once.sh
3. Set environmental variables:
source src/.tcshrc
4. Activate conda environment
sh src/activate_env.sh 
4. Then to run the whole project, run the following in the shell:
sh src/run_project.sh
