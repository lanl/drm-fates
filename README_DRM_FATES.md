About
--------------------------------------------------------------------------------

This project will create and run point (e.g. each HPU) simulations of ELM-FATES in parallel and extract outputs.

Pre-requisites
--------------------------------------------------------------------------------
New Users:

1. To be able to access this repo, ask github-register@lanl.gov to add you to https://github.com/lanl
2. Ask Rutuja Chitra-Tarak (rutuja@lanl.gov) to add you to git@github.com:lanl/DRM team.
3. Ask Chonggang Xu (cxu@lanl.gov) to add you to the CESM project space on HPC. (After adding, takes a day for activation)
4. Ask Eunmo Koo (koo_e@lanl.gov) to add you to the HIGRAD project space on HPC.

Workflow
--------------------------------------------------------------------------------

1. First set variables in config.yaml (Pre-configured for a toy model run, so no changes required, but if you are running the same simulations again, to avoid overwriting existing simulations you will need be prompted for inputs. To avoid, change the case (simulation) name by changing the TAG name in config.yaml)

2. Set environmental variables. Also, re-run this if you get disconnected from HPC at any point in the workflow:

salloc -N 1 -t 00:20:00 --qos=interactive

source tools/.tcshrc

3. Run once to create and activate a conda environment (takes ~10-15 min) and load ELM (and FATES) branches. Also, re-run this if you get disconnected from HPC at any point in the workflow after #2:

sh tools/run_once.sh

4. Activate conda environment

! sh tools/activate_env.sh 
! Somehow the conda activate command throws an error (needs fixing), so run it outside:

conda activate elm_env

5. Then to run DRM:
sh tools/mod_run_drm.sh # modifies sbatch commands based on config.yaml
sbatch run_drm.sh
