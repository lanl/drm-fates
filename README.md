About
--------------------------------------------------------------------------------

This project will create and run point (e.g. each HPU) simulations of ELM-FATES in parallel and extract outputs.

Pre-requisites
--------------------------------------------------------------------------------
New Users:

1. To be able to access this repo, ask github-register@lanl.gov with your Github username to add you to https://github.com/lanl
2. First time users of HPC need to create a UNIX group in their name (e.g. LANL moniker) to be able to be added to other UNIX groups.
3. Ask Rutuja Chitra-Tarak (rutuja@lanl.gov) with your Z# to add you to git@github.com:lanl/DRM team.
4. Ask Chonggang Xu (cxu@lanl.gov) with your Z# to add you to the CESM project space on HPC. (After adding, takes a day for activation)
5. Ask Eunmo Koo (koo_e@lanl.gov) with your Z# to add you to the HIGRAD project space on HPC.
6. Ask Adam Atchley to add you to w22_fire HPC charging account.
Workflow
--------------------------------------------------------------------------------
1. To set-up this project on LANL HPC, do this and follow instructions in README.md:

ssh user@wtrw.lanl.gov

ssh ch-fe

! If you get disconnected from HPC, you wont loose the screen; and wont have load the environment if you use screen utility.
! While on screen, you could go up and down the screen with ctrl+A esc; and then esc to go back to writing mode.
! To detach and attach to a screen, start with:

screen -d -r

cd /usr/projects/higrad/$user

mkdir -p E3SM_cases

cd E3SM_cases

! Before cloning repo, you may need to ssh auntheticate, if you havent done that already: 
! Follow instructions at: https://docs.github.com/en/authentication/connecting-to-github-with-ssh

mkdir -p proj1

cd proj1

git clone git@github.com:lanl/drm-fates.git .

git checkout dev

!vim README.md

2. Copy LANL internal software at root

cp SOURCE/7.QUICFIRE-MODEL.zip /usr/projects/higrad/$user/proj1

unzip 7.QUICFIRE-MODEL.zip

sed -i '/7.QUICFIRE-MODEL\/projects\/ftFiles/c"'"$PWD/7.QUICFIRE-MODEL/projects/ftFiles/"'"' 7.QUICFIRE-MODEL/projects/Tester/QUIC_fire.inp

3. First set variables in config.yaml (Pre-configured for a toy model run). For example, you could change experiment name (TAG), WALL_TIME & N_NODE to reserve on HPC back-end, area of each FATES simulation or grid-cell size (FATES_RES), No. of FATES grid cells (SIM_END), FATES parameter files to use, turn on sensitivity analysis etc. 

4. Set environmental variables. Note: If you did not use screen utility, if you get disconnected from HPC at any point in the workflow, re-run #2:

source tools/.tcshrc

5. Run once to create and activate a conda environment (takes ~10-15 min) and load ELM (and FATES) branches. Note: If you did not use screen utility, if you get disconnected from HPC at any point in the workflow, re-run #2 & # 3:

sh tools/run_once.sh

6. Activate conda environment

! sh tools/activate_env.sh 
! Somehow the conda activate command throws an error (needs fixing), so run it outside:

conda activate elm_env

7. Rebuild QF

cd 5.TREES-QUICFIRE

make clean

make

cd ..

8. Then to run DRM:

! Next sbatch command calls DRM_framework_coupling.py on the back-end.
! Edit DRM_framework_coupling.py. Under ---main--- set VDM (FATES/LLM); VDM spin-up years (nyears), no of fires/loops (ncycle), and VDM run duration i.e. fire-interval (ncycyear)

!To view progress during short runs on the front-end, run (will need to continue Quickfire run): 

python DRM_framework_coupling.py

! If running on back-end, next command modifies sbatch commands based on config.yaml

sh tools/mod_run_drm.sh

! Next sbatch command calls DRM_framework_coupling.py on the back-end.
sbatch run_drm.sh 

! On an interactive mode, salloc -N 1 -t 00:20:00 --qos=interactive, then re-run source tools/.tcshrc and conda activate elm_env. Then sh run_drm.sh. (ssh does not work so, cant run run_once.sh)
