About
--------------------------------------------------------------------------------

This project will create and run point (e.g. each HPU) simulations of ELM-FATES in parallel and extract outputs.

Workflow
--------------------------------------------------------------------------------

0. To set-up this project on LANL HPC, do this and follow instructions in README.md:
ssh user@wtrw.lanl.gov

ssh ba-fe

cd /usr/projects/cimmid/users/$user

mkdir -p ELM_cases

cd ELM_cases

! Before cloning repo, you may need to ssh auntheticate, if you haven't done that already:
! Follow instructions at: https://docs.github.com/en/authentication/connecting-to-github-with-ssh

mkdir -p proj1

cd proj1

git clone git@github.com:lanl/cimmid-e3sm.git .

git checkout cimmid

!vim README.md

1. First set variables in config.yaml (Pre-configured for a toy model run, so no changes required.)

2. Set environmental variables. Also, re-run this if you get disconnected from HPC at any point in the workflow:

source src/.tcshrc

3. Run once to create and activate a conda environment (takes ~10-15 min) and load ELM (and FATES) branches:

sh src/run_once.sh

4. Activate conda environment

! sh src/activate_env.sh 
! Somehow the conda activate command throws an error (needs fixing), so run it outside:

conda activate elm_env

5. Then to set-up and run ELM simulations in parallel on the back node, run: 

sh src/run_elm_parallel.sh

6. Run this script until you get "sbatch has finished running" & "a simulation ran successfully":

sh src/run_success.sh

7. To find which simulations (cases) are complete (output/Filter.txt) and which are not (output/Missing.txt), run:

python src/create.filter.py

8.  If all are not complete, try submitting the job again (cases that were already complete, will be skipped):

sbatch src/run_elm.sh

9. Then to extract ELM outputs, run:

sh src/extract_outputs.sh

Acknowledgements
--------------------------------------------------------------------------------

This workflow was developed by Rutuja Chitra-Tarak (rutuja@lanl.gov). Some python scripts, developed by Chonggang Xu (cxu@lanl.gov), were repurposed for use in this workflow. Those two should be acknowledged in publications. 

The Energy Exascale Earth System Model (E3SM) Project should be acknowledged in publications as the origin of the model using 
[these guidelines](https://e3sm.org/resources/policies/acknowledge-e3sm/).

In addition, the software should be cited. For your convenience, the following BibTeX entry is provided.

@misc{e3sm-model,
	title = {{Energy Exascale Earth System Model (E3SM)}},
	author = {{E3SM Project}},
	abstractNote = {{E3SM} is a state-of-the-art fully coupled model of the {E}arth's 
		climate including important biogeochemical and cryospheric processes.},
	howpublished = {[Computer Software] \url{https://dx.doi.org/10.11578/E3SM/dc.20180418.36}},
	url = {https://dx.doi.org/10.11578/E3SM/dc.20180418.36},
	doi = {10.11578/E3SM/dc.20180418.36},
	year = 2018,
	month = apr,
}

License
--------------------------------------------------------------------------------

The E3SM model became open development at the time of first model and data release. Please see [LICENSE](https://github.com/rutujact/E3SM/blob/master/LICENSE) for details.
