All codes in this repository are subject to the following 
© 2025. Triad National Security, LLC. All rights reserved.

This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S. Department of Energy/National Nuclear Security Administration. All rights in the program are reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear Security Administration. The Government is granted for itself and others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare. derivative works, distribute copies to the public, perform publicly and display publicly, and to permit others to do so.





# FATES-DRM
### Step 1: Create the FATES case 
setup and run ./create_FATES_DRM.sh
You will need to set: 
user: your moniker
TAG: A name of your case
SITE_BASE_DIR: Where are you building your case
ex: "/lustre/scratch5/.mdt1/${user}/E3SM_cases/Eglin"   # Where is the site folder located? (SITE_NAME)
DATA_DIR: Where is this repository (where the internal data is stored)
ex: "/lustre/scratch5/.mdt1/zjrobbins/FATES-DRM" 

And the time of simulation at line 100-109


### Step 2: Run FATES_DRM.csh this will call the additional files
This will need you to set the case the following variables: 
   E3SM_CASE_CLONE: What is the name of your case. 
   CLONE_ROOT: Where are the FATES runs stored. 
   SRC: The location of these scripts
   RUN_ROOT: Where is the root where FATES is run
   STOP_N: How often to stop and run Quicfire

Results: Outputs will be stored in the following places: 
## Fates outputs
$RUN_ROOT/$E3SM_CASE_CLONE/run/$E3SM_CASE_CLONE.elm.h0*
$RUN_ROOT/$E3SM_CASE_CLONE/run/$E3SM_CASE_CLONE.elm.r*
## Fire Outputs
$RUN_ROOT/drm-fates/1.FATES-MODEL/VDM2FM/Pre_fire*
$RUN_ROOT/drm-fates/1.FATES-MODEL/FM2VDM/Afterfire*
$RUN_ROOT/drm-fates/1.FATES-MODEL/VDM2FM/fuels*
