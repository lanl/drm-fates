## 1. Make sure you are on the desired ELM branch. 
 
cd $ACME_ROOT/
git checkout rutuja/ELM_FATES_HYDRO_DFLT_pedotrf_HKSAT_ADJ
# <!-- git clone git@github.com:rutujact/E3SM.git --branch rutuja/ELM_FATES_HYDRO_default_pedotransfer -->
## If a submodule is not found, update them

# git submodule update --init --recursive

## 2. Make sure you are on the desired fates branch. 

cd $FATES_ROOT
git branch 
git checkout xuchongang/coastal_veg_ready_to_merge_July_2019_hybrid_static

##  - If not, do:
  
# cd /turquoise/usr/projects/veg/rutuja/ACME/components/clm/src/external_models

# git clone git@github.com:rutujact/fates.git --branch xuchongang/coastal_veg_ready_to_merge_July_2019_hybrid_static --single-branch

## 3. Add conda path.

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
if ( -f "/usr/projects/hpcsoft/common/x86_64/anaconda/2020.07-python-3.8/etc/profile.d/conda.csh" ) then
    source "/usr/projects/hpcsoft/common/x86_64/anaconda/2020.07-python-3.8/etc/profile.d/conda.csh"
else
    setenv PATH "/usr/projects/hpcsoft/common/x86_64/anaconda/2020.07-python-3.8/bin:$PATH"
endif
# <<< conda initialize <<<

# 4. Activate conda env
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
cd ../$SCRIPT_DIR
conda activate conda_env
