# 3. Activate conda env
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
PROJECT_ROOT=`realpath "$SCRIPT_DIR/.."`
cd $PROJECT_ROOT
conda activate elm_env
