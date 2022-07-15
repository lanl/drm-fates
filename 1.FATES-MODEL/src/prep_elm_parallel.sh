#!/bin/sh
# =======================================================================================
#
# This script will create, setup and build and run a multi-site ELM-FATES simulation  
# based on config.yaml.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================
# 1. READ YAML VARIABLES
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )" #Locate the directory of this script no matter where it is called from
source `realpath "$SCRIPT_DIR/../../tools/yaml.sh"`
create_variables "$SCRIPT_DIR/../../config.yaml"

# 2.  Unzip sample climate data:

unzip -u data/ConvertMetCSVtoCLM/NCOut/eglin_0.1x0.1_met.v2pio.zip -d data/eglin_0.1x0.1_v4.0i/

# 3. If MUTATE is true, then if sensitivity analysis is true, build or use exisitng parameter table:

if [ ${MUTATE} == TRUE ] || [ ${MUTATE} == true ]; then
        if [ ${SENSITIVITY} == TRUE ] || [ ${SENSITIVITY} == true ]; then
                if [ ${BUILD_PARAM_TABLE} == TRUE ] || [ ${BUILD_PARAM_TABLE} == true ]; then
                        echo "Building Parameter Sensitivity Table"
                        python src/generate.param.table.py
                else
                        echo "Parameter Sensitivity Table given, Not building"
                fi
        fi

        # 4. To generate multiple parameter files based on paramter ensembles in Param.Table.csv, run:

        python src/generate.inputs.py
fi

# 5. To generate a base case of ELM, run:

python src/create.basecase.py
 
# 6. To generate ELM clone cases each associated with an ensemble member, run:

python src/create.ELM.ensemble.py
 
