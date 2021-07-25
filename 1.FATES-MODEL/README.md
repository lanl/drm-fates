1. To create and activate a conda environment with all the python and R modules needed
./create.conda.env.sh

2. To Generate parameter files based on paramter ensembles in HPU.Table.csv, run:
python generate.inputs.py 

3. To generate a base case of ELM via a shell script, run:
python create.basecase.py

3.To generate ELM clone cases via a shell script, run:
python generate_Sen_Cases.py 

3.To run the toy model, run the following command:
python run_elm.py

