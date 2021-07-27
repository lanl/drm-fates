1. To create and activate a conda environment with all the python and R modules needed, run this in a shell:
```
./create.conda.env.sh
```

2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv, run:
```
python generate.inputs.py 
```

3. To generate a base case of ELM, run:
```
python create.basecase.py
```

4. To generate ELM clone cases each associated with an ensemble member, run:
```
python create.ELM.ensembles.py 
```

5. To run the parallel emsemble simulations as per elm.py, run the following command:
```
sbatch run_elm.sh
```

6. To find which cases are complete (OutputExtract/Filter.txt) and which are not (OutputExtract/Missing.txt), run:
```
python create.filter.py
```

