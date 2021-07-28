1. To create and activate a conda environment with all the python and R modules needed, run this in a shell:
```
./create.conda.env.sh
```
   Alternatively, you could also create the environment by running
```
conda env create -f environment.yml
```

2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv, run:
```
python src/generate.inputs.py 
```

3. To generate a base case of ELM, run:
```
python src/create.basecase.py
```

4. To generate ELM clone cases each associated with an ensemble member, run:
```
python src/create.ELM.ensembles.py 
```

5. To run the parallel emsemble simulations as per elm.py, run the following command:
```
sbatch src/run_elm.sh
```

6. To find which cases are complete (OutputExtract/Filter.txt) and which are not (OutputExtract/Missing.txt), run:
```
python src/create.filter.py
```

7. To extract outputs from ELM ensembles, run:
```
python src/extract.output.py
``` 
