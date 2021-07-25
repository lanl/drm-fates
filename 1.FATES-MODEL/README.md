1. To create and activate a conda environment with all the python and R modules needed, run this in a shell
```
./create.conda.env.sh
```

2. To generate multiple parameter files based on paramter ensembles in HPU.Table.csv via an R script by the same name, run:
```
python generate.inputs.py 
```

3. To generate a base case of ELM via a shell script by the same name, run:
```
python create.basecase.py
```

4. To generate ELM clone cases via a shell script by the same name, run:
```
python generate_Sen_Cases.py 
```

5. To run the toy model, run the following command:
```
python run_elm.py
```

6. To find which cases are complete (OutputExtract/Filter.txt) and which are not (OutputExtract/Missing.txt), run:
```
python create.filter.py
```

