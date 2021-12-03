#!/bin/sh

# 7. Make sure cases have successfully finished, else re-run sbatch:
# Typically no harm if sbatch is submitted again. It will skip those cases that are completed:

if grep -w "Finished on"  slurm*.out; then
   echo "sbatch has finished running"
   # Then test if at least one case has finished successfully.
   if grep -w "CASE.RUN HAS FINISHED" slurm*.out; then
        echo "At least one simulation ran successfully."
   else
       	echo "None of the simulation ran successfully. Rerunning batch";
	sbatch src/run_elm.sh
   fi
else
echo "sbatch has not yet finished running"
fi
