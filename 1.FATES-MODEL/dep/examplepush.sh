#!/bin/sh
export testsample= "lustre/scratch5/.mdt1/zjrobbins/"
nano testsample/testsample/ 
.xmlchange --file env_run.xml --id CONTINUE_RUN --val FALSE

for case_i in $(seq 1 1000)
	echo "this is run $case_i"
	cd  testsample/runnumber$case_i/
	./runthis 
