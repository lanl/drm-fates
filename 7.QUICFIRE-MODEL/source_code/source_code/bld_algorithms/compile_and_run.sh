#!/bin/sh
#!/bin/bash

#testcase="EglandCase"
testcase="LineFire_Mac"
compile_code=1

clear 

if [ $compile_code -eq 1 ]; then
	rm *.o 2> /dev/null
	rm *.mod 2> /dev/null
	if [ -d compiled_files ]; then
		echo "Clean compiled files"
		cd compiled_files
		rm quicfire 2> /dev/null
    	cd ..
	fi
fi

# Get parent directory
parent=$(dirname $PWD)

if [ $compile_code -eq 1 ]; then
#   if ifort -V; then
#       echo Using IFORT
#       FORTRANCOMPILER=ifort
#       COMPILEFLAGS=-c
#       LINKERFLAGS="-O3 -fopenmp -free -static-intel -qopenmp-link=static -parallel -fpe:0 -heap-arrays "
##       LINKERFLAGS="-O3 -fopenmp -free -static-intel -qopenmp-link=static -parallel -warn all -fpe:0 -heap-arrays"
#   else
	    echo Using GFORTRAN
	    FORTRANCOMPILER=gfortran
	    COMPILEFLAGS=-c
      # Debug flags
       LINKERFLAGS="-fopenmp -O3 -ffree-form -finit-real=snan -finit-integer=-999 -ffpe-trap=zero,overflow -Wall -fbounds-check -Wno-tabs"

		# Release flags
       #LINKERFLAGS="-fopenmp -O3 -ffree-form -finit-real=snan -finit-integer=-999 -ffpe-trap=zero"
#   fi

   $FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/file_handling_module_file.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/mersenne.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/atmos_mod.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/datamodule.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/datamodule_fire.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/utilities.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/FireCA_LDRD.f90 $LINKERFLAGS	
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/FireCA_common_subs.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/plantinit.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/bisect.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/poisson.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/buoyant_plume.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/datefunctions.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/read_hotmac_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/read_ittmm5_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/building_damage.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/read_quic_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/building_parameterizations.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/regress.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/interpolatewinds.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/sensorinit.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/diffusion.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/sor3d.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/divergence.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/sort.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/euler.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/surface_coords.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/init.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/interpolation_for_fire_grid.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/turbulence_model.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/utmll.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/outfile.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/wallbc.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/GenerateQUWinds.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/update_canopy.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/cleanup.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/disturbedflow.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS $parent/source_code/source_code/main.f90 $LINKERFLAGS

	$FORTRANCOMPILER $LINKERFLAGS -o compiled_files/quicfire *.o
else
	echo "------------------------------------------"
	echo "WARNING: Source code was not recompiled"	
	echo "------------------------------------------"
fi


if [ -f compiled_files/quicfire ];then
    #Copy the exe to the project folder
	cp compiled_files/quicfire "${parent}/projects/${testcase}/quicfire"
	cp drawfire.py "${parent}/projects/${testcase}/drawfire.py"
	cd "${parent}/projects//${testcase}"
	rm *.bin 2> /dev/null
	rm *.jpg 2> /dev/null
	rm *.log 2> /dev/null
	rm Plots 2> /dev/null
	./quicfire
	
	if [ -s "timelog.log" ];then
		
		echo "Drawing"
		if ! [ -z "$1" ]; then
			gen_vtk=$1
		else
			gen_vtk=0
		fi
		echo "Generate VTK (0=no, 1=yes): ", $gen_vtk		

		if ! [ -z "$2" ]; then
			#"/Users/sbrambilla/VirtualEnvPython/NISAC_Water/VEnv/bin/activate"
			echo "Using virtual env at: "
			echo $2
			source $2
		fi
		# python3 -m pip install --upgrade pip
		# python3 -m pip install numpy
		# python3 -m pip install matplotlib
		# python3 -m pip install scipy
		python3 drawfire.py $gen_vtk
		if ! [ -z "$2" ]; then
			deactivate
		fi
   else
		echo "***************************************************************"
   	echo "QUIC-FIRE did not complete the calculations. Execution aborted."
		echo "***************************************************************"
	fi
else
	echo "********************************"
	echo "Failed compilation"
	echo "********************************"
fi

