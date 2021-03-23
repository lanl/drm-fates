#!/bin/sh
#!/bin/bash

clear 

# Example
# compile=1 run=1 fireca='ldrd' compiler='intel' compiler_mode='release' testcase='Linefire/Linefire' exename="quicfire.exe" python_env="/Users/sbrambilla/VirtualEnvPython/QUIC_Fire/VEnv/bin/activate" ./adv_compile_and_run.sh
# compile=1 run=0 fireca='fst' compiler='intel' compiler_mode='release' exename="quicfst_MACI.exe" ./adv_compile_and_run.sh

if [ -z "$fireca" ]; then
   fireca="fst"
else
   if [ "$fireca" != "ldrd" ] && [ "$fireca" != "fst" ]; then
      echo "********************************"
      echo "Invalid fireca file [ldrd/fst]: "$fireca
      echo "Execution will be interrupted"
      echo "********************************"
      exit   
   fi
fi
echo "FireCA file: "$fireca
if [ "$fireca" == "ldrd" ]; then
   filename="FireCA_LDRD.f90"
else
   filename="FireCA_FST.f90"
fi

if [ -z "$compiler" ]; then
   compiler="gfortran"
else
   if [ "$compiler" == "intel" ] || [ "$compiler" == "gfortran" ]; then
      echo "Chosen compiler: "$compiler
   else
      echo "********************************"
      echo "Invalid compiler [gfortran/intel]: "$compiler
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   fi
fi

if [ -z "$compiler_mode" ]; then
   compiler_mode="debug"
else
   if [ "$compiler_mode" == "debug" ] || [ "$compiler_mode" == "release" ]; then
      echo "Chosen compiler mode: "$compiler_mode
   else
      echo "********************************"   
      echo "Invalid compiler mode [debug/release]: "$compiler_mode
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   fi
fi

if [ -z "$compile" ]; then
   compile=1
else
   if [ "$compile" -lt 0 ] || [ "$compile" -gt 1 ]; then
      echo "********************************"
      echo "Invalid compile option [0-1]: "$compile
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   fi
fi

if [ -z "$exename" ]; then
   exename="quicfire.exe"
fi

if [ -z "$run" ]; then
   run=1
else
   if [ "$run" -lt 0 ] || [ "$run" -gt 1 ]; then
      echo "********************************"
      echo "Invalid run option [0-1]: "$run
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   fi
fi

if [ "$run" -eq 1 ]; then
   # Only check test case if run is 1
   if [ -z "$testcase" ]; then
      echo "********************************"
      echo "Please provide case"
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   else
      echo "Test case to run: "$testcase
   fi
fi

if [ "$compile" -eq 0 ]; then
   echo "Script will NOT compile the code"
elif [ "$compile" -eq 1 ]; then
   echo "Script will compile code"
fi
if [ "$run" -eq 0 ]; then
   echo "Script will NOT run the code"
elif [ "$run" -eq 1 ]; then
   echo "Script will run code"
fi

if ! [ -z "gen_vtk" ]; then
   gen_vtk=0
else
   if [ "$gen_vtk" -lt 0 ] || [ "$gen_vtk" -gt 1 ]; then
      echo "********************************"
      echo "Invalid vtk option [0-1]: "$gen_vtk
      echo "Execution will be interrupted"
      echo "********************************"
      exit
   fi
fi
if [ "$gen_vtk" -eq 0 ]; then
   echo "VTK files will NOT be generated"
else
   echo "VTK files will be generated"
fi

# Get parent directory
parent=$(dirname $PWD)

# Remove old code
if [ $compile -eq 1 ]; then
   echo "Remove old compiled code"
	rm *.o 2> /dev/null
	rm *.mod 2> /dev/null
   # If the folder compiled_files exist
	if [ -d compiled_files ]; then		
		cd compiled_files
		rm *.exe 2> /dev/null
    	cd ..
	fi
fi

# Compile
COMPILEFLAGS=-c
if [ $compile -eq 1 ]; then
  if [ "$compiler" == "intel" ]; then      
      FORTRANCOMPILER=ifort      
      if [ "$compiler_mode" == "debug" ]; then
         LINKERFLAGS="-O0 -fopenmp -free -static-intel -qopenmp-link=static -fpe:0 -heap-arrays -warn all"
      else
         LINKERFLAGS="-O3 -fopenmp -free -static-intel -qopenmp-link=static -parallel -heap-arrays "
      fi
  else      
      FORTRANCOMPILER=gfortran      
      if [ "$compiler_mode" == "debug" ]; then
         # Debug flags
         LINKERFLAGS="-fopenmp -O3 -ffree-form -finit-real=snan -finit-integer=-999 -ffpe-trap=zero,overflow -Wall -fbounds-check -Wno-tabs"
      else
		   # Release flags
         LINKERFLAGS="-fopenmp -O3 -ffree-form -finit-real=snan -finit-integer=-999 -ffpe-trap=zero"
      fi
  fi

   $FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/file_handling_module_file.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/mersenne.f90 $LINKERFLAGS	
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/datamodule.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/datamodule_fire.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/utilities.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/${filename} $LINKERFLAGS	
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/FireCA_common_subs.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/plantinit.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/bisect.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/poisson.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/buoyant_plume.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/datefunctions.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/read_hotmac_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/read_ittmm5_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/building_damage.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/read_quic_met.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/building_parameterizations.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/regress.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/interpolatewinds.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/sensorinit.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/diffusion.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/sor3d.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/divergence.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/sort.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/euler.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/surface_coords.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/init.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/interpolation_for_fire_grid.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/turbulence_model.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/utmll.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/outfile.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/wallbc.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/GenerateQUWinds.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/update_canopy.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/cleanup.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/disturbedflow.f90 $LINKERFLAGS
	$FORTRANCOMPILER $COMPILEFLAGS ${parent}/source_code/source_code/main.f90 $LINKERFLAGS

	$FORTRANCOMPILER $LINKERFLAGS -o compiled_files/${exename} *.o

   echo "******************************************"
	echo "Code compiled successfully"	
	echo "******************************************"

else
	echo "******************************************"
	echo "WARNING: Source code was not recompiled"	
	echo "******************************************"
fi

# Run case
if [ -f compiled_files/${exename} ]; then
   if [ $run -eq 1 ]; then
      #Copy the exe to the project folder
      cp compiled_files/${exename} "${parent}/projects/${testcase}/${exename}"
      cp drawfire.py "${parent}/projects/${testcase}/drawfire.py"
      cd "${parent}/projects//${testcase}"
      # Remove existing results
      rm *.bin 2> /dev/null
      rm *.jpg 2> /dev/null
      rm *.log 2> /dev/null
      rm Plots 2> /dev/null

      # Run code
      ./${exename}
      
      # Check if code has run
      if [ -s "timelog.log" ];then
         
         echo "Drawing"		
         if ! [ -z "$python_env" ]; then
            #"/Users/sbrambilla/VirtualEnvPython/QUIC_Fire/VEnv/bin/activate"
            echo "Using virtual env at: "
            echo $python_env
            source $python_env
         fi
         # python3 -m pip install --upgrade pip
         # python3 -m pip install numpy
         # python3 -m pip install matplotlib
         # python3 -m pip install scipy

         python3 drawfire.py $gen_vtk
         if ! [ -z "$python_env" ]; then
            deactivate
         fi
      else
         echo "***************************************************************"
         echo "QUIC-FIRE did not complete the calculations. Execution aborted."
         echo "***************************************************************"
         exit
      fi
   fi
else
	echo "********************************"
	echo "Executable not found"
	echo "********************************"
   exit
fi
