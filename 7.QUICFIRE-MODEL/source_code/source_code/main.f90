	!/Qprofile-loops:all
	!/Qprofile-functions
	!                                 Notice
	!  This program was prepared by the University of California (University)
	!  under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
	!  All rights in the program are reserved by DOE on behalf of the Government
	!  and the University pursuant to the contract. You are authorized to use
	!  this program for Government purposes but it is not to be released or
	!  distributed to the public.
	!  NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
	!  NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
	!  MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
	!  OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
	!  ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
	!  THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
	!
	program main
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! QWIC-URB is a multi-building flow solver using 3 dimensional
	! successive overrelaxation method solver and is based on the
	! work of Rockle (1990) and Kaplan and Dinar (1996).
	!
	! p1 and p2 are cell centered quantites   and called Lagrange
	! Multipliers.
	! u,v,w and uo,vo,wo are cell face quantites
	! uo,vo, and wo are the inital non-mass conserved velocities
	! u,v,w is the final mass conserved velocity field
	!
	! bcsetup.f does much of the work in the code, here the ground
	! and buildings are defined, the empircal flow parameterizations
	! input and the boundary conditions for the SOR solver set up.
	!
	!  - the final velocity field is written out in euler.f
	!  - the Lagrangde multiplier field is writtein out in sor3d.f
	!
	! UPDATES
	! Eric Pardyjak
	! QWIC-URBv1.00Beta September 2002
	!  September 2002 Update:
	!  This version has an updated input file. The height of zref
	!  for the power inlet velocity profile has been added.
	! QWIC-URBv1.00 October 2002
	! This version has multiple updates that include
	!  1. new output format of output files & new output files
	!  2. input.dat has a new line for a rooftop recirculation
	!  3. new coordinate system
	!  4. fixes in array out of bounds in sor3d.f
	!
	! QWIC-URBv1.00f90 January 2003
	!  1. Winds out of the North and South can be handle (ERP 12/17)
	!  2. A bug in the street canyon routine was fixed (MDW 1/8)
	!  3. Allocatable arrays are now deallocated allowing the qwicurb
	!        to be run multiple time in the GUI. (TWH 1/8)
	!
	! Note if this version of the code is being used with a Matlab GUI
	! and Visual Compaq Fortran Compiler 6.6, it will be necessary to do the
	! following:
	!  1. Download the 6.6B update (VFRUN66BI.exe). It can be found at
	!      the following location:
	!        http://h18009.www1.hp.com/fortran/visual/redist.html
	!  2. Once patch has been installed move the file dformd.dll from
	!     the C:\Windows\system32 on XP machines or
	!        C:\WINNT\system32 on Windows2000 machines folder to
	!        C:\MATLAB6p5\bin\win32
	!  3. Be sure to run mex -setup so that the new df66opts.bat is
	!      updated.
	!
	! erp 1/29/04
	! This verson of QUICURB has Petra Kastner-Klein's street canyon
	! models available. See subroutine init.f90 and bcsetup.f90
	!
	! ERP 10/04/04
	! This version of QUICURB allows for multiple runs to simulate
	! "quasi-steady" time dependent flows
	! ERP 3/8/05
	! - This version of the code contains the new building sorting logic to sort
	! building from small to tall so that the algorithms are applied in that order
	! - THis version also uses a local velocity to determine the orientation of a
	! street canyon vortex
	! ERP 6/17/05
	! This version of QUICURB has the basic data assimilation algorithms based
	! on Barnes objective mapping as implemented by Tom Booth
	! new subroutines include: sensorinit
	!
	!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	use rnd_num_vars
	use fireca_module
	use omp_lib
	use time_module
	use constants
	use file_handling_module
	use sor_module
	use canopy_module
	use wind_profile_module
	use flags_module
	use bld_module
	use diffusion_module
	use landuse_module
	use plume_module
	use ignitions_module
	use grid_module
	use damage_module
	use winds_module
	use time_module

	implicit none
	
	integer ::							&
		qu_it,							& ! N/A, time step iterator for quic winds
		fca_it,							& ! N/A, time step iterator for fire
		i,									&
		qu_updates,						&
		qu_updates_avewinds,			&
		fireca_updates,				&
		fireca_updates_emiss,		&
		cpu_start,cpu_end,count_rate

	integer :: 									&
		max_num_cpu,							& ! N/A, maximum number of cpu used
		cpu_start_subs,cpu_end_subs		  ! N/A, cpu time counters

	integer, dimension(4) :: tot_cpu_time_subs

	!--------------------------------------------------------------------------
	call system_clock(cpu_start,count_rate)
	cpu_start = cpu_start / count_rate

	tot_cpu_time_subs = 0

	 IF(IARGC() .gt. 0)THEN
      CALL GETARG(1, workingDirectory)
      DO i=1,256
         IF(workingDirectory(i:i) .eq. '/')THEN
            fileSeparator='/'
            EXIT
         ELSEIF(workingDirectory(i:i) .eq. '\')THEN
            fileSeparator='\'
            EXIT
         ENDIF
      ENDDO
      IF(i .gt. 256)THEN
         workingDirectory=' '
         fileSeparator=' '
      ENDIF
   ELSE
      workingDirectory=' '
      fileSeparator=' '
   ENDIF
	flag%errorWrite = 0

	call read_cpu(max_num_cpu)
	numThreads = 1
	!$OMP PARALLEL
	numThreads = min(OMP_GET_NUM_THREADS()	, max_num_cpu)
	!$OMP END PARALLEL
	print *, "Number of threads=",numThreads

	!$ CALL OMP_SET_NUM_THREADS(numThreads)
	!$ print *, "OpenMP Enabled.  Allocating", numThreads, "threads"


	open(ID_FILE_TIMELOG,file=TRIM(workingDirectory)//TRIM(fileSeparator)//'timelog.log',status='replace')
	
	! Initialize world
	call init()

	! Import fire parameters
	call readQuicFireFile()
	
	! Initialize time steps
	call init_times(qu_updates, qu_updates_avewinds, fireca_updates,	&
		fireca_updates_emiss, fca_it)
	
	! sort buildings by height small to tall
	if(bld%number .eq. 0 &
			.and. windProfile%num_sites .eq. 1 &
			.and. canopy%number .eq. 0 &
			.and. flag%isfire .eq. 0 &
			.and. (landuse%flag .eq. 0 .or. &
			(landuse%veg_flag .eq. 0 .and. landuse%urb_flag .eq. 0)))then
		diff%flag = 0
		sor%itermax = 0
	endif
	if(bld%number > 0) then
		call sort
		if(damage%flag .eq. 1) call building_damage
	endif
	
	! QUIC time steps
	do qu_it = 1,qutime%nsteps
		qutime%current_iter = qu_it
		
		! Reset SOR
		sor%p1 = 0.
		
		call output_time()
		
		! Read met input file for each time step
		call sensorinit()
		
		! Call boundary condition set up routine
		if(flag%isfire == 1)then
			call building_parameterizations_fire()
		else
			call building_parameterizations_qu()
		endif

		if(flag%isfire == 0) then
			! Generate the wind field for this QU time
			call GenerateQUWinds()
		else
			call ResetSORMCVars(1)
			call GenerateQUWinds()			
			
			! Done each time that there is a new QU wind
			if (bld%number > 0 .and. flag%bld_to_fuel == 1) call ConvertBldgToCanopy()
			if(flag%fuel_to_canopy == 1) call ConvertFuelToCanopy()
			
			if(qu_it == 1) then
				call PrintQUOutputs(0)
				if(flag%quwinds_ave_file == 1) then
					quwinds%u_ave = quwinds%u
					quwinds%v_ave = quwinds%v
					quwinds%w_ave = quwinds%w					
					call PrintQUOutputsAveraged(0,1.)
				endif
			endif

			if(fire_ignition%unix_time < qutime%local(qutime%current_iter)) then
				call save_mc_winds()

				!Intrepolate wind from QUIC grid (Face centered) to fire grid (body centered)
				if(sum(fire%energy_to_atmos) > 0) then
					call buoyant_plume(fire%energy_to_atmos)
					quwinds%uo = quwinds%u_mc
					quwinds%vo = quwinds%v_mc
					quwinds%wo = quwinds%w_mc + quwinds%wplume
					call GenerateQUWinds()
					fire%energy_to_atmos = 0.
				endif

				call interpolate_for_fire_grid()

				! keep running the fire until either the time is > of the fire sim time or a new QU field is available
				do while(fcatime%current_int < qutime%local(qutime%current_iter) .and.  &
					fcatime%current_int < fcatime%len_simulation)

					do while(fca_it < qutime%firesteps_to_update .and.								&
						fcatime%current_int < qutime%local(qutime%current_iter) .and.		&
						fcatime%current_int < fcatime%len_simulation) ! do fire steps until it's time to update QU

						call increment_time(fca_it, fireca_updates, fireca_updates_emiss)

						! Ignore the first time step
						if(fcatime%current_int > fcatime%dt_int .and.					&
							(fire_ignition%flag == 5 .or.  fire_ignition%flag == 7)) then
							call update_ignitions_pattern()
						endif

						call system_clock(cpu_start_subs,count_rate)
						cpu_start_subs = cpu_start_subs / count_rate
						call update_fire(fb%flag, qutime%firesteps_to_update*fcatime%dt_int)
						call system_clock(cpu_end_subs,count_rate)
						cpu_end_subs = cpu_end_subs / count_rate
						tot_cpu_time_subs(1) = tot_cpu_time_subs(1) + cpu_end_subs - cpu_start_subs

						! Output Fireca variables
						if(fireca_updates == flag%num_fireca_updates) then
							fireca_updates = 0
							
							!if(fcatime%current_int >= 200) 
							call PrintFireCAOutputs(fcatime%current_int)
						endif		

						! Output emissions
						if(fireca_updates_emiss == flag%num_emission_and_rad_updates .and. &
							(flag%emissions_file == 1 .or. flag%emissions_file == 2 .or. &
							flag%emissions_file == 3 .or. flag%thermalrad_file == 1)) then
							call PrintFireCAEmissionsAndRadiationOutputs(fcatime%current_int)
							fireca_updates_emiss = 0
						endif
					enddo

					! They are not equal in the case QU needs to update its wind
					if(fca_it == qutime%firesteps_to_update) then						
						fca_it = 0
						qu_updates = qu_updates + 1
						qu_updates_avewinds = qu_updates_avewinds + 1

						! Generate plumes & redo mass consistency
						quwinds%wplume = 0.
						fire%energy_to_atmos = fire%energy_to_atmos / real(qutime%firesteps_to_update)
						
						!fire%energy_to_atmos = 0
						!fire%energy_to_atmos(25, 50, 1) = 400
						!fire%energy_to_atmos(25, 60, 1) = 400
						
						call system_clock(cpu_start_subs,count_rate)
						cpu_start_subs = cpu_start_subs / count_rate
						call buoyant_plume(fire%energy_to_atmos)
						!if(mod(fcatime%current_int, flag%num_fireca_updates) == 0) call output_en2atm_w(fcatime%current_int)
						call system_clock(cpu_end_subs,count_rate)
						cpu_end_subs = cpu_end_subs / count_rate
						tot_cpu_time_subs(2) = tot_cpu_time_subs(2) + cpu_end_subs - cpu_start_subs
						
						quwinds%uo = quwinds%u_mc
						quwinds%vo = quwinds%v_mc
						quwinds%wo = quwinds%w_mc + quwinds%wplume
						
						! Update winds with new fuel density
						if(canopy%UPDATE_WINDS == 1)then
							call update_canopy()
						endif
						
						call ResetSORMCVars(0)
						! Impose mass consistency
						call system_clock(cpu_start_subs,count_rate)
						cpu_start_subs = cpu_start_subs / count_rate
						! Make it an average over the time steps
						call GenerateQUWinds()
						call system_clock(cpu_end_subs,count_rate)
						cpu_end_subs = cpu_end_subs / count_rate
						tot_cpu_time_subs(3) = tot_cpu_time_subs(3) + cpu_end_subs - cpu_start_subs

						! Calculate QU average winds
						if(flag%quwinds_ave_file == 1)then
							quwinds%u_ave = quwinds%u_ave + quwinds%u
							quwinds%v_ave = quwinds%v_ave + quwinds%v
							quwinds%w_ave = quwinds%w_ave + quwinds%w
						endif

						!Intrepolate wind from QUIC grid (Face centered) to fire grid (body centered)
						call system_clock(cpu_start_subs,count_rate)
						cpu_start_subs = cpu_start_subs / count_rate						
						call interpolate_for_fire_grid()
						call system_clock(cpu_end_subs,count_rate)
						cpu_end_subs = cpu_end_subs / count_rate
						tot_cpu_time_subs(4) = tot_cpu_time_subs(4) + cpu_end_subs - cpu_start_subs

					endif

					! Output QU ave winds
					if(qu_updates_avewinds == flag%num_qu_updates_avewinds .and. flag%quwinds_ave_file == 1) then
						call PrintQUOutputsAveraged(fcatime%current_int, real(flag%num_qu_updates))
						qu_updates_avewinds = 0						
					endif
					
					! Output QU instantaneous winds
					if(qu_updates == flag%num_qu_updates) then
						!if(fcatime%current_int >= 200) 
						call PrintQUOutputs(fcatime%current_int)
						qu_updates = 0
					endif

					if(fca_it == 0) then
						! Winds just updated
						fire%energy_to_atmos = 0.						
					endif
				enddo
			endif
		endif

		call disturbedflow()
		call outfile()

	enddo

	close(DBGFILE)

	! deallocate allocatable arrays
	call CLEANUP()

	if(flag%isfire == 1) deallocate(qutime%local)

	call system_clock(cpu_end,count_rate)
	cpu_end = cpu_end / count_rate

	write(ID_FILE_TIMELOG,*)'Simulation time = ', cpu_end-cpu_start, 's'
	do i = 1,4
		write(ID_FILE_TIMELOG,*)'Sub time #', i,': ', tot_cpu_time_subs(i)
	enddo
	!write(ID_FILE_TIMELOG,*)'!=========================================='
	!write(ID_FILE_TIMELOG,*)'Init plumes: ', time_bp(1)
	!write(ID_FILE_TIMELOG,*)'Move: ', time_bp(2)
	!write(ID_FILE_TIMELOG,*)'Update prop: ', time_bp(3)
	!write(ID_FILE_TIMELOG,*)'Merge: ', time_bp(4)
	!write(ID_FILE_TIMELOG,*)'WWCalc: ', time_bp(5)
	!write(ID_FILE_TIMELOG,*)'PlumeFinalization: ', time_bp(6)
	!write(ID_FILE_TIMELOG,*)'SOR MC: ', time_bp(7)
	!write(ID_FILE_TIMELOG,*)'Alpha2: ', time_bp(8)
	write(ID_FILE_TIMELOG,*)'!=========================================='
	write(ID_FILE_TIMELOG,*)'Maximum plume height [m]: ',plume_diagnostics%max_height
	write(ID_FILE_TIMELOG,*)'Maximum plume w [m/s]: ', plume_diagnostics%max_w
	write(ID_FILE_TIMELOG,*)'Maximum domain height [m]: ', qugrid%Lz

	close(ID_FILE_TIMELOG)
	write(*,*)'Simulation time = ',cpu_end-cpu_start,' s'
	
	end
 !=================================================================================================
 !=================================================================================================
	subroutine initializeSeed(fix_random_seed,seed_val)

	!use constants
	use rnd_num_vars

	implicit none
	integer,intent(IN) :: fix_random_seed,seed_val
	integer,dimension(8) :: values
	integer :: TID

	if (fix_random_seed == 1) then
		iseed = seed_val
	else
		call date_and_time(values=values)
		iseed = values(8)
	endif
	allocate(stat(numThreads))
	do TID = 1,numThreads
		call random_setseed(iseed + TID,stat(TID))
	enddo

	end subroutine
 !=================================================================================================
 !=================================================================================================
	subroutine PrintQUOutputs(fire_time)

	use flags_module

	implicit none

	integer,intent(IN) :: fire_time

	if(flag%quwinds_inst_file == 1)		then
		print*,'Writing QU output files'
		call output_wind_qu_inst(fire_time)
	endif

	end
!===================================================================================
!===================================================================================
	subroutine PrintQUOutputsAveraged(simtime, num_updates)

   use constants
   use grid_module
   use file_handling_module
   use winds_module
	use flags_module

	implicit none

	integer,intent(IN) :: simtime
	real, intent(IN) :: num_updates
	character(len=35) :: filename
	integer i,j,k

	print*,'Writing QU output files (averaged)'
	quwinds%u_ave = quwinds%u_ave / num_updates
	quwinds%v_ave = quwinds%v_ave / num_updates
	quwinds%w_ave = quwinds%w_ave / num_updates

	
	write(filename, "(A12,I0.5,A4)") "qu_windu_ave", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_U, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_U) (((0.5*(quwinds%u_ave(i,j,k)+quwinds%u_ave(i+1,j,k)),  &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
   close(ID_FILE_WIND_QU_U)

	write(filename, "(A12,I0.5,A4)") "qu_windv_ave", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_V, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_V) (((0.5*(quwinds%v_ave(i,j,k)+quwinds%v_ave(i,j+1,k)),  &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
   close(ID_FILE_WIND_QU_V)

	write(filename, "(A12,I0.5,A4)") "qu_windw_ave", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_W, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_W) (((0.5*(quwinds%w_ave(i,j,k)+quwinds%w_ave(i,j,k+1)),  &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
	close(ID_FILE_WIND_QU_W)
	
	quwinds%u_ave = 0
	quwinds%v_ave = 0
	quwinds%w_ave = 0
						
	return

1101 write(msgoutfile,*)'Error while opening file '//TRIM(filename)//'.'
	 call TerminateProgram()
	

	end
!===================================================================================
!===================================================================================
	subroutine PrintFireCAOutputs(fire_time)

	use flags_module

	implicit none

	integer,intent(IN) :: fire_time

	print*,'Writing FireCA output files'
	if(flag%reactrate_file == 1)							call output_react_rate(fire_time)
	if(flag%fuelmass_file == 1 .or. fire_time == 0)	call output_fuels(fire_time)
	if(flag%qfwinds_file == 1)								call output_wind(fire_time)
	if(flag%moisture_file == 1)							call output_moisture(fire_time)
	if(flag%int_massburnt_file == 1)						call output_massburnt(fire_time)
	if(flag%en2atm_file == 1)								call output_en2atm(fire_time)

	end
!===================================================================================
!===================================================================================
	SUBROUTINE PrintFireCAEmissionsAndRadiationOutputs(fire_time)

	use flags_module
	use fireca_module

	implicit none

	integer,intent(IN) :: fire_time

	
	print*,'Writing emission and radiation output files'

	if(flag%emissions_file > 0)then
		call output_emissions(fire_time)
		
		! Reset variables
		fuels%mu_soot = 0
		fuels%sigma_soot = 0
	endif

	if(flag%thermalrad_file == 1)then
		call output_thermal_radiation(fire_time)
		fuels%conv_human = 0.
	endif

	end
!===================================================================================
!===================================================================================
	SUBROUTINE output_time()
	
	use time_module
	use grid_module
	
	implicit none
	
	integer*8 :: idate(6)
	
	call unix2c(qutime%unix(qutime%current_iter)+qugrid%utc_offset*3600,idate)
	write(*,14001)idate(2),idate(3),idate(1),idate(4),idate(5),idate(6)
14001	format(' Calculating winds for ',i2.2,'/',i2.2,'/',i4,' ',i2.2,':',i2.2,':',i2.2)
	
	end
!===================================================================================
!===================================================================================
	SUBROUTINE save_mc_winds()
	
	use winds_module
	
	implicit none

	quwinds%u_mc = quwinds%u
	quwinds%v_mc = quwinds%v
	quwinds%w_mc = quwinds%w
			
	end
!===================================================================================
!===================================================================================
	SUBROUTINE increment_time(fca_it, fireca_updates, fireca_updates_emiss)

	use time_module
	
	implicit none
	
	integer :: fca_it, fireca_updates, fireca_updates_emiss
	
	fca_it = fca_it + 1
	fcatime%current_int = fcatime%current_int + fcatime%dt_int
	fcatime%current_real = real(fcatime%current_int)
	fireca_updates = fireca_updates + 1
	fireca_updates_emiss = fireca_updates_emiss + 1
	
	write(*,10)fca_it,fcatime%current_int
10			format(' Fire time step #: ',i3,'; Fire time [s]: ',i6)
	
	end
!===================================================================================
!===================================================================================
	SUBROUTINE PrintFuelHeight()

	use file_handling_module
	use fireca_module
	
	implicit none
	
	open(unit=ID_FILE_FUELHEIGHT, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'groundfuelheight.bin',form='unformatted')
	write (ID_FILE_FUELHEIGHT) fuels%actual_height
	close(ID_FILE_FUELHEIGHT)

	end
!===================================================================================
!===================================================================================
	SUBROUTINE read_cpu(max_num_cpu)

	use constants
	use file_handling_module

	implicit none

	integer, intent(OUT) :: max_num_cpu
	character(STR_LEN) :: fname

	fname = 'Runtime_Advanced_User_Inputs.inp'
	open(ID_FILE_CPU, file=TRIM(workingDirectory)//TRIM(fileSeparator)//fname, &
		status='old',err=1101)
	read(ID_FILE_CPU, *,err=1103,end=1102) max_num_cpu
	if(max_num_cpu <= 0) then
		write (msgoutfile,*) &
			'Invalid maximum number of CPU. Must be >=1.'
		call TerminateProgram()
	endif
	if(max_num_cpu > 8)then
		max_num_cpu = 8
		write (msgoutfile,*) 'Maximum number of CPU limited to 8.'
	endif
	
	return

1101 print*,'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
