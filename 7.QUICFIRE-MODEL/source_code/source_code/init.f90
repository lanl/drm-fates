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
	subroutine init
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! subroutine to initialize main array variables
	! reads in data from the input.dat file
	! init.f90 subroutine is called by main.f90
	! init.f90 calls zwake.f90 function
	! ERP 2001
	! Variable information:
	! uo,vo,wo - are the initial velocity that are prescribed prior to mass
	!            conservation.
	! u,v,w - final velocity field
	! xfo,yfo - denote the x and y locations that represent the
	!           left side center of the building
	! nx,ny, and nz are the number of cells in the x,y and z directions
	! theta - wind angle in standard meteorological format. 0 degrees is
	!    out of the north, 90 degrees is a wind out of the east, etc
	! bld%number - number of buildings in building array
	!
	! Building type designations
	!  bld%btype = 1 regular building (rectangular parrelpiped)
	!  bld%btype = 2 cylindrical/elliptical
	!  bld%btype = 3 pentagon shape
	!  bld%btype  = 9 vegetation, requires an attenuation coefficient (see Cionco, 1965)
	!
	!
	! * note that the velocity field is initialized at the end of the
	!   subroutine.
	! erp 6/2/03 modifications to allow for data input velocity profiles. For
	!     example wind direction can be varied as a function of height
	! erp 6/2/03 modifications to allow for variable grid resolutions to
	!     be entered via the ddx variable. ddx is specified in meters
	! erp 6/5/03 added empirical parameterization subdomain. A subdomain box
	!     may be defined in which the empirical parameterizations are
	!     applied. Outside of this domain conservation of mass is applied only
	!     to an incoming boundary layer profile.
	! erp 7/25/03   This version of qwicurb has the added array zfo(ibuild)
	!     which allows buildings of different sizes to be stacked on one another.
	! erp 8/14/03 This version of the code incorporates Tom Booth's upwind urban
	!     boundary layer work. As such it calls the function zwake.
	! erp 2/11/04 This version has been modified for Nilesh's upwind vortex
	! erp 6/08/04 Modifications to canyon parameterization
	! erp 10/05/04  This version removes the meteorological input information
	!     and puts it in the subroutine met_init.f90 to allow for multi-run
	!     capability.
	! erp 6/30/05 This version adds variable dx,dy and dz capability based on the
	!     work of Tau.
	!
	! cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	
	use constants
	use file_handling_module
	use sor_module
	use diffusion_module
	use canopy_module
	use bld_module
	use grid_module
	use flags_module

	implicit none

	integer nx1,ny1,nz1,nx1ny1

	flag%domain_initialized = 0
	flag%output_initialized = 0
	bld%flag%blending = 3

	call Read_QU_simparams(nx1ny1,nx1,ny1,nz1)
	
	call Read_QU_bld(nx1, nx1ny1)
	
	call Read_QU_fileopt()

	call OpenOutputFiles()
	
	call InitWindVariables()

	call CheckBuildingsToVegetation()

	call InitCanopyVariables()
	
	call Read_QU_landuse()
	
	call Read_QU_parameterRange()
	
	call InitShips()

	call ReadWindsInit()
	
	call InitSorVariables()
	
	return
	
	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadWindsInit
	
	use file_handling_module
	
	implicit none
	
	integer :: met_input_flag
	
	open(unit=ID_FILE_QU_METPARAMS,file=TRIM(workingDirectory)// &
   	TRIM(fileSeparator)//'QU_metparams.inp',status='old')

   read(ID_FILE_QU_METPARAMS,*) ! QUIC version header line
   read(ID_FILE_QU_METPARAMS,*)met_input_flag
   if(met_input_flag .le. 1)then
      call read_quic_met
   elseif(met_input_flag .eq. 2)then
      call read_ittmm5_met
   else
      call read_hotmac_met
   endif
	close(ID_FILE_QU_METPARAMS)
	
	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE InitCanopyVariables()
	
	use canopy_module
	use grid_module
	
	implicit none
	
	allocate(																&
		canopy%ustar(qugrid%nx-1,qugrid%ny-1),						&
		canopy%atten(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
	canopy%atten = 0.
	canopy%ustar = 0.
	if(canopy%number > 0)then
		allocate(															&
			canopy%ktop(qugrid%nx-1,qugrid%ny-1),					&
			canopy%top(qugrid%nx-1,qugrid%ny-1),					&
			canopy%zo(qugrid%nx-1,qugrid%ny-1),						&
			canopy%d(qugrid%nx-1,qugrid%ny-1))
	endif
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE InitSorVariables()

	use sor_module
	use diffusion_module
	use landuse_module
	use bld_module
	use flags_module
	use canopy_module
	use wind_profile_module
	use grid_module
	
	implicit none

	
	sor%alpha1 = 1.			!coefficient for u,v
	sor%alpha2 = 1.			!coefficient for w

	sor%eta = (sor%alpha1/sor%alpha2)**2

	sor%ovalph1 = 0.5 / sor%alpha1**2
	sor%ovalph2 = 0.5 / sor%alpha2**2

	sor%alpha1sq = sor%alpha1**2
	sor%alpha2sq = sor%alpha2**2
	
	sor%ALPHA2_FIRE_VAL = 10.

	sor%omegarelax = 1.78

	sor%one_minus_omegarelax = (1. - sor%omegarelax)	
	sor%res_red_factor = 1./(10.**sor%residual_reduction)
	
	sor%invomegarelax = 1./sor%omegarelax
	
	allocate(																	&
		sor%p1(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),					&
		sor%p2(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),					&
		sor%r(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),					&
		sor%denom0(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),			&
		sor%denom(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
	
	sor%p1 = 0.  ! initialized here because it might not be initialized in the sor depending on SOR_OPTION

	allocate(																	&
		sor%bc%e(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),				&
		sor%bc%f(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),				&
		sor%bc%g(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),				&
		sor%bc%h(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),				&
		sor%bc%m(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),				&
		sor%bc%n(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
	
	if(diff%flag .gt. 0) sor%itermax = sor%itermax/diff%step
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadInitFlags(fname)

	use constants	
	use file_handling_module
	use rnd_num_vars

	implicit none

	character(STR_LEN),intent(IN) :: fname

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)rnd_opt
	if(rnd_opt < 0 .and. rnd_opt /= -1) then
		write (msgoutfile,*)'Invalid random number generator option. Must be > 0 or -1.'
		call TerminateProgram()
	endif
	if(rnd_opt > 0) then
		call initializeSeed(1,rnd_opt)
	else
		call initializeSeed(0,1)
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFireTimes(fname, wind_time_step)

	use constants	
	use time_module
	use file_handling_module
	use ignitions_module
	use flags_module

	implicit none

	character(STR_LEN),intent(IN) :: fname
	integer :: wind_time_step

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)! FIRE TIMES

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%unix_time ! when the fire is ignited
	if(fire_ignition%unix_time < 0) then
		write (msgoutfile,*)'Invalid ignition time. Must be > 0 seconds.'
		call TerminateProgram()
	elseif(fire_ignition%unix_time < qutime%unix(1)) then
		write (msgoutfile,*)'The fire cannot start before the first wind field'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fcatime%len_simulation
	if(fcatime%len_simulation <= 0) then
		write (msgoutfile,*)'Invalid fire simulation time. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fcatime%dt_int
	if(fcatime%dt_int <= 0) then
		write (msgoutfile,*)'Invalid fire simulation time step. Must be > 0.'
		call TerminateProgram()
	endif
	fcatime%dt_real = real(fcatime%dt_int)

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)qutime%firesteps_to_update
	if(qutime%firesteps_to_update < 1) then
		write (msgoutfile,*)'Invalid number of fire steps before updating the winds. Must be >= 1.'
		call TerminateProgram()
	endif
	wind_time_step = fcatime%dt_int * qutime%firesteps_to_update

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%num_fireca_updates
	if(flag%num_fireca_updates < 1) then
		write (msgoutfile,*)'Invalid number of fireca time-steps before dumping on file. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%num_qu_updates
	if(flag%num_qu_updates < 1) then
		write (msgoutfile,*)'Invalid number of wind updates before dumping on file. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%num_emission_and_rad_updates
	if(flag%num_emission_and_rad_updates < 1) then
		write (msgoutfile,*)'Invalid number of fireca time-steps before dumping the emissions and radiation on file. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%num_qu_updates_avewinds
	if(flag%num_qu_updates < 1) then
		write (msgoutfile,*)'Invalid number of wind updates before dumping the averaged winds on file. Must be >= 1.'
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE init_times(qu_updates, qu_updates_avewinds, fireca_updates,	&
		fireca_updates_emiss, fca_it)

	use flags_module
	use time_module
	use ignitions_module
	
	implicit none
	
	integer ::							&
		qu_updates,						&
		qu_updates_avewinds,			&
		fireca_updates,				&
		fireca_updates_emiss,		&
		fca_it
	
	if (flag%isfire == 1) then
		if(qutime%nsteps == 1) then
			allocate(qutime%local(1))
			qutime%local = fcatime%len_simulation*2			
			fire_ignition%unix_time = 0
		else
			allocate(qutime%local(qutime%nsteps))
			qutime%local(1:qutime%nsteps-1) = qutime%unix(2:qutime%nsteps) - qutime%unix(1)
			qutime%local(qutime%nsteps) = fcatime%len_simulation*2
			fire_ignition%unix_time = fire_ignition%unix_time - qutime%unix(1)			
		endif
		qutime%current_iter = 1
		fcatime%current_int = 0
		fcatime%current_real = 0.
		qu_updates = 0
		qu_updates_avewinds = 0
		fireca_updates = 0
		fireca_updates_emiss = 0
		fca_it = 0
		call initialize_fire()  ! Initialize fire and fuel properties		
		call PrintFireCAOutputs(0)
		call PrintFireCAEmissionsAndRadiationOutputs(0)
		call PrintFuelHeight()

	else
		call initializeSeed(0,1)
	endif
		
	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFireGrid(fname)

	use constants
	use fireca_constants_module
	use fireca_module
	use fireca_constants_module
	use file_handling_module
	use plume_const_module
	use grid_module
	use interpolation_module

	implicit none

	character(STR_LEN),intent(IN) :: fname
	real :: temp, domain_height
	integer :: i,k, stretchflag, nextra

	
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! FIRE GRID

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)firegrid%nz
	if(firegrid%nz < 1) then
		write (msgoutfile,*)'Invalid number cells in the vertical direction. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)xratio_int
	if(xratio_int /= 1) then
		!write (msgoutfile,*)'Invalid grid ratio in the x-direction. Must be >= 1.'
		write (msgoutfile,*)'Invalid grid ratio in the x-direction. Must be = 1.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)yratio_int
	if(yratio_int /= 1) then
		!write (msgoutfile,*)'Invalid grid ratio in the y-direction. Must be >= 1.'
		write (msgoutfile,*)'Invalid grid ratio in the y-direction. Must be = 1.'
		call TerminateProgram()
	endif

	if(xratio_int == 1 .and. yratio_int == 1)then
		matching_grids = 1
	else
		matching_grids = 0
	endif

	xratio_real = real(xratio_int)
	yratio_real = real(yratio_int)

	firegrid%dx = qugrid%dx / xratio_real
	firegrid%dy = qugrid%dy / yratio_real
	if(firegrid%dx >= RESOLUTION_SWITCH .or. firegrid%dy >= RESOLUTION_SWITCH) then		
		! Large cells
		print*,'**** Large cell physics turned on'
		call CheckVer(2)
		call SetParamsFireSim()
	else
		! Small cells
		print*,' **** Small cell physics turned on'
		call CheckVer(1)
		call SetParamsFireSim()
	endif

	firegrid%dxi = 1./firegrid%dx
	firegrid%dyi = 1./firegrid%dy
	firegrid%reciprocal = (/firegrid%dxi,firegrid%dyi/)
	firegrid%avg_cellsize = (firegrid%dx+firegrid%dy)*0.5
	firegrid%cellarea = firegrid%dx*firegrid%dy

	firegrid%nx = (qugrid%nx - 1)*xratio_int
	firegrid%ny = (qugrid%ny - 1)*yratio_int
	firegrid%Lx = qugrid%Lx
	firegrid%Ly = qugrid%Ly
	RADIATON_LOSS_FRACT_COEFF = 0.2 + (0.5 - 0.2)/(20. - 2.) * (firegrid%avg_cellsize - 2.)
	RADIATON_LOSS_FRACT_COEFF = max(0.2, RADIATON_LOSS_FRACT_COEFF)
	RADIATON_LOSS_FRACT_COEFF = min(0.5, RADIATON_LOSS_FRACT_COEFF)

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)stretchflag
	if(stretchflag /= 0 .and. stretchflag /= 1) then
		write (msgoutfile,*)'Invalid fire grid stretching flag. Must be [0 1].'
		call TerminateProgram()
	endif

	allocate(												&
		firegrid%dz_array(firegrid%nz+1),			&
		firegrid%cellvol(firegrid%nz+1),				&
		firegrid%delta(firegrid%nz+1),				&
		firegrid%z(firegrid%nz+1),						&
		firegrid%dzmi(firegrid%nz),					&
		firegrid%zm(firegrid%nz+1))

	if(stretchflag == 0) then
		read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)temp
		if(temp <= 0) then
			write (msgoutfile,*)'Invalid fire cell dimension in the z-direction. Must be > 0.'
			call TerminateProgram()
		endif
		firegrid%dz_array = temp
	else
		do i = 1,firegrid%nz
			read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)temp
			if(temp <= 0) then
				write (msgoutfile,*)'Invalid fire cell dimension in the z-direction. Must be > 0.'
				call TerminateProgram()
			endif
			firegrid%dz_array(i) = temp
		enddo
		firegrid%dz_array(firegrid%nz+1) = firegrid%dz_array(firegrid%nz)
	endif

	firegrid%z(1) = 0.
	do i = 1,firegrid%nz
		firegrid%z(i+1) = firegrid%z(i) + firegrid%dz_array(i)
	enddo
	firegrid%Lz = firegrid%z(firegrid%nz+1)
	if(firegrid%Lz > qugrid%Lz) then
		write (msgoutfile,*)'The fire domain is taller than the QUIC domain.'
		call TerminateProgram()
	endif

	do i = 1,firegrid%nz
		firegrid%zm(i) = firegrid%z(i) + 0.5*firegrid%dz_array(i)
	enddo
	firegrid%zm(firegrid%nz+1) = firegrid%z(firegrid%nz+1) + 0.5*firegrid%dz_array(firegrid%nz)

	do k = 1,firegrid%nz
		firegrid%dzmi(k) = 2./(firegrid%dz_array(k)+firegrid%dz_array(k+1))
	enddo

	firegrid%cellvol = firegrid%dz_array * firegrid%cellarea
	firegrid%delta = 4.*firegrid%cellvol**ONE_THIRD

	firegrid%utmx = qugrid%utmx
	firegrid%utmy = qugrid%utmy
	
	! nz cells for num_ignitions_not_absorbed and en_to_atmos: use an additional MAX_FLAME_HEIGHT over firegrid%Lz
	domain_height = min(firegrid%Lz + MAX_FLAME_HEIGHT, qugrid%Lz) - firegrid%Lz
	nextra = max(ceiling(domain_height / firegrid%dz_array(firegrid%nz)), 1)	
	firegrid%nz_en2atmos = firegrid%nz + nextra
	allocate(																	&
		firegrid%z_en2atmos(firegrid%nz_en2atmos + 1),				&
		firegrid%zm_en2atmos(firegrid%nz_en2atmos + 1),				&
		firegrid%dz_array_en2atmos(firegrid%nz_en2atmos + 1),		&
		firegrid%cellvol_en2atmos(firegrid%nz_en2atmos + 1))
	
	firegrid%z_en2atmos(1:firegrid%nz+1) = firegrid%z(1:firegrid%nz+1)
	firegrid%zm_en2atmos(1:firegrid%nz+1) = firegrid%zm(1:firegrid%nz+1)
	
	firegrid%dz_array_en2atmos(1:firegrid%nz+1) = firegrid%dz_array(1:firegrid%nz+1)
	firegrid%dz_array_en2atmos(min(firegrid%nz+2, firegrid%nz_en2atmos+1): firegrid%nz_en2atmos+1) = &
		firegrid%dz_array(firegrid%nz)	
	do i = firegrid%nz+1, firegrid%nz_en2atmos
		firegrid%z_en2atmos(i+1) = firegrid%z_en2atmos(i) + firegrid%dz_array_en2atmos(i)
		firegrid%zm_en2atmos(i+1) = firegrid%z_en2atmos(i+1) + firegrid%dz_array_en2atmos(i+1) * 0.5
	enddo
	! Check last cell => cannot be higher than qugrid%Lz
	if(firegrid%z_en2atmos(firegrid%nz_en2atmos+1) > qugrid%Lz) then
		firegrid%z_en2atmos(firegrid%nz_en2atmos+1) = qugrid%Lz
		firegrid%dz_array_en2atmos(firegrid%nz_en2atmos+1) =	&
			firegrid%z_en2atmos(firegrid%nz_en2atmos+1) - firegrid%z_en2atmos(firegrid%nz_en2atmos)
	endif
	firegrid%cellvol_en2atmos = firegrid%dz_array_en2atmos * firegrid%cellarea
	
	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFireFuel(fname)

	use constants
	use ignitions_module
	use file_handling_module
	use fireca_module

	implicit none

	character(STR_LEN),intent(IN) :: fname

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! FUEL

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%density_flag
	if(fuels%density_flag /= 1 .and. fuels%density_flag /= 2 .and. &
		fuels%density_flag /= 3 .and. fuels%density_flag /= 4 .and. &
		fuels%density_flag /= 66) then
		write (msgoutfile,*)'Invalid fuel density flag. Accepted values [1 4].'
		call TerminateProgram()
	endif

	if(fuels%density_flag == 1) then !uniform
		read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%input_density
		if(fuels%input_density <= 0) then
			write (msgoutfile,*)'Invalid fuel density. Must be > 0.'
			call TerminateProgram()
		endif
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%moisture_flag
	if(fuels%moisture_flag /= 1 .and. fuels%moisture_flag /= 2 .and. &
		fuels%moisture_flag /= 3 .and. fuels%moisture_flag /= 4 .and. &
		fuels%moisture_flag /= 66) then
		write (msgoutfile,*)'Invalid fuel moisture flag. Accepted values [1 4].'
		call TerminateProgram()
	endif
	if(fuels%moisture_flag == 1) then !uniform
		read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%input_moisture
		if(fuels%input_moisture < 0 .or. fuels%input_moisture > 1) then
			write (msgoutfile,*)'Invalid fuel moisture. Must be [0 1].'
			call TerminateProgram()
		endif
	endif

	if(fuels%density_flag == 1) then
		read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%height_flag
		if(fuels%height_flag /= 1 .and. fuels%height_flag /= 2 .and. &
			fuels%height_flag /= 3 .and. fuels%height_flag /= 4) then
			write (msgoutfile,*)'Invalid fuel height flag. Accepted values [1 4].'
			call TerminateProgram()
		endif

		if(fuels%height_flag == 1)then !uniform
			read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fuels%input_height
			if(fuels%input_height <= 0 ) then
				write (msgoutfile,*)'Invalid fuel height. Must be > 0.'
				call TerminateProgram()
			endif
		endif
	endif

	! -- Need to read in the Firetech grid info
	if(fuels%height_flag == 4 .or. fuels%density_flag == 4. .or. &
	   fuels%density_flag == 66. .or. &
		fuels%moisture_flag == 4 .or. fire_ignition%flag == 6) then

		call ImportFTGrid()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ImportLineFire(fname)

	use constants
	use file_handling_module
	use ignitions_module
	use grid_module
	
	implicit none

	character(STR_LEN),intent(IN) :: fname

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xo
	if(fire_ignition%xo > qugrid%Lx .or. fire_ignition%xo < 0) then
		write (msgoutfile,*)'Fire south-west x-coordinate is outside the domain.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%yo
	if(fire_ignition%yo > qugrid%Ly .or. fire_ignition%yo < 0) then
		write (msgoutfile,*)'Fire south-west y-coordinate is outside the domain.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xlen
	if(fire_ignition%xlen <= 0) then
		write (msgoutfile,*)'Fire length in the x-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%xlen + fire_ignition%xo > qugrid%Lx) then
		write (msgoutfile,*)'Fire in the x-direction exits the domain.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%ylen
	if(fire_ignition%ylen <= 0) then
		write (msgoutfile,*)'Fire length in the y-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%ylen + fire_ignition%yo > qugrid%Ly) then
		write (msgoutfile,*)'Fire in the y-direction exits the domain.'
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ImportSquareRing(fname)
	
	use constants
	use ignitions_module
	use file_handling_module
	use grid_module
	
	implicit none

	character(STR_LEN),intent(IN) :: fname

	call ImportLineFire(fname)

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xwidth

	if(fire_ignition%xwidth <= 0) then
		write (msgoutfile,*)'Fire ring width in the x-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%xlen - fire_ignition%xwidth * 2. < qugrid%dx) then
		write (msgoutfile,*)'Fire ring width in the x-direction must be < ',	&
			fire_ignition%xlen-qugrid%dx
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%ywidth

	if(fire_ignition%ywidth <= 0) then
		write (msgoutfile,*)'Fire ring width in the y-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%ylen - fire_ignition%ywidth * 2. < qugrid%dy) then
		write (msgoutfile,*)'Fire ring width in the y-direction must be < ',	&
			fire_ignition%ylen-qugrid%dy
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ImportCircularRing(fname)

	use constants
	use file_handling_module
	use ignitions_module
	use grid_module

	implicit none

	character(STR_LEN),intent(IN) :: fname

	! circular ring
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xo
	if(fire_ignition%xo > qugrid%Lx .or. fire_ignition%xo < 0) then
		write (msgoutfile,*)'Fire south-west x-coordinate is outside the domain.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%yo
	if(fire_ignition%yo > qugrid%Ly .or. fire_ignition%yo < 0) then
		write (msgoutfile,*)'Fire south-west y-coordinate is outside the domain.'
		call TerminateProgram()
	endif
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xlen
	if(fire_ignition%xlen <= 0) then
		write (msgoutfile,*)'Fire length in the x-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%xlen + fire_ignition%xo > qugrid%Lx) then
		write (msgoutfile,*)'Fire in the x-direction exits the domain.'
		call TerminateProgram()
	endif
	fire_ignition%ylen = fire_ignition%xlen

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%xwidth
	fire_ignition%ywidth = fire_ignition%xwidth
	if(fire_ignition%xwidth <= 0) then
		write (msgoutfile,*)'Fire ring width in the x-direction must be > 0.'
		call TerminateProgram()
	elseif(fire_ignition%xlen - fire_ignition%xwidth * 2. < qugrid%dx) then
		write (msgoutfile,*)'Fire ring width in the x-direction must be < ',fire_ignition%xlen-qugrid%dx
		call TerminateProgram()
	elseif(fire_ignition%xlen - fire_ignition%xwidth * 2. < qugrid%dy) then
		write (msgoutfile,*)'Fire ring width in the y-direction must be < ',fire_ignition%ylen-qugrid%dy
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFireIgnitions(fname)

	use constants
	use ignitions_module
	use file_handling_module

	implicit none

	character(STR_LEN),intent(IN) :: fname

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! IGNITION TXT

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%flag
	if(fire_ignition%flag /= 1 .and. fire_ignition%flag /= 2 .and. fire_ignition%flag /= 3 .and. &
		fire_ignition%flag /= 4 .and. fire_ignition%flag /= 5 .and. fire_ignition%flag /= 6 .and. &
		fire_ignition%flag /= 7) then
		write (msgoutfile,*)'Invalid fire source flag. Must be [1 7].'
		call TerminateProgram()
	endif

	if(fire_ignition%flag == 1)then
		call ImportLineFire(fname)

	elseif(fire_ignition%flag == 2)then
		call ImportSquareRing(fname)

	elseif(fire_ignition%flag == 3)then
		call ImportCircularRing(fname)

	! elseif(fire_ignition%flag == 4)then
	! 	call ReadIgnitionPattern()

	elseif(fire_ignition%flag == 5)then
		call ReadIgnitionAerial()

	elseif(fire_ignition%flag == 7)then
		call ReadIgnitionAerial()
	endif

	if(fire_ignition%flag /= 4 .and. fire_ignition%flag /= 5)then
		read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fire_ignition%num_ign_percell
		if(fire_ignition%num_ign_percell < 0) then
			write (msgoutfile,*)'Invalid number of ignitions. Must be > 0.'
			call TerminateProgram()
		endif
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFireFlags(fname)

	use constants
	use grid_module
	use file_handling_module
	use flags_module
	use fireca_module

	implicit none

	character(STR_LEN),intent(IN) :: fname

	!==============================================================================
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! FIREBRANDS
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)fb%flag
	if(fb%flag /= 0 .and. fb%flag /= 1) then
		write (msgoutfile,*)'Invalid firebrand flag. Must be [0 1].'
		call TerminateProgram()
	endif
	call OpenFBOutputFiles()

	!==============================================================================
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)  ! FILES FLAGS
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%en2atm_file
	if(flag%en2atm_file /= 0 .and. flag%en2atm_file /= 1) then
		write (msgoutfile,*)'Invalid energy to atmosphere file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%reactrate_file
	if(flag%reactrate_file /= 0 .and. flag%reactrate_file /= 1) then
		write (msgoutfile,*)'Invalid reaction rate file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%fuelmass_file
	if(flag%fuelmass_file /= 0 .and. flag%fuelmass_file /= 1) then
		write (msgoutfile,*)'Invalid fuel mass file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%qfwinds_file
	if(flag%qfwinds_file /= 0 .and. flag%qfwinds_file /= 1) then
		write (msgoutfile,*)'Invalid fire winds file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%quwinds_inst_file
	if(flag%quwinds_inst_file /= 0 .and. flag%quwinds_inst_file /= 1) then
		write (msgoutfile,*)'Invalid QUIC-URB fire-determined winds (instantaneous) file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%quwinds_ave_file
	if(flag%quwinds_ave_file /= 0 .and. flag%quwinds_ave_file /= 1) then
		write (msgoutfile,*)'Invalid QUIC-URB fire-determined winds (averaged) file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%plume_traj_file
	if(flag%plume_traj_file /= 0 .and. flag%plume_traj_file /= 1 .and. flag%plume_traj_file /= 2) then
		write (msgoutfile,*)'Invalid plume trajectories file flag. Must be [0 2].'
		call TerminateProgram()
	endif
	call OpenTrajectoryFiles()

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%moisture_file
	if(flag%moisture_file /= 0 .and. flag%moisture_file /= 1) then
		write (msgoutfile,*)'Invalid moisture file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%int_massburnt_file
	if(flag%int_massburnt_file /= 0 .and. flag%int_massburnt_file /= 1) then
		write (msgoutfile,*) &
			'Invalid vertically-integrated % of mass burnt file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%plume_loc_file
	if(flag%plume_loc_file /= 0 .and. flag%plume_loc_file /= 1) then
		write (msgoutfile,*) &
			'Invalid plume location file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%emissions_file
	if(flag%emissions_file < 0 .or. flag%emissions_file > 4) then
		write (msgoutfile,*) &
			'Invalid emissions file flag. Must be [0 4].'
		call TerminateProgram()
	endif
	if(flag%emissions_file == 4)then
		print*,'Fuel density flag is set to 4 to use the emission library approach'
		flag%fuelmass_file = 1	
	endif

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%thermalrad_file
	if(flag%thermalrad_file /= 0 .and. flag%thermalrad_file /= 1) then
		write (msgoutfile,*) &
			'Invalid thermal radiation file flag. Must be [0 1].'
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE OpenFBOutputFiles()
	
	use file_handling_module
	use fireca_module

	implicit none

	character(STR_LEN) :: fname

	if(fb%flag == 1) then
		fname = 'firebrands.bin'
		open(ID_FILE_FB_OUT,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
			form='unformatted',status='replace',err = 1101)

		!fname = 'firebrands_lost.bin'
		!open(ID_FILE_FB_LOST,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
		!	form='unformatted',status='replace',err = 1101)
	endif

	return

1101 print*,'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE OpenTrajectoryFiles()

	use file_handling_module
	use constants
	use flags_module

	implicit none

	if(flag%plume_traj_file == 2)then
		! Open to delete previous results
		open(ID_FILE_TRAJ,file=TRIM(workingDirectory)//TRIM(fileSeparator)//&
			'plume_trajectory.bin',status='replace', form='unformatted')
		open(ID_FILE_TRAJ_MERGE,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
			'plume_mergetrajectory.bin',status='replace', form='unformatted')
		
	elseif(flag%plume_traj_file == 1)then
		open(ID_FILE_TRAJ,file=TRIM(workingDirectory)//TRIM(fileSeparator)//&
			'plume_trajectory.csv',status='replace')
		open(ID_FILE_TRAJ_MERGE,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
			'plume_mergetrajectory.csv',status='replace')
				
		write(ID_FILE_TRAJ,'(a)') 'id [-], fire time [s], plume time [s], '// &
			'x [m], y [m], z [m], uc [m/s], vc [m/s], wc [m/s], plume_radius h [m], ehx, ehy, ehz, '// &
			'plume_radius n [m], enx, eny, enz, '
	endif
	
	return

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE CanopyConstantsInitialization

	use canopy_module

	implicit none

	canopy%FUEL_ATTEN_COEFF	= 2.
	canopy%FUEL_CANOPY_Z0	= 0.01
	canopy%BLDG_ATTEN_COEFF = 2.
	canopy%BLDG_CANOPY_Z0	= 0.01
	canopy%BLDG_FUEL_DENS	= 0.5
	canopy%UPDATE_WINDS		= 1

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE FirebrandsConstantsInitialization()

	use fireca_constants_module
	use fireca_module

	implicit none 

	fb%time_step = 1
	fb%FRACTION_LAUNCHED = 0.05
	fb%c_s = 40.
	fb%num_deposited = 500.
	fb%FRACTION_LAUNCHED_to_RT_ratio = 10.
	fb%frac_of_max_size = 0.75
	fb%min_b_value_coef = 50.
	fb%min_theta_value = GREEK_PI/6.
	fb%w_mult = 5.				
	fb%LAUNCH_TIME	= 10
	fb%min_number_of_ignitions = 50
	fb%max_number_of_ignitions = 100
	fb%germination_delay = 300

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE readQuicFireFile()

	use constants
	use file_handling_module
	use fireca_module
	use sor_module
	use plume_const_module
	use plume_module
	use bld_module
	use interpolation_module
	use grid_module
	use flags_module
	use winds_module

	implicit none

	character(STR_LEN) :: fname, fnameign
	integer :: i, wind_time_step

	call PlumeConstantsInitialization()
	call CanopyConstantsInitialization()
	call FirebrandsConstantsInitialization()
	
	!==============================================================================
	fnameign = 'ignite_selected.dat'
	open(ID_FILE_IGNITE_SEL,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fnameign), &
		form='unformatted',status='replace')

	! read in input files
	fname = 'QUIC_fire.inp'
	open(ID_FILE_QFIRE_INP,file=TRIM(workingDirectory)//TRIM(fileSeparator)//TRIM(fname), &
		status='old',err = 1101)

	! initialization in case of error before they are defined
	msgoutfile = STDOUT !DBGFILE
	open(DBGFILE,file=TRIM(workingDirectory)//TRIM(fileSeparator)//'DebugMsg.txt',status='replace')

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)flag%isfire
	if(flag%isfire /= 0 .and. flag%isfire /= 1) then
		write (msgoutfile,*)'Invalid fire flag. Must be [0 1].'
		call TerminateProgram()
	endif
	if (flag%isfire == 0)then
		close(ID_FILE_QFIRE_INP)
		return
	endif

	! THESE OPERATIONS ARE COMPLETED ONLY IF flag%isfire == 1
	call ReadInitFlags(fname)

	call ReadFireTimes(fname, wind_time_step)

	call ReadFireGrid(fname)
	call output_vert_grid_qu()
	
	call ReadFiretecFilePath(fname)

	call ReadFireFuel(fname)
	
	call ReadFireTopo(fname)

	call ReadFireIgnitions(fname)

	call ReadFireFlags(fname)

	close(ID_FILE_QFIRE_INP)
	!==============================================================================

	call ReadAdvancedPlumeFile(wind_time_step)

	call ReadAdvancedBldgFile()
	!==============================================================================

	! Index of the start and end cell in the fire grid for each vertical cell of the quic domain
	allocate(fca_2_qu_kstart(qugrid%nz),fca_2_qu_kend(qugrid%nz))
	allocate(quwinds%wplume(qugrid%nx,qugrid%ny,qugrid%nz))

	allocate(fcawinds%u(firegrid%nx,firegrid%ny,firegrid%nz+1), &
		fcawinds%v(firegrid%nx,firegrid%ny,firegrid%nz+1), &
		fcawinds%w(firegrid%nx,firegrid%ny,firegrid%nz+1))
	allocate(firegrid%xcenters(firegrid%nx),firegrid%ycenters(firegrid%ny))
	allocate(quwinds%u_mc(qugrid%nx,qugrid%ny,qugrid%nz),	&
		quwinds%v_mc(qugrid%nx,qugrid%ny,qugrid%nz),quwinds%w_mc(qugrid%nx,qugrid%ny,qugrid%nz))
	fcawinds%u = 0.
	fcawinds%v = 0.
	fcawinds%w = 0.

	allocate(plume(plume_const%MAX_NUM_PLUMES_TIMESTEP))
	allocate(to_delete(plume_const%MAX_NUM_PLUMES_TIMESTEP))

	!$OMP parallel do private(i)	
	do i = 1,plume_const%MAX_NUM_PLUMES_TIMESTEP
		plume(i)%id = 0
	enddo
	!$OMP end parallel do 
	max_id_plume = 0

	do i = 1,firegrid%nx
		firegrid%xcenters(i) = (real(i)-0.5)*firegrid%dx
	enddo

	do i = 1,firegrid%ny
		firegrid%ycenters(i) = (real(i)-0.5)*firegrid%dy
	enddo

	total_fires = 0

	sor%mc_buffer = sor%itermax+1

	firegrid%cellvol = firegrid%dx * firegrid%dy * firegrid%dz_array

	return

1101 print*,'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	END
	!============================================================================
	!============================================================================
	SUBROUTINE ReadFireTopo(fname)
	
	use constants
	use flags_module
	use topo_module
	
	implicit none

	character(STR_LEN),intent(IN) :: fname	  ! N/A, file name

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! header
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) topo%flag
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) topo%file_path

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()


	END
	!============================================================================
	!============================================================================
	SUBROUTINE ReadAerialIgnitionPattern(fname)

	use constants
	use fireca_module
	use file_handling_module
	use ignitions_module
	use grid_module

	implicit none

	character(STR_LEN),intent(IN) :: fname	  ! N/A, file name

	integer ::					&
		i,							& ! N/A, loop index
		ip,jp,					& ! N/A, cell ignited
		count_ign				  ! N/A, number of ignitions in the FireCA domain

	real ::						&
		temp_time,				& ! s, used to read in time
		temp_x,temp_y			  ! m, used to compute the position of the ignition

	character(STR_LEN) :: tchar1  ! N/A, used to read in a namelist file

	real,allocatable,dimension(:) ::			&
		temp_x_ign,temp_y_ign,					& ! m, temporary arrays for ignition location
		temp_time_ign								  ! s, temporary arrays for ignition time


	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,'(a)')tchar1
	i = scan(tchar1,"=")
	read(tchar1(i+1:),*)fire_ignition%npoints
	write (*,*) 'number of ignition points',fire_ignition%npoints
	allocate(											&
		temp_x_ign(fire_ignition%npoints),			&
		temp_y_ign(fire_ignition%npoints),			&
		temp_time_ign(fire_ignition%npoints))
	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,*)
	count_ign = 0
	do i = 1,fire_ignition%npoints
		read(ID_FILE_IGNITE,*)ip,jp,temp_time
		write (*,*) i, ip,jp,temp_time

		! Input check
		if(ip < 1 .or. ip > ft%nx) then
			write (msgoutfile,*)'Invalid x-coordinate (out of Firetec domain) in ',fname,' at line ',i
			call TerminateProgram()
		endif
		if(jp < 1 .or. jp > ft%ny) then
			write (msgoutfile,*)'Invalid y-coordinate (out of Firetec domain) in ',fname,' at line ',i
			call TerminateProgram()
		endif
		if(temp_time < 0) then
			write (msgoutfile,*)'Invalid time < 0 in ',fname,' at line ',i
			call TerminateProgram()
		endif

		! Check if ignition is in the FireCA domain
		temp_x = (real(ip)-0.5)*ft%dx + ft%utmx - qugrid%utmx
		temp_y = (real(jp)-0.5)*ft%dy + ft%utmy - qugrid%utmy

		if(temp_x >= 0. .and. temp_x <= qugrid%Lx .and. &
			temp_y >= 0. .and. temp_y <= qugrid%Ly) then

			count_ign = count_ign + 1
			temp_x_ign(count_ign) = temp_x
			temp_y_ign(count_ign) = temp_y
			temp_time_ign(count_ign) = temp_time
		endif
	enddo

	if(count_ign == 0)then
		write (msgoutfile,*)'No ignitions present in the FireCA domain'
		call TerminateProgram()
	endif

	allocate(										&
		fire_ignition%x(count_ign),			&
		fire_ignition%y(count_ign),			&
		fire_ignition%z(count_ign),			&
		fire_ignition%time(count_ign),		&
		fire_ignition%radius(count_ign),		&
		fire_ignition%new_num(count_ign))

	fire_ignition%z = 0.
	fire_ignition%radius = min(ft%dx,ft%dy)
	fire_ignition%new_num = 100
	fire_ignition%x = temp_x_ign(1:count_ign)
	fire_ignition%y = temp_y_ign(1:count_ign)
	fire_ignition%time = temp_time_ign(1:count_ign)
	fire_ignition%npoints = count_ign

	deallocate(temp_x_ign,temp_y_ign,temp_time_ign)

	END
	!============================================================================
	!============================================================================
	SUBROUTINE ATVIgnitionPattern(fname)

	use constants
	use fireca_module
	use time_module
	use file_handling_module
	use ignitions_module
	use grid_module

	implicit none

	integer ::					&
		i,j,						& ! N/A, loop index
		nlines,					& ! N/A, number of ATV lines
		count_ign				  ! N/A, number of ignitions in the FireCA domain

	real ::						&
		dtt,						& ! s, time step
		dxx,dyy,					& ! m, space covered over one time step in x and y directions
		x,y,						& ! m, ignition location
		t,							& ! s, ignition time
		ireal

	character(STR_LEN) :: fname	  ! N/A, file name
	character(STR_LEN) :: tchar1  ! N/A, used to read in a namelist file


	integer,allocatable,dimension(:) ::		&
		atv_np										  ! N/A, number of points per line

	real,allocatable,dimension(:,:) ::		&
		atv_x,atv_y,								& ! m, position of the ignitions deposited by ATV
		atv_time									     ! s, time of the ignitions deposited by ATV

	real,dimension(:),allocatable ::			&
		temp_time,									& ! s, used to read in time
		temp_x,temp_y								  ! m, used to compute the position of the ignition


	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,'(a)')tchar1
	i = scan(tchar1,"=")
	read(tchar1(i+1:),*)nlines
	allocate(atv_x(nlines,2),atv_y(nlines,2),atv_time(nlines,2),atv_np(nlines))
	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,*)
	read(ID_FILE_IGNITE,*)

	fire_ignition%npoints = 0
	do i = 1,nlines
		read(ID_FILE_IGNITE,*)					&
			atv_x(i,1),atv_y(i,1),				&
			atv_x(i,2),atv_y(i,2),				&
			atv_time(i,1),atv_time(i,2)

		do j = 1,2
			if(atv_time(i,j) < 0) then
				write (msgoutfile,*)'Invalid time < 0 in ',fname,' at line ',i
				call TerminateProgram()
			endif
			atv_x(i,j) = atv_x(i,j) + ft%utmx - qugrid%utmx
			atv_y(i,j) = atv_y(i,j) + ft%utmy - qugrid%utmy
		enddo
		
		atv_np(i) =	max(ceiling( abs(atv_x(i,1)-atv_x(i,2)) / ft%dx) , &
			ceiling( abs(atv_y(i,1)-atv_y(i,2)) / ft%dy) )

	enddo
	fire_ignition%npoints	= sum(atv_np)

	allocate(								&
		temp_x(fire_ignition%npoints),		&
		temp_y(fire_ignition%npoints),		&
		temp_time(fire_ignition%npoints))

	count_ign = 0
	do j = 1,nlines
		dtt = (atv_time(j,2) - atv_time(j,1)) / real(atv_np(j))
		dxx = (atv_x(j,2) - atv_x(j,1)) / real(atv_np(j))
		dyy = (atv_y(j,2) - atv_y(j,1)) / real(atv_np(j))
		do i = 1,atv_np(j)
			ireal = real(i-1)
			x = atv_x(j,1) + dxx * ireal
			y = atv_y(j,1) + dyy * ireal
			t = atv_time(j,1) + dtt * ireal

			if(x >= 0. .and. x <= qugrid%Lx .and. y >= 0. .and. y <= qugrid%Ly .and. t <= fcatime%len_simulation) then
				count_ign = count_ign + 1

				temp_time(count_ign) = t
				temp_x(count_ign) = x
				temp_y(count_ign) = y
			endif
		enddo
	enddo

	if(count_ign == 0)then
		write (msgoutfile,*)'No ignitions present in the FireCA domain'
		call TerminateProgram()
	endif

	allocate(										&
		fire_ignition%x(count_ign),			&
		fire_ignition%y(count_ign),			&
		fire_ignition%z(count_ign),			&
		fire_ignition%time(count_ign),		&
		fire_ignition%radius(count_ign),		&
		fire_ignition%new_num(count_ign))

	fire_ignition%x = temp_x(1:count_ign)
	fire_ignition%y = temp_y(1:count_ign)
	fire_ignition%z = 0.
	fire_ignition%time = temp_time(1:count_ign)
	fire_ignition%radius = max(ft%dx,ft%dy)	
	fire_ignition%new_num = 100
	fire_ignition%npoints = count_ign

	deallocate(atv_x,atv_y,temp_x,temp_y,temp_time)


	END
	!============================================================================
	!============================================================================
	SUBROUTINE ReadIgnitionAerial()

	use constants
	use fireca_module
	use file_handling_module
	use ignitions_module
	use fireca_module
	use grid_module

	implicit none

	integer ::			&
		i,					&
		ignition_flag	! N/A, what time of ignite file it is

	character(STR_LEN) :: fname	  ! N/A, file name
	character(STR_LEN) :: tchar1  ! N/A, used to read in a namelist file


	! If it has not already been read in
	if(fuels%height_flag /= 4 .and. fuels%density_flag /= 4 .and. &
	     fuels%density_flag /= 66 .and. &
		fuels%moisture_flag /= 4 .and. fire_ignition%flag /= 6) then
		call ImportFTGrid()
	endif

	fname = 'ignite.dat'
	open(ID_FILE_IGNITE,file=TRIM(ft%file_path)//TRIM(fileSeparator)// &
		trim(fname),status='old',err = 1101)
	read(ID_FILE_IGNITE,'(a)')tchar1
	i = scan(tchar1,"=")
	read(tchar1(i+1:),*)ignition_flag
	if(ignition_flag /= 4 .and. ignition_flag /= 5)then
		write (msgoutfile,*)'Invalid aerial flag at line 1. Accepted values [4 5].'
		call TerminateProgram()
	endif

	if(ignition_flag == 4) then
		write (*,*) 'reading aerial ignition pattern'
		call ReadAerialIgnitionPattern(fname)

	elseif(ignition_flag == 5) then
		call ATVIgnitionPattern(fname)

	endif

	return

1101 write(msgoutfile,*)'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE PlumeConstantsInitialization()
	
	use plume_const_module
	use plume_module	
	use constants
	use grid_module

	implicit none
	
	real, parameter ::				&
		cp = 1.004,						& ! kJ/kg/K, air heat capacity
		air_temp = 300					  ! K, air temperature
	integer :: i, j, k
	
	plume_const%en2atmos_mult = GRAVITY / (air_temp * cp * RHO_AIR)		! Fb = g/T *(H/rho/cp)  => Fb = mult * H * volume   (dvol used to convert units)
	plume_const%max_fraction_traj = 0.7	
	plume_const%max_angle = 30. * PI	/ 180.										! rad, default 60 deg
	plume_const%ALPHA_ENTR = 0.11
	plume_const%BETA_ENTR = 0.6											! N/A, entrainment coefficient
	plume_const%BETA_ENTR2 = plume_const%BETA_ENTR**2				! entrainment coefficient squared			
	plume_const%WC_EXPONENT = 3											! N/A, exponent to elevate w to compute overlap
	plume_const%INV_WC_EXPONENT = 1./plume_const%WC_EXPONENT		! N/A, exponent to elevate w to compute overlap
	plume_const%SPEEDS_RATIO = 0.1										! N/A, if wc < SPEEDS_RATIO*wind_speed, the plume is deleted
	plume_const%WC_MAX = 200												! m/s, maximum updraft velocity
	plume_const%WC_MIN = 0.1												! m/s, if wc < WC_MIN, the plume is deleted
	plume_const%time_step_flag = 0
	
	if(plume_const%time_step_flag == 1)	then	
		plume_const%time_step = 1e8
		plume_const%min_time_step = 0.1
	else
		plume_const%time_step = 1.
	endif
	
	plume_const%do_w_fill = 1
	plume_const%MAX_NUM_PLUMES_TIMESTEP = 85e3	
	plume_const%BRUNT_VAISALA_FREQ2 = 0
	plume_const%zvirt0 = sqrt(qugrid%dx*qugrid%dy / PI) / plume_const%BETA_ENTR
	plume_const%init_w_en_fract = 0.
	
	plume_diagnostics%max_height = 0
	plume_diagnostics%max_w = 0
	
	!allocate(plume_lock(qugrid%nx, qugrid%ny, qugrid%nz))
	!do k = 1, qugrid%nz
	!	do j = 1, qugrid%ny
	!		do i = 1, qugrid%nx
	!			CALL omp_init_lock (plume_lock(i, j, k))
	!		enddo
	!	enddo
	!enddo
	
	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadIgnitionPattern()
	
	use constants
	use file_handling_module
	use ignitions_module
	use grid_module
	
	implicit none

	integer :: i
	character(STR_LEN) :: fname

	fname = 'QF_IgnitionPattern.inp'
	open(ID_FILE_IGNPATTERN,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(fname),status='old',err = 1101)
	read(ID_FILE_IGNPATTERN,*) fire_ignition%npoints
	read(ID_FILE_IGNPATTERN,*) ! header
	allocate(														&
		fire_ignition%x(fire_ignition%npoints),			&
		fire_ignition%y(fire_ignition%npoints),			&
		fire_ignition%z(fire_ignition%npoints),			&
		fire_ignition%time(fire_ignition%npoints),		&
		fire_ignition%radius(fire_ignition%npoints),		&
		fire_ignition%new_num(fire_ignition%npoints))
	do i = 1,fire_ignition%npoints
		read(ID_FILE_IGNPATTERN,*,err=1102,end=1103)		&
			fire_ignition%time(i),								&
			fire_ignition%x(i),									&
			fire_ignition%y(i),									&
			fire_ignition%z(i),									&
			fire_ignition%radius(i),							&
			fire_ignition%new_num(i)
		if(fire_ignition%time(i) < 0) then
			write (msgoutfile,*)'Invalid time < 0 in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
		if(fire_ignition%x(i) < 0 .or. fire_ignition%x(i) > qugrid%Lx) then
			write (msgoutfile,*)'Invalid x-coordinate (out of domain) in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
		if(fire_ignition%y(i) < 0 .or. fire_ignition%y(i) > qugrid%Ly) then
			write (msgoutfile,*)'Invalid y-coordinate (out of domain) in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
		if(fire_ignition%z(i) < 0 .or. fire_ignition%z(i) > firegrid%Lz) then
			write (msgoutfile,*)'Invalid z-coordinate (out of domain) in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
		if(fire_ignition%radius(i) <= 0) then
			write (msgoutfile,*)'Invalid radius (< 0) in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
		if(fire_ignition%new_num(i) <= 0)then
			write (msgoutfile,*)'Invalid number of new ignitions (<= 0) in QF_IgnitionPattern.inp at line ',i
			call TerminateProgram()
		endif
	enddo
	close(ID_FILE_IGNPATTERN)

	return

1101 write(msgoutfile,*)'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadAdvancedPlumeFile(wind_time_step)

	use constants
	use fireca_module
	use file_handling_module
	use sor_module
	use plume_const_module
	use plume_module	

	implicit none

	character(STR_LEN) :: fname
	integer, intent(IN) :: wind_time_step

	sor%option = 0	
	
	! Read advanced input file
	fname = 'QFire_Plume_Advanced_User_Inputs.inp'
	open(ID_FILE_ADV_PLUME,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
		status='old',err = 1104) ! 1104 = continue

	read(ID_FILE_ADV_PLUME,*,err=1102,end=1103)plume_const%MAX_NUM_PLUMES_TIMESTEP
	if(plume_const%MAX_NUM_PLUMES_TIMESTEP <= 0 .or. plume_const%MAX_NUM_PLUMES_TIMESTEP > 500e3) then
		write (msgoutfile,*) &
			'Invalid maximum number of plumes per timestep. Must be (0 500k].'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1102,end=1103)plume_const%WC_MIN
	if(plume_const%WC_MIN <= 0) then
		write (msgoutfile,*) 'Invalid minimum plume vertical velocity. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1102,end=1103)plume_const%WC_MAX
	if(plume_const%WC_MAX <= plume_const%WC_MIN) then
		write (msgoutfile,*) &
			'Invalid maximum plume vertical velocity. Must be > of minimum vertical velocity.'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1102,end=1103)plume_const%SPEEDS_RATIO
	if(plume_const%SPEEDS_RATIO <= 0 .or. plume_const%SPEEDS_RATIO > 1) then
		write (msgoutfile,*) 'Invalid minimum plume speed ratio. Must be (0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%BRUNT_VAISALA_FREQ2
	if(plume_const%BRUNT_VAISALA_FREQ2 < 0) then
		write (msgoutfile,*) 'Invalid Brunt-Vaisala frequency. Must be >= 0.'
		call TerminateProgram()
	endif
	if(plume_const%BRUNT_VAISALA_FREQ2 > 0)then
		plume_const%max_time = PI / sqrt(plume_const%BRUNT_VAISALA_FREQ2)
	else
		plume_const%max_time = -1
	endif

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)CREEPING_FLAG
	if(CREEPING_FLAG /= 0 .and. CREEPING_FLAG /= 1) then
		write (msgoutfile,*) 'Invalid creeping flag. Must be [0 1].'
		call TerminateProgram()
	endif
	
	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%time_step_flag
	if(plume_const%time_step_flag /= 0 .and. plume_const%time_step_flag /= 1) then
		write (msgoutfile,*) 'Invalid plume time step. Must be [0 1].'
		call TerminateProgram()
	endif	

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%min_time_step
	if(plume_const%min_time_step < 0.05) then
		write (msgoutfile,*) 'Invalid plume time step. Must be >= 0.05 s'
		call TerminateProgram()
	endif
	if(plume_const%time_step_flag == 0) plume_const%time_step = plume_const%min_time_step

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)sor%option
	if(sor%option /= 0 .and. sor%option /= 1) then
		write (msgoutfile,*) 'Invalid SOR option. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)sor%ALPHA2_FIRE_VAL
	if(sor%ALPHA2_FIRE_VAL <= 0 .or. sor%ALPHA2_FIRE_VAL > 20.) then
		write (msgoutfile,*) 'Invalid alpha 2 value. Must be (0 20].'
		call TerminateProgram()
	endif

	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%init_w_en_fract
	if(plume_const%init_w_en_fract < 0) then
		write (msgoutfile,*) 'Invalid fraction of ignition energy that is released to the atmosphere. Must be >= 0.'
		call TerminateProgram()
	endif
		
	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%max_angle
	if(plume_const%max_angle <= 0) then
		write (msgoutfile,*) 'Invalid maximum angle for merging. Must be > 0.'
		call TerminateProgram()
	endif
	plume_const%max_angle = plume_const%max_angle * PI / 180.
	
	read(ID_FILE_ADV_PLUME,*,err=1104,end=1104)plume_const%max_fraction_traj	
	if(plume_const%max_fraction_traj <= 0 .or. plume_const%max_fraction_traj > 100.) then
		write (msgoutfile,*) 'Invalid maximum trajectory fraction for merging. Must be 0 < max_fraction_traj <= 1.'
		call TerminateProgram()
	endif
		
	close(ID_FILE_ADV_PLUME)

1104 continue

	plume_const%time_step = min(plume_const%time_step,real(wind_time_step))	
	n_plume_steps = ceiling(real(wind_time_step) / real(plume_const%time_step))
	! Recompute to make sure they are consisten
	plume_const%time_step = real(wind_time_step) / real(n_plume_steps)
	print*,'Plume time step [s]: ', plume_const%time_step
	

	minw_plume = -1
	idx_minw_plume = -1
	
	return

1102 write(msgoutfile,*)'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	end
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadAdvancedBldgFile

	use constants
	use file_handling_module
	use canopy_module
	use flags_module

	implicit none

	character(STR_LEN) :: fname

	fname = 'QFire_Bldg_Advanced_User_Inputs.inp'
	open(ID_FILE_BLDG_ADVANCED, file = TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
		status='old',err = 1101)

	read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) flag%bld_to_fuel
	if(flag%bld_to_fuel /= 0 .and. flag%bld_to_fuel /= 1) then
		write (msgoutfile,*) 'Invalid building-to-fuel flag. Must be [0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) canopy%BLDG_FUEL_DENS
	if(canopy%BLDG_FUEL_DENS <= 0) then
		write (msgoutfile,*) 'Invalid building-to-fuel density. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) canopy%BLDG_ATTEN_COEFF
	if(canopy%BLDG_ATTEN_COEFF <= 0) then
		write (msgoutfile,*) 'Invalid building-to-fuel attenuation coefficient. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) canopy%BLDG_CANOPY_Z0
	if(canopy%BLDG_CANOPY_Z0 <= 0) then
		write (msgoutfile,*) 'Invalid building-to-fuel z0. Must be > 0.'
		call TerminateProgram()
	endif
	
	read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) flag%fuel_to_canopy
	if(flag%fuel_to_canopy /= 0 .and. flag%fuel_to_canopy /= 1) then
		write (msgoutfile,*) 'Invalid fuel-to-canopy flag. Must be [0 1].'
		call TerminateProgram()
	endif

	if(flag%fuel_to_canopy == 1)then
		read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103)canopy%UPDATE_WINDS
		if(canopy%UPDATE_WINDS /= 0 .and. canopy%UPDATE_WINDS /= 1) then
			write (msgoutfile,*) 'Invalid flag to update canopy winds. Must be [0 1].'
			call TerminateProgram()
		endif

		read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) canopy%FUEL_ATTEN_COEFF
		if(canopy%FUEL_ATTEN_COEFF <= 0) then
			write (msgoutfile,*) 'Invalid fuel attenuation coefficient. Must be > 0.'
			call TerminateProgram()
		endif

		read(ID_FILE_BLDG_ADVANCED,*,err=1102,end=1103) canopy%FUEL_CANOPY_Z0
		if(canopy%FUEL_CANOPY_Z0 <= 0) then
			write (msgoutfile,*) 'Invalid fuel z0. Must be > 0.'
			call TerminateProgram()
		endif
	else
		canopy%UPDATE_WINDS = 0
	endif


	close(ID_FILE_BLDG_ADVANCED)

	return

1101 print*,'Error in opening file '//TRIM(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE ReadFiretecFilePath(fname)

	use file_handling_module
	use grid_module	
	use constants
	
	implicit none

	character(STR_LEN),intent(IN) :: fname	

	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102) ! header
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)ft%file_path
	read(ID_FILE_QFIRE_INP,*,err=1103,end=1102)ft%fuel_file_type
	if(ft%fuel_file_type /= 1 .and. ft%fuel_file_type /= 2) then
		write (msgoutfile,*)'Invalid fuel flag. Accepted values [1 2].'
		call TerminateProgram()
	endif

	return

1102 write(msgoutfile,*)'End of file encountered while reading file '//TRIM(fname)//'.'
	call TerminateProgram()

1103 print*,'Error in reading file '//TRIM(fname)//'.'
	call TerminateProgram()

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_simparams(nx1ny1,nx1,ny1,nz1)
	
	use file_handling_module
	use sor_module
	use bld_module
	use diffusion_module
	use damage_module
	use grid_module
	use time_module
	
	implicit none
	
	integer :: stretchgridflag, i, j, k, nx1ny1, nx1, ny1, nz1, qcfd_flag

	open(unit=ID_FILE_QU_SIMPARAMS,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QU_simparams.inp',status='old')

	! read input data from QU_domain.inp
	read(ID_FILE_QU_SIMPARAMS,*) ! QUIC version header line
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%nx      !nx defined in input file
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%ny      !ny defined in input file
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%nz      !nz defined in input file
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%dx      !dx defined in input file
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%dy      !dy defined in input file
	
	qugrid%dxy=min(qugrid%dx,qugrid%dy)
	qugrid%dxi=1./qugrid%dx
	qugrid%dyi=1./qugrid%dy
	! man 1/14/05 account for difference in grid cell definitions
	qugrid%nx=qugrid%nx+1
	qugrid%ny=qugrid%ny+1
	qugrid%nz=qugrid%nz+2
	nx1=qugrid%nx-1
	ny1=qugrid%ny-1
	nz1=qugrid%nz-1
	nx1ny1 = nx1*ny1
	! calculate domain Length width and height erp 1/30/2003
	qugrid%Lx = real(qugrid%nx-1)*qugrid%dx
	qugrid%Ly = real(qugrid%ny-1)*qugrid%dy
	qugrid%one_over_nxnynz = 1./real((qugrid%nx-1)*(qugrid%ny-1)*(qugrid%nz-1))

	! end man 1/14/05 account for difference in grid cell definitions
	! MAN 07/25/2008 stretched vertical grid
	allocate(qugrid%z(qugrid%nz),qugrid%zm(qugrid%nz),	&
		qugrid%dz_array(qugrid%nz),qugrid%dzmi(qugrid%nz))
	read(ID_FILE_QU_SIMPARAMS,*)stretchgridflag   !Stretched grid flag (0= dz constant with z)
	qugrid%z=0.
	qugrid%zm=0.
	qugrid%dz_array=0.
	select case(stretchgridflag)
	case(0) !uniform
		read(ID_FILE_QU_SIMPARAMS,*)qugrid%dz      !dz defined in input file
		qugrid%dz_array = qugrid%dz
	case(1) !custom
		read(ID_FILE_QU_SIMPARAMS,*)
		do k=2,qugrid%nz-1
			read(ID_FILE_QU_SIMPARAMS,*)qugrid%dz_array(k)
		enddo
	case(2,3,4) !parabolic or exponential
		read(ID_FILE_QU_SIMPARAMS,*)
		read(ID_FILE_QU_SIMPARAMS,*)
		read(ID_FILE_QU_SIMPARAMS,*)
		do k=2,qugrid%nz-1
			read(ID_FILE_QU_SIMPARAMS,*)qugrid%dz_array(k)
		enddo
	endselect
	qugrid%dz_array(1) = qugrid%dz_array(2)
	qugrid%dz_array(qugrid%nz) = qugrid%dz_array(qugrid%nz-1)
	qugrid%zm(1) = -0.5 * qugrid%dz_array(1)
	qugrid%z(1) = 0.0
	do k=2,qugrid%nz
		qugrid%z(k) = qugrid%z(k-1) + qugrid%dz_array(k)
		qugrid%zm(k) = qugrid%z(k) - 0.5*qugrid%dz_array(k)
	enddo
	do k=2,qugrid%nz-1
		qugrid%dzmi(k) = 2./(qugrid%dz_array(k)+qugrid%dz_array(k-1))
	enddo

	qugrid%dz = minval(qugrid%dz_array)
	! MAN 07/25/2008 stretched vertical grid
	qugrid%Lz = qugrid%z(qugrid%nz-1)
	! MAN 09/02/2008 time information	
	read(ID_FILE_QU_SIMPARAMS,*)qutime%nsteps    !total time increments
	allocate(qutime%unix(qutime%nsteps))	
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%utc_offset ! UTC conversion
	read(ID_FILE_QU_SIMPARAMS,*) ! header line
	do i = 1,qutime%nsteps
		read(ID_FILE_QU_SIMPARAMS,*)qutime%unix(i)
	enddo
	! building parameterization flags
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%roof   ! rooftop recirc flag
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%upwind ! upwind cavity flag
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%streetcanyon    ! street canyon initialization method, added PKK 05/12/03
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%intersection    ! MAN 7/11/2006
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%wake ! MAN 06/29/2007 added wake flag to QU_simparams.inp
	read(ID_FILE_QU_SIMPARAMS,*)bld%flag%sidewall ! MAN 04/24/2013 added sidewall flag to QU_simparams.inp
	read(ID_FILE_QU_SIMPARAMS,*)
	read(ID_FILE_QU_SIMPARAMS,*)
	
	! MAN 7/10/2006 convergence criteria
	read(ID_FILE_QU_SIMPARAMS,*)sor%itermax ! max number of iterations
	read(ID_FILE_QU_SIMPARAMS,*)sor%residual_reduction     ! MAN 09/26/2006 added residual reduction to input
	! AAG 08/25/2006 turbulent diffusion parameters
	read(ID_FILE_QU_SIMPARAMS,*)diff%flag   ! turns on diffusion
	read(ID_FILE_QU_SIMPARAMS,*)diff%step      ! diffusion iterations
	! MAN 02/05/2007 Geo-referencing parameters
	!!! Scot - look here
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%domain_rotation
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%utmx
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%utmy
	read(ID_FILE_QU_SIMPARAMS,*)qugrid%utmzone
	read(ID_FILE_QU_SIMPARAMS,*)
	! MAN 05-24-2007 QUIC-CFD flag turns off the SOR solver
	read(ID_FILE_QU_SIMPARAMS,*)qcfd_flag
	read(ID_FILE_QU_SIMPARAMS,*)damage%flag
	
	close(ID_FILE_QU_SIMPARAMS)

	if(diff%flag > 0) &
		allocate(												&
			diff%Fxd(qugrid%nx,qugrid%ny,qugrid%nz),	&
			diff%Fyd(qugrid%nx,qugrid%ny,qugrid%nz),	&
			diff%Fzd(qugrid%nx,qugrid%ny,qugrid%nz),	&
			diff%visc(qugrid%nx,qugrid%ny,qugrid%nz))
			
	allocate(													&
		qugrid%xcenters(qugrid%nx),					&
		qugrid%ycenters(qugrid%ny),					&
		qugrid%xedge(qugrid%nx),						&
		qugrid%yedge(qugrid%ny))
	!$omp parallel do private(i)
	do i=1,qugrid%nx
		qugrid%xcenters(i)=(real(i)-0.5)*qugrid%dx
		qugrid%xedge(i) = real(i-1)*qugrid%dx
	enddo
	!$omp end parallel do
	!$omp parallel do private(j)
	do j=1,qugrid%ny
		qugrid%ycenters(j)=(real(j)-0.5)*qugrid%dy
		qugrid%yedge(j) = real(j-1)*qugrid%dy
	enddo
	!$omp end parallel do
	
	if(qcfd_flag .gt. 0)then
		sor%itermax = 0
		diff%flag = 0
	endif

	allocate(bld%uProfile(qugrid%nz),bld%vProfile(qugrid%nz),bld%speedProfile(qugrid%nz))
	
		
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_bld(nx1, nx1ny1)
	
	use file_handling_module	
	use constants
	use bld_module
	use interpolation_module
	use grid_module
	use flags_module
	
	implicit none
	
	integer, intent(IN) :: nx1, nx1ny1
	integer :: i,j,k,numnodes,idx,inumpolygon
	
	open(unit=ID_FILE_QU_BLD,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QU_buildings.inp',status='old')

	! input from QU_buildings.inp
	read(ID_FILE_QU_BLD,*) ! QUIC version header line
	read(ID_FILE_QU_BLD,*)bld%zo			 ! MAN 8-19-2005 Updated input output file structures
	read(ID_FILE_QU_BLD,*)bld%number  ! number of buildings
	read(ID_FILE_QU_BLD,*)inumpolygon ! number of polygon nodes

	if(bld%number > 0) then
		allocate(											&
			bld%Ht(bld%number),							&
			bld%Wti(bld%number),							&
			bld%Lti(bld%number),							&
			bld%num(bld%number),							&
			bld%btype(bld%number),						&
			bld%group_id(bld%number),					&
			bld%invnum(bld%number),						&
			bld%atten(bld%number),						&
			bld%rooftop_flag(bld%number),				&
			bld%zfo_actual(bld%number),				&
			bld%geometry(bld%number),					&
			bld%numpolygons(bld%number),				&
			bld%startidx(bld%number),					&
			bld%stopidx(bld%number),					&
			bld%roof(bld%number),						&
			bld%wall(bld%number),						&
			bld%cx(bld%number),							&
			bld%cy(bld%number),							&
			bld%xfo(bld%number),							&
			bld%yfo(bld%number),							&
			bld%zfo(bld%number),							&
			bld%gamma(bld%number),						&
			bld%aa(bld%number),							&
			bld%bb(bld%number),							&		
			bld%damage(bld%number))
		
		allocate(											&
			bld%x(inumpolygon),							&
			bld%y(inumpolygon),							&
			bld%LrNode(inumpolygon),					&
			bld%LrFace(inumpolygon),					&
			bld%FaceRelWindDir(inumpolygon),				&
			bld%FaceWakeDir(inumpolygon))
	
		do i=1,bld%number
			read(ID_FILE_QU_BLD,*)
			bld%num(i)=i
			read(ID_FILE_QU_BLD,*)bld%group_id(i)
			read(ID_FILE_QU_BLD,*)bld%geometry(i)
			read(ID_FILE_QU_BLD,*)bld%btype(i)
			if(bld%btype(i) .eq. 2)read(ID_FILE_QU_BLD,*)bld%atten(i)
			if(bld%geometry(i) .eq. 4 .or. bld%geometry(i) .eq. 5)then
				read(ID_FILE_QU_BLD,*)bld%wall(i)
				read(ID_FILE_QU_BLD,*)bld%roof(i)
			endif
			read(ID_FILE_QU_BLD,*)bld%Ht(i)
			read(ID_FILE_QU_BLD,*)bld%zfo(i)
			read(ID_FILE_QU_BLD,*)bld%cx(i)
			read(ID_FILE_QU_BLD,*)bld%cy(i)
			if(bld%geometry(i) .eq. 6)then
				read(ID_FILE_QU_BLD,*)bld%numpolygons(i)
				do j=1,bld%numpolygons(i)
					read(ID_FILE_QU_BLD,*)
					read(ID_FILE_QU_BLD,*)numnodes
					read(ID_FILE_QU_BLD,*)
					if(j .eq. 1)bld%startidx(i)=idx+1
					do k=1,numnodes
						idx=idx+1
						read(ID_FILE_QU_BLD,*)bld%x(idx),bld%y(idx)
					enddo
					read(ID_FILE_QU_BLD,*)
				enddo
				bld%stopidx(i)=idx
			else
				read(ID_FILE_QU_BLD,*)bld%xfo(i)
				read(ID_FILE_QU_BLD,*)bld%yfo(i)
				read(ID_FILE_QU_BLD,*)bld%Lti(i)
				read(ID_FILE_QU_BLD,*)bld%Wti(i)
				read(ID_FILE_QU_BLD,*)bld%gamma(i)
				if(bld%geometry(i) .eq. 2 .or. bld%geometry(i) .eq. 5)then    !if the building is a cylinder/ellipse
					bld%bb(i)=bld%Wti(i)/2.         !set minor axis to input Width
					bld%aa(i)=bld%Lti(i)/2.         !set major axis to input Lenth
				endif
				if(bld%geometry(i) .eq. 3)then      ! if the building is a Pentagon
					bld%bb(i)=bld%Wti(i)/2.                 ! Radius Pentagon is inscribed in
					bld%xfo(i)=bld%xfo(i)-bld%bb(i)
				endif
			endif
			read(ID_FILE_QU_BLD,*)
			bld%Ht(i)=bld%Ht(i)+bld%zfo(i)
			bld%gamma(i)=bld%gamma(i)*pi/180.  !erp 7/25/03
		enddo
		
		allocate(bld%ufarwake(qugrid%nx,qugrid%ny),bld%vfarwake(qugrid%nx,qugrid%ny))
		
		bld%rooftop_flag = 0 ! AAG 09/13/06  intialized rooftop flag to 1
		bld%damage = 0 ! MAN 08/17/2009
		bld%roof = 0
		bld%Wti = 0.
		bld%Lti = 0.
		
		call Read_QU_bldarray(nx1, nx1ny1)
		call Read_QU_bldflag(nx1, nx1ny1)		
	endif

	close(ID_FILE_QU_BLD)

	allocate(																&
		bld%icellflag(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),	&
		bld%flag%isbld(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
	if(flag%isfire == 0) then
		allocate(bld%icellwake(qugrid%nx-1,qugrid%ny-1))
	endif
	
	if(bld%number == 0) then
		bld%icellflag = 1
		bld%icellflag(:,:,1) = 0
		bld%flag%isbld = 0
	endif
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_fileopt()
	
	use file_handling_module	
	use flags_module
	
	implicit none
	
	open(unit=ID_FILE_QU_FILEOPT,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QU_fileoptions.inp',status='old')

	! read file read writing option data from QU_options.dat
	read(ID_FILE_QU_FILEOPT,*) ! QUIC version header line
	read(ID_FILE_QU_FILEOPT,*)flag%frm				!output format flag 1=ascii,2=binary,3=both
	read(ID_FILE_QU_FILEOPT,*)flag%uofield			!write out uofield.dat 1=yes, 0=no
	read(ID_FILE_QU_FILEOPT,*)flag%uosensor		!write out flag for sensor velocities 1=yes, 0=no
	read(ID_FILE_QU_FILEOPT,*)flag%staggered		!write out flag for staggered velocities 1=yes, 0=no

	close(ID_FILE_QU_FILEOPT)

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE OpenOutputFiles()
	
	use file_handling_module
	use flags_module
	
	implicit none
	
	open(unit=ID_FILE_QU_BUILDOUT,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QP_buildout.inp',status='unknown')
	
	if(flag%errorWrite .gt. 0)then
		open(unit=ID_FILE_QU_ERROR,file=TRIM(workingDirectory)// &
			TRIM(fileSeparator)//'QU_error.dat',status='unknown')
	endif
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE InitWindVariables()
		
	use grid_module
	use winds_module
	
	implicit none
	
	allocate(																	&
		quwinds%uo(qugrid%nx,qugrid%ny,qugrid%nz),					&
		quwinds%vo(qugrid%nx,qugrid%ny,qugrid%nz),					&
		quwinds%wo(qugrid%nx,qugrid%ny,qugrid%nz),					&
		quwinds%u(qugrid%nx,qugrid%ny,qugrid%nz),						&
		quwinds%v(qugrid%nx,qugrid%ny,qugrid%nz),						&
		quwinds%w(qugrid%nx,qugrid%ny,qugrid%nz),						&
		quwinds%undisturbed(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),&
		quwinds%uint(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),			&
		quwinds%vint(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))

	quwinds%uo = 0.
	quwinds%vo = 0.
	quwinds%wo = 0.
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_bldarray(nx1, nx1ny1)
	
	use file_handling_module	
	use bld_module
	use grid_module
	
	implicit none
	
	integer, intent(IN) :: nx1, nx1ny1
	integer :: i,j,k,idx
	integer, allocatable:: bin_int_read(:)
	
	allocate(bin_int_read((qugrid%nx-1)*(qugrid%ny-1)*(qugrid%nz-1)))
	
	open(unit=ID_FILE_QU_BLDARRAY,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//"QU_buildarray.inp",form="unformatted",status="old")
	read(ID_FILE_QU_BLDARRAY)bin_int_read
	idx=0
	!$omp parallel do private(i,j,k,idx)
	do k=1,qugrid%nz-1
		do j=1,qugrid%ny-1
			do i=1,qugrid%nx-1
				idx=nx1ny1*(k-1)+nx1*(j-1)+i
				bld%icellflag(i,j,k)=bin_int_read(idx)
			enddo
		enddo
	enddo
	!$omp end parallel do
	close(ID_FILE_QU_BLDARRAY)
	deallocate(bin_int_read)
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_bldflag(nx1, nx1ny1)
	
	use file_handling_module
	USE bld_module
	use grid_module
	
	implicit none
	
	integer, intent(IN) :: nx1, nx1ny1
	integer :: i,j,k,idx
	integer, dimension(:), allocatable :: bin_int_read

	allocate(bin_int_read((qugrid%nx-1)*(qugrid%ny-1)*(qugrid%nz-1)))
	
	open(unit=ID_FILE_QU_BLDFLAG,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//"QU_buildflag.inp",form="unformatted",status="old")
	read(ID_FILE_QU_BLDFLAG)bin_int_read
	idx=0
	!$omp parallel do private(i,j,k,idx)
	do k=1,qugrid%nz-1
		do j=1,qugrid%ny-1
			do i=1,qugrid%nx-1
				idx=nx1ny1*(k-1)+nx1*(j-1)+i
				bld%flag%isbld(i,j,k)=bin_int_read(idx)
			enddo
		enddo
	enddo
	!$omp end parallel do
	close(ID_FILE_QU_BLDFLAG)
	deallocate(bin_int_read)
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_landuse()
	
	use file_handling_module
	use canopy_module
	use landuse_module
	use grid_module
	
	implicit none
	
	integer :: i,j,idx
	integer landuse_params(3)
	integer, allocatable, dimension(:) :: bin_int_read
	real, allocatable, dimension(:) :: bin_real_read
	
	open(unit=ID_FILE_QU_LANDUSE,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QU_landuse.inp',form='unformatted',status='old')
		
	read(ID_FILE_QU_LANDUSE)landuse_params
	landuse%flag = landuse_params(1)
	landuse%veg_flag = landuse_params(2)
	landuse%urb_flag = landuse_params(3)
	if(landuse%flag .eq. 1)then
		allocate(landuse%val(qugrid%nx-1,qugrid%ny-1),bin_int_read((qugrid%nx-1)*(qugrid%ny-1)))
		read(ID_FILE_QU_LANDUSE)bin_int_read
		idx = 1
		do i=1,qugrid%nx-1
			do j=1,qugrid%ny-1
				landuse%val(i,j)=bin_int_read(idx)
				idx=idx+1
			enddo
		enddo
		deallocate(bin_int_read)
		if(landuse%veg_flag .eq. 1 .or. landuse%urb_flag .eq. 1)then
			landuse%canopy_flag=1
			allocate(landuse%height(qugrid%nx-1,qugrid%ny-1),	&
				landuse%atten(qugrid%nx-1,qugrid%ny-1),bin_real_read((qugrid%nx-1)*(qugrid%ny-1)))
			read(ID_FILE_QU_LANDUSE)bin_real_read
			idx = 1
			do i=1,qugrid%nx-1
				do j=1,qugrid%ny-1
					landuse%height(i,j)=bin_real_read(idx)
					idx=idx+1
				enddo
			enddo
			read(ID_FILE_QU_LANDUSE)bin_real_read
			idx = 1
			do i=1,qugrid%nx-1
				do j=1,qugrid%ny-1
					landuse%atten(i,j)=bin_real_read(idx)
					idx=idx+1
				enddo
			enddo
			deallocate(bin_real_read)
			if(canopy%number == 0)then
				allocate(											&
					canopy%ktop(qugrid%nx-1,qugrid%ny-1),	&
					canopy%top(qugrid%nx-1,qugrid%ny-1),	&
					canopy%zo(qugrid%nx-1,qugrid%ny-1),		&
					canopy%d(qugrid%nx-1,qugrid%ny-1))
			endif
		endif
	else
		landuse%veg_flag=0
		landuse%urb_flag=0
		landuse%canopy_flag=0
	endif
	close(ID_FILE_QU_LANDUSE)

	
	if(canopy%number .gt. 0 .or. landuse%canopy_flag .gt. 0)then
		canopy%flag = 1
	else
		canopy%flag = 0
	endif

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE Read_QU_parameterRange()
	
	use file_handling_module	
	use bld_module
	
	implicit none
	
	logical :: there

	INQUIRE( FILE=TRIM(workingDirectory)//TRIM(fileSeparator)//"QU_parameterRange.inp", EXIST=THERE )
	IF(THERE)then
		open(unit=ID_FILE_QU_PARAMRANGE,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
			"QU_parameterRange.inp",status="old") ! material or atmos. cell
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%upwindAngleRange
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%sidewallAngleRange
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%upwindfraction
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%downwindfraction
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%canyonwidthfactor
		read(ID_FILE_QU_PARAMRANGE,*)bld%params%vortexHeightFactor
		close(ID_FILE_QU_PARAMRANGE)
	else
		bld%params%upwindAngleRange=50
		bld%params%sidewallAngleRange=10.0
		bld%params%upwindfraction=0.6
		bld%params%downwindfraction=0.4
		bld%params%canyonwidthfactor=1.
		bld%params%vortexHeightFactor=1.5
	endif

	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE InitShips()
	
	use ship_module
	use file_handling_module
	use time_module
	
	implicit none
	
	logical :: there
	integer :: idx, i
	
	allocate(												&
		ship%speed(qutime%nsteps),						&
		ship%bearing(qutime%nsteps),					&
		ship%currentSpeed(qutime%nsteps),			&
		ship%currentDirection(qutime%nsteps))
	ship%movingCoordsFlag = 0
	ship%relativeBearing = 0.0
	ship%speed = 0.0
	ship%bearing = 0.0
	ship%currentSpeed = 0.0
	ship%currentDirection = 0.0

	INQUIRE( FILE=TRIM(workingDirectory)//TRIM(fileSeparator)// &
		"QU_movingcoords.inp", EXIST=THERE )
	IF(THERE)then
		open(unit=ID_FILE_QU_MOVINGCOORD,file=TRIM(workingDirectory)// &
			TRIM(fileSeparator)//"QU_movingcoords.inp",status="old")
		read(ID_FILE_QU_MOVINGCOORD,*) !skip version header
		read(ID_FILE_QU_MOVINGCOORD,*)ship%movingCoordsFlag
		read(ID_FILE_QU_MOVINGCOORD,*)ship%relativeBearing
		read(ID_FILE_QU_MOVINGCOORD,*) !skip header
		do i = 1,qutime%nsteps
			read(ID_FILE_QU_MOVINGCOORD,*)				&
				idx,												&
				ship%speed(i),									&
				ship%bearing(i),								&
				ship%currentSpeed(i),						&
				ship%currentDirection(i)
		enddo
		close(ID_FILE_QU_MOVINGCOORD)
	endif
	
	END
!===========================================================================================
!===========================================================================================
	SUBROUTINE CheckBuildingsToVegetation()
	
	use canopy_module
	use bld_module
	
	implicit none
	
	integer :: i
	
	canopy%number = 0
	bld%numberneg = 0
	bld%number_garage = 0
	do i = 1,bld%number
		if(bld%btype(i) .eq. 2)then
			canopy%number = canopy%number + 1			! total number of vegative canopies
		endif
		if(bld%btype(i) .eq. 3)then
			bld%number_garage = bld%number_garage+1							! total number of garage buildings
		endif
		if(bld%btype(i) .eq. 0)then
			bld%numberneg=bld%numberneg+1					! total number of negative buildings
		endif
	enddo

	END
!===========================================================================================
!===========================================================================================