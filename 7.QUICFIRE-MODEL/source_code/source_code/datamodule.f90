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
	module constants

	implicit none

	real,parameter ::					&
		vk = 0.4,						& !Von Karman's constant
		GRAVITY = 9.81,				& ! m/s^2
		RHO_AIR = 1.225,				& ! kg/m3, density of air at sea level at 15 degrees celsius
		PI = 4.*atan(1.0)

	integer ::							&
		STDOUT = 6,						&
		DBGFILE = 90,					&
		msgoutfile

	end module constants
	!======================================================================================
	!======================================================================================
	module sor_module
	
	implicit none
	
	TYPE :: BoundaryConditions_Type
		real, dimension(:,:,:), allocatable :: &
			e,f,g,h,m,n
		
	END TYPE BoundaryConditions_Type
	
	TYPE :: sor_type
		TYPE(BoundaryConditions_Type) :: bc
		
		integer ::										&
			itermax,										&
			option,										&
			mc_is,mc_ie,								&
			mc_js,mc_je,								&
			mc_ks,mc_ke, 								& ! N/A, index in the mass consitency grid where to do mass consistency
			mc_buffer								 	  ! N/A, number of cells around the plume where to compute mass consistency

		real ::											&
			alpha1,										&
			alpha2,										&
			eta,											&
			ovalph1,										&
			ovalph2,										&
			alpha1sq,									&
			alpha2sq,									&
			omegarelax,									&
			one_minus_omegarelax,					&
			res_red_factor,							&
			residual_reduction,						&
			ALPHA2_FIRE_VAL,							&
			invomegarelax
			
		real,dimension(:,:,:),allocatable ::	&
			alpha2_fire,								&
			p1,											&
			p2,											&
			denom,										&
			denom0,										&
			r		
	END TYPE sor_type
	
	TYPE(sor_type) :: sor
	
	end module sor_module
	!======================================================================================
	!======================================================================================
	module canopy_module
	
	implicit none
	
	TYPE :: canopy_type
		integer ::											&
			number,											& !N/A, number of vegetative canopies in domain
			flag,												&
			UPDATE_WINDS				  					  ! N/A, 0 = do not call update_canopy, 1 = call update_canopy

		real, dimension(:,:), allocatable ::		&
			top,												&
			zo,												&
			ustar,											&
			d
		real ::												&
			FUEL_ATTEN_COEFF,								&
			FUEL_CANOPY_Z0,								&
			BLDG_ATTEN_COEFF,								&
			BLDG_CANOPY_Z0,								&
			BLDG_FUEL_DENS  								 ! kg/m3, fuel density of buildings when treated like canopy in FireCA

	
		real, dimension(:,:,:), allocatable ::	atten
		integer, dimension(:,:), allocatable :: ktop
	END TYPE canopy_type
		
	TYPE(canopy_type) :: canopy
	
	end module canopy_module
	!======================================================================================
	!======================================================================================
	MODULE ship_module
	
	implicit none
	
	TYPE :: ShipType
		real, dimension(:), allocatable ::		&
			speed,										&
			bearing,										&
			currentSpeed,								&
			currentDirection
		real ::											&
			relativeBearing,							&
			shipU,										&
			shipV,										&
			currentU,									&
			currentV	
		integer ::movingCoordsFlag
	END TYPE ShipType
	
	TYPE(ShipType) :: ship
	
	END MODULE ship_module
	!======================================================================================
	!======================================================================================
	MODULE wind_profile_module
	
	implicit none
	
	TYPE :: site_type
		real, allocatable, dimension(:) ::					&
			uoint, voint,											& !Interpolated U and V velocities
			xcoord, ycoord											 ! x,y coordinates of each sensor site
		real, allocatable, dimension(:,:,:) ::				&
			z_data,													&
			ws_data,													& ! wind speed
			wd_data,													& ! wind direction
			u_data,													&
			v_data,													&
			wm,wms					  								  !the "Weight" of each site on each grid point
	integer, allocatable, dimension(:,:) ::				&
		nz_data,														&
		blayer_flag													  !number of points in the z-direction and the boundary layer flag
	real, allocatable, dimension(:,:) ::					&
		pp,															& ! exp/zo for each site
		H,																&
		ac,															&
		rL,															&
		u_prof, v_prof												  !U and V vertical velocity profiles
	END TYPE site_type

	TYPE :: windProfile_type
		TYPE(site_type) :: site
		real, allocatable, dimension(:) ::					&
			ws_data,													&
			wd_data,													&
			z_data
		integer ::													&
			num_sites,												& !number of sensor sites
			num_vert_points										  !MAN 4/5/2007 number of data points to allocate in data points profiles
	END TYPE windProfile_type
	
	TYPE(windProfile_type) :: windProfile	
	
	END MODULE wind_profile_module
	!======================================================================================
	!======================================================================================
	MODULE flags_module
	
	implicit none
	
	TYPE :: flag_type
		integer ::										&
			errorWrite,									&
			frm,											& ! 1 = ASCII only, 2 = binary only, 3 = both
			uofield,										&
			uosensor,									&
			staggered,									& ! write out flag for staggered velocities 1=yes, 0=no
			emissions_file,							& ! N/A, 0 = no, 1 = yes, print emissions in file
			thermalrad_file,							& ! N/A, 0 = no, 1 = yes, print thermal radiation in file
			en2atm_file,								& ! N/A, 0 = no, 1 = yes, print energy-to-atmosphere in file
			reactrate_file,							& ! N/A, 0 = no, 1 = yes, print reaction rate in file
			fuelmass_file,								& ! N/A, 0 = no, 1 = yes, print fuel mass in file
			qfwinds_file,								& ! N/A, 0 = no, 1 = yes, print FIRECA winds in file
			quwinds_inst_file,						& ! N/A, 0 = no, 1 = yes, print QUIC-URB winds in file
			quwinds_ave_file,							& ! N/A, 0 = no, 1 = yes, print QUIC-URB winds in file
			plume_traj_file,							& ! N/A, 0 = no, 1 = yes, print plume trajectory in file
			moisture_file,								& ! N/A, 0 = no, 1 = yes, print fuel moisture in file
			int_massburnt_file,						& ! N/A, 0 = no, 1 = yes, print vertically-integrated % of mass burnt in file
			plume_loc_file,							& ! N/A, 0 = no, 1 = yes, print file with locations of the plumes
			num_qu_updates,							& ! N/A, after how many quic-urb updates, the code prints out the output files
			num_qu_updates_avewinds,				& ! N/A, after how many quic-urb updates, the code prints out the averaged output files
			num_fireca_updates,						& ! N/A, after how many fireca time steps, the code prints out the output files
			num_emission_and_rad_updates, 		& ! N/A, after how many fireca time steps, the code prints out the emission and radiation output files
			isfire,										& ! N/A, 0 = no fire, 1 = yes fire
			fuel_to_canopy,							& ! N/A, convert fuel to canopy for fire
			bld_to_fuel,								& ! N/A, convert buildings to fuel for fire
			domain_initialized,						& 			
			output_initialized

	END TYPE flag_type
	
	TYPE(flag_type) :: flag
	
	END MODULE flags_module
	!======================================================================================
	!======================================================================================
	MODULE topo_module
		
	use file_handling_module
	
	implicit none
	
	TYPE :: topo_type
		integer :: flag									 ! N/A, 0 = no topo, 1 = read in topo
		character(STR_LEN) :: file_path
	END TYPE topo_type
	
	TYPE(topo_type) :: topo
		
	END MODULE topo_module
	!======================================================================================
	!======================================================================================	
	MODULE damage_module
	
	implicit none
	
	TYPE :: damage_type
		
		integer :: flag
		real ::												&
			explosionx,explosiony,						&
			hemass,											&
			Rdestroyed,										&
			Rdamaged,										&
			Rbuild
			
	END TYPE damage_type
	
	TYPE(damage_type) :: damage
	
	END MODULE damage_module
	!======================================================================================
	!=====================================================================================+
	MODULE bld_module
	
	implicit none
	
	TYPE :: flag_type
		integer, allocatable, dimension(:,:,:) :: isbld
		
		integer ::											&
			roof,												&
			upwind,											&
			intersection,									&
			wake,												&
			sidewall,										&
			streetcanyon,									&
			blending
	END TYPE flag_type
	
	TYPE :: params_type
		real ::												&
			upwindAngleRange,								&
			sidewallAngleRange,							&
			upwindfraction,								&
			downwindfraction,								&
			canyonwidthfactor,							&
			vortexHeightFactor
		
	END TYPE params_type
	
	TYPE :: bld_type
		TYPE(flag_type) :: flag
		TYPE(params_type) :: params
		
		integer ::											&
			number,											&
			numberneg,										&
			number_garage
		
		real :: zo
		integer, dimension(:), allocatable ::		&
			btype,											&
			geometry,										&
			invnum,											&
			istart, jstart, kstart,						&
			iend, jend, kend,								&
			startidx,										&
			stopidx,											&
			numpolygons,									&
			roof,												&
			num,												&
			group_id,										&
			rooftop_flag,									&
			damage
		
		real, dimension(:), allocatable ::			&
			Lf,												&
			Lr,												&
			Lt,												&
			Weff,												&
			Leff,												&
			Rscale,											&
			Rcx,												&
			Wt,												&
			x, y,												&
			cx,cy,											&
			wall,												&
			Ht,												&
			Wti,												&
			Lti,												&
			aa, bb,											&
			xfo,yfo,zfo,									&
			gamma,											&
			atten,											&
			zfo_actual,										&
			LrNode,											&
			LrFace,											&
			FaceWakeDir,									&
			FaceRelWindDir,								&
			uProfile,										&
			vProfile,										&
			speedProfile
		integer, allocatable, dimension(:,:) :: 	&
			icellwake,										&
			ufarwake,										&
			vfarwake
		integer, allocatable, dimension(:,:,:) :: icellflag		
		
	END TYPE bld_type
	
	TYPE(bld_type) :: bld
	
	END MODULE bld_module
	!======================================================================================
	!======================================================================================
	MODULE diffusion_module
	
	implicit none
	
	TYPE :: diffusion_type
		real, allocatable, dimension(:,:,:) ::					&
			visc,															&
			Fxd, Fyd, Fzd
		
		integer :: flag,step
	END TYPE diffusion_type
	
	TYPE(diffusion_type) :: diff
	
	END MODULE diffusion_module
	!======================================================================================
	!======================================================================================
	MODULE landuse_module
	
	implicit none
	
	TYPE :: landuse_type
		integer ::						&
			flag,							&
			veg_flag,					&
			urb_flag,					&
			canopy_flag
		integer, allocatable, dimension(:,:) :: val
		real, allocatable, dimension(:,:) :: height,atten
	END TYPE landuse_type
	
	TYPE(landuse_type) :: landuse
	
	END MODULE landuse_module
	!======================================================================================
	!======================================================================================
	module winds_module

	implicit none

	TYPE :: winds_type
		real, allocatable, dimension(:,:,:) :: 		&
			uint, vint,											&
			undisturbed,										&
			uo, vo, wo,											&
			u, v, w,												&
			u_ave, v_ave, w_ave,								&
			wplume,												&
			uo_roof, vo_roof,									&
			u_mc, v_mc, w_mc									  ! m/s, QUIC-URB wind fields	

		real, allocatable, dimension(:) :: 				&
			uo_before_params, vo_before_params,			&
			sigma
		
		real ::													&
			max_velmag											 ! m/s
			
	END TYPE winds_type

	TYPE(winds_type) :: quwinds, fcawinds

	end module winds_module
	!======================================================================================
	!======================================================================================
	module ignitions_module
	
	implicit none
	
	TYPE :: ign_type
		real, allocatable,dimension(:) ::			&
			x,y,z,											& ! m, x,y,z coordinates of the ignition pattern
			time,												& ! s, time of the ignition for the ignition pattern
			radius											  ! m, radius over which to distribute ignitions for the ignition pattern
		real :: 												&
			xo,yo,											& ! m, fire source south-west corner in the QUIC-URB domain
			xlen,	ylen,										& ! m, fire source dimensions in the x and y directions
			xwidth,ywidth 									  ! m, for the rings, this is the width of the ring itself

		integer, allocatable,dimension(:) ::		&
			new_num											  ! N/A, number of new ignition for an ignition pattern
		integer ::											&
			num_ign_percell,								& ! N/A, number of ignitions per cell
			npoints,											& ! N/A, number of waypoints in the ignition pattern
			flag,												& ! N/A, 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp)
			unix_time										  ! s, when is the fire ignited in Unix Epoch Time
					
	END TYPE ign_type
	
	TYPE(ign_type) :: fire_ignition
	
	END module ignitions_module
	!======================================================================================
	!======================================================================================
	module rnd_num_vars
	
	use mersenne_twister

	implicit none

	!! RANDOM NUMBER VARIABLES
	integer :: envErr,numThreads,iseed
	character*40 :: head1
	type (random_stat), allocatable :: stat(:) 	! state variables for mersenne PRNG
	real,dimension(1) :: harvest
	real,dimension(2) :: harvest2
	real,dimension(3) :: harvest3
	integer,dimension(1) :: iharvest
	integer :: rnd_opt									 ! N/A, 0 = random number seed is time/date, > 0 = random number seed

	end module rnd_num_vars
	!======================================================================================
	!======================================================================================
	module time_module

	implicit none

	TYPE :: time_type
		integer :: 								&
			firesteps_to_update,				& ! N/A, after how many fireca timesteps, the plumes are calculated
			current_iter,						& ! N/A, iteration count
			len_simulation,					& ! s, length of the simulation
			current_int,						& ! s, current time
			dt_int,								& ! s, time step
			nsteps 								  ! N/A, number of time steps

		real ::									&
			current_real,						& ! s, current time
			dt_real								  ! s, time step

		integer, allocatable, dimension(:) ::	&
			unix, 										& ! s, time values in unix time
			local

	END TYPE time_type

	TYPE(time_type) :: qutime, fcatime

	end module time_module
	!======================================================================================
	!======================================================================================
	MODULE interpolation_module

	implicit none

	integer, allocatable,dimension(:) ::		&		
		en2atm_2_qu_kstart,en2atm_2_qu_kend,	& ! N/A, indexes in the QUIC-URB grid where there is fuel		
		fca_2_qu_kstart,fca_2_qu_kend,			& ! N/A, indexes in the QUIC-URB grid where there is fuel		
		qu_2fca_kuv,qu_2fca_w,						& ! N/A, indexes to interpolate QU winds to fireca winds
		kmap_start,kmap_end							  ! N/A, in which QU cells is a certain FireCA cell (vertical dir)

	integer ::											&
			istart_interp,iend_interp,				&
			jstart_interp,jend_interp,				&
			xratio_int, yratio_int, 				& ! N/A, horizontal grid spacing ratio with QU
			matching_grids								  ! N/A, 0 = no, 1 = yes, the quic-urb and fireca grids use the same horizontal cell size	

	real :: xratio_real, yratio_real

	TYPE :: FT_GridConversion
		real, allocatable,dimension(:) :: length			! m, part of the Firetec cell that is in a FireCA cell
		integer, allocatable,dimension(:) :: index		! N/A, index of the Firetec cell that is in a FireCA cell
		integer :: nelem											! N/A, how many  Firetec cell that are in a FireCA cell
	END TYPE FT_GridConversion

	TYPE(FT_GridConversion),dimension(:),allocatable ::	&
		ft_2_fca_conv_x,												&
		ft_2_fca_conv_y,												&
		ft_2_fca_conv_z

	TYPE :: EnRatio
		real,allocatable,dimension(:) :: volRatio ! [0 1], ratio of the fireca cell that is in a quic-urb cell
	END TYPE EnRatio
	TYPE(EnRatio),allocatable,dimension(:) :: fca2quic

	TYPE :: GridConversion
		real, allocatable,dimension(:) :: en_ratio
		integer, allocatable,dimension(:) :: index
		integer :: nelem
	END TYPE GridConversion

	TYPE(GridConversion),dimension(:,:,:),allocatable :: grid_conv

	END MODULE interpolation_module
	!======================================================================================
	!======================================================================================
	MODULE plume_const_module
		
	implicit none
	
	TYPE :: plume_const_type
				
		integer ::										&
			do_w_fill,									&
			MAX_NUM_PLUMES_TIMESTEP					  ! N/A, maximum number of plume in a time step
		
		real ::											&
			en2atmos_mult,								& ! GRAVITY / (air_temp * cp * RHO_AIR)  ! Fb = g/T *(H/rho/cp)  => Fb = mult * H * volume   (dvol used to convert units)
			max_fraction_traj,						& ! N/A, fraction of plume ip1 trajectory length that a point on ip2 can be for merging
			max_angle,									& ! rad, max angle between plumes to merge (default 30 deg)
			init_w_en_fract,							& ! N/A, fraction of ignition energy that is released to the atmosphere
			zvirt0,									  	& ! m, virtual source height if the fire was on the ground
			BRUNT_VAISALA_FREQ2,						& ! 1/s2, N^
			min_time_step,								& ! s, flag == 0: time step to advance the plumes; flag == 1: min time step
			time_step,									& ! s, time step to advance the plumes
			total_time,									& ! s, total time in a time step
			max_time,									& ! s, for stable conditions, the eqn are valid only during the rise phase
			ALPHA_ENTR,									& 
			BETA_ENTR,									& ! N/A, entrainment coefficient
			BETA_ENTR2,									& ! entrainment coefficient squared
			WC_EXPONENT,								& ! N/A, exponent to elevate w to compute overlap
			INV_WC_EXPONENT,							& ! N/A, exponent to elevate w to compute overlap
			SPEEDS_RATIO,								& ! N/A, if wc < SPEEDS_RATIO*wind_speed, the plume is deleted
			WC_MAX = 10,								& ! m/s, maximum updraft velocity
			WC_MIN = 0.1								  ! m/s, if wc < WC_MIN, the plume is deleted
		
		integer :: time_step_flag					  ! N/A, 0 = fixed, 1 = adaptable
	END TYPE plume_const_type
	
	TYPE(plume_const_type) :: plume_const

	END MODULE plume_const_module
	!======================================================================================
	!======================================================================================
	MODULE plume_module
	
	use omp_lib
	
	implicit none
		
	! --------------------------------------------------
	TYPE :: plume_diagnostics_type
		real ::											&
			max_height,									& ! m, plume maximum height
			max_w											  ! m/s, plume maximum w
	END TYPE plume_diagnostics_type
	TYPE(plume_diagnostics_type) :: plume_diagnostics
	
	! --------------------------------------------------
	
	TYPE :: loopidx_type
		integer,dimension(2) ::						&
			increm,										& ! N/A, direction of the loops (-1 or +1)
			istart,iend,								& ! N/A, start and end index for the loops to initialize plumes (x-dir)
			jstart,jend,								& ! N/A, start and end index for the loops to initialize plumes (y-dir)
			kstart,kend									  ! N/A, start and end index for the loops to initialize plumes (k-dir)
		
		integer :: list_count					! N/A, [1 8] how to run the loops to initialize the plumes
		
		integer,dimension(8,3) :: list  ! N/A, combinations of the loops orders to initialize plumes
	ENDTYPE loopidx_type
	TYPE(loopidx_type) :: loopidx
	! --------------------------------------------------
	integer, allocatable,dimension(:) ::	&
		to_delete									  ! N/A, 0 = no, 1 = yes
	
	integer :: 										&	
		n_plume_steps,								& ! N/A, number of plume time steps for each QU time step
		total_fires,								& ! N/A, current number of plumes
		idx_minw_plume,							& ! N/A, index in the arrays of the plume with the lowest w
		max_id_plume							 	  ! N/A, max id used to identify plumes

	real :: minw_plume							  ! m/s, min value of w
	
	TYPE :: PlumeType
		real, dimension(3) ::					&
			traj_coord,								& ! m, x, y, z coordinates of the trajectory
			traj_vel,								& ! m/s, plume centerline speed components
			traj_coord_prev,						& ! m, x, y, z coordinates of the trajectory at the previous time step
			en, eh,									& ! m, unit vectors
			s											  ! m, trajectory direction vector (not unit)	
		real ::										&
			perimeter,								& ! m, ellipse perimeter
			perimeter_over_rph,					& ! m, ellipse perimeter / rph
			traj_len,								& ! m, length of the trajectory between time t and t+dt
			gamma,									& ! N/A, ellipse semi-axes ratio
			zmax,										& ! m, maximum height that can be reached in stable conditions
			heat_of_fire,							& ! kW, heat released by the fire in a QU cell
			wc_prev,									& ! m/s,	updraft at the previous point on the trajectory
			ws,										& ! m/s, wind speed at the plume centerline
			start_time,								& ! s, when the plume was born
			total_time,								& ! s, plume total time (used for stable conditions)
			radius_h,								& ! m, plume radius
			radius_n,								& ! m, plume radius			
			zvirt,									& ! m, plume virtual source
			zfire_orig,								& ! m, height where the plume originated
			total_z									 ! m, vertical distance covered by the plume because of its buoyancy used to compute the updraft

		integer ::									&
			i,j,k,									& ! N/A, cell of the plume centerline
			i0,j0,k0,								& ! N/A, cell of the plume centerline at the origin
			istart,iend,							& ! N/A, to compute the quwinds%wplume
			jstart,jend,							& ! N/A, to compute the quwinds%wplume
			kstart,kend,							& ! N/A, to compute the quwinds%wplume
			id											  ! N/A, unique plume segment ID for plotting
	END TYPE PlumeType

	TYPE(PlumeType),dimension(:),allocatable :: plume
	
	integer(omp_lock_kind), dimension(:,:,:), allocatable :: plume_lock
	
	END MODULE plume_module
!======================================================================================
!======================================================================================
	MODULE grid_module
	
	use file_handling_module

	implicit none 

	TYPE :: grid_type
		integer ::									&
			nx, ny, nz,								&									
			num_fuel_cells							 ! N/A, number of cells in the domain that contain fuel
		
		real :: 										&			
			domain_rotation,						& ! deg, domain rotation
			dx, dy,									& ! m, horizontal grid spacing
			dz,										& ! m, vertical min grid spacing			
			dxi, dyi,								& ! 1/m, inverse of horizontal grid spacing			
			utmx, utmy,								& ! m, coordinates of southwest corner
			Lx,										& ! m, extent of the fire domain in the x-direction
			Ly,										& ! m, extent of the fire domain in the y-direction
			Lz 										  ! m, top of the fire domain			

		real, allocatable, dimension(:) ::	&
			xcenters, ycenters,					& ! m, cell centers
			dz_array,								& ! m, vertical grid spacing						
			z,											& ! m, cell bottom or top (z(1) = 0)
			zm,										& ! m, cell center (zm(1) > z(1))						
			dzmi										  ! 1/m, inverse of the spacing between the centers of the cells
		
		integer,dimension(:),allocatable :: cell_index			! N/A, index of the cell in the compressed array scheme
		integer,dimension(:,:),allocatable :: ijk_cell_index	! N/A, i,j,k valies of the index of the cell in the com
		
	END TYPE grid_type
	
	TYPE, EXTENDS(grid_type) :: qugrid_type
		integer ::									&
			utmzone,									& ! N/A, utm zone number
			kfire_top_en2atmos,					&
			kfire_top								  ! N/A max cell in the qu-domain that has a bit of a fireca cell in it

		integer*8 :: utc_offset
		
		real :: 										&
			dxy,										& ! m2, dx*dy
			one_over_nxnynz						  ! 1/m3, 1/(nx*ny*nz)
		
			real, allocatable, dimension(:) ::	&	
			xedge, yedge 							  ! m, QUIC-URB cell edges
			
	END TYPE qugrid_type

	TYPE, EXTENDS(grid_type) :: firegrid_type
		integer, dimension(:), allocatable :: idx

		integer ::									&
			num_active_fuel_cells,				&
			nz_en2atmos							     ! N/A, max nz for the arrays num_ignitions_not_absorbed and en_to_atmos
		
		real :: 										&
			avg_cellsize, 	 				  		& ! m, 0.5*(dx+dy)
			cellarea									  ! m2, dx*dy

		real, allocatable, dimension(:) ::	&	
			delta,									& ! m, delta = 4.*curr_cellvol**ONE_THIRD
			dz_array_en2atmos,					& 
			z_en2atmos,								& ! m, cell bottom or top in the extended vertical grid
			zm_en2atmos,							& ! m, cell center in the extended vertical grid
			cellvol,									& ! m3, cell volumes
			cellvol_en2atmos 						  ! m3, cell volumes

		real, dimension(2) :: reciprocal

	END TYPE firegrid_type

	TYPE, EXTENDS(grid_type) :: ft_type
		integer ::									&
			num_ignitions,							& ! N/A, used for the ft domain
			n_fuel_types,							&
			n_grid_top,								&
			fuel_file_type 						  ! N/A, 1 = all fuels together, otherwise = separate files
		
			character(STR_LEN) :: file_path
	END TYPE ft_type

	TYPE(qugrid_type) :: qugrid
	TYPE(firegrid_type) :: firegrid
	TYPE(ft_type) :: ft 

	END MODULE grid_module
!======================================================================================
!======================================================================================
