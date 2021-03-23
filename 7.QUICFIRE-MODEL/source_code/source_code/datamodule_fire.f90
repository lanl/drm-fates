	module fireca_constants_module

	implicit none

	real,parameter ::														&
		MAX_FLAME_HEIGHT = 20,											& ! m, maximum height of the flames to compute num_ignitions_not_absorbed and en_to_atmos
		BURNCENTER_INIT = -100.,										& ! N/A, initial value of fire%burncenter
		ONE_THIRD = 1./3.,												& ! N/A
		ROOTTHREE = sqrt(3.),											& ! N/A
		FOUR_THIRDS = 4./3.,												& ! N/A
		ONE_OVER_ROOTTHREE = 1./ROOTTHREE,							& ! N/A
		GREEK_PI = 4.*atan(1.0),										& ! N/A
		MIN_FUEL_DENSITY = 1e-3,										& ! kg/m3, minimum fuel density, below the fuel is considered completely burned
		MIN_MOISTURE = 1e-3,												& ! kg water/kg fuel, minimum moisture content, below the fuel is considered dry		
		X_O2MAX = 0.21,													& ! %vol, maximum oxygen fraction in air
		X_N2MAX = 0.79,													& ! %vol, maximum nitrogen fraction in air
		SMAG_COEFF = 0.2,													& ! N/A
		PARTICAL_BURN_VEL = 2.,											& ! 2 new value
		HEAT_OF_COMBUSTION = 18620.,									& ! KJ/kg, wood heat of combustion
		ENERGY_PER_IGNITION = 50.0,									& ! kW/ignition, energy of each ignition
		BURNOUT_TIME = 30.,												& ! s, time for fine fuel element to burn
		IGNITION_REACTION_RATE = 0.1,									& ! kg/m3/s, used to initialize the reaction rate where there are ignitions
		HORIZONTAL_VARIATION = 0.5 * GREEK_PI,						& ! N/A
		RHO_FB = 200.,														& ! kg/m^3, material density of fire brands
		DRAG_COEFF = 1.17,												& ! N/A, value for disk assuming angle of attack=90 degrees
		MIN_FB_THICKNESS_AT_BURNOUT = 5e-6,							& ! m, minimum thickness of the firebrand below which it is considered burned out
		FB_MIN_THICKNESS = MIN_FB_THICKNESS_AT_BURNOUT*2.,		& ! m, minimum firebrand thickness to decide to launch it
		LATENT_HEAT_WATER_T0 = 2260.,									& ! KJ/Kg, water latent heat
		PRESSURE = 101325,												& ! Pa
		RGAS = 8.3144598,                                     & ! m3 Pa/K mol
		CP_WATER = 4.184,													& ! KJ/kg/K, water specific heat
		CP_WOOD = 1.7,														& ! KJ/kg/K, wood specific heat
		AMBIENT_TEMPERATURE = 300.,									& ! K, ambient temperature
		INITIAL_O2DENSITY = 0.257,										& ! kg/m3, ambient oxygen concentration  P*MW/(R*T) * 0.21
		RESOLUTION_SWITCH = 10.,										& ! m
		NET_ENERGY_WOOD_KG = HEAT_OF_COMBUSTION - (500.-AMBIENT_TEMPERATURE)*CP_WOOD,			& ! kJ/kg
		LATENT_HEAT_WATER_KG = LATENT_HEAT_WATER_T0 + CP_WATER*(373.-AMBIENT_TEMPERATURE)	  ! kJ/kg

	integer, parameter ::												&
		FB_TIME_DEFAULT = 1e7,											& ! s, time assigned by defaul to the delay time
		FB_NCELLS_DISTRIB = 5,											& ! N/A, cells over which distributing the firebrands total = (FB_NCELLS_DISTRIB*2+1)^2
		dp = kind(1.0d0)

	real :: RADIATON_LOSS_FRACT_COEFF						        ! N/A, coefficient connecting the reaction rate to energy lost by radiation
	! val	= 0.2 + (0.5 - 0.2)/(20 - 2) * (mean(dx,dy) - 2)

	integer :: CELL_SIZE_MODE
	integer,parameter ::													&
		CELL_SIZE_MODE_LARGE = 1,										&
		CELL_SIZE_MODE_SMALL = 0
		
	real ::																	&
		moist_depl_var_init,												&
		init_var_length													  ! m
	
	end module fireca_constants_module
!===========================================================================
!===========================================================================
	module coord_module

	implicit none

	! Defined to avoid temporary arrays when only one cell is passed to a subroutine
	! It was found that avoiding temporary arrays sped up the code a couple of times
	type :: coord
		real,dimension(:),allocatable :: val
		real,dimension(:,:),allocatable :: val2
	end type coord

	end module coord_module
!===========================================================================
!===========================================================================
	module fuel_module
	
	use coord_module
	
	implicit none
	
	type:: fuelgrid
		type(coord),allocatable ::				&
			moist_depl_var(:),					&
			moist_depl_center(:),				& ! N/A, position of the moisture depletion center normalized
														  !     over the cell size [0 1], 2nd index = x,y,z
			depl_center(:),						& ! N/A, position of the mass depletion center normalized over the cell size [0 1], 2nd index = x,y,z
			depl_var(:)

		integer, allocatable, dimension(:) :: &
			is_urban									  ! N/A, 0 = no, 1 = yes
		real, allocatable ::						&
			density(:),								& ! kg/m3, fuel density
			density_initial(:),					& ! kg/m3, fuel density	before ignition
			moisture(:),							& ! kg water/kg fuel, moisture content
			moisture_initial(:),					& ! kg water/kg fuel, moisture content before ignition
			height_initial(:,:),					& ! m, initial total fuel height (nx,ny)
			actual_height(:,:),					& ! m, fuel height in the first layer in the vertical (nx,ny)
			mu_soot(:),								& ! m, mu in the log normal pdf for soot diameter
			sigma_soot(:),							& ! m, sigma in the log normal pdf for soot diameter
			pm_emissions(:),                 & ! g, soot produced
			co_emissions(:),                 & ! g, co produced
			conv_human(:),							& ! W/m^2 skin, convective heat flux to a human
			therm_dose(:),							& ! (W/m^2)^(4/3)*s, the thermal dosage, a cummulative quantity equal to (heat flux)**4/3*delta_t
			windmag(:),								& ! m/s, magnitude of the wind in a fuel cell
			o2density(:),							& ! kg/m3, oxygen concentration in air (reflects depeletion due to combustion)
			l_gap(:),								& ! m, length scale for fuel gaps
			l_patch(:)								  ! m, length scale for fuel patches

		integer :: 									&
			height_flag,						   & ! N/A, 1 = uniform, 2 = user provides QF_FuelHeight.inp,
														  !      3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)
			density_flag,							& ! N/A, 1 = uniform, 2 = user provides QF_FuelDensity.inp,
														  !      3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)
			moisture_flag							  ! N/A, 1 = uniform, 2 = user provides QF_FuelMoisture.inp,
														  !		3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)		
	
		real ::										&
			input_density,							& ! kg/m3, uniform fuel density value
			input_moisture,						& ! [0-1], uniform fuel moisture value
			input_height							  ! m, uniform fuel height value
	
	end type fuelgrid

	end module fuel_module
!===========================================================================
!===========================================================================
	module fire_module
	
	use coord_module
	
	implicit none
	
	type:: fire_type
		type(coord),allocatable ::				&
			burncenter(:),							& ! N/A, centroid of the ignitions coming into a cell normalized by the cell size, in x,y,z direction
			burncenter_info(:),					& ! to keep the history of the ignitions coming into a cell
			spat_var(:)								  ! N/A, standard deviation of the centroid of the ignitions normalized by the cell size.
														  !     Assuming a uniform distribution, the length in x and y dir containing the ignitions
														  !	  is 2*sqrt(3)*spat_var and the area is 12*spat_var(1)*spat_var(2)

		integer, allocatable ::					&
			ignitions(:),							& ! N/A, number of ignitions			
			num_ignitions_not_absorbed(:,:,:)  ! N/A, ignitions that are launched but not absorbed and thus can be given to atmosphere

		real, allocatable ::						&
			reaction_rate(:),						& ! kg/m3/s, reaction rate
			energy_to_atmos(:,:,:),				& ! kW/m3, energy released to the atmosphere
			timedecay(:),							& ! N/A, fraction of the maximum burning at each time step of the burning
			energy_used_to_evaporate_water(:)

	end type fire_type
	
	TYPE fire_cell_type
		real,dimension(2) ::						&
			sum_pos,									& ! N/A, sum of the positions of new ignitions (used to update the centroid)
			sum_pos_squared						  ! N/A, sum of the square of the positions of new ignitions
		
		real,dimension(2) ::						&
			sum_dcw_dcxy,							& ! summed over time for calculation of standard deviation of depl
			sum_dcw_dcxy2							  ! summed over time for calculation of standard deviation
		
	end TYPE fire_cell_type

	end module fire_module
!===========================================================================
!===========================================================================
	module firebrand_module
	
	use fireca_constants_module
	
	implicit none
	
	real, parameter :: tree_height = 15.	  ! m
	
	type :: FB_Type
		integer ::													&
			flag,														& ! N/A, 0 = no firebrands, 1 = yes firebrands
			time,														& ! s, total time
			LAUNCH_TIME	= 10,										& ! s, how often to launch firebrands
			min_number_of_ignitions = 50,						&
			max_number_of_ignitions = 100,					&
			germination_delay = 300

		real ::														&			
			time_step,												& ! s, firebrand transport time step
			FRACTION_LAUNCHED,									& ! N/A, fraction of the cells on fire that will launch a firebrand
			c_s,														& ! N/A, constant
			num_deposited,											&
			FRACTION_LAUNCHED_to_RT_ratio,					&
			frac_of_max_size,										&
			min_b_value_coef,										&
			min_theta_value,					 					&
			w_mult 				  								 	  ! N/A, w-enhancement for firebrands
		integer, dimension(:),allocatable ::				&		
			time_delay,												&
			num_ignitions
		
	end type FB_Type
	
	end module firebrand_module
!===========================================================================
!===========================================================================
	module fireca_module
	
	use constants
	use fuel_module
	use fire_module	
	use fireca_constants_module
	use firebrand_module
	use omp_lib

	implicit none
	
	integer ::										&
		CREEPING_FLAG = 0							  ! N/A, 0 = creeping off, 1 = creeping on

	! Working variables
	real, dimension(:),allocatable ::		&
		w_rr_o2_turb,								& ! added to track o2 and turbulence contribution to reaction rate
		w_rr,											& ! kg/m3/s, current reaction rate
		w_mass,										& ! kg, current fuel mass
		w_moisture,									& ! kg water/kg fuel, current moisture
		w_fueldensity,								& ! kg/m3, current fuel density
		w_q,											& ! kW/m3, (reaction rate) * (wood combustion entalphy) = energy released per cubic meter
		n_new_starts,								&
		w_depl,										& ! added to track weights for depletion center calc
		init_dens_canopy

	integer, dimension(:),allocatable ::	&
		w_ignitions

	real, dimension(:,:),allocatable ::		&
		mass_int_0									  ! kg, vertically integrated initial mass

	integer, allocatable,dimension(:,:,:) ::	&
		grid_where ! N/A, it tells where a cell i,j,k is in the compressed array scheme. It is -1 if not present

	integer ::										&
		nbtime,										& ! N/A, = int(BURNOUT_TIME/real(fcatime%dt_int))
		max_init_ign								  ! N/A, maximum number of initial ignitions

	real ::											&
		nbtime_real,								& ! N/A, = BURNOUT_TIME/real(fcatime%dt_int)
		sumdecay,									& ! N/A, cumulative sum over fire%timedecay
		partical_burn_vel_dt						  ! partical_burn_vel*dt

	real,dimension(:),allocatable :: timedecay_ratio  ! N/A, ratio between consecutive values of timedecay


	type(fire_type) :: fire
	type(fire_cell_type),dimension(:),allocatable :: firecell
	type(fuelgrid) :: fuels	
	type(FB_Type) :: fb
	integer(omp_lock_kind), dimension(:,:,:), allocatable :: lock
	integer(omp_lock_kind), dimension(:), allocatable :: fire_spread_lock
	
	end module fireca_module
!======================================================================================
!======================================================================================
	
	

