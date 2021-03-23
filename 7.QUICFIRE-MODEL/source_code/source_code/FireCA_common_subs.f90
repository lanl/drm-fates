	subroutine ComputeNewFireSpread()

	use fireca_constants_module
	use fireca_module
	use time_module
	use omp_lib	
	use file_handling_module
	use flags_module
	use grid_module
	use winds_module

	implicit none

	integer,dimension(firegrid%num_fuel_cells) :: &
		new_ignitions,						& ! N/A, how many ignitions are actually generated, assuming that
												  !      only the first 1/3 most recent ignited fuel is in flaming combustion
		supported_ignitions,				& ! N/A, how many ignitions could be generated with this much energy
												  !      to the fuel with extra_energy assigned stocastically
		new_smold_ignitions				  ! N/A, how many ignitions are actually generated, assuming that
												  !      only final 2/3 or burn time contributes to smoldering

	integer ::								&
		index,								&
		idx,									&
		TID,									&
		ew,ns,ud

	real :: 									&
		w_fraction_to_atmos,				& ! N/A, fraction of energy to the atmosphere
		net_energy_rate,					& ! KW/m3, energy released by the combustion
		energy_to_fuels,					& ! kW/m3, energy to the fuel
		lengthscale_area,					& ! m2, is the minimum between the horizontal area of a cell and the surface area that is on fire
		avg_u, avg_v, avg_w,				& ! m/s, average u,v,w when the quic and fire grid are not consistent; quic value otherwise
		avg_rr,								& ! kg/m3/s, average reaction rate when the quic and fire grid are not consistent;
		!	previous fireca value otherwise
		avg_q,								& ! (kW/m3)**(1/3), average heat release in the combustion, to be used to calculate a wprime
		wprime,								& ! m/s, perturbation to the average vertical velocity calculated by quic
		local_w,								& ! m/s, local vertical velocity = average + perturbation (wprime)
		hwindmag,							& ! m/s, horizontal wind magnitude (sqrt(u**2+v**2))
		wind2d_sq,							& ! (m/s)^2, u^2 + v^2
		wind3d,								& ! m/s, sqrt(u^2 + v^2	+ w^2)
		psi,									& ! N/A, fraction of mass in a cell burning because of the ignition energy
		k,										&
		mixing								  ! m2/s, turbulent viscosity = 0.1 * (turbulence scale) * sqrt(turbulent kinetic energy)

	real, dimension(firegrid%num_fuel_cells) ::		&
		lengthscale,						& ! m, flame length
		k2sqrt,								&
		wstar,								& ! m/s
		mass_burned						     ! g, mass burned
	
	mass_burned = 0
	new_ignitions = 0
	supported_ignitions = 0
	new_smold_ignitions = 0
	
	!$OMP parallel do private(idx,index,TID,lengthscale_area,wind2d_sq,wind3d,			&
	!$OMP avg_u,avg_v,avg_w,avg_rr,avg_q,ew,ns,ud,wprime,local_w,							&
	!$OMP hwindmag,k,mixing,psi,w_fraction_to_atmos,net_energy_rate,energy_to_fuels)
	do idx = 1,firegrid%num_active_fuel_cells
		index = firegrid%idx(idx)

		TID = 1

		!!!!!!! DO NOT REMOVE THE FOLLOWING LINE ABOUT TID
		!$ TID = TID + OMP_GET_THREAD_NUM()

		ew = firegrid%ijk_cell_index(index,1)
		ns = firegrid%ijk_cell_index(index,2)
		ud = firegrid%ijk_cell_index(index,3)

		! -- lengthscale_area is the minimum between the horizontal area of a cell and the surface area that is on fire, m2
		! Assuming a uniform distribution for the ignitions centroid,
		! the normalized surface area that is on fire =
		! (2*sqrt(3)*fire%spat_var(index,1))*(2*sqrt(3)*fire%spat_var(index,2)) =
		! 12*fire%spat_var(index,1)*fire%spat_var(index,2)
		call lengthscale_area_calc(fire%spat_var(index)%val, firegrid%cellarea, lengthscale_area)			
		
		! -- flame length based on empirical jet correlation by Turns, m
		wind2d_sq = fcawinds%u(ew,ns,ud)**2 + fcawinds%v(ew,ns,ud)**2
		wind3d = sqrt(wind2d_sq + fcawinds%w(ew,ns,ud)**2)
		call fire_length_scale_calc(wind2d_sq,wind3d,									&
			lengthscale_area,fire%ignitions(index),										&
			fuels%o2density(index),fuels%height_initial(ew,ns),						&
			firegrid%cellarea,firegrid%dz_array(ud),TID, lengthscale(index))			 

		! -- average velocities and reaction rate when grids are not consistent
		call ComputeAverages(ew,ns,ud,index,lengthscale(index),						&
			avg_u,avg_v,avg_w,avg_rr,avg_q)

		call ComputeWindComponents(w_q(index),avg_q,avg_u,avg_v,avg_w,				&
			firegrid%zm(ud),firegrid%cellvol(ud),wprime,local_w,hwindmag,			&
			fuels%windmag(index))

		call ComputeTurbulence(ew,ns,ud,fuels%windmag(index),wprime,				&
			lengthscale(index),firegrid%delta(ud),fcawinds%sigma(index),k,			&
			mixing,k2sqrt(index))

		call ComputeReactionVariables(firegrid%cellarea,lengthscale_area,       &
			firegrid%cellvol(ud),fuels%windmag(index),avg_rr,local_w,mixing,		&
			fuels%density(index),fuels%density_initial(index),							&
			fire%ignitions(index),wstar(index),psi,w_fraction_to_atmos,				&
			net_energy_rate,w_rr(index),w_rr_o2_turb(index),							&
			w_fueldensity(index),w_mass(index),energy_to_fuels,						&
			lengthscale(index),fuels%o2density(index),mass_burned(index),			&
			fuels%height_initial(ew,ns),firegrid%dz_array(ud),TID)

		call UpdateBurnCenter(nbtime,fire%timedecay,										&
			w_rr(index),w_rr_o2_turb(index),sumdecay,w_depl(index),					&
			firecell(index)%sum_dcw_dcxy,firecell(index)%sum_dcw_dcxy2,				&
			fire%burncenter_info(index)%val2,fire%burncenter(index)%val,			&
			fire%spat_var(index)%val,fuels%depl_var(index)%val,						&
			fuels%depl_center(index)%val)

		if(flag%thermalrad_file == 1) then
			call ComputeConvectiveHeatTransfer(fuels%windmag(index),					&
				k2sqrt(index),psi,fuels%conv_human(index),fuels%density(index))
					
			fuels%therm_dose(index) = fuels%therm_dose(index) +						&
				fuels%conv_human(index)**FOUR_THIRDS*fcatime%dt_real

		endif

		if (w_fueldensity(index) > MIN_FUEL_DENSITY .and. w_rr(index) > 0.0) then

			call GenerateNewIgnitions(TID,nbtime,fire%ignitions(index),				&
				firegrid%cellvol(ud),energy_to_fuels,w_rr_o2_turb(index),			&
				w_rr(index),fire%burncenter_info(index)%val2,w_ignitions(index),	&
				new_ignitions(index),new_smold_ignitions(index),						&
				supported_ignitions(index))

		else
			call NoFuel(fire%energy_to_atmos(ew,ns,ud),w_mass(index),				&
				w_moisture(index),w_ignitions(index))
		endif		
	enddo
	!$OMP end parallel do
	
	call LoopOverIgnitions(new_ignitions,1,										&
		lengthscale,k2sqrt,wstar,supported_ignitions,mass_burned)

	if(CREEPING_FLAG == 1) then
		call LoopOverIgnitions(new_smold_ignitions,2,							&
			lengthscale,k2sqrt,wstar,supported_ignitions,mass_burned)
	endif

	fire%energy_to_atmos = fire%energy_to_atmos +		&
		real(fire%num_ignitions_not_absorbed) 	  *		&
		ENERGY_PER_IGNITION	
	!$OMP parallel do private(ud)
	do ud = 1, firegrid%nz_en2atmos
		fire%energy_to_atmos(:,:,ud) = fire%energy_to_atmos(:,:,ud) / firegrid%cellvol_en2atmos(ud)			
	enddo
	!$OMP end parallel do

	end
!======================================================================================
!======================================================================================
	SUBROUTINE LoopOverIgnitions(new_ignitions, function_switch,	&
		lengthscale, k2sqrt, wstar, supported_ignitions, mass_burned)

	use time_module
	use fireca_module
	use grid_module
	use flags_module
	use omp_lib	
	use rnd_num_vars
	use winds_module

	implicit none

	integer, intent(IN) ::								&			
		function_switch									  ! N/A, to know to call the main wind propagation (1) or the upwind propagation routine (2)

	integer,intent(IN), dimension(firegrid%num_fuel_cells) :: & 
		supported_ignitions,								&
		new_ignitions										  ! N/A, number of ignitions that have been generated
	
	real,intent(IN), dimension(firegrid%num_fuel_cells) :: & 
		lengthscale,										&
		k2sqrt,												&
		wstar,												& ! m/s
		mass_burned
	
	integer ::												&
		ud_destination,									&
		TID,													& ! N/A, thread counter
		index,												& ! N/A, index of the current cell
		ncell,												& ! N/A, linear index
		it,													& ! N/A, time counter
		ig,													& ! N/A, ignition counter
		ew_temp,ns_temp,ud_temp,						& ! N/A, cell indexes
		w_ew,w_ns,w_ud,									& ! N/A, cell indexes
		nosubseqent_jump									  ! N/A, flag to say that an ignition has or has not generated a new fire

	real,dimension(3) ::									&
		start_coord,										&
		landlocation,										&  ! N/A, where the energy packet lands [0,1]
		launchlocation,									&  !N/A, where energy packet originates from [0,1]
		flame_components									   ! distance in x,y,z that energy packet is moved (m)

	real ::													&
		current_soot,										&
		cell_soot_sum,										& ! g, sum of the soot in one cell [er time step
		flame_counter,										& ! N/A, total number of flames in one cell per time step
		mean_mu,												& ! m, mean of mu_soot
		mean_sigma,											& ! m, mean of sigma_soot
		d,														& ! m, flame length
		single_flame_soot,								& ! N/A, single flame value of the soot fraction
		single_flame_mu,									& ! m, single flame value of mu
		single_flame_sigma 								  ! m, single flame valu of sigma	
	
	integer, dimension(firegrid%nx, firegrid%ny, firegrid%nz_en2atmos, numThreads) :: num_ignitions_not_absorbed
	
	num_ignitions_not_absorbed = 0
	
	!$OMP parallel do private(index, cell_soot_sum, flame_counter, mean_sigma, ig,				&
	!$OMP ew_temp, ns_temp, ud_temp, start_coord, nosubseqent_jump, it, TID,						&
	!$OMP w_ew, w_ns, w_ud, landlocation, d, ncell,															&
	!$OMP single_flame_soot, single_flame_mu, single_flame_sigma, current_soot, mean_mu,		&
	!$OMP launchlocation, ud_destination, flame_components)
	do index = 1,firegrid%num_fuel_cells
		!initialize the soot variables
		cell_soot_sum = 0
		flame_counter = 0
		mean_mu = 0
		mean_sigma = 0

		TID = 1
		!$ TID = TID + OMP_GET_THREAD_NUM()
		
		do ig = 1, new_ignitions(index)

			ew_temp = firegrid%ijk_cell_index(index,1)
			ns_temp = firegrid%ijk_cell_index(index,2)
			ud_temp = firegrid%ijk_cell_index(index,3)

			start_coord = (/fire%burncenter(index)%val(1),fire%burncenter(index)%val(2),0.5/)

			nosubseqent_jump = 0

			it = 0
			do while(it < fcatime%dt_int .and. nosubseqent_jump == 0)
				it = it + 1

				if(function_switch == 1) then
					call FirePropagationMainWind(TID,it,ew_temp,ns_temp,ud_temp,					&
						lengthscale(index),k2sqrt(index),wstar(index),start_coord,					&
						fire%spat_var(index)%val,w_ew,w_ns,w_ud,d,landlocation,						&
						launchlocation, flame_components, fuels%actual_height(ew_temp,ns_temp))
				else
					call FirePropagationCreeping(TID,it,ew_temp,ns_temp,ud_temp,					&
						partical_burn_vel_dt,start_coord,fire%spat_var(index)%val,					&
						fcawinds%u(ew_temp,ns_temp,ud_temp),fcawinds%v(ew_temp,ns_temp,ud_temp),&
						w_ew,w_ns,w_ud, landlocation, launchlocation, flame_components)
				endif


				if(flag%emissions_file >= 1 .and. flag%emissions_file <= 3) then
					! Get soot fraction and soot distribution paramters for the cell, one flame at a time
					if (function_switch == 2) then
						call ComputeSootYield(0.1,0.25,fuels%o2density(index),						&
							single_flame_soot,single_flame_mu,single_flame_sigma)
					else
						call ComputeSootYield(d,fuels%windmag(index),fuels%o2density(index),		&
							single_flame_soot,single_flame_mu,single_flame_sigma)
					endif
					cell_soot_sum = cell_soot_sum + single_flame_soot
					flame_counter = flame_counter + 1 !total flames in the cell
					mean_mu = mean_mu + single_flame_soot*single_flame_mu !add up the mu values to be normalized later
					mean_sigma = mean_sigma + single_flame_soot*single_flame_sigma !add up the sigma values to be normalized later
				endif
		
				if(w_ew .gt. 0 .and. w_ew .lt. firegrid%nx .and. &
					w_ns .gt. 0 .and. w_ns .lt. firegrid%ny .and. &
					w_ud .gt. 0 .and. w_ud .le. firegrid%nz) then
					
					! Is there fuel in this cell					
					ncell = grid_where(w_ew,w_ns,w_ud)
					
					if(ncell > 0)then
												
						if(fuels%density(ncell) > MIN_FUEL_DENSITY)then
							CALL omp_set_lock(fire_spread_lock(ncell))
							call SpreadFire(w_ew,w_ns,w_ud,TID,landlocation,launchlocation,							&
								fuels%actual_height(ew_temp,ns_temp), fuels%actual_height(w_ew,w_ns),				&
								flame_components,																						&
								fuels%l_gap(ncell), fuels%l_patch(ncell),firegrid%cellvol(w_ud),						&
								fuels%density_initial(ncell),fuels%density(ncell),											&
								fuels%moisture_initial(ncell),fuels%moisture(ncell),										&
								fuels%depl_var(ncell)%val,fuels%depl_center(ncell)%val,n_new_starts(ncell),		&
								w_ignitions(ncell),nosubseqent_jump,ew_temp,ns_temp,ud_temp,							&
								firecell(ncell)%sum_pos,firecell(ncell)%sum_pos_squared,									&
								fuels%moist_depl_var(ncell)%val,fuels%moist_depl_center(ncell)%val,start_coord,	&
								fire%num_ignitions_not_absorbed(w_ew,w_ns,w_ud),											&
								fire%energy_used_to_evaporate_water(ncell),													&
								fire%burncenter(ncell)%val,fire%spat_var(ncell)%val,function_switch,					&
								supported_ignitions(ncell), fcawinds%w(w_ew,w_ns,w_ud), fuels%windmag(ncell),		&
								w_depl(ncell))
							CALL omp_unset_lock(fire_spread_lock(ncell))
						else							
							nosubseqent_jump = 1							
							num_ignitions_not_absorbed(w_ew,w_ns,w_ud,TID) =												&
								num_ignitions_not_absorbed(w_ew,w_ns,w_ud,TID) + 1							
						endif
					else
						num_ignitions_not_absorbed(w_ew,w_ns,w_ud,TID) =													&
								num_ignitions_not_absorbed(w_ew,w_ns,w_ud,TID) + 1						
						nosubseqent_jump = 1
					endif
					
				else
					if (w_ew .gt. 0 .and. w_ew .lt. firegrid%nx .and. &
							w_ns .gt. 0 .and. w_ns .lt. firegrid%ny) then						
						ud_destination = min(w_ud, firegrid%nz_en2atmos)						
						num_ignitions_not_absorbed(w_ew,w_ns,ud_destination,TID) =										&
								num_ignitions_not_absorbed(w_ew,w_ns,ud_destination,TID) + 1						
					endif
					nosubseqent_jump = 1
				endif
			enddo
		enddo
	
		if(flag%emissions_file >= 1 .and. flag%emissions_file <= 3) then
			if(cell_soot_sum > 0) then
				current_soot = cell_soot_sum/real(flame_counter)*mass_burned(index)
				! Get the average soot fraction from all flames (even non-sooting)
				fuels%pm_emissions(index) = fuels%pm_emissions(index) + current_soot
				! Update mean_mu to normalize it by the total soot produced in this time step
				mean_mu = mean_mu/cell_soot_sum
				! Update fuels%mu_soot to normalize it by the total soot produced since the last printout.
				! This is done in several steps: 1) The previous weight, non-normalized sum is recoverd
				! from prior time steps by multiplying the normalized sum by the normalization factor (
				! the previous fuels%pm_emissions value). On the first time step, this results in a negative
				! normalization factor, but fuels%mu_soot is 0 on the first time step, so this is inconse-
				! quential. 2) The current mean_mu is weighted by the current amout of soot and added to
				! the weighted, non-normalized sum. 3) The weighted sum is renormalized by the cumulative
				! sum of soot since the last printout. The same process is applied to fuels%sigma_soot.
				fuels%mu_soot(index) = (fuels%mu_soot(index)*						&
					(fuels%pm_emissions(index)-current_soot) +						&
					mean_mu*current_soot) / fuels%pm_emissions(index)
				! Update mean_sigma to normalize it by the total soot produced in this time step
				mean_sigma = mean_sigma/cell_soot_sum
				fuels%sigma_soot(index) = (fuels%sigma_soot(index)*				&
					(fuels%pm_emissions(index)-current_soot)+							&
					mean_sigma*current_soot) / fuels%pm_emissions(index)
				fuels%co_emissions(index) = fuels%co_emissions(index) +			&
					current_soot*28.01*PRESSURE/(RGAS*13.9*AMBIENT_TEMPERATURE)
			endif
		endif
	enddo
	!$OMP end parallel do
	
	
	fire%num_ignitions_not_absorbed = fire%num_ignitions_not_absorbed + sum(num_ignitions_not_absorbed, DIM=4)
	
	end
!======================================================================================
!======================================================================================
	subroutine output_thermal_radiation(simtime)

	use fireca_module
	use file_handling_module	

	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A17,I0.5,A4)") "thermalradiation-", simtime, ".bin"
	open (unit=ID_FILE_CONVHUMAN, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_CONVHUMAN) fuels%conv_human
	close(ID_FILE_CONVHUMAN)


	write(filename, "(A12,I0.5,A4)") "thermaldose-", simtime, ".bin"
	open (unit=ID_FILE_THERMDOSE, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_THERMDOSE) fuels%therm_dose
	close(ID_FILE_THERMDOSE)

	end
!======================================================================================
!======================================================================================
	subroutine output_emissions(simtime)

	use fireca_module
	use file_handling_module
	
	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A13,I0.5,A4)") "pm_emissions-", simtime, ".bin"
	open (unit=ID_FILE_PM_EMISSIONS, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_PM_EMISSIONS) fuels%pm_emissions
	close(ID_FILE_PM_EMISSIONS)

	write(filename, "(A23,I0.5,A4)") "emissions_distribution-", simtime, ".bin"
	open (unit=ID_FILE_PM_EMISSIONS_DISTR, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_PM_EMISSIONS_DISTR) fuels%mu_soot
	write (ID_FILE_PM_EMISSIONS_DISTR) fuels%sigma_soot
	close(ID_FILE_PM_EMISSIONS_DISTR)

	write(filename, "(A13,I0.5,A4)") "co_emissions-", simtime, ".bin"
	open (unit=ID_FILE_CO_EMISSIONS, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_CO_EMISSIONS) fuels%co_emissions
	close(ID_FILE_CO_EMISSIONS)

	end
!======================================================================================
!======================================================================================
	subroutine output_react_rate(simtime)
	
	use fireca_module
	use file_handling_module
	
	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A19,I0.5,A4)") "fire-reaction_rate-", simtime, ".bin"
	open (unit=ID_FILE_REACTRATE, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_REACTRATE) fire%reaction_rate
	close(ID_FILE_REACTRATE)

	end
!======================================================================================
!======================================================================================
	subroutine output_en2atm(simtime)
	
	use fireca_module
	use file_handling_module

	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A21,I0.5,A4)") "fire-energy_to_atmos-", simtime, ".bin"
	open (unit=ID_FILE_EN2ATM, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_EN2ATM) fire%energy_to_atmos	
	close(ID_FILE_EN2ATM)

	end
!======================================================================================
!======================================================================================
	subroutine output_massburnt(simtime)

	use fireca_module
	use file_handling_module
	use grid_module

	implicit none

	integer,intent(IN) :: simtime
	character(len=32) :: filename
	integer :: i,j,k
	real,dimension(firegrid%nx,firegrid%ny) :: mass_int

	write(filename, "(A13,I0.5,A4)") "mburnt_integ-", simtime, ".bin"
	open (unit=ID_FILE_MASSBURNT, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	mass_int = 0.
	
	!$OMP parallel do private(k)
	do k = 1,firegrid%num_fuel_cells
		mass_int(firegrid%ijk_cell_index(k,1),firegrid%ijk_cell_index(k,2)) =		&
			mass_int(firegrid%ijk_cell_index(k,1),firegrid%ijk_cell_index(k,2)) +	&
			fuels%density(k) * firegrid%dz_array(firegrid%ijk_cell_index(k,3))
	enddo
	!$OMP end parallel do
	
	!$OMP parallel do private(i,j)
	do j = 1,firegrid%ny
		do i = 1,firegrid%nx
			if(mass_int_0(i,j) > 0) then
				mass_int(i,j) = (mass_int_0(i,j) - mass_int(i,j)) / mass_int_0(i,j) * 100.
			else
				mass_int(i,j) = 0.
			endif
		enddo
	enddo
	!$OMP end parallel do

	write (ID_FILE_MASSBURNT) mass_int
	close(ID_FILE_MASSBURNT)

	end
!======================================================================================
!======================================================================================
	subroutine output_fuels(simtime)

	use fireca_module
	use file_handling_module
	
	implicit none

	character(len=32) :: filename
	integer,intent(IN) :: simtime

	write(filename, "(A11,I0.5,A4)") "fuels-dens-", simtime, ".bin"
	open (unit=ID_FILE_FUEL_OUT, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_FUEL_OUT) fuels%density
	close(ID_FILE_FUEL_OUT)

	end
!======================================================================================
!======================================================================================
	subroutine output_wind(simtime)

	use fireca_module
	use file_handling_module
	use winds_module

	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A5,I0.5,A4)") "windu", simtime, ".bin"
	open (unit=ID_FILE_WIND_FIRE_U, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_WIND_FIRE_U) fcawinds%u
	close(ID_FILE_WIND_FIRE_U)

	write(filename, "(A5,I0.5,A4)") "windv", simtime, ".bin"
	open (unit=ID_FILE_WIND_FIRE_V, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_WIND_FIRE_V) fcawinds%v
	close(ID_FILE_WIND_FIRE_V)

	write(filename, "(A5,I0.5,A4)") "windw", simtime, ".bin" 
	open (unit=ID_FILE_WIND_FIRE_W, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_WIND_FIRE_W) fcawinds%w
	close(ID_FILE_WIND_FIRE_W)

	return

	end
!======================================================================================
!======================================================================================
	subroutine output_moisture(simtime)

	use fireca_module
	use file_handling_module

	implicit none

	integer, intent(in) :: simtime
	character(len=32) :: filename

	write(filename, "(A12,I0.5,A4)") "fuels-moist-", simtime, ".bin"
	open (unit=ID_FILE_MOISTURE_OUT, file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(filename),form='unformatted')
	write (ID_FILE_MOISTURE_OUT) fuels%moisture
	close(ID_FILE_MOISTURE_OUT)

	end
!===================================================================================
!===================================================================================
	subroutine IndexCalc(ew,ns,ud,index)

	! Calculation of the linear index from the location of the cell in the
	! 3D volume

	use grid_module

	implicit none

	integer,intent(IN) :: ew,ns,ud   ! N/A, indexes in a 3D array
	integer :: index						! N/A, linear unique index

	index = ew + (ns-1) * firegrid%nx + (ud-1) * firegrid%nx * firegrid%ny

	end
!===================================================================================
!===================================================================================
	subroutine AllocateArrays()

	! Allocate and initialize some of the arrays used in the fire simulation
	use omp_lib
	use fireca_module
	use grid_module
	use constants
	use fireca_constants_module
	use flags_module
	use winds_module

	implicit none

	integer ::			&
		ew,ns,ud,		& ! N/A, indexes to loop over the x,y,z directions
		ncell,			& ! N/A, index to loop only over the cells with fuel	
		ii					  ! N/A, counter
	real ::				&
		ii_real,			& ! N/A, real(ii)
		spat_var_val_x,spat_var_val_y


	! allocate arrays only with cells with fuel
	allocate(grid_where(firegrid%nx,firegrid%ny,firegrid%nz))
	grid_where = 0
	! Calculate indexes for cells with fuel
	!$OMP parallel do private(ncell, ew, ns, ud)
	do ncell = 1,firegrid%num_fuel_cells
		ew = firegrid%ijk_cell_index(ncell,1)
		ns = firegrid%ijk_cell_index(ncell,2)
		ud = firegrid%ijk_cell_index(ncell,3)
		grid_where(ew,ns,ud) = ncell
	enddo
	!$OMP end parallel do

	allocate (fcawinds%sigma(firegrid%num_fuel_cells))

	allocate(																	&
		fuels%density(firegrid%num_fuel_cells),						&
		fuels%moisture(firegrid%num_fuel_cells),						&
		fuels%moisture_initial(firegrid%num_fuel_cells),			&		
		fuels%density_initial(firegrid%num_fuel_cells),				&
		fuels%depl_center(firegrid%num_fuel_cells),					&
		fuels%depl_var(firegrid%num_fuel_cells),						&
		fuels%moist_depl_center(firegrid%num_fuel_cells),			&
		fuels%moist_depl_var(firegrid%num_fuel_cells),				&		
		fuels%windmag(firegrid%num_fuel_cells),						&
		fuels%l_gap(firegrid%num_fuel_cells), 							&
		fuels%l_patch(firegrid%num_fuel_cells),						&
		fuels%o2density(firegrid%num_fuel_cells))
	
	fuels%o2density = INITIAL_O2DENSITY
	
	allocate(fire_spread_lock(firegrid%num_fuel_cells))
	!$OMP parallel do private(ud)
	do ud = 1, firegrid%num_fuel_cells
		call omp_init_lock(fire_spread_lock(ud))
	enddo
	!$OMP end parallel do 
	
	allocate(firecell(firegrid%num_fuel_cells))

	if(flag%emissions_file >= 1 .and. flag%emissions_file <= 3) then
		allocate(																&
			fuels%pm_emissions(firegrid%num_fuel_cells),				&
			fuels%co_emissions(firegrid%num_fuel_cells),				&
			fuels%mu_soot(firegrid%num_fuel_cells),					&
			fuels%sigma_soot(firegrid%num_fuel_cells))
	endif

	if(flag%thermalrad_file == 1) then
		allocate(fuels%conv_human(firegrid%num_fuel_cells),		&
			fuels%therm_dose(firegrid%num_fuel_cells))
	endif

	allocate(																	&
		fire%ignitions(firegrid%num_fuel_cells),						&		
		fire%reaction_rate(firegrid%num_fuel_cells),					&
		fire%burncenter(firegrid%num_fuel_cells),						&
		fire%burncenter_info(firegrid%num_fuel_cells),				&
		fire%spat_var(firegrid%num_fuel_cells))

	allocate(																	&
		fb%time_delay(firegrid%num_fuel_cells),						&
		fb%num_ignitions(firegrid%num_fuel_cells))

	! Initialize value
	if(flag%emissions_file >= 1 .and. flag%emissions_file <= 3) then
		fuels%pm_emissions = 0
		fuels%co_emissions = 0
		fuels%mu_soot = 0
		fuels%sigma_soot = 0
	endif
	if(flag%thermalrad_file == 1) then
		fuels%conv_human = 0.
		fuels%therm_dose = 0.
	endif
	fb%time_delay = FB_TIME_DEFAULT
	fb%num_ignitions = 0
	fire%ignitions = 0
	fire%reaction_rate = 0.0	

	spat_var_val_x = init_var_length*ONE_OVER_ROOTTHREE * min(1., firegrid%dxi)
	spat_var_val_y = init_var_length*ONE_OVER_ROOTTHREE * min(1., firegrid%dyi)

	do ii = 1,firegrid%num_fuel_cells
		allocate(												&
			fuels%moist_depl_var(ii)%val(2),				&
			fuels%moist_depl_center(ii)%val(2),			&
			fuels%depl_var(ii)%val(2),						&
			fuels%depl_center(ii)%val(2),					&
			fire%spat_var(ii)%val(2),						&
			fire%burncenter(ii)%val(2),					&
			fire%burncenter_info(ii)%val2(nbtime,3))


		fuels%moist_depl_var(ii)%val = moist_depl_var_init
		fuels%moist_depl_center(ii)%val = 0.0
		fuels%depl_var(ii)%val = 0.0
		fuels%depl_center(ii)%val = 0.0

		fire%spat_var(ii)%val(1) = spat_var_val_x
		fire%spat_var(ii)%val(2) = spat_var_val_y

		fire%burncenter(ii)%val = BURNCENTER_INIT
		fire%burncenter_info(ii)%val2 = 0
	enddo

	allocate(																&
		w_rr_o2_turb(firegrid%num_fuel_cells),						&
		w_rr(firegrid%num_fuel_cells),								&
		w_mass(firegrid%num_fuel_cells),								&
		w_moisture(firegrid%num_fuel_cells),						&
		w_fueldensity(firegrid%num_fuel_cells),					&
		w_q(firegrid%num_fuel_cells),									&
		n_new_starts(firegrid%num_fuel_cells),						&
		w_depl(firegrid%num_fuel_cells),								&
		w_ignitions(firegrid%num_fuel_cells))

	allocate(fire%energy_to_atmos(firegrid%nx,firegrid%ny,firegrid%nz_en2atmos))
	fire%energy_to_atmos = 0.
	allocate(fire%num_ignitions_not_absorbed(firegrid%nx,firegrid%ny,firegrid%nz_en2atmos))
	allocate(fire%energy_used_to_evaporate_water(firegrid%num_fuel_cells))

	allocate(fire%timedecay(nbtime),timedecay_ratio(nbtime))
	sumdecay = 0.0
	do ii = 1,nbtime
		ii_real = float(ii-1)
		fire%timedecay(ii) = (nbtime_real-ii_real)/nbtime_real
		sumdecay = sumdecay + fire%timedecay(ii)
	enddo

	timedecay_ratio(1) = 1.
	do ii = 2,nbtime
		timedecay_ratio(ii) = fire%timedecay(ii)/fire%timedecay(ii-1)
	enddo
	w_depl = -1

	end
!===================================================================================
!===================================================================================
	subroutine ReadFBFile()

	use constants
	use fireca_module
	use file_handling_module
	
	implicit none

	open(ID_FILE_FB_IN,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//'QFire_Advanced_User_Inputs.inp', &
		status='old',err = 1106)

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%FRACTION_LAUNCHED
	if(fb%FRACTION_LAUNCHED <= 0 .or. fb%FRACTION_LAUNCHED > 1) then
		write (msgoutfile,*)'Invalid franction of firebrands launched. Must be (0 1].'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%c_s
	if(fb%c_s <= 0) then
		write (msgoutfile,*) &
			'Invalid scaling factor of the radius represented by the firebrands launched. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%time_step
	if(fb%time_step <= 0) then
		write (msgoutfile,*) &
			'Invalid timestep for the firebrands trajectory calculation. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%LAUNCH_TIME
	if(fb%LAUNCH_TIME <= 0) then
		write (msgoutfile,*)'Invalid firebrands launch interval. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%num_deposited
	if(fb%num_deposited <= 0) then
		write (msgoutfile,*) &
			'Invalid number of firebrands distributed over the landing area. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%FRACTION_LAUNCHED_to_RT_ratio
	if(fb%FRACTION_LAUNCHED_to_RT_ratio <= 0) then
		write (msgoutfile,*) &
			'Invalid value for fb%FRACTION_LAUNCHED_to_RT_ratio. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107) fb%min_b_value_coef
	if(fb%min_b_value_coef <= 0) then
		write (msgoutfile,*) &
			'Invalid value for min_b_value_coef. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%frac_of_max_size
	if(fb%frac_of_max_size <= 0) then
		write (msgoutfile,*) &
			'Invalid fraction of max size. Must be > 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%germination_delay
	if(fb%germination_delay < 0) then
		write (msgoutfile,*) 'Invalid germination delay. Must be >= 0.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%w_mult
	if(fb%w_mult < 1) then
		write (msgoutfile,*) 'Invalid fraction of cell on fire. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%min_number_of_ignitions
	if(fb%min_number_of_ignitions <= 0) then
		write (msgoutfile,*) 'Invalid minimum number of ignitions via firebrands. Must be >= 1.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%max_number_of_ignitions
	if(fb%max_number_of_ignitions <= 0) then
		write (msgoutfile,*) &
			'Invalid maximum number of ignitions via firebrands. Must be >= 1.'
		call TerminateProgram()
	endif

	if(fb%max_number_of_ignitions <= fb%min_number_of_ignitions) then
		write (msgoutfile,*) &
			'Minimum and maximum number of ignitions via firebrands are inconsistent.'
		call TerminateProgram()
	endif

	read(ID_FILE_FB_IN,*,err=1107,end=1107)	fb%min_theta_value
	if(fb%min_theta_value < 0 .or. fb%min_theta_value > 0.5*GREEK_PI) then
		write (msgoutfile,*) &
			'Invalid min_theta_value. Must be [0 PI/2].'
		call TerminateProgram()
	endif

	close(ID_FILE_FB_IN)

1106 return

1107 write(msgoutfile,*)'Error in reading file QFire_Advanced_User_Inputs.inp.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	SUBROUTINE ReadGridlist()

	use fireca_module
	use file_handling_module
	use grid_module
	
	implicit none

	character(STR_LEN) :: fname								! N/A, file name
	character(STR_LEN) :: tt									! N/A, used to read in file lines
	integer,dimension(10) :: is_read
	integer :: do_stop,iequal,ispace,i,istring_start
	real :: aa1,aa3,zb

	fname = 'gridlist'
	open(ID_FILE_GRIDLIST,file=TRIM(ft%file_path)// &
		TRIM(fileSeparator)//trim(fname),err=1101,status='old')

	is_read = 0

	do while(.true.)
		read(ID_FILE_GRIDLIST,'(a)',end=1345,err=1102)tt
		tt = adjustl(tt)

		do_stop = 0
		istring_start = 1
		do while(do_stop == 0)
			iequal = scan(tt,"=")
			if(iequal == 0) then
				do_stop = 1
			else
				ispace = scan(tt(iequal+1:)," ") + iequal

				if( tt(istring_start:iequal-1) == 'n')then
					read(tt(iequal+1:ispace-1),*)ft%nx
					is_read(1) = 1
				elseif( tt(istring_start:iequal-1) == 'm')then
					read(tt(iequal+1:ispace-1),*)ft%ny
					is_read(2) = 1
				elseif( tt(istring_start:iequal-1) == 'l')then
					read(tt(iequal+1:ispace-1),*)ft%nz
					is_read(3) = 1
				elseif( tt(istring_start:iequal-1) == 'dx')then
					read(tt(iequal+1:ispace-1),*)ft%dx
					is_read(4) = 1
				elseif( tt(istring_start:iequal-1) == 'dy')then
					read(tt(iequal+1:ispace-1),*)ft%dy
					is_read(5) = 1
				elseif( tt(istring_start:iequal-1) == 'dz')then
					read(tt(iequal+1:ispace-1),*)zb
					is_read(6) = 1
				elseif( tt(istring_start:iequal-1) == 'aa1')then
					read(tt(iequal+1:ispace-1),*)aa1
					is_read(7) = 1
				elseif( tt(istring_start:iequal-1) == 'l_bld_top')then
					read(tt(iequal+1:ispace-1),*)ft%n_grid_top
					is_read(8) = 1
				elseif( tt(istring_start:iequal-1) == 'nfueltypes')then
					read(tt(iequal+1:ispace-1),*)ft%n_fuel_types
					is_read(9) = 1
				elseif( tt(istring_start:iequal-1) == 'numign')then
					read(tt(iequal+1:ispace-1),*)ft%num_ignitions
					is_read(10) = 1
				endif
				tt = adjustl(tt(ispace+1:))
			endif
		enddo
	enddo
1345 continue
	close(ID_FILE_GRIDLIST)

	if(is_read(1)== 0) then
		print*,'Unspecified number of cells in the x-direction in gridlist.'
		call TerminateProgram()
	endif

	if(is_read(2)== 0) then
		print*,'Unspecified number of cells in the y-direction in gridlist.'
		call TerminateProgram()
	endif

	if(is_read(3)== 0) then
		print*,'Unspecified number of cells in the z-direction in gridlist.'
		call TerminateProgram()
	endif

	if(is_read(4)== 0) then
		print*,'Unspecified cell size in the x-direction in gridlist.'
		call TerminateProgram()
	endif

	if(is_read(5)== 0) then
		print*,'Unspecified cell size in the y-direction in gridlist.'
		call TerminateProgram()
	endif

	if(is_read(6)== 0) then
		print*,'zb not specified in gridlist. Invalid input file.'
		call TerminateProgram()
	endif

	if(is_read(7)== 0) then
		print*,'aa1 not specified in gridlist. Invalid input file.'
		call TerminateProgram()
	endif

	if(is_read(8)== 0) then
		ft%n_grid_top = ft%nz
	endif

	if(is_read(9)== 0) then
		ft%n_fuel_types = 1
	endif
	
	ft%dxi = 1. / ft%dx
	ft%dyi = 1. / ft%dy

	allocate(ft%zm(ft%nz),ft%z(ft%nz+1),ft%dz_array(ft%nz))
	aa3 = (1.-aa1)/(zb * real(ft%nz))**2	
	do i = 1,ft%nz
		ft%zm(i) = (real(i) - 0.5)*zb
		ft%zm(i) = aa3*ft%zm(i)**3 + aa1*ft%zm(i)
		ft%z(i) = (real(i) - 1.)*zb
		ft%z(i) = aa3*ft%z(i)**3 + aa1*ft%z(i)
		write (*,*) 'FIRETEC method',i,ft%z(i),ft%zm(i)
	enddo

	ft%z(1) = 0
	ft%z(2) = 2. * ft%zm(1)
	do i = 2,ft%nz
		ft%z(i+1) = (ft%zm(i) - ft%z(i)) * 2. + ft%z(i)
		write (*,*) 'midpoint method',i,ft%z(i)
	enddo
	
	do i = 1,ft%nz
		ft%dz_array(i) = ft%z(i+1) - ft%z(i)
	enddo

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()


	END
!===================================================================================
!===================================================================================
	SUBROUTINE ReadRasterOrigin()

	use fireca_module
	use file_handling_module
	use grid_module
	use interpolation_module
	
	implicit none

	real,parameter :: FIRETECH_GRID_ALIGNMENT_MARGIN = 1.! m, maximum difference between the southwest corners of the firetech and quic grids to consider them the same
	character(STR_LEN) :: fname								! N/A, file name

	fname = 'rasterorigin.txt'
	open(ID_FILE_RASTERORIGIN,file=TRIM(ft%file_path)// &
		TRIM(fileSeparator)//trim(fname),err=1101,status='old')
	read(ID_FILE_RASTERORIGIN,*,err=1102,end=1102) ft%utmx
	read(ID_FILE_RASTERORIGIN,*,err=1102,end=1102) ft%utmy
	close(ID_FILE_RASTERORIGIN)

	! Move SW corner if QUIC and firetech ones are close enough
	if(abs(ft%utmx - qugrid%utmx) <= FIRETECH_GRID_ALIGNMENT_MARGIN) ft%utmx = qugrid%utmx
	if(abs(ft%utmy - qugrid%utmy) <= FIRETECH_GRID_ALIGNMENT_MARGIN) ft%utmy = qugrid%utmy

	! Check if QUIC domain fits within the Firetech domain
	if(qugrid%utmx < ft%utmx)then
		write(msgoutfile,*) ' QUIC domain x-SW corner is outside the firetech domain'
		call TerminateProgram()
	endif
	if(qugrid%utmy < ft%utmy)then
		write(msgoutfile,*) ' QUIC domain y-SW corner is outside the firetech domain'
		call TerminateProgram()
	endif
	if(qugrid%utmx + real(qugrid%nx-1)*qugrid%dx > ft%utmx + real(ft%nx)*ft%dx)then
		write(msgoutfile,*) ' QUIC domain extends beyond firetech domain in the x-direction'		
		call TerminateProgram()
	endif
	if(qugrid%utmy + real(qugrid%ny-1)*qugrid%dy > ft%utmy + real(ft%ny)*ft%dy)then
		write(msgoutfile,*) ' QUIC domain extends beyond firetech domain in the y-direction'
		call TerminateProgram()
	endif

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()

	END
!===================================================================================
!===================================================================================
	subroutine ImportFTGrid()

	use grid_module

	implicit none

	call ReadGridlist()

	call ReadRasterOrigin()
	
	call DefineInterpolationArrays()

	end
!===================================================================================
!===================================================================================
	subroutine DetermineFuelTop(fuel_height_array)

	use constants
	use fireca_constants_module
	use grid_module

	implicit none

	real,intent(IN),dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, array of fuels heights
	real :: maxh
	integer :: k

	qugrid%kfire_top = -1
	qugrid%kfire_top_en2atmos = -1

	maxh = maxval(fuel_height_array)

	k = 1
	do while(k < qugrid%nz-1 .and. qugrid%kfire_top == -1)
		k = k + 1
		if(qugrid%z(k) >= maxh)then
			qugrid%kfire_top = k
		endif
	enddo

	if(qugrid%kfire_top == -1) then
		write(msgoutfile,*) 'Error in evaluating qugrid%kfire_top'
		call TerminateProgram()
	endif
	
	k = 1
	do while(k < qugrid%nz-1 .and. qugrid%kfire_top_en2atmos == -1)
		k = k + 1
		if(qugrid%z(k) >= firegrid%Lz + MAX_FLAME_HEIGHT)then
			qugrid%kfire_top_en2atmos = k
		endif
	enddo

	if(qugrid%kfire_top_en2atmos == -1) then
		write(msgoutfile,*) 'Error in evaluating qugrid%kfire_top_en2atmos'
		call TerminateProgram()
	endif

	end
!===================================================================================
!===================================================================================
	subroutine DetermineGridConversionParams()

	use constants
	use fireca_module
	use grid_module
	use interpolation_module

	implicit none

	integer,dimension(:),allocatable :: itemp			! N/A, used to reallocate arrays
	real,dimension(:),allocatable ::	rtemp				! N/A, used to reallocate
	real :: z1,z2, &
		one_over_xr,one_over_yr								! N/A, reciprocal: 1./xgrid_ratio_real
	integer :: i,j,k,ii,jj,kk,kqf,count,icount,ncell

	fca_2_qu_kstart = 0
	fca_2_qu_kend = 0

	allocate(fca2quic(qugrid%kfire_top_en2atmos)) ! qugrid%kfire_top =  max cell in the qu-domain that has a bit of a fireca cell in it
	allocate(kmap_start(firegrid%nz_en2atmos),kmap_end(firegrid%nz_en2atmos))
	kmap_start = qugrid%nz+10
	kmap_end = 0
	do k = 2,qugrid%kfire_top_en2atmos
		! Indexes of the FireCA grid cells in this QUIC-URB cell
		kqf = 1
		do while(qugrid%z(k-1) >= firegrid%z_en2atmos(kqf) .and. kqf < firegrid%nz_en2atmos)
			kqf = kqf + 1
		enddo
		kqf = kqf - 1
		fca_2_qu_kstart(k) = max(kqf,1)

		kqf = 1
		do while(qugrid%z(k) > firegrid%z_en2atmos(kqf) .and. kqf < firegrid%nz_en2atmos)
			kqf = kqf + 1
		enddo
		if(qugrid%z(k) > firegrid%z_en2atmos(firegrid%nz_en2atmos)) then
			fca_2_qu_kend(k) = firegrid%nz_en2atmos
		else
			kqf = max(kqf - 1,1)
			fca_2_qu_kend(k) = min(kqf,firegrid%nz_en2atmos)
		endif

		count = fca_2_qu_kend(k) - fca_2_qu_kstart(k) + 1
		allocate(fca2quic(k)%volRatio(count))
		do icount = 1,count
			! Portion of the FireCA cell in this QUIC-URB cell
			z1 = max(firegrid%z_en2atmos(fca_2_qu_kstart(k)+icount-1), qugrid%z(k-1))
			z2 = min(firegrid%z_en2atmos(fca_2_qu_kstart(k)+icount), qugrid%z(k))
			fca2quic(k)%volRatio(icount) = (z2-z1)/firegrid%dz_array_en2atmos(fca_2_qu_kstart(k)+icount-1)

			! In which QU cells is a certain FireCA cell
			kmap_start(fca_2_qu_kstart(k)+icount-1) = min(k,kmap_start(fca_2_qu_kstart(k)+icount-1))
			kmap_end(fca_2_qu_kstart(k)+icount-1) = max(k,kmap_end(fca_2_qu_kstart(k)+icount-1))
		enddo
	enddo

	do k = 1,firegrid%nz_en2atmos
		if(kmap_start(k) > kmap_end(k) .and. kmap_end(k) > 0) then
			write(msgoutfile,*) 'Error in evaluating kmap at level ',k
			call TerminateProgram()
		endif
	enddo

	! --- grid_conv

	allocate(grid_conv(qugrid%nx-1,qugrid%ny-1,qugrid%kfire_top))
	!$OMP parallel do private(i,j,k)
	do k = 1,qugrid%kfire_top
		do j = 1,qugrid%ny-1
			do i = 1,qugrid%nx-1
				grid_conv(i,j,k)%nelem = 0
			enddo
		enddo
	enddo
	!$OMP end parallel do
	one_over_xr = 1./xratio_real
	one_over_yr = 1./yratio_real
	do count = 1,firegrid%num_fuel_cells
		! Indexes in the FireCA grid
		ii = firegrid%ijk_cell_index(count,1)
		jj = firegrid%ijk_cell_index(count,2)
		kk = firegrid%ijk_cell_index(count,3)
		! Indexes in the QUIC-URB grid
		i = ceiling(real(ii)*one_over_xr)
		j = ceiling(real(jj)*one_over_yr)
		if(kmap_start(kk) <= qugrid%kfire_top)then
			
			do k = kmap_start(kk),min(kmap_end(kk),qugrid%kfire_top)   ! in which QU cells is a certain FireCA cell (vertical dir)
				if(grid_conv(i,j,k)%nelem == 0)then
					ncell = (fca_2_qu_kend(k) - fca_2_qu_kstart(k) + 1) * xratio_int * yratio_int
					allocate(grid_conv(i,j,k)%en_ratio(ncell),grid_conv(i,j,k)%index(ncell))
				endif
				grid_conv(i,j,k)%nelem = grid_conv(i,j,k)%nelem + 1
				grid_conv(i,j,k)%index(grid_conv(i,j,k)%nelem) = count  !index of the fireca cell in the QU cell
				z1 = max(firegrid%z(kk), qugrid%z(k-1))
				z2 = min(firegrid%z(kk+1), qugrid%z(k))
				grid_conv(i,j,k)%en_ratio(grid_conv(i,j,k)%nelem) = (z2-z1)/firegrid%dz_array(kk)
			enddo
		endif
	enddo

	allocate(rtemp(qugrid%nz),itemp(qugrid%nz))
	do k = 2,qugrid%kfire_top
		count = (fca_2_qu_kend(k) - fca_2_qu_kstart(k) + 1) * xratio_int * yratio_int
		do j = 1,qugrid%ny-1
			do i = 1,qugrid%nx-1
				if(grid_conv(i,j,k)%nelem < count .and. grid_conv(i,j,k)%nelem	> 0)then
					rtemp(1:grid_conv(i,j,k)%nelem) = grid_conv(i,j,k)%en_ratio(1:grid_conv(i,j,k)%nelem)
					itemp(1:grid_conv(i,j,k)%nelem) = grid_conv(i,j,k)%index(1:grid_conv(i,j,k)%nelem)
					deallocate(grid_conv(i,j,k)%en_ratio,grid_conv(i,j,k)%index)
					allocate(grid_conv(i,j,k)%en_ratio(grid_conv(i,j,k)%nelem),grid_conv(i,j,k)%index(grid_conv(i,j,k)%nelem))
					grid_conv(i,j,k)%en_ratio = rtemp(1:grid_conv(i,j,k)%nelem)
					grid_conv(i,j,k)%index = itemp(1:grid_conv(i,j,k)%nelem)
				endif
			enddo
		enddo
	enddo
	deallocate(rtemp,itemp)

	! Fuel in the QUIC world
	call InitFuelQUIC()

	
	! To map energy to atmosphere
	allocate(en2atm_2_qu_kstart(qugrid%nz), en2atm_2_qu_kend(qugrid%nz))
	en2atm_2_qu_kstart = 0
	en2atm_2_qu_kend = 0
	
	do k = 2,qugrid%kfire_top_en2atmos
		! Start
		kqf = 1		
		do while(.not.(firegrid%zm_en2atmos(kqf) < qugrid%z(k) .and. &
			firegrid%zm_en2atmos(kqf) > qugrid%z(k-1) ) .and. kqf < firegrid%nz_en2atmos)
			kqf = kqf + 1
		enddo
		if(kqf > firegrid%nz_en2atmos)then
			print*, 'Error while computing en2atm_2_qu_kstart'
			call TerminateProgram()
		endif
		en2atm_2_qu_kstart(k) = kqf
		! End
		kqf = firegrid%nz_en2atmos
		do while(.not.(firegrid%zm_en2atmos(kqf) < qugrid%z(k) .and. &
			firegrid%zm_en2atmos(kqf) > qugrid%z(k-1) ) .and. kqf > 1)
			kqf = kqf - 1
		enddo
		if(kqf < 1)then
			print*, 'Error while computing en2atm_2_qu_kend'
			call TerminateProgram()
		endif
		en2atm_2_qu_kend(k) = kqf

	enddo

	end
!===================================================================================
!===================================================================================
	subroutine InitLoopIndex()

	use plume_module
	use grid_module

	implicit none

	loopidx%istart = (/1,qugrid%nx-1/)
	loopidx%jstart = (/1,qugrid%ny-1/)
	loopidx%kstart = (/2,qugrid%kfire_top_en2atmos/)
	loopidx%iend = (/qugrid%nx-1,1/)
	loopidx%jend = (/qugrid%ny-1,1/)
	loopidx%kend = (/qugrid%kfire_top_en2atmos,2/)

	loopidx%increm = (/1,-1/)

	loopidx%list(1,:) = (/1,1,1/)
	loopidx%list(2,:) = (/2,1,1/)
	loopidx%list(3,:) = (/1,2,1/)
	loopidx%list(4,:) = (/1,1,2/)
	loopidx%list(5,:) = (/2,2,1/)
	loopidx%list(6,:) = (/1,2,2/)
	loopidx%list(7,:) = (/2,1,2/)
	loopidx%list(8,:) = (/2,2,2/)

	end
!===================================================================================
!===================================================================================
	subroutine InitInterpolationVars

	use constants
	use interpolation_module
	use grid_module

	implicit none

	integer :: found,kuv,kw,k, &
		extra_x,extra_y		! N/A, indexes to interpolate between the quic and fireca grids


	! Variables for interpolating the winds
	! Interpolate on the terrain following grid
	allocate(qu_2fca_kuv(firegrid%nz+1),qu_2fca_w(firegrid%nz+1))

	do k = 1,firegrid%nz+1

		!  - u,v (face centered, middle of the cell)
		found = 0
		kuv = 0
		do while(found == 0 .and. kuv < qugrid%nz-1)
			kuv = kuv + 1
			if(qugrid%zm(kuv) >= firegrid%zm(k)) then
				found = 1
			endif
		enddo
		if(found == 0)then
			write(msgoutfile,*)'Error in the calculation of interpolation indexes for the wind field (u,v).'
			call TerminateProgram()
		endif
		qu_2fca_kuv(k) = kuv

		! - w (face centered, bottom of the cell)
		found = 0
		kw = 2  ! w(1) => -z(1); w(2) => 0   => kw = 1 or 2 are not valid
		do while(found == 0 .and. kw < qugrid%nz)
			kw = kw + 1
			if(qugrid%z(kw-1)	>= firegrid%zm(k)) then
				found = 1
			endif
		enddo

		if(found == 0)then
			write(msgoutfile,*)'Error in the calculation of interpolation indexes for the wind field (w).'
			call TerminateProgram()
		endif
		qu_2fca_w(k) = kw
	enddo

	if(mod(xratio_real,2.) == 0)then
		extra_x = 1
	else
		extra_x = 0
	endif
	if(mod(yratio_real,2.) == 0)then
		extra_y = 1
	else
		extra_y = 0
	endif

	istart_interp = int(ceiling(xratio_real * 0.5)) + extra_x
	jstart_interp = int(ceiling(yratio_real * 0.5)) + extra_y
		
	iend_interp = firegrid%nx - istart_interp + 1
	jend_interp = firegrid%ny - jstart_interp + 1

	end
!===================================================================================
!===================================================================================
	subroutine InitIntegratedFuelMass

	use fireca_module
	use grid_module
	use flags_module

	implicit none

	integer :: k

	if(flag%int_massburnt_file == 1)then
		allocate(mass_int_0(firegrid%nx,firegrid%ny))
		mass_int_0 = 0.
		!$OMP parallel do private(k)
		do k = 1,firegrid%num_fuel_cells
			mass_int_0(firegrid%ijk_cell_index(k,1),firegrid%ijk_cell_index(k,2)) =		&
				mass_int_0(firegrid%ijk_cell_index(k,1),firegrid%ijk_cell_index(k,2)) +	&
				fuels%density_initial(k) * firegrid%dz_array(firegrid%ijk_cell_index(k,3))
		enddo
		!$OMP end parallel do
	endif

	end
!===================================================================================
!===================================================================================
	subroutine initialize_fire()

	use rnd_num_vars
	use constants
	use fireca_module
	use time_module
	use fuel_module
	use plume_module
	use sor_module	
	use grid_module
	
	implicit none

	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, array of fuels heights

	print*,'Initializing fire variables'

	! -------- Constants
	nbtime = floor(BURNOUT_TIME/fcatime%dt_real)
	nbtime_real = real(nbtime)
	partical_burn_vel_dt = PARTICAL_BURN_VEL * fcatime%dt_real

	! -------- Firebrands
	fb%time = 0.
	
	! -------- Initialize fuel and ignitions
	if(fb%flag == 1)then
		call ReadFBFile()
	endif

	print *, "specifying Fuel "
	! -- Specify density and moisture
	call SpecifyFuel(fuel_height_array)

	print *, "Fuel specified"
	fuels%l_patch = 0.05    !rrl temporary.  These two lines should be moved once we figure out the details of how they will be put in.
	fuels%l_gap = 0.0			! rrl as stated above

	! -- Specify ignitions
	call SpecifyIgnitions(fuel_height_array)

	! -- Fuel top
	call DetermineFuelTop(fuel_height_array)

	! -- Grid conversion params
	call DetermineGridConversionParams()

	allocate(sor%alpha2_fire(qugrid%nx,qugrid%ny,qugrid%nz))
	sor%alpha2_fire = sor%alpha2 ! Initialized to the default, then adjusted where there is fire

	call random_index(8,iharvest)
	loopidx%list_count = iharvest(1)

	call InitLoopIndex()

	call InitInterpolationVars()

	call InitIntegratedFuelMass()
	
	call CheckFuelArrays()

	return

	end
!===================================================================================
!===================================================================================
	subroutine InitVarFire(i,j,k,n_ignitions)

	use file_handling_module
	use constants
	use grid_module
	use fireca_module
	use plume_const_module

	implicit none

	integer,intent(IN) ::		&
		i,j,k,						& ! N/A, locations in a 3D grid
		n_ignitions					  ! N/A, number of ignitions
	integer :: ncell				  ! N/A, position in the linear arrays of the i,j,k cell

	ncell = grid_where(i,j,k)

	if(ncell > 0) then
		fire%ignitions(ncell) = n_ignitions
		fuels%moisture(ncell) = 0.0
		fire%burncenter_info(ncell)%val2(1,1:2) = -1.
		fire%burncenter(ncell)%val = 0.5
		!this should tie the spatial variance to an actual physical size for a rectangle with a side of init_var_length (m)
		fire%spat_var(ncell)%val(1) = init_var_length*ONE_OVER_ROOTTHREE * min(1., firegrid%dxi)
		fire%spat_var(ncell)%val(2) = init_var_length*ONE_OVER_ROOTTHREE * min(1., firegrid%dyi)
		fire%reaction_rate(ncell) = IGNITION_REACTION_RATE
		
		! Add energy to cells with ignitions to generate w
		fire%energy_to_atmos(i,j,k) = real(n_ignitions) * ENERGY_PER_IGNITION * &
			plume_const%init_w_en_fract / firegrid%cellvol(k) ! kW/m3
		
		write(ID_FILE_IGNITE_SEL)i,j,k
	endif

	end
!===================================================================================
!===================================================================================
	subroutine LineFireInition(fire_source_xo,fire_source_yo, &
		fire_source_xlen,fire_source_ylen,fuel_height_array,n_ignitions)

	use grid_module
	
	implicit none

	real,intent(IN), dimension(firegrid%nx,firegrid%ny) :: &
		fuel_height_array								  ! m
	real,intent(IN) ::								&
		fire_source_xo,fire_source_yo,			& ! m, line source sw corner
		fire_source_xlen,fire_source_ylen		  ! m, length of the fire line in the x and y directions
	integer,intent(IN) :: n_ignitions			  ! N/A, number of ignition per cell
	integer :: i,j,k,is,js,ie,je

	is = ceiling(fire_source_xo * firegrid%dxi)
	if(mod(fire_source_xo,firegrid%dx) == 0) is = is + 1
	ie = ceiling((fire_source_xo+fire_source_xlen) * firegrid%dxi)
	js = ceiling(fire_source_yo * firegrid%dyi)
	if(mod(fire_source_yo,firegrid%dy) == 0) js = js + 1
	je = ceiling((fire_source_yo+fire_source_ylen) * firegrid%dyi)

	do j = js,je
		do i = is,ie
			k = 1
			do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
				call InitVarFire(i,j,k,n_ignitions)
				k = k + 1
			enddo
		enddo
	enddo

	end
!===================================================================================
!===================================================================================
	subroutine SquareRingIgnition(fire_source_xo,fire_source_yo,	&
		fire_source_xlen,fire_source_ylen,									&
		fire_source_xwidth,fire_source_ywidth,								&
		fuel_height_array,n_ignitions)

	use grid_module

	implicit none

	real,intent(IN), dimension(firegrid%nx,firegrid%ny) :: &
		fuel_height_array								     ! m
	real,intent(IN) ::									&
		fire_source_xo,fire_source_yo,				& ! m, ring sw corner
		fire_source_xwidth,fire_source_ywidth,		& ! m, ring width
		fire_source_xlen,fire_source_ylen		     ! m, length of the fire line in the x and y directions
	integer,intent(IN) :: n_ignitions			     ! N/A, number of ignition per cell
	integer ::	&
		idelta,jdelta,			& ! N/A, number of cell spanned by the ring in the x and y directions
		i,j,k,is,js,ie,je

	idelta = ceiling(fire_source_xwidth * firegrid%dxi)
	jdelta = ceiling(fire_source_ywidth * firegrid%dyi)
	is = ceiling(fire_source_xo * firegrid%dxi)
	if(mod(fire_source_xo,firegrid%dx) == 0) is = is + 1
	ie = ceiling((fire_source_xo+fire_source_xlen) * firegrid%dxi)
	js = ceiling(fire_source_yo * firegrid%dyi)
	if(mod(fire_source_yo,firegrid%dy) == 0) js = js + 1
	je = ceiling((fire_source_yo+fire_source_ylen) * firegrid%dyi)

	do i = is,ie
		! bottom
		do j = js,js+jdelta-1
			k = 1
			do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
				call InitVarFire(i,j,k,n_ignitions)
				k = k + 1
			enddo
		enddo
		! top
		do j = je-jdelta+1,je
			k = 1
			do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
				call InitVarFire(i,j,k,n_ignitions)
				k = k + 1
			enddo
		enddo
	enddo

	do j = js,je
		! right
		do i = is,is+idelta-1
			k = 1
			do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
				call InitVarFire(i,j,k,n_ignitions)
				k = k + 1
			enddo
		enddo
		! left
		do i = ie-idelta+1,ie
			k = 1
			do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
				call InitVarFire(i,j,k,n_ignitions)
				k = k + 1
			enddo
		enddo
	enddo

	end
!===================================================================================
!===================================================================================
	subroutine CircularRingIgnition(fire_source_xo,fire_source_yo,	&
		fire_source_xlen,fire_source_ylen,fire_source_xwidth,			&
		fuel_height_array,n_ignitions)

	use grid_module

	implicit none

	real,intent(IN), dimension(firegrid%nx,firegrid%ny) :: &
		fuel_height_array								     ! m
	real,intent(IN) ::									&
		fire_source_xo,fire_source_yo,				& ! m, ring sw corner
		fire_source_xwidth,								& ! m, ring width
		fire_source_xlen,fire_source_ylen			  ! m, length of the fire line in the x and y directions
	real ::													&
		radius,											& ! m, radius of the ring ignition pattern
		xringcenter,yringcenter,					& ! m, coordinates of the center of the ring ignition pattern
		x,y,												& ! m, cell center coordinates
		dist												  ! m, distance between a cell center and the center of the ring ignition pattern
	integer,intent(IN) :: n_ignitions			     ! N/A, number of ignition per cell
	integer :: i,j,k,is,js,ie,je

	is = ceiling(fire_source_xo * firegrid%dxi)
	if(mod(fire_source_xo,firegrid%dx) == 0) is = is + 1
	ie = ceiling((fire_source_xo+fire_source_xlen)*firegrid%dxi)
	js = ceiling(fire_source_yo * firegrid%dyi)
	if(mod(fire_source_yo,firegrid%dy) == 0) js = js + 1
	je = ceiling((fire_source_yo+fire_source_ylen)*firegrid%dyi)

	radius = fire_source_xlen * 0.5
	xringcenter = fire_source_xo + radius
	yringcenter = fire_source_yo + radius
	do j = js,je
		y = firegrid%ycenters(j)
		do i = is,ie
			x = firegrid%xcenters(i)
			dist = sqrt( (x-xringcenter)**2 + (y-yringcenter)**2 )
			if(dist <= radius .and. dist >= radius - fire_source_xwidth) then
				k = 1
				do while(fuel_height_array(i,j) > firegrid%z(k) .and. k <= firegrid%nz)
					call InitVarFire(i,j,k,n_ignitions)
					k = k + 1
				enddo
			endif
		enddo
	enddo

	end
!===================================================================================
!===================================================================================
	subroutine UseQFIgnitions()

	use constants
	use grid_module
	use file_handling_module
	
	implicit none

	character(STR_LEN) :: fname							  ! N/A, file name
	integer,dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: ignitions_array	! N/A, number of ignitions per cell
	integer :: index,i,j,k


	fname = 'QF_Ignitions.inp'
	open(ID_FILE_QF_IGNITION,file=TRIM(workingDirectory)// &
		TRIM(fileSeparator)//trim(fname), &
		status = 'old',form='unformatted',err=1101)
	read(ID_FILE_QF_IGNITION,err=1102,end=1102) (((ignitions_array(i,j,k),i=1,firegrid%nx),j=1,firegrid%ny),k=1,firegrid%nz)
	close(ID_FILE_QF_IGNITION)

	do index = 1,firegrid%num_fuel_cells
		i = firegrid%ijk_cell_index(index,1)
		j = firegrid%ijk_cell_index(index,2)
		k = firegrid%ijk_cell_index(index,3)
		if(ignitions_array(i,j,k) > 0)then
			call InitVarFire(i,j,k,ignitions_array(i,j,k))
		endif
	enddo

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	subroutine FiretechIgnition()

	use rnd_num_vars
	use fireca_module
	use file_handling_module
	use ignitions_module
	use grid_module

	implicit none

	character(STR_LEN) :: fname						  ! N/A, file name
	real,allocatable,dimension(:,:) ::	prb_ign_thin ! N/A, [0,1] probability of ignition
	real ::					&
		dxx,dyy,				& ! m, difference in the sw corner of firetech and quic grids
		loc
	integer ::				&
		i,j,k,				&
		ii,jj,				&
		index,				&
		is,ie,js,je,		&
		ntot_active			  ! N/A, total number of ignitions in the Firetech input file


	allocate(prb_ign_thin(5,ft%num_ignitions))

	! ignition pattern
	fname = 'ignite.dat'   ! columns are (1) I, (2) J, (3) K, (4) prob ignition thin fuel, (5) prob ignition thick fuel
	open(ID_FILE_IGNITE,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
		status = 'old',form='unformatted', access='stream',err=1101)
	read(ID_FILE_IGNITE) prb_ign_thin
	close(ID_FILE_IGNITE)

	dxx = qugrid%utmx - ft%utmx
	dyy = qugrid%utmy - ft%utmy
	ntot_active = 0
	do i = 1,ft%num_ignitions
		call random_number(harvest,stat(1))
		if(harvest(1) + prb_ign_thin(4,i) > 1.) then
			loc = ((prb_ign_thin(1,i)-1)*ft%dx-dxx)
			is = ceiling( loc * firegrid%dxi)
			if(mod(loc,firegrid%dx) == 0.) is = is + 1
			is = max(is,1)
			loc = prb_ign_thin(1,i)*ft%dx-dxx
			ie = ceiling( loc * firegrid%dxi)
			ie = min(ie,firegrid%nx)
			loc =  ((prb_ign_thin(2,i)-1)*ft%dy-dyy)
			js = ceiling( loc * firegrid%dyi)
			if(mod(loc,firegrid%dy) == 0.) js = js + 1
			js = max(js,1)
			loc = prb_ign_thin(2,i)*ft%dy-dyy
			je = ceiling( loc * firegrid%dyi)
			je = min(je,firegrid%ny)

			j = int(prb_ign_thin(3,i))
			do k = 1,firegrid%nz
				if(ft%z(j+1) - firegrid%z(k) > 0.01  .and. firegrid%z(k+1) - ft%z(j) > 0.01)then
					do ii = is,ie
						do jj = js,je
							index = grid_where(ii,jj,k)
							if(index > 0) then
								if(fire%ignitions(index) <= 0)then
									ntot_active = 1
									call InitVarFire(ii,jj,k,fire_ignition%num_ign_percell)									
								endif
							endif
						enddo
					enddo
				endif
			enddo
		endif
	enddo
	
	if(ntot_active == 0) then
		write(msgoutfile,*)'No ignitions in the QUIC-Fire domain'
		call TerminateProgram()
	else
		print*,'Number of cells with ignitions = ',count(fire%ignitions > 0)
	endif

	deallocate(prb_ign_thin)

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	subroutine InitFuelAndFireVars()

	use fireca_module
	use grid_module

	implicit none

	real ::											&
		center_to_grid_dist,						&  ! m, max(fire%burncenter(kk,1),1.-fire%burncenter(kk,1))
		x_val_p1,x_val_m1,y_val_p1,y_val_m1    ! N/A, used to sample around the cell
	integer ::		&
		kk,n1,n2,	&		! N/A, used to sample around the cell
		i,j,k

	!$OMP parallel do private(i,j,k,kk,n1,n2,x_val_p1,x_val_m1,center_to_grid_dist, &
	!$OMP y_val_p1,y_val_m1)
	do kk = 1,firegrid%num_fuel_cells
	if(fire%ignitions(kk) > 0)then

		i = firegrid%ijk_cell_index(kk,1)
		j = firegrid%ijk_cell_index(kk,2)
		k = firegrid%ijk_cell_index(kk,3)

		fuels%moisture(kk) = 0.
		fire%burncenter_info(kk)%val2(1,1:2) = -1  !signifying initial burn location

		n1 = grid_where(i+1,j,k)
		n2 = grid_where(i-1,j,k)
		if(n1 > 0)then
			x_val_p1 = real(fire%ignitions(n1))
		else
			x_val_p1 = 0.
		endif
		if(n2 > 0)then
			x_val_m1 = real(fire%ignitions(n2))
		else
			x_val_m1 = 0.
		endif

		fire%burncenter(kk)%val(1) = (0.5+x_val_p1)/ &
			(1.+x_val_p1+x_val_m1)

		center_to_grid_dist = min(fire%burncenter(kk)%val(1),1.-fire%burncenter(kk)%val(1))

		fire%spat_var(kk)%val(1) = min(min(ONE_OVER_ROOTTHREE*center_to_grid_dist,  &
			max(0.5, 0.5*(x_val_p1+x_val_m1))), &
			ONE_OVER_ROOTTHREE*firegrid%dxi)

		n1 = grid_where(i,j+1,k)
		n2 = grid_where(i,j-1,k)
		if(n1 > 0)then
			y_val_p1 = real(fire%ignitions(n1))
		else
			y_val_p1 = 0.
		endif
		if(n2 > 0)then
			y_val_m1 = real(fire%ignitions(n2))
		else
			y_val_m1 = 0.
		endif

		fire%burncenter(kk)%val(2) = (0.5+y_val_p1)/ &
			(1.+y_val_p1+y_val_m1)

		center_to_grid_dist = min(fire%burncenter(kk)%val(2),1-fire%burncenter(kk)%val(2))

		fire%spat_var(kk)%val(2) = min(min(ONE_OVER_ROOTTHREE*center_to_grid_dist, max(.5, &
			0.5*(y_val_p1+y_val_m1))), &
			ONE_OVER_ROOTTHREE*firegrid%dyi)

		fire%reaction_rate(kk) = IGNITION_REACTION_RATE
	endif
	enddo
	!$OMP end parallel do

	max_init_ign = maxval(fire%ignitions)

	end
!===================================================================================
!===================================================================================
	subroutine SpecifyIgnitions(fuel_height_array)

	use grid_module
	use ignitions_module
	use time_module

	implicit none
	
	real,dimension(firegrid%nx,firegrid%ny),intent(IN) :: fuel_height_array		! m, height of the fuel at a certain (x,y) location

	if(fire_ignition%flag == 1)then
		call LineFireInition(fire_ignition%xo,fire_ignition%yo, 			&
		fire_ignition%xlen,fire_ignition%ylen,fuel_height_array,			&
		fire_ignition%num_ign_percell)

	elseif(fire_ignition%flag == 2)then
		call SquareRingIgnition(fire_ignition%xo,fire_ignition%yo, 	&
			fire_ignition%xlen,fire_ignition%ylen,								&
			fire_ignition%xwidth,fire_ignition%ywidth,						&
			fuel_height_array,fire_ignition%num_ign_percell)

	elseif(fire_ignition%flag == 3) then
		call CircularRingIgnition(fire_ignition%xo,fire_ignition%yo,	&
			fire_ignition%xlen,fire_ignition%ylen,fire_ignition%xwidth,	&
			fuel_height_array,fire_ignition%num_ign_percell)

	elseif(fire_ignition%flag == 4) then
		call UseQFIgnitions()

	elseif(fire_ignition%flag == 5) then
		call update_ignitions_pattern()

	elseif(fire_ignition%flag == 6) then
		call FiretechIgnition()

	endif

	call InitFuelAndFireVars()

	return

	end
!===================================================================================
!===================================================================================
	subroutine ReadFuelHeight(fuel_density_array,fuel_height_array)

	use grid_module
	use fireca_module
	use constants
	use file_handling_module

	implicit none

	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: fuel_density_array	! kg/m3, density of the fuel at each cell
	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array					! m, height of the fuel at a certain (x,y) location
	character(STR_LEN) :: fname																! N/A, file name
	real,dimension(:,:),allocatable :: temp_h												! m, height of the fuel at a certain (x,y) location in the firetech grid
	real,dimension(:,:,:,:),allocatable :: temp4											! m, used to read in the Firetec files
	real,dimension(:,:,:),allocatable :: temp3											! m, used to read in the Firetec files
	integer :: i,j,k
	
	allocate(fuels%actual_height(firegrid%nx,firegrid%ny))
	
	fuel_height_array = 0.
	
	! Fuel height
	if(fuels%height_flag == 1)then
		! Constant fuel height (the borders don't have fuel)
		fuel_height_array(2:firegrid%nx-1, 2:firegrid%ny) = fuels%input_height
		fuels%actual_height = min(fuels%input_height, firegrid%dz_array(1))

	elseif(fuels%height_flag == 2)then
		fname = 'QF_FuelHeight.inp'
		open(ID_FILE_QF_FUELHEIGHT,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
			status = 'old',form='unformatted',err=1101)
		read(ID_FILE_QF_FUELHEIGHT,err=1102,end=1102) ((fuel_height_array(i,j),i=1,firegrid%nx),j=1,firegrid%ny)
		close(ID_FILE_QF_FUELHEIGHT)
		
		call GetGroundFuelHeight(fuel_height_array)
		
	elseif(fuels%height_flag == 3)then
		fname = 'treesfueldepth.dat'
		open(ID_FILE_TREESFUELDEPTH,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
			status = 'old',form='unformatted', access='stream',err=1101)
		read(ID_FILE_TREESFUELDEPTH) fuel_height_array
		close(ID_FILE_TREESFUELDEPTH)
		
		call GetGroundFuelHeight(fuel_height_array)
		
	else
		! fuel_height_array is determined from the fuel density
		if(ft%fuel_file_type == 1)then
			allocate(temp4(ft%n_fuel_types,ft%nx,ft%ny,ft%n_grid_top), temp_h(ft%nx,ft%ny))

			fname = 'treesfueldepth.dat'
			open(ID_FILE_TREESFUELDEPTH,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESFUELDEPTH,err=1102,end=1102) temp4
			close(ID_FILE_TREESFUELDEPTH)

			! Interpolate ground layer, only thin fuel
			temp_h = temp4(1,:,:,1)
			call AverageFiretechFuelDepth(fuels%actual_height,temp_h)
			deallocate(temp4, temp_h)
		else
			allocate(temp3(ft%nx,ft%ny,ft%n_grid_top))

			fname = 'treesfueldepth_1.dat'
			open(ID_FILE_TREESFUELDEPTH,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESFUELDEPTH,err=1102,end=1102) temp3
			close(ID_FILE_TREESFUELDEPTH)
		
			! Interpolate ground layer, only thin fuel			
			call AverageFiretechFuelDepth(fuels%actual_height, temp3(:,:,1))
			deallocate(temp3)
		endif
			
	endif

	if(fuels%density_flag == 1) then
		!$OMP parallel do private(i,j,k)
		do j = 1,firegrid%ny
			do i = 1,firegrid%nx
				do k = 1,firegrid%nz
					if(fuel_height_array(i,j) > firegrid%z(k+1)) then
						fuel_density_array(i,j,k) = fuels%input_density		! kg/m^3
					elseif(fuel_height_array(i,j) > firegrid%z(k) .and. fuel_height_array(i,j) <= firegrid%z(k+1)) then
						! Reduce bulk density to compute the correct mass
						fuel_density_array(i,j,k) = fuels%input_density * &
							(fuel_height_array(i,j) - firegrid%z(k)) / firegrid%dz_array(k)		! kg/m^3
					endif
				enddo
			enddo
		enddo
		!$OMP end parallel do 
		fuel_density_array=3.0*fuel_density_array
	endif
	
	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	subroutine GetGroundFuelHeight(fuel_height_array)
	
	use grid_module
	use fireca_constants_module
	use fireca_module
	
	implicit none
	
	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, height of the fuel at a certain (x,y) location
	integer :: i,j
	
	!$OMP parallel do private(i, j)
	do j = 1, firegrid%ny
		do i = 1, firegrid%nx
			fuels%actual_height(i,j) = min(fuel_height_array(i,j), firegrid%dz_array(1))			
		enddo
	enddo
	!$OMP end parallel do 
	
	end
!===================================================================================
!===================================================================================	
	subroutine ReadFuelDensity(fuel_density_array, fuel_density_array_thick, fuel_height_array)
	
	use fireca_module
	use file_handling_module
	use grid_module
	
	implicit none

	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz) ::	&
		fuel_density_array,												& ! kg/m3, density of the fuel at each cell
		fuel_density_array_thick										  ! kg/m3, density of the thick fuel at each cell
	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, height of the fuel at a certain (x,y) location
	character(STR_LEN) :: fname			! N/A, file name
	real,dimension(:,:,:),allocatable :: temp_dens, temp_ft, temp3
	real,dimension(:,:,:,:),allocatable :: temp4	
	integer :: i,j,k,ii

	! The fuel height is needed only in case the density is constant throughout the domain
	if(fuels%density_flag == 1) then
		call ReadFuelHeight(fuel_density_array,fuel_height_array)

	elseif(fuels%density_flag == 2)then
		fname = 'QF_FuelDensity.inp'
		open(ID_FILE_QF_FUELDENSITY,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
			status = 'old',form='unformatted',err=1101)
		read(ID_FILE_QF_FUELDENSITY,err=1102,end=1102) &
			(((fuel_density_array(i,j,k),i=1,firegrid%nx),j=1,firegrid%ny),k=1,firegrid%nz)
		close(ID_FILE_QF_FUELDENSITY)

	elseif(fuels%density_flag == 3)then
		allocate(temp4(1,firegrid%nx,firegrid%ny,firegrid%nz))

		fname = 'treesrhof.dat'
		open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
			status = 'old',form='unformatted', access='stream',err=1101)
		read(ID_FILE_TREESRHOF,err=1102,end=1102) temp4
		close(ID_FILE_TREESRHOF)
		
		fuel_density_array = temp4(1,:,:,:)

		deallocate(temp4)
	elseif(fuels%density_flag == 66)then
		call ReadFuelHeight(fuel_density_array, fuel_height_array)
		
		if(ft%fuel_file_type == 1)then
			allocate(temp4(ft%n_fuel_types,ft%nx,ft%ny,ft%n_grid_top), &
				temp_ft(ft%nx,ft%ny,ft%nz), &
				temp_dens(firegrid%nx,firegrid%ny,firegrid%nz))

			fname = 'treesrhof.dat'
			open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESRHOF,err=1102,end=1102) temp4 ! kg/m3
			close(ID_FILE_TREESRHOF)
			write (*,*) 'before loop'
			do ii = 1,1 !ft%n_fuel_types
				write (*,*) 'before fuel density'
				temp4(ii,:,:,2:13)=temp4(ii,:,:,2:13)/10. !RRL /10. 66 for estimated correction of Spring Hill Creek data provided to LANL 8/2020 from TT
				temp4(ii,:,:,1:1)=temp4(ii,:,:,1:1)*.666    !  RRL test of half surface
				!temp4(ii,550:560,:,:)=1.e-6   !RRL road cuts
				!temp4(ii,40:50,:,:)=1.e-6     !RRL road cuts
				!temp4(ii,:,550:560,:)=1.e-6   !RRL road cuts
				!temp4(ii,:,40:50,:)=1.e-6   !RRL road cuts
				write (*,*) 'before fuel density'
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(ii,:,:,:)
				call InterpolateFiretechFile(temp_dens, temp_ft)

				if(ii == 1)then
					fuel_density_array = temp_dens
				endif
			enddo
			
			if(ft%n_fuel_types > 1)then
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(3,:,:,:)
				call InterpolateFiretechFile(fuel_density_array_thick, temp_ft)				
			endif
			
			deallocate(temp4,temp_dens,temp_ft)
		else
			allocate(temp3(ft%nx,ft%ny,ft%n_grid_top), temp_ft(ft%nx,ft%ny,ft%nz))

			fname = 'treesrhof_1.dat'
			open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESRHOF,err=1102,end=1102) temp3 ! kg/m3
			close(ID_FILE_TREESRHOF)
			
			temp_ft = 0.
			temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
			call InterpolateFiretechFile(fuel_density_array, temp_ft)
			
			if(ft%n_fuel_types > 1)then
				fname = 'treesrhof_3.dat'
				open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESRHOF,err=1102,end=1102) temp3 ! kg/m3
				close(ID_FILE_TREESRHOF)
			
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
				call InterpolateFiretechFile(fuel_density_array_thick, temp_ft)			
			endif
			
			deallocate(temp3, temp_ft)
		endif
	else		
		
		call ReadFuelHeight(fuel_density_array, fuel_height_array)
		
		if(ft%fuel_file_type == 1)then
			allocate(temp4(ft%n_fuel_types,ft%nx,ft%ny,ft%n_grid_top), &
				temp_ft(ft%nx,ft%ny,ft%nz), &
				temp_dens(firegrid%nx,firegrid%ny,firegrid%nz))

			fname = 'treesrhof.dat'
			open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESRHOF,err=1102,end=1102) temp4 ! kg/m3
			print *,"shape rhof",shape(temp4)
			close(ID_FILE_TREESRHOF)

			do ii = 1,1 !ft%n_fuel_types
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(ii,:,:,:)
				call InterpolateFiretechFile(temp_dens, temp_ft)

				if(ii == 1)then
					fuel_density_array = temp_dens
				endif
			enddo
			
			if(ft%n_fuel_types > 1)then
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(3,:,:,:)
				call InterpolateFiretechFile(fuel_density_array_thick, temp_ft)				
			endif
			
			deallocate(temp4,temp_dens,temp_ft)
		else
			allocate(temp3(ft%nx,ft%ny,ft%n_grid_top), temp_ft(ft%nx,ft%ny,ft%nz))

			fname = 'treesrhof_1.dat'
			open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESRHOF,err=1102,end=1102) temp3 ! kg/m3
			close(ID_FILE_TREESRHOF)
			
			temp_ft = 0.
			temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
			call InterpolateFiretechFile(fuel_density_array, temp_ft)
			
			if(ft%n_fuel_types > 1)then
				fname = 'treesrhof_3.dat'
				open(ID_FILE_TREESRHOF,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESRHOF,err=1102,end=1102) temp3 ! kg/m3
				close(ID_FILE_TREESRHOF)
			
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
				call InterpolateFiretechFile(fuel_density_array_thick, temp_ft)			
			endif
			
			deallocate(temp3, temp_ft)
		endif
		
	endif
	
	! Eliminate borders
	fuel_density_array(1,:,:) = 0
	fuel_density_array(firegrid%nx,:,:) = 0
	fuel_density_array(:,1,:) = 0
	fuel_density_array(:,firegrid%ny,:) = 0

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	SUBROUTINE DefineUrban(fuel_dens_thick)
	
	use grid_module
	use fireca_module
	
	implicit none
	
	real,intent(IN), dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: fuel_dens_thick
	integer :: index, i, j, k
	
	allocate(fuels%is_urban(firegrid%num_fuel_cells))
	fuels%is_urban = 0
	!$OMP parallel do private(index, i, j, k)
	do index = 1, firegrid%num_fuel_cells
		i = firegrid%ijk_cell_index(index,1)
		j = firegrid%ijk_cell_index(index,2)
		k = firegrid%ijk_cell_index(index,3)
		
		if(fuel_dens_thick(i,j,k) > MIN_FUEL_DENSITY) fuels%is_urban(index) = 1
	enddo
	!$OMP end parallel do 
	
	end
!===================================================================================
!===================================================================================
	subroutine UpdateFuelHeight(fuel_height_array,fuel_density_array)

	use constants
	use grid_module
	use fireca_module
	use bld_module

	implicit none

	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: fuel_density_array	! kg/m3, density of the fuel at each cell
	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, height of the fuel at a certain (x,y) location
	integer :: i,j,k

	if(fuels%density_flag == 4 .or. fuels%density_flag == 66) then
		! Need to compute the fuel_height_array from density
		fuel_height_array = 0.
		!$OMP parallel do private(i,j,k)
		do j = 2,firegrid%ny-1
			do i = 2,firegrid%nx-1
				k = firegrid%nz
				do while(k > 1 .and. fuel_density_array(i,j,k) < MIN_FUEL_DENSITY)
					k = k - 1
				enddo
				if(k == 1 .and. fuel_density_array(i,j,k) < MIN_FUEL_DENSITY) then
					fuel_height_array(i,j) = 0.
					fuels%actual_height(i,j) = 0.
				else
					fuel_height_array(i,j) = firegrid%z(k+1)
				endif
			enddo
		enddo
		!$OMP end parallel do
	endif
	
	! Check that fuels%density, fuel_height_array, and fuels%actual_height are consistent
	!$OMP parallel do private(i,j,k)
	do k = 1,firegrid%nz
		do j = 2,firegrid%ny-1
			do i = 2,firegrid%nx-1
				if(fuel_density_array(i,j,k) < MIN_FUEL_DENSITY)then
					fuel_density_array(i,j,k) = 0.										
				elseif((k == 1 .and. fuels%actual_height(i,j) == 0.) .or. fuel_height_array(i,j) == 0) then
					fuel_density_array(i,j,k) = 0.
				endif
			enddo
		enddo
	enddo
	!$OMP end parallel do
	
	!$OMP parallel do private(i,j)
	do j = 2,firegrid%ny-1
		do i = 2,firegrid%nx-1
			! No density in the vertical column
			if(sum( fuel_density_array(i,j,:)) == 0) then
				fuels%actual_height(i,j) = 0.
				fuel_height_array(i,j) = 0.
			endif
			! No density in the first cell above the ground
			if(fuel_density_array(i,j,1) == 0) then
				fuels%actual_height(i,j) = 0.
			endif
		enddo
	enddo
	!$OMP end parallel do
	
	if(sum(fuel_density_array) == 0) then
		write(msgoutfile,*) 'No fuel in the domain'
		call TerminateProgram()
	endif

	allocate(fuels%height_initial(firegrid%nx,firegrid%ny))
	fuels%height_initial = fuel_height_array

	end
!===================================================================================
!===================================================================================
	SUBROUTINE ComputeLinearIndex(fuel_density_array)

	use file_handling_module
	use grid_module
	use fireca_module

	implicit none

	integer :: ew,ns,ud,ncell,index
	real,intent(IN),dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: fuel_density_array

	! Count number of cells with fuel
	firegrid%num_fuel_cells = 0
	do ud = 1,firegrid%nz
		do ns = 2,firegrid%ny-1
			do ew = 2,firegrid%nx-1
				if(fuel_density_array(ew,ns,ud) > MIN_FUEL_DENSITY) firegrid%num_fuel_cells = firegrid%num_fuel_cells + 1
			enddo
		enddo
	enddo

	firegrid%num_active_fuel_cells = firegrid%num_fuel_cells
	print *,firegrid%num_active_fuel_cells, "num fuel cells"
	allocate(firegrid%idx(firegrid%num_fuel_cells))
	
	print*,'Percentage of cells with fuel (fire domain): ', &		
		real(firegrid%num_fuel_cells)/real(firegrid%nx*firegrid%ny*firegrid%nz)*100.,'%'
	write(ID_FILE_TIMELOG,*)'Percentage of cells with fuel (fire domain): ', &		
		real(firegrid%num_fuel_cells)/real(firegrid%nx*firegrid%ny*firegrid%nz)*100.,'%'

	allocate(firegrid%cell_index(firegrid%num_fuel_cells), &
		firegrid%ijk_cell_index(firegrid%num_fuel_cells,3))
	! Calculate indexes for cells with fuel
	ncell = 0
	print *,"firegrid dims",firegrid%nx,firegrid%ny,firegrid%nz
	print *, shape(firegrid%ijk_cell_index),shape(firegrid%idx)
	print *, shape(fuel_density_array)
	do ud = 1,firegrid%nz
		do ns = 2,firegrid%ny-1
			do ew = 2,firegrid%nx-1
				if(fuel_density_array(ew,ns,ud) > MIN_FUEL_DENSITY)then
					ncell = ncell + 1
					!print *, ncell, shape(firegrid%idx),fuel_density_array(ew,ns,ud)
					call IndexCalc(ew,ns,ud,index)
					firegrid%cell_index(ncell) = index
					firegrid%ijk_cell_index(ncell,1) = ew
					firegrid%ijk_cell_index(ncell,2) = ns
					firegrid%ijk_cell_index(ncell,3) = ud
					firegrid%idx(ncell) = index
				endif
			enddo
		enddo
	enddo
	print *, "cleared"


	open(ID_FILE_FIREINDEX,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
		'fire_indexes.bin',status='replace',form='unformatted')
	write(ID_FILE_FIREINDEX) firegrid%num_fuel_cells
	write(ID_FILE_FIREINDEX) maxval(firegrid%ijk_cell_index(:,3))
	write(ID_FILE_FIREINDEX) (firegrid%cell_index(ew),ew=1,firegrid%num_fuel_cells)
	write(ID_FILE_FIREINDEX) ((firegrid%ijk_cell_index(ew,ns),ew= 1,firegrid%num_fuel_cells),ns=1,3)
	close(ID_FILE_FIREINDEX)

	end
!======================================================================================
!======================================================================================
	subroutine InitializeFCAdensity(fuel_density_array)

	use fireca_module
	use grid_module

	implicit none

	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz),intent(IN) :: fuel_density_array	! kg/m3, density of the fuel at each cell
	integer :: index,i,j,k

	! Initialize fuels%density
	!$OMP parallel do private(index, i, j, k)
	do index = 1,firegrid%num_fuel_cells
		i = firegrid%ijk_cell_index(index,1)
		j = firegrid%ijk_cell_index(index,2)
		k = firegrid%ijk_cell_index(index,3)

		fuels%density(index) = fuel_density_array(i,j,k)
	enddo
	!$OMP end parallel do 
	fuels%density_initial = fuels%density

	end
!===================================================================================
!===================================================================================
	subroutine ReadFuelMoisture()

	use grid_module
	use fireca_module
	use file_handling_module

	implicit none

	real,dimension(:,:,:),allocatable :: fuel_moisture_array
	character(STR_LEN) :: fname			! N/A, file name
	real,dimension(:,:,:,:),allocatable :: temp4
	real,dimension(:,:,:),allocatable :: temp_ft, temp3
	integer :: index,i,j,k

	if(fuels%moisture_flag == 1)then
		! Constant moisture
		fuels%moisture = fuels%input_moisture		! fraction: mass of water/mass of dry fuel

	else
		allocate(fuel_moisture_array(firegrid%nx,firegrid%ny,firegrid%nz))

		if(fuels%moisture_flag == 2)then
			fname = 'QF_FuelMoisture.inp'
			open(ID_FILE_QF_FUELMOIST,file=TRIM(workingDirectory)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted',err=1101)
			read(ID_FILE_QF_FUELMOIST,err=1102,end=1102) &
				(((fuel_moisture_array(i,j,k),i=1,firegrid%nx),j=1,firegrid%ny),k=1,firegrid%nz)
			close(ID_FILE_QF_FUELMOIST)

		elseif(fuels%moisture_flag == 3)then
			allocate(temp4(1,ft%nx,ft%ny,ft%n_grid_top))

			fname = 'treesmoist.dat'
			fuel_moisture_array = 0.
			open(ID_FILE_TREESMOIST,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
				status = 'old',form='unformatted', access='stream',err=1101)
			read(ID_FILE_TREESMOIST,err=1102,end=1102) temp4
			fuel_moisture_array = temp4(1,:,:,:)
			close(ID_FILE_TREESMOIST)

			deallocate(temp4)

		elseif(fuels%moisture_flag == 66)then
			if(ft%fuel_file_type == 1)then
				allocate(temp4(ft%n_fuel_types,ft%nx,ft%ny,ft%n_grid_top), &
					temp_ft(ft%nx,ft%ny,ft%nz))

				fname = 'treesmoist.dat'
				open(ID_FILE_TREESMOIST,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESMOIST,err=1102,end=1102) temp4
				close(ID_FILE_TREESMOIST)
				write (*,*) 'before moisture'
				temp4(1,:,:,2:13)=min(temp4(1,:,:,2:13),1.) !RRL 66 for estimated correction of Spring Hill Creek data provided to LANL 8/2020 from TT
				!temp4(1,:,:,:)=temp4(1,:,:,:)*1.5
				write (*,*) 'after moisture'
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(1,:,:,:)
				call InterpolateFiretechFile(fuel_moisture_array,temp_ft)

				deallocate(temp4,temp_ft)
			else
			
				allocate(temp3(ft%nx,ft%ny,ft%n_grid_top), temp_ft(ft%nx,ft%ny,ft%nz))

				fname = 'treesmoist_1.dat'
				open(ID_FILE_TREESMOIST,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESMOIST,err=1102,end=1102) temp3
				close(ID_FILE_TREESMOIST)
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
				call InterpolateFiretechFile(fuel_moisture_array, temp_ft)

				deallocate(temp3,temp_ft)
			endif
		else
			
			if(ft%fuel_file_type == 1)then
				allocate(temp4(ft%n_fuel_types,ft%nx,ft%ny,ft%n_grid_top), &
					temp_ft(ft%nx,ft%ny,ft%nz))

				fname = 'treesmoist.dat'
				open(ID_FILE_TREESMOIST,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESMOIST,err=1102,end=1102) temp4
				close(ID_FILE_TREESMOIST)
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp4(1,:,:,:)
				call InterpolateFiretechFile(fuel_moisture_array,temp_ft)

				deallocate(temp4,temp_ft)
			else
			
				allocate(temp3(ft%nx,ft%ny,ft%n_grid_top), temp_ft(ft%nx,ft%ny,ft%nz))

				fname = 'treesmoist_1.dat'
				open(ID_FILE_TREESMOIST,file=TRIM(ft%file_path)//TRIM(fileSeparator)//trim(fname), &
					status = 'old',form='unformatted', access='stream',err=1101)
				read(ID_FILE_TREESMOIST,err=1102,end=1102) temp3
				close(ID_FILE_TREESMOIST)
				temp_ft = 0.
				temp_ft(:,:,1:ft%n_grid_top) = temp3(:,:,:)
				call InterpolateFiretechFile(fuel_moisture_array, temp_ft)

				deallocate(temp3,temp_ft)
			endif
				
		endif

		where(fuel_moisture_array < MIN_MOISTURE)
			fuel_moisture_array = 0
		end where

		fuels%moisture = 0.
		!$OMP parallel do private(index, i, j, k)
		do index = 1,firegrid%num_fuel_cells
			i = firegrid%ijk_cell_index(index,1)
			j = firegrid%ijk_cell_index(index,2)
			k = firegrid%ijk_cell_index(index,3)
			fuels%moisture(index) = fuel_moisture_array(i,j,k)
		enddo
		!$OMP end parallel do 

		deallocate(fuel_moisture_array)
	endif

	fuels%moisture_initial = fuels%moisture

	return

1101 write(msgoutfile,*)'Error in opening input file: '//trim(fname)//'.'
	call TerminateProgram()

1102 write(msgoutfile,*)'Error in reading file: '//trim(fname)//'.'
	call TerminateProgram()

	end
!===================================================================================
!===================================================================================
	subroutine SpecifyFuel(fuel_height_array)
	
	use constants
	use fireca_constants_module
	use fireca_module
	use grid_module
	use bld_module
	use flags_module
	use interpolation_module
	
	implicit none

	real,dimension(firegrid%nx,firegrid%ny) :: fuel_height_array	! m, height of the fuel at a certain (x,y) location
	real,dimension(:,:,:),allocatable ::		&
		fuel_density_array,							& ! kg/m3, density of the fuel at each cell
		fuel_density_array_thick					  ! kg/m3, density of the thick fuel at each cell
	integer :: i
	
	allocate(fuel_density_array(firegrid%nx,firegrid%ny,firegrid%nz), &
		fuel_density_array_thick(firegrid%nx,firegrid%ny,firegrid%nz))

	fuel_density_array = 0.
	fuel_density_array_thick = 0.
	! --- Fuel density
	print *, "Fuel density"
	call ReadFuelDensity(fuel_density_array, fuel_density_array_thick, fuel_height_array)

	! --- Convert buildings into fuel
	if(bld%number > 0 .and. flag%bld_to_fuel == 1) call ConvertBld2Fuel(fuel_density_array)

	! --- Update fuel_height_array with the info about buildings or
	! for fuel density fields provided by the user or fuel height
	print *, "Fuel height"
	call UpdateFuelHeight(fuel_height_array,fuel_density_array)

	! --- Compute the linear
	print *, "linear index"
	call ComputeLinearIndex(fuel_density_array)

	! -- Now that we have the density field, we can allocate the arrays
	print *, "allocate arrays"
	call AllocateArrays()

	! --- Initialize fuels%density
	print *, "fca density"
	call InitializeFCAdensity(fuel_density_array)

	! --- Fuel moisture
	print *, "moisture"
	call ReadFuelMoisture()
	
	if(ft%n_fuel_types > 1) call DefineUrban(fuel_density_array_thick)
	
	deallocate(fuel_density_array, fuel_density_array_thick)
	
	return

	end
!===================================================================================
!===================================================================================
	subroutine InitFuelQUIC

	use fireca_module
	use flags_module
	use bld_module
	use canopy_module
	use grid_module
	use interpolation_module
	use winds_module
	
	implicit none

	integer,dimension(:,:,:),allocatable :: has_fuel
	integer :: cont,index,i,j,k

	allocate(has_fuel(qugrid%nx,qugrid%ny,qugrid%kfire_top))

	has_fuel = 0
	qugrid%num_fuel_cells = 0
	
	do k = 2,qugrid%kfire_top
		do j = 1,qugrid%ny-1
			do i = 1,qugrid%nx-1
				if(grid_conv(i,j,k)%nelem > 0) then
					if(sum(fuels%density(grid_conv(i,j,k)%index)) > 0)then
						has_fuel(i,j,k) = 1
						qugrid%num_fuel_cells = qugrid%num_fuel_cells + 1
					endif
				endif
			enddo
		enddo
	enddo
	print*,'Percentage of cells with fuel (QU domain): ', &
		real(qugrid%num_fuel_cells) / real((qugrid%nx-1)*(qugrid%ny-1)*(qugrid%nz-2)) * 100.
	write(ID_FILE_TIMELOG,*)'Percentage of cells with fuel (QU domain): ', &
		real(qugrid%num_fuel_cells) / real((qugrid%nx-1)*(qugrid%ny-1)*(qugrid%nz-2)) * 100.,'%'

	allocate(qugrid%ijk_cell_index(qugrid%num_fuel_cells,3))
	allocate(quwinds%uo_before_params(qugrid%num_fuel_cells),quwinds%vo_before_params(qugrid%num_fuel_cells))
	if(flag%quwinds_ave_file == 1)then
		allocate(quwinds%u_ave(qugrid%nx,qugrid%ny,qugrid%nz),	&
			quwinds%v_ave(qugrid%nx,qugrid%ny,qugrid%nz),quwinds%w_ave(qugrid%nx,qugrid%ny,qugrid%nz))
		quwinds%u_ave = 0.
		quwinds%v_ave = 0.
		quwinds%w_ave = 0.
	endif

	cont = 0
	do k = 2,qugrid%kfire_top
		do j = 1,qugrid%ny-1
			do i = 1,qugrid%nx-1
				if(has_fuel(i,j,k) == 1) then
					cont = cont + 1
					qugrid%ijk_cell_index(cont,1) = i
					qugrid%ijk_cell_index(cont,2) = j
					qugrid%ijk_cell_index(cont,3) = k
					! not already accounted for in the buildings
					if(bld%icellflag(i,j,k) /= 8 .and. bld%icellflag(i,j,k) /= 0)then
						quwinds%uo_before_params(cont) = quwinds%uo(i,j,k)
						quwinds%vo_before_params(cont) = quwinds%vo(i,j,k)
					endif
				endif
			enddo
		enddo
	enddo

	if(canopy%UPDATE_WINDS == 1)then
		allocate(init_dens_canopy(qugrid%num_fuel_cells))
		!$OMP parallel do private(index, i, j, k)
		do index = 1,qugrid%num_fuel_cells
			i = qugrid%ijk_cell_index(index,1)
			j = qugrid%ijk_cell_index(index,2)
			k = qugrid%ijk_cell_index(index,3)

			init_dens_canopy(index) = sum(fuels%density_initial(grid_conv(i,j,k)%index) * grid_conv(i,j,k)%en_ratio)
		enddo
		!$OMP end parallel do 
	endif

	deallocate(has_fuel)


	end
!===================================================================================
!===================================================================================
	subroutine ConvertBld2Fuel(fuel_density_array)

	use bld_module
	use canopy_module
	use grid_module
	use interpolation_module

	implicit none

	integer :: i,j,k,iii,jjj,kkk,kqf,kf_start,kf_end		! N/A, counter
	real :: fuel_density_array(firegrid%nx,firegrid%ny,firegrid%nz)	! kg/m3, fuel density array

	!$OMP parallel do private(i,j,k,iii,jjj,kkk,kqf,kf_start,kf_end)
	do k = 2,qugrid%nz-1
		! check if bottom of the quic cell is within the QUIC-FIRE domain
		if(qugrid%z(k-1) < firegrid%z(firegrid%nz)+firegrid%dz_array(firegrid%nz)) then

			kqf = 1
			do while(qugrid%z(k-1) >= firegrid%z(kqf) .and. kqf < firegrid%nz)
				kqf = kqf + 1
			enddo
			kqf = kqf - 1
			kf_start = max(kqf,1)

			kqf = 1
			do while(qugrid%z(k) > firegrid%z(kqf) .and. kqf < firegrid%nz+1)
				kqf = kqf + 1
			enddo
			kqf = max(kqf - 1,1)
			kf_end = min(kqf,firegrid%nz)

			do j = 1,qugrid%ny-1
				do i = 1,qugrid%nx-1
					if(bld%icellflag(i,j,k) == 0) then
						do kkk = kf_start,kf_end
							do jjj = 1 + (j-1) * yratio_int, j * yratio_int
								do iii = 1 + (i-1) * xratio_int , i * xratio_int
									if(fuel_density_array(iii,jjj,kkk) == 0) then
										fuel_density_array(iii,jjj,kkk) = canopy%BLDG_FUEL_DENS
									endif
								enddo
							enddo
						enddo
					endif
				enddo
			enddo
		endif
	enddo
	!$OMP end parallel do

	end
!===================================================================================
!===================================================================================
	subroutine update_ignitions_pattern()

	use grid_module
	use fireca_module
	use ignitions_module
	use time_module	

	implicit none
	
	integer :: i,j,k,m,									& ! N/A, counter
		is,ie,js,je,										& ! N/A, indexes of starting and ending cell around the point of ignition
		ncell													  ! N/A, index of the cell
	real ::													&
		tprec,tcurr,										& ! s, current and current time+dt
		x,y,													& ! m, cell centers
		dist													  ! m, distance between the ignition point and the cell center

	tprec = real(fcatime%current_int - fcatime%dt_int)
	tcurr = fcatime%current_real

	do m = 1,fire_ignition%npoints
		if(fire_ignition%time(m) >= tprec .and. fire_ignition%time(m) < tcurr) then
			! add ignitions
			is = max(ceiling((fire_ignition%x(m)-fire_ignition%radius(m))*firegrid%dxi) , 1)
			ie = min(ceiling((fire_ignition%x(m)+fire_ignition%radius(m))*firegrid%dxi) , firegrid%nx)
			js = max(ceiling((fire_ignition%y(m)-fire_ignition%radius(m))*firegrid%dyi) , 1)
			je = min(ceiling((fire_ignition%y(m)+fire_ignition%radius(m))*firegrid%dyi) , firegrid%ny)
			k = 1
			do while(fire_ignition%z(m) >= firegrid%z(k))
				k = k + 1
			enddo
			k = k - 1

			do j = js,je
				y = (real(j)-0.5)*firegrid%dy
				do i = is,ie
					x = (real(i)-0.5)*firegrid%dx
					dist = sqrt((fire_ignition%x(m)-x)**2 + (fire_ignition%y(m)-y)**2)
					if(dist <= fire_ignition%radius(m))then
						ncell = grid_where(i,j,k)
						if(ncell > 0) then
							if(fire%ignitions(ncell) <= 0) then
								! set cell on fire
								call InitVarFire(i,j,k,fire_ignition%new_num(m))
							elseif(fire%ignitions(ncell) >0) then
								! cell is already on fire, therefore, only add ignitions
								fire%ignitions(ncell) =	fire%ignitions(ncell) + fire_ignition%new_num(m)
								!else
								!do nothing, not enough fuel
							endif
						endif
					endif
				enddo
			enddo
		endif
	enddo

	call CheckFuelArrays()
	
	end
!======================================================================================
!======================================================================================
	subroutine InitWorkArrays()

	use fireca_module
	use grid_module
	use time_module
	use winds_module

	implicit none

	integer :: ii,ew,ns,ud

	fcawinds%sigma = 0.
	w_rr = 0.0   !reaction rate
	w_q = fire%reaction_rate * NET_ENERGY_WOOD_KG  ! heat release rate
	
	!$OMP parallel do private(ii, ew, ns, ud)
	do ii = 1,firegrid%num_fuel_cells
		ew = firegrid%ijk_cell_index(ii,1)
		ns = firegrid%ijk_cell_index(ii,2)
		ud = firegrid%ijk_cell_index(ii,3)
		! Initialize the wind magnitude field
		fuels%windmag(ii) = sqrt(fcawinds%u(ew,ns,ud)**2+ &
			fcawinds%w(ew,ns,ud)**2+fcawinds%v(ew,ns,ud)**2)

		firecell(ii)%sum_pos = 0
		firecell(ii)%sum_pos_squared = 0
		
		w_mass(ii) = fuels%density(ii)* firegrid%cellvol(ud)    ! mass of fuel in cell
	enddo
	!$OMP end parallel do 
	if(CELL_SIZE_MODE == CELL_SIZE_MODE_SMALL) then
		w_depl = fuels%density_initial - fuels%density
	else
		w_depl = -1
	endif
	
	w_moisture = fuels%moisture		! mass of water/mass of dry fuel
	w_fueldensity = fuels%density		! mass of fuel/ volume of cell
	w_ignitions = 0
	n_new_starts = 0.0

	fire%num_ignitions_not_absorbed = 0
	fire%energy_used_to_evaporate_water = 0.


	end
!======================================================================================
!======================================================================================
	SUBROUTINE AddFBIgnitions(coord_fb,cell_fb,tot_fb_time, &
		num_ignitions_from_firebrands,RT,ufb_x,ufb_y,ufb_z,Lx,Ly)

	use rnd_num_vars
	use grid_module
	use time_module
	use fireca_module

	implicit none

	real,intent(IN) ::				&
		Lx,Ly,							&
		ufb_x,ufb_y,ufb_z,			& ! m/s, firebrand velocity components
		RT									  ! m, radius of the envelope this firebrand corresponds to
	real,intent(IN),dimension(3) :: coord_fb		! m, firebrand landing position
	integer,intent(IN),dimension(3) :: cell_fb	! N/A, firebrand landing cell
	integer,intent(IN) ::					&
		tot_fb_time,							& ! s, firebrand time of flight
		num_ignitions_from_firebrands		  ! N/A, number of ignitions associated to a firebrand
	real,dimension(3) :: coord_ig			  ! m, location where ignitions are deposited
	integer,dimension(3) :: cell_ig		  ! N/A, cell where ignitions are deposited
	integer ::				&
		ncell,				& ! N/A, cell in the linear array
		ir,igamma,			& ! N/A, counters
		n_r, 					& ! N/A, number of concentric ellipses
		number_of_ignitions_via_firebrands_at_point ! N/A, number of ignitions per point
	real ::					&
		arg,					&
		gamma,				& ! rad, angle in the horizontal plane
		alpha,				& ! rad, current horizontal angle over which ignitions are distributed
		theta,				& ! rad, angle between the horizontal plane and teh z axis
		ufbhoriz,			& ! m/s, wind speed in the horizontal
		RT_s,					& ! m, major ellipses radius
		rad,					& ! m, current ellipses radius
		xprime,yprime,		& ! m, distance from the landing point
		n_r_real,n_gamma_real

	integer,parameter :: n_gamma = 8 ! N/A, number of radial positions

	call CalcCellFb(coord_fb,cell_fb)
	ufbhoriz = sqrt(ufb_x**2+ufb_y**2)

	gamma = atan2(ufb_y,ufb_x)

	if(ufb_z == 0.) then
		theta = GREEK_PI * 0.5
		n_r = 1
	else
		arg = min(max(ufbhoriz/sqrt(ufb_x**2+ufb_y**2+ufb_z**2), -0.999999), 0.999999)
		theta = max(fb%min_theta_value, acos(arg))	
		n_r = min(2,int(RT/sqrt(firegrid%dx**2+firegrid%dy**2)))
	endif
	RT_s = RT/sin(theta)

	if (n_r < 1 .or. num_ignitions_from_firebrands == 1) then

		ncell = grid_where(cell_fb(1),cell_fb(2),cell_fb(3))
		if (ncell > 0)then
			if (fuels%density(ncell) > MIN_FUEL_DENSITY)then
				fb%time_delay(ncell) = min(tot_fb_time, fb%time_delay(ncell))
				fb%num_ignitions(ncell) = fb%num_ignitions(ncell) + num_ignitions_from_firebrands
			endif
		endif

	else
		n_r_real = real(n_r)
		n_gamma_real = real(n_gamma)

		number_of_ignitions_via_firebrands_at_point =  &
			num_ignitions_from_firebrands/(n_gamma*(n_r+1))
		number_of_ignitions_via_firebrands_at_point = min(max( &
			fb%min_number_of_ignitions, &
			number_of_ignitions_via_firebrands_at_point), &
			fb%max_number_of_ignitions) ! 100  acceptable but fast

		if(number_of_ignitions_via_firebrands_at_point > 0) then

			ncell = grid_where(cell_fb(1),cell_fb(2),cell_fb(3))

			if(ncell > 0)then
				if (fuels%density(ncell) > MIN_FUEL_DENSITY)then
					fb%time_delay(ncell) = min(tot_fb_time, fb%time_delay(ncell))
					fb%num_ignitions(ncell) = fb%num_ignitions(ncell)+ &
						n_gamma*number_of_ignitions_via_firebrands_at_point
				endif
			endif
			call random_number(harvest,stat(1))
			do ir = 1,n_r
				do igamma = 1,n_gamma
					alpha = (harvest(1)+float(igamma)/n_gamma_real)*2.*PI
					rad = float(ir)/n_r_real

					xprime = rad*(RT_s*cos(alpha)*cos(gamma) - RT*sin(alpha)*sin(gamma))
					yprime = rad*(RT_s*cos(alpha)*sin(gamma) + RT*sin(alpha)*cos(gamma))

					coord_ig(1) = coord_fb(1) + xprime
					coord_ig(2) = coord_fb(2) + yprime
					coord_ig(3) = coord_fb(3)
					if(coord_ig(1) > 0 .and. coord_ig(1) < Lx .and. &
						coord_ig(2) > 0 .and. coord_ig(2) < Ly)then
						call CalcCellFb(coord_ig,cell_ig)
						ncell = grid_where(cell_ig(1),cell_ig(2),cell_ig(3))
						if (ncell > 0)then
							if (fuels%density(ncell) > MIN_FUEL_DENSITY)then

								fb%time_delay(ncell) = min(tot_fb_time, fb%time_delay(ncell))
								fb%num_ignitions(ncell) = fb%num_ignitions(ncell) + &
									number_of_ignitions_via_firebrands_at_point

							endif
						endif
					endif
				enddo
			enddo
		endif
	endif

	end
!======================================================================================
!======================================================================================
	SUBROUTINE PlaceFirebrandTimeDelay(cell_fb_loc,TID,num_cell_ignited,num_fb_ignitions)

	use time_module
	use fireca_module
	use grid_module
	use file_handling_module
	
	implicit none

	integer :: num_cell_ignited,num_fb_ignitions
	integer,intent(IN) :: TID
	integer,intent(IN),dimension(3) :: cell_fb_loc
	integer :: ncell_loc	 ! N/A, index of the cell where the firebrand is added

	ncell_loc = grid_where(cell_fb_loc(1),cell_fb_loc(2),cell_fb_loc(3))
	if(ncell_loc > 0)then
		if(fuels%density(ncell_loc) > MIN_FUEL_DENSITY .and. w_ignitions(ncell_loc) < 2.*max_init_ign)then

			write(ID_FILE_FB_OUT)fcatime%current_int,cell_fb_loc,1
			num_cell_ignited = num_cell_ignited + 1

			call SpreadFireSpotting(firegrid%cellvol(cell_fb_loc(3)),fuels%density(ncell_loc),		&
				fuels%moisture(ncell_loc),fuels%moisture_initial(ncell_loc),								&
				fuels%moist_depl_var(ncell_loc)%val,fuels%moist_depl_center(ncell_loc)%val,			&
				w_depl(ncell_loc),w_ignitions(ncell_loc),firecell(ncell_loc)%sum_pos,					&
				firecell(ncell_loc)%sum_pos_squared,fuels%density_initial(ncell_loc),					&
				TID,n_new_starts(ncell_loc),num_fb_ignitions)
		else
			write(ID_FILE_FB_OUT)fcatime%current_int,cell_fb_loc,-1  ! -1 = cell without fuel
		endif
	endif

	END
!======================================================================================
!======================================================================================
	SUBROUTINE CalcCellFb(coord_fb,cell_fb)

	use grid_module

	implicit none

	real,intent(IN),dimension(3) :: coord_fb
	integer,dimension(3) :: cell_fb

	cell_fb(1) = ceiling(coord_fb(1) * firegrid%dxi)
	cell_fb(2) = ceiling(coord_fb(2) * firegrid%dyi)
	cell_fb(3) = 1
	do while(firegrid%z(cell_fb(3)) < coord_fb(3) .and. cell_fb(3) < firegrid%nz)
		cell_fb(3) = cell_fb(3) + 1
	enddo
	cell_fb(3) = max(cell_fb(3)-1,1)

	end
!======================================================================================
!======================================================================================
	real function erfinv(x0)
	
	implicit none
	
	integer, parameter :: dp = kind(1.0d0)
	real, intent(in) :: x0

	real(kind=dp) :: w, p, x
	x = real(x0, dp)
	if (dabs(x) .lt. 1.0d0) then
		w = -dlog(1.0d0-x*x)
	else
		w = -dlog(1.0d0-0.9999d0*0.9999d0)
	endif
	if ( w < 5.000000 ) then
		w = w - 2.500000
		p = 2.81022636e-08
		p = 3.43273939e-07 + p*w
		p = -3.5233877e-06 + p*w
		p = -4.39150654e-06 + p*w
		p = 0.00021858087 + p*w
		p = -0.00125372503 + p*w
		p = -0.00417768164 + p*w
		p = 0.246640727 + p*w
		p = 1.50140941 + p*w
	else
		w = sqrt(w) - 3.000000
		p = -0.000200214257
		p = 0.000100950558 + p*w
		p = 0.00134934322 + p*w
		p = -0.00367342844 + p*w
		p = 0.00573950773 + p*w
		p = -0.0076224613 + p*w
		p = 0.00943887047 + p*w
		p = 1.00167406 + p*w
		p = 2.83297682 + p*w
	end if
	erfinv = real(p*x)
	end function erfinv
!======================================================================================
!======================================================================================
	subroutine FirebrandTransport(coord_fb0 , depth_fb)
	
	use fireca_constants_module
	use rnd_num_vars
	use time_module
	use fireca_module
	use grid_module
	use constants

	implicit none

	real,parameter ::									&
		viscosity_air = 1.47755e-5					  ! m2/s, air kinematic viscosity =  1.81e-5 / RHO_AIR

	real,intent(IN) :: depth_fb					  ! m, firebrand thickness
	real,dimension(3),intent(IN) :: coord_fb0	  ! m, initial firebrand coordinates

	integer :: is_flying, inside_domain,		& ! N/A, conditions to keep the firebrand flaying
		tot_fb_time,									& ! s, time of flight
		num_ignitions_from_firebrands				  ! N/A, number of ignitions associated with a firebrand

	real :: shear,										& ! dummy
		beta,b_fb,										& ! firebrand burnout variables
		ws,												& ! m/s, wind speed
		radius_fb,										& ! m, firebrand radius
		wrel,												& ! m/s, firebrand relative velocity
		wfb,												& ! m/s, firebrand absolute velocity
		RT,												& ! m, firebrand envelope radius
		urpurp_barsqrt,								& ! m/s
		sfb,												& ! m, length scale
		min_b_value_over_r2,							& !
		n_ign_real

	real,dimension(3)	::		&
		coord_fb,				& ! m, firebrand coordinates
		ua						     ! m/s, wind components at the firebrand location

	integer,dimension(3) :: cell_fb

	! Initialization
	radius_fb = 4.5 * depth_fb ! m, firebrand radius
	RT = sqrt(firegrid%cellarea / (fb%FRACTION_LAUNCHED_to_RT_ratio*fb%FRACTION_LAUNCHED)) ! m, initial envelope size
	sfb = fb%c_s * RT
	coord_fb = coord_fb0

	!rrl check check .52 below, maybe should be .104
	wrel = sqrt(2.*RHO_FB/RHO_AIR*depth_fb*GRAVITY/DRAG_COEFF)		!relative velocity
	beta = 1.04	* RHO_AIR / RHO_FB * sqrt(viscosity_air)		! m2/s
	b_fb = beta*sqrt(wrel/radius_fb**3)								! 1/s, radial burning

	! Conditions for flight
	is_flying = 1
	inside_domain = 1
	tot_fb_time = fcatime%current_int + fb%germination_delay

	! KJ         1      m3     kg      m2
	! ---     * ---- * ---- * ----  = ----
	! s * ign	 m      kg     kJ       s
	min_b_value_over_r2 = fb%min_b_value_coef * ENERGY_PER_IGNITION / fb%num_deposited / &
		(GREEK_PI*depth_fb*RHO_FB*HEAT_OF_COMBUSTION)

	do while(is_flying == 1 .and. b_fb > min_b_value_over_r2/radius_fb**2 .and. &
		radius_fb > MIN_FB_THICKNESS_AT_BURNOUT .and. inside_domain == 1)

		tot_fb_time = tot_fb_time + int(fb%time_step)

		call GetWinds(coord_fb,ua,shear,1,fb%w_mult)
		ws = sqrt(sum(ua**2))		! wind speed
		wfb = ua(3) - wrel			! speed of the firebrand in the vertical direction

		! Move firebrand
		coord_fb(1) = coord_fb(1) + fb%time_step * ua(1)
		coord_fb(2) = coord_fb(2) + fb%time_step * ua(2)
		coord_fb(3) = coord_fb(3) + fb%time_step * wfb

		call CalcCellFb(coord_fb,cell_fb)

		! Firebrand combustion
		b_fb = beta * sqrt(wrel/radius_fb**3)
		radius_fb = radius_fb - b_fb*radius_fb*0.5 * fb%time_step
		urpurp_barsqrt = 0.05 * max(wrel , ws)
		RT = RT + fb%time_step * (0.3*sfb*urpurp_barsqrt / RT)

		! Check if the firebrand is outside the domain
		if(coord_fb(1) < 0 .or. coord_fb(1) > qugrid%Lx .or. &
			coord_fb(2) < 0 .or. coord_fb(2) > qugrid%Ly .or. &
			coord_fb(3) > qugrid%Lz) then

			inside_domain = 0

			! Check if firebrand is on the ground
		elseif(coord_fb(3) <=  0.   .or.															&
			(coord_fb(3) <=  fuels%height_initial(cell_fb(1),cell_fb(2)) .and.	&
			fuels%height_initial(cell_fb(1),cell_fb(2)) >= tree_height)) then

			if(coord_fb(3) <= 0) coord_fb(3) = firegrid%zm(1)

			! Number of ignition	associated to this firebrand, considering its burning rate
			n_ign_real = fb%num_deposited*b_fb*GREEK_PI*radius_fb**2*									&
				depth_fb*RHO_FB*HEAT_OF_COMBUSTION/ENERGY_PER_IGNITION
			num_ignitions_from_firebrands = int(n_ign_real)
			if(num_ignitions_from_firebrands <= 0) then
				call random_number(harvest,stat(1))
				num_ignitions_from_firebrands = int(n_ign_real + harvest(1))
			endif
			if(num_ignitions_from_firebrands >= 1) then
				call AddFBIgnitions(coord_fb,cell_fb,tot_fb_time,							&
					num_ignitions_from_firebrands,RT,ua(1),ua(2),wfb,qugrid%Lx,qugrid%Ly)
				is_flying = 0
			endif
		endif
	enddo

	end 
!======================================================================================
!======================================================================================
	subroutine RollOverWorkingArrays()

	use fireca_module

	implicit none

	fire%reaction_rate = w_rr
	fire%ignitions = w_ignitions
	fuels%density = w_fueldensity

	end 
!======================================================================================
!======================================================================================
	subroutine ComputeFirebrand()

	use fireca_module
	use time_module
	use fireca_constants_module
	use grid_module
	use fire_module

	implicit none

	integer ::								&
		TID,									&
		index,								&
		ncell,								&
		ig,ew,ns,ud

	integer,dimension(3) :: cell_fb

	! Generate new ignitions
	TID = 1
	ig = 0
	ncell = 0
	do index = 1,firegrid%num_fuel_cells
		if(w_ignitions(index) > 0)then
			ig = ig + 1
			ew = firegrid%ijk_cell_index(index,1)
			ns = firegrid%ijk_cell_index(index,2)
			ud = firegrid%ijk_cell_index(index,3)

			call CheckFBLaunch(ew,ns,ud,TID,ncell)
		endif
	enddo

	ncell = 0
	! Put them in cells
	do index = 1,firegrid%num_fuel_cells
		if(fb%time_delay(index) <= fcatime%current_int) then
			cell_fb = (/firegrid%ijk_cell_index(index,1),firegrid%ijk_cell_index(index,2), &
				firegrid%ijk_cell_index(index,3)/)
			call PlaceFirebrandTimeDelay(cell_fb,TID,ncell, &
				fb%num_ignitions(index))
			fb%time_delay(index) = FB_TIME_DEFAULT
			fb%num_ignitions(index) = 0.0
		endif
	enddo

	end 
!======================================================================================
!======================================================================================
	subroutine CheckFBLaunch(ew,ns,ud,TID,num_fb_launched)

	use rnd_num_vars
	use grid_module
	use fireca_module

	implicit none

	integer,intent(IN) ::	&
		ew,ns,ud,				& ! N/A, current cell
		TID						  ! N/A, used for random number generation

	integer,intent(INOUT) :: num_fb_launched

	real ::						&
		disk_thickness,		& ! m, thickness of the firebrand, assumed to be a disk
		dummy						  ! shear, not used here

	real,dimension(3) ::		&
		winds,					& ! m/s, 3D wind components
		coord_fb					  ! m, coordinates of the firebrand launch

	! Firebrand is launched from the center of the cell and its top
	coord_fb = (/firegrid%xcenters(ew),firegrid%ycenters(ns),firegrid%z(ud+1)/)
	call GetWinds(coord_fb,winds,dummy,1,fb%w_mult)

	if(winds(3) > 0)then
		disk_thickness = fb%frac_of_max_size * DRAG_COEFF * 0.5 * &
			RHO_AIR/RHO_FB * winds(3)**2 / GRAVITY

		if(disk_thickness >= FB_MIN_THICKNESS)then
			call random_number(harvest,stat(TID))
			if(harvest(1) <= fb%FRACTION_LAUNCHED) then
				num_fb_launched = num_fb_launched  + 1
				call FirebrandTransport(coord_fb,disk_thickness)
			endif
		endif
	endif

	end 
!======================================================================================
!======================================================================================
	subroutine ComputeAverages(ew,ns,ud,lin_index,lengthscale, &
		avg_u,avg_v,avg_w,avg_rr,avg_q)

	use fireca_constants_module
	use grid_module
	use fireca_module
	use winds_module
	use interpolation_module

	implicit none

	real,intent(IN) :: lengthscale			  ! m, flame length
	integer,intent(IN) ::						&
		lin_index,									& ! N/A, linear indexes of the current cell
		ew,ns,ud										  ! N/A, indexes of the current cell

	real,intent(OUT) ::							&
		avg_u,avg_v,avg_w,						& ! m/s, average u,v,w when the fireca and quic grids are not matching
		avg_rr,										& ! kg/m3/s, average reaction rate when the fireca and quic grids are not matching
		avg_q											  ! kg water/kg air, average moisture when the fireca and quic grids are not matching

	real ::											&
		wgt_spatial_average,						& ! N/A, function to compute the average of 3D variables
		wgt_spatial_average_3D,					& ! N/A, function to compute the average of compressed array variables
		len_cell										  ! 1/m, lengthscale/firegrid_avg_cellsize

	! Average calculated over the flame length
	if (matching_grids .eq. 0) then
		len_cell = lengthscale / firegrid%avg_cellsize
		avg_u = wgt_spatial_average_3D(fcawinds%u(:,:,ud),ns,ew,len_cell,1.0)
		avg_v = wgt_spatial_average_3D(fcawinds%v(:,:,ud),ns,ew,len_cell,1.0)
		avg_w = wgt_spatial_average_3D(fcawinds%w(:,:,ud),ns,ew,len_cell,1.0)

		avg_rr = wgt_spatial_average(fire%reaction_rate,ew,ns,ud,len_cell,1.0)
		! average heat release rate to 1/3 power, (kW/m3)**(1/3)
		avg_q = wgt_spatial_average(w_q,ew,ns,ud,lengthscale*firegrid%zm(ud),ONE_THIRD)

	else
		avg_u = fcawinds%u(ew,ns,ud)
		avg_v = fcawinds%v(ew,ns,ud)
		avg_w = fcawinds%w(ew,ns,ud)
		avg_rr = fire%reaction_rate(lin_index)
		avg_q = w_q(lin_index)**ONE_THIRD  ! (kW/m3)**(1/3)
	endif

	end 
!======================================================================================
!======================================================================================
	real function wgt_spatial_average_3D(x,ns,ew,l,p)

	use grid_module

	implicit none

	integer, intent(in) :: ew,ns ! N/A, current cell in the x,y directions


	real, intent(in), dimension(firegrid%nx,firegrid%ny) :: x  ! variable to average
	real, intent(in) ::						&
		l,											& ! N/A, flame length / average cell size (dx+dy)*0.5
		p											  ! N/A, exponenent used in the averaging

	real :: numerator, denominator, wgt
	integer :: i,j,is,ie,js,je,			& ! N/A, counters
		int_l									     ! N/A, calculate int(l) just once

	numerator = 0.
	denominator = 0.
	int_l = int(l)
	is = max(ew-int_l,1)
	ie = min(ew+int_l,firegrid%nx)
	js = max(ns-int_l,1)
	je = min(ns+int_l,firegrid%ny)	
	do j = js,je
		do i = is,ie
			wgt = 1./max( float((i-ew)**2+(j-ns)**2), 0.0001)
			numerator = numerator + (x(i,j)**p) * wgt
			denominator = denominator + wgt
		enddo
	enddo
	wgt_spatial_average_3D = numerator/denominator

	end function wgt_spatial_average_3D
!======================================================================================
!======================================================================================
	real function wgt_spatial_average(x,ew,ns,ud,l,p)

	use grid_module
	use fireca_module

	implicit none

	integer, intent(in) :: ns,ew,ud ! N/A, current cell in the x,y,z directions

	real, intent(in), dimension(firegrid%num_fuel_cells) :: x  ! variable to average
	real, intent(in) ::						&
		l,											& ! N/A, flame length / average cell size (dx+dy)*0.5
		p											  ! N/A, exponenent used in the averaging

	real :: numerator, denominator, wgt
	integer :: i,j,is,ie,js,je,			& ! N/A counters
		int_l,									& ! N/A, calculate int(l) just once
		ncell

	numerator = 0.
	denominator = 0.
	int_l = int(l)
	is = max(ew-int_l,1)
	ie = min(ew+int_l,firegrid%nx)
	js = max(ns-int_l,1)
	je = min(ns+int_l,firegrid%ny)
	do j = js,je
		do i = is,ie
			wgt = 1./max( float((i-ew)**2+(j-ns)**2), 0.0001)
			if(grid_where(i,j,ud) > 0) then
				ncell = grid_where(i,j,ud)
				numerator = numerator + (x(ncell)**p) * wgt
			endif
			denominator = denominator + wgt
		enddo
	enddo
	wgt_spatial_average = numerator/denominator

	end function wgt_spatial_average
!======================================================================================
!======================================================================================
	real function calc_shear(u,v,w,ew,ns,ud)

	use grid_module

	implicit none

	real, intent(in), dimension(firegrid%nx,firegrid%ny,firegrid%nz+1) :: u,v,w		! m/s, wind components in the x,y,z direction
	integer, intent(in) :: ew,ns,ud						! N/A, current cell indexes
	real :: dudx, dvdy, dwdz, shear, dzz

	dudx = (abs(u(ew+1,ns,ud)-u(ew,ns,ud))+ &
		abs(u(ew,ns,ud)-u(ew-1,ns,ud)))*.5*firegrid%dxi
	dvdy = (abs(v(ew,ns+1,ud)-v(ew,ns,ud))+ &
		abs(v(ew,ns,ud)-v(ew,ns-1,ud)))*.5*firegrid%dyi

	if (ud.eq.1) then
		dwdz = (abs(w(ew,ns,ud+1)-w(ew,ns,ud))+ &
			abs(2*w(ew,ns,ud)))*.5*firegrid%dzmi(ud)
	else
		dwdz = abs(w(ew,ns,ud+1)-w(ew,ns,ud))*firegrid%dzmi(ud)+ &
			abs(w(ew,ns,ud)-w(ew,ns,ud-1))*firegrid%dzmi(ud-1)  !rrl_flag
	endif

	shear = 2.0*(dudx**2+dvdy**2+ dwdz**2)
	!lines below are taking averages of 4 values on the 4 lateral sides.
	shear = shear + 0.25*( &
		((firegrid%dyi*(u(ew,ns+1,ud)-u(ew,ns,ud)))+ &
		(.25* firegrid%dxi*(v(ew+1,ns,ud)-v(ew-1,ns,ud)+  &
		v(ew+1,ns+1,ud)-v(ew-1,ns+1,ud))))**2+ &
		((firegrid%dyi*(u(ew,ns,ud)-u(ew,ns-1,ud)))+ &
		(.25* firegrid%dxi*(v(ew+1,ns,ud)-v(ew-1,ns,ud)+  &
		v(ew+1,ns-1,ud)-v(ew-1,ns-1,ud))))**2+ &
		((firegrid%dxi*(v(ew,ns+1,ud)-v(ew,ns,ud)))+ &
		(.25* firegrid%dyi*(u(ew+1,ns+1,ud)-u(ew+1,ns-1,ud)+  &
		u(ew,ns+1,ud)-u(ew,ns-1,ud))))**2+ &
		((firegrid%dxi*(v(ew,ns,ud)-v(ew,ns-1,ud)))+ &
		(.25* firegrid%dyi*(u(ew,ns+1,ud)-u(ew,ns-1,ud)+  &
		u(ew-1,ns+1,ud)-u(ew-1,ns-1,ud))))**2)
	if (ud.gt.1) then
		!lines below are taking averages of 4 values on the xz lateral sides.
		dzz = 1./(firegrid%dz_array(ud)+0.5*firegrid%dz_array(ud-1)+0.5*firegrid%dz_array(ud+1))
		shear = shear + 0.25*( &
			( firegrid%dzmi(ud)  *(u(ew,ns,ud+1)-u(ew,ns,ud))  + &
			0.25* firegrid%dxi*(w(ew+1,ns,ud)  - w(ew-1,ns,ud)  + w(ew+1,ns,ud+1) - w(ew-1,ns,ud+1))   )**2 + &
			( firegrid%dzmi(ud-1)*(u(ew,ns,ud)  -u(ew,ns,ud-1))+ &
			0.25* firegrid%dxi*(w(ew+1,ns,ud)  - w(ew-1,ns,ud)  + w(ew+1,ns,ud-1) - w(ew-1,ns,ud-1))   )**2 + &
			( firegrid%dxi      *(w(ew+1,ns,ud)-w(ew,ns,ud))  + &
			0.5 * dzz			 *(u(ew+1,ns,ud+1)-u(ew+1,ns,ud-1) + u(ew,ns,ud+1)	  - u(ew,ns,ud-1))     )**2 + &
			( firegrid%dxi      *(w(ew,ns,ud)-w(ew-1,ns,ud))  + &
			0.5 * dzz         *(u(ew,ns,ud+1)  -u(ew,ns,ud-1)   + u(ew-1,ns,ud+1) - u(ew-1,ns,ud-1))	  )**2)
		!lines below are taking averages of 4 values on the yz sides.
		shear = shear + 0.25*( &
			( firegrid%dyi*(w(ew,ns+1,ud)-w(ew,ns,ud))      + &
			0.25* dzz         * (v(ew,ns,ud+1)   - v(ew,ns,ud-1)   + v(ew,ns+1,ud+1) - v(ew,ns+1,ud-1))   )**2 + &
			( firegrid%dyi*(w(ew,ns,ud)-w(ew,ns-1,ud))      + &
			0.25* dzz         * (v(ew,ns,ud+1)   - v(ew,ns,ud-1)   + v(ew,ns-1,ud+1) - v(ew,ns-1,ud-1))   )**2 + &
			( firegrid%dzmi(ud)*(v(ew,ns,ud+1)-v(ew,ns,ud))  + &
			0.25* firegrid%dyi* (w(ew,ns+1,ud+1) - w(ew,ns-1,ud+1) + w(ew,ns+1,ud)   - w(ew,ns-1,ud))     )**2 + &
			( firegrid%dzmi(ud-1)*(v(ew,ns,ud)-v(ew,ns,ud-1))+ &
			0.25* firegrid%dyi* (w(ew,ns+1,ud)   - w(ew,ns-1,ud)   + w(ew,ns+1,ud-1) - w(ew,ns-1,ud-1))   )**2)
	else
		!lines below are taking averages of 4 values on the xz lateral sides.****neglects vertical shear at the ground
		shear = shear + 0.25*( &
			( firegrid%dzmi(ud)*(u(ew,ns,ud+1)-u(ew,ns,ud))+ &
			0.25* firegrid%dxi*(w(ew+1,ns,ud)-w(ew-1,ns,ud)+w(ew+1,ns,ud+1)-w(ew-1,ns,ud+1))   )**2 + &
			(                                              0.25* firegrid%dxi*2.*(w(ew+1,ns,ud)-w(ew-1,ns,ud)))**2 + &
			( firegrid%dxi*(w(ew+1,ns,ud)-w(ew,ns,ud)))**2 + &
			( firegrid%dxi*(w(ew,ns,ud)  -w(ew-1,ns,ud)))**2)

		!lines below are taking averages of 4 values on the yz sides.****neglects vertical shear at the ground
		shear = shear + 0.25*( &
			( firegrid%dyi*(w(ew,ns+1,ud)-w(ew,ns,ud))														)**2 + &
			( firegrid%dyi*(w(ew,ns,ud)-w(ew,ns-1,ud))														)**2 + &
			( firegrid%dzmi(ud)*(v(ew,ns,ud+1)-v(ew,ns,ud)) + &
			0.25* firegrid%dyi*(w(ew,ns+1,ud+1)-w(ew,ns-1,ud+1)+ w(ew,ns+1,ud)-w(ew,ns-1,ud))	)**2 + &
			( firegrid%dzmi(ud)*2.*v(ew,ns,ud)				  + &
			0.25* firegrid%dyi*(-w(ew,ns-1,ud))																	)**2)
		shear=shear*0.5
	endif
	calc_shear = shear

	end function calc_shear
!======================================================================================
!======================================================================================
	subroutine NoFuel(fire_energy_to_atmos,w_mass,w_moisture,w_ignitions)

	use fireca_constants_module

	implicit none

	real,intent(INOUT) :: fire_energy_to_atmos

	real,intent(OUT) :: w_mass,w_moisture
	integer,intent(OUT) :: w_ignitions

	w_mass = 0.
	w_moisture = 0.
	w_ignitions = 0
	! decay of energy to atmos
	fire_energy_to_atmos = min(ENERGY_PER_IGNITION,fire_energy_to_atmos)

	end 
!======================================================================================
!======================================================================================
	subroutine SlideBurnCenterInfo()

	use fireca_module
	use grid_module

	implicit none

	integer ::				&
		icount,				&
		index

	!$OMP parallel do private(index,icount)
	do index = 1,firegrid%num_fuel_cells
		if(fire%burncenter_info(index)%val2(1,1) > 0.) then
			do icount = nbtime,3,-1
				fire%burncenter_info(index)%val2(icount,1) = fire%burncenter_info(index)%val2(icount-1,1)
				fire%burncenter_info(index)%val2(icount,2) = fire%burncenter_info(index)%val2(icount-1,2)
				fire%burncenter_info(index)%val2(icount,3) = &
					fire%burncenter_info(index)%val2(icount-1,3) * timedecay_ratio(icount)
			enddo

			icount = 2
			fire%burncenter_info(index)%val2(icount,1) = max(0.0,min(1.,fire%burncenter_info(index)%val2(icount-1,1)))
			fire%burncenter_info(index)%val2(icount,2) = max(0.0,min(1.,fire%burncenter_info(index)%val2(icount-1,2)))
			fire%burncenter_info(index)%val2(icount,3) = &
				fire%burncenter_info(index)%val2(icount-1,3) * timedecay_ratio(icount)
		endif
	enddo
	!$OMP end parallel do

	end 
!======================================================================================
!======================================================================================
	subroutine ComputeSootYield(flame_length,gas_velocity,O2density,emissions,mu,sigma)

	! This subroutine takes a number of fitting parameters, a gas velocity, normalized O2
	! depletion, and the predicted flame length to predict soot production from a give flame
	! or packet of energy. O2d_o is the initial mass density of O2 at room temperature, and
	! is used to normalize the current o2density. FL is the flame length predicted elsewhere
	! in QUICFIRE, and gas_velocity is the mean wind magnitude (which includes) buoyancy. Tau is effectively
	! a residence time that indicates whether or not the soot has time to be completely consumed
	! within the flame envelope. Finally, the rational polynomial fit of soot data from more detailed
	! simulations is good but not perfect. A small error is occasionally introduced when soot values
	! exceed physical limitations in extreme case. This is corrected by if statements that enforce
	! the physical limits of soot production.

	!S_yield is the fraction of the dry mass consumed that went to soot.

	implicit none

	real,dimension(7), parameter ::	& ! These are rational polynomial parameters used as a surrogate for the detailed soot model
		SOOT_PARAMS = (/					&
		20091.385,							&
		28606.449,							&
		34751.318,							&
		1618.158,							&
		-67747.090,							&
		39516.816,							&
		-634.157/)

	real,dimension(5), parameter ::	& ! These are polynomial parameters as a surrogate for the detailed model
		SIGMA_PARAMS = (/					&
		1.10281176,							&
		0.0324021,							&
		-0.920322,							&
		0.51317,								&
		0.023716/)

	real,dimension(4), parameter ::	& !These are polynomial parameters as a surrogate for the detailed model
	MU_PARAMS = (/						&
		-14.1689,							&
		0.890782,							&
		-0.478186,							&
		-0.261891/)

	real,parameter ::						&
		O2D_O = 0.257,						& ! kg/m3, density of O2 in air at ambient temperature and see level
		SOOT_POTENTIAL = 0.53			  ! N/A, sooting potential (if all heavy pyrolysates go to soot) for Douglas Fir or similar

	real,intent(IN) ::					&
		flame_length,						& ! m, flame length
		gas_velocity,						& ! m/s, gas velocity
		O2density							  ! kg/m3, O2 density

	real,intent(OUT) ::					&
		emissions,							& !
		mu,									& ! m
		sigma									  ! m

	real ::									&
		ofrac

	! Check the residence time cutoff
	if(flame_length / gas_velocity > 0.5) then
		emissions = 0.0
		mu = 0.
		sigma = 0.
	else
		! Compute soot yield fraction, and the first and second moments of the log-normal distribution
		ofrac = min(O2density/O2D_O, 1.)

		emissions = SOOT_POTENTIAL - (SOOT_PARAMS(1) * flame_length) / &
			(SOOT_PARAMS(2) + SOOT_PARAMS(3) * flame_length + SOOT_PARAMS(4) * gas_velocity + &
			SOOT_PARAMS(5) * ofrac + SOOT_PARAMS(6) * ofrac**2 + SOOT_PARAMS(7) * gas_velocity * ofrac)

		emissions = min(emissions, 1.)
		emissions = max(emissions, 0.)

		mu = MU_PARAMS(1) + MU_PARAMS(2) * ofrac + MU_PARAMS(3) * ofrac**2 + &
			MU_PARAMS(4) * ofrac * flame_length

		sigma = SIGMA_PARAMS(1) + SIGMA_PARAMS(2) * flame_length + &
			SIGMA_PARAMS(3) * ofrac + SIGMA_PARAMS(4) * ofrac**2 + SIGMA_PARAMS(5) * flame_length * ofrac
	endif

	end 
!======================================================================================
!======================================================================================
	Subroutine Psi_to_Temp(psi,temps)

	! This function accepts a value of Psi and returns a mean temperature. Other temperatures in
	! the cell can be established via the  pdf with known mean and variance.
	implicit none

	real, parameter ::						&
		c1 = 0.5,								& !     c1, c2, and c3 are constants to shift from the mean of 0 and std=0.5 of the error function
		c2 = 0.0058,							&
		c3 = 1.0,								&
		tcrit = 600.0,							& !     tcrit is the critical temperature at which "burning" begins, and tfstep is the calibrated
		tfstep = 310.0,						& !     a parameter that (in conjunction with c1, c2, and c3) fit the temperature pdf
		p1_2 = 0.283013,						& !     All p1_* and p2_* are coefficients for the piecewise polynomial that maps the inverse error function to psi
		p1_1 = 0.807656,						&
		p1_0 = 0.00420903,					&
		p2_4 = 241.333,						&
		p2_3 = -740.04132,					&
		p2_2 = 848.515,						&
		p2_1 = -429.210,						&
		p2_0 = 81.2893

	real,intent(IN) ::						&
		psi										  !The fraction of the cell that is burning


	real,intent(OUT) ::						&
		temps										  ! K, The mean solid temperature of the cell

	real ::										&
		er_inv,									& !The result of the surrogate error inverse function
		rev_psi,									& !The transformation of psi to an appropriate input for the error inverse function
		temp_i									  !An intermediated temperature in the calculation

	rev_psi = psi/c1-c3 !This reverse engineers psi to the value that the inverse error function needs as an input.

	! The following series of if statements mimics the inverse error function with an acceptable level
	! of error in the final predicted mean solid temperature. Usually this means to within aobut 0.5K,
	! though the extreme ends of the distribution become progressively more uncertain. In particular,
	! if the reverse engineered psi is greater than 0.99 or less than -0.99, we are in the extreme ends
	! of the temperature curve, and the actual temperature value could be anything above 900 K or below
	! 298 K. In practical terms, the lower limit is ambient temperature and the upper limit is the adiabatic
	! flame temperature (although we are usually quite a bit cooler than that).

	!Compute the piecewise polynomial fit of the inverse error function and the resulting temperature.
	if (abs(rev_psi) >= 0.99) then
		if (rev_psi > 0) then
			temps = 900.  !This value really means that we are at 900 or more,with high uncertainty
		else
			temps = 298.  !This value really means we are at ambient temperature
		endif

	elseif (abs(rev_psi) < 0.605) then
		er_inv = p1_2*abs(rev_psi)**2 + p1_1*abs(rev_psi) + p1_0
		temp_i = er_inv/c2
		if (rev_psi < 0) then
			temp_i = -temp_i
		endif
		temps = temp_i + tcrit

	else
		er_inv = p2_4*abs(rev_psi)**4 + p2_3*abs(rev_psi)**3 +  &
			p2_2*abs(rev_psi)**2 + p2_1*abs(rev_psi)+p2_0
		temp_i = er_inv/c2
		if (rev_psi < 0) then
			temp_i = -temp_i
		endif
		temps = temp_i+tcrit
	endif

	end 
!======================================================================================
!======================================================================================
	SUBROUTINE ComputeConvectiveHeatTransfer(windmag,k2sqrt,psi,conv_human,current_density)

	use fireca_constants_module

	implicit none

	real,parameter ::									&
		ECON1 = 0.193,									& ! N/A, empirical constant
		ECON2 = 0.618,									& ! N/A, second empirical constant, Equation 7.52 Incropera & DeWitt
		ECONK = 30.49e-3,								& ! W/mK thermal conductivity of air at 360K
		PRANDTL = 15.89/22.5,						& ! N/A, Prandtl number
		HUMAN_ARM_DIAMETER = 0.1,					& ! m, diameter of a human arm
		HUMAN_CORE_TEMPERATURE = 310.15,			& ! K, human core temperature
		air_viscosity=0.00005099					  ! Kinematic viscosity of air at 600K

	real, intent(IN) ::								&
		windmag,											&
		psi,												&
		k2sqrt,											&
		current_density							      ! The current density of the fuel in the cell (kg/m^3)

	real,intent(OUT) ::								&
		conv_human										  ! W/m2, convective heat flux per human m2

	real ::												&
		turb_wind,										& ! m/s, wind magnitude enhanced by turbulent mixing
		h_conv,											& ! J/m^2/K/s, convective heat transfer coefficient
		gas_temperature,								& ! K, mean gas temperature
		solid_temperature,							& ! K, mean solid temperature
		decoupling_coeff								  ! A coefficient that decreases the coupling
	!between the solid and gas temperature as the density of the solid decreases
	turb_wind = windmag + k2sqrt !get the wind magnitude enhanced by turbulent mixing
	call Psi_to_temp(psi,solid_temperature) !get the mean solid temperature based on Psi

	!heat transfer parameters

	!Get gas temperature from the solid temperature
	decoupling_coeff=exp(-current_density*40.) !This factor should begin to deviate notably from 0
	!For fuel density values of less than 0.1 kg/m^3 and be 96% decoupled by the minimum fuel density
	gas_temperature = decoupling_coeff*(AMBIENT_TEMPERATURE - solid_temperature)+solid_temperature

	h_conv = ECON1 * ((((turb_wind * HUMAN_ARM_DIAMETER) / air_viscosity) ** ECON2) &
		* PRANDTL) * (Econk / HUMAN_ARM_DIAMETER)            ! convective heat transfer coefficient
	conv_human = h_conv * (gas_temperature - HUMAN_CORE_TEMPERATURE) !convective heat flux in watts/m^2 (for human)
	conv_human = max(0.0,conv_human) !because ambient temperatures fluxuate, and we don't wantto through off our heat flux score with negligible cooling

	end 
!======================================================================================
!======================================================================================
	subroutine update_fire(firebrands_flag,quic_update_time)

	use omp_lib
	use constants
	use mersenne_twister
	use rnd_num_vars
	use grid_module
	use time_module
	use fireca_module

	implicit none

	integer, intent(IN) ::			&		
		firebrands_flag,				&
		quic_update_time

	integer :: num_fires

	! Initialize work arrays
	call InitWorkArrays()

	! Spread fire
	call ComputeNewFireSpread()

	! Firebrands
	if(firebrands_flag == 1) then
		fb%time = fb%time + fcatime%dt_int
		if(fb%time > max(fb%LAUNCH_TIME , quic_update_time) .and. mod(fb%time,FB%LAUNCH_TIME) == 0) then
			call ComputeFirebrand()
		endif
	endif

	! Slide the burncenter info
	call SlideBurnCenterInfo()

	! Moisture variables
	call UpdateMoisture()

	! Roll over working arrays
	call RollOverWorkingArrays()

	num_fires = count(fire%ignitions > 0)
	
	! Check if there is any fire
	if(num_fires > 0) then
		!print*,'   - number of cells on fire after fire spread = ',num_fires
	else
		write (msgoutfile,*)'******** Fire is out ********'
	endif

	call CheckFuelArrays()
	
	return

	end
!======================================================================================
!======================================================================================
	SUBROUTINE CheckFuelArrays()
	
	use fireca_module
	use fireca_constants_module
	use grid_module
	
	implicit none
	
	integer :: index
		
	firegrid%num_active_fuel_cells = 0
	do index = 1, firegrid%num_fuel_cells
		if (fuels%density(index) >= MIN_FUEL_DENSITY .and. fire%ignitions(index) > 0) then
			firegrid%num_active_fuel_cells = firegrid%num_active_fuel_cells + 1
			firegrid%idx(firegrid%num_active_fuel_cells) = index !firegrid%cell_index(index)
		endif
	enddo
	
	end
!======================================================================================
!======================================================================================
	subroutine ComputeWindComponents(ww_q,avg_q,avg_u,avg_v,avg_w,zb,curr_cellvol, &
		wprime,local_w,hwindmag,windmag)

	use constants
	use fireca_constants_module
	use interpolation_module

	implicit none

	real,parameter ::									&
		p_t0 = 300.,									& ! K, approximate ambient temperature
		p_w = 1.0,										& ! N/A, parameter controlling magnitude of w_prime
		c_cp = 1.004									  ! KJ/Kg/K, specific heat of air

	real,intent(IN) ::								&
		ww_q,												& ! kW/m3, (reaction rate) * (wood combustion entalphy) = energy released per cubic meter
		avg_q,											& ! (kW/m3)**(1/3), average w_qq
		avg_u,avg_v,avg_w,							& ! m/s, average u,v,w components of the mean wind
		zb,												& ! m, cell center
		curr_cellvol									  ! m3, cell volume

	real,intent(OUT) ::								&
		wprime,											& ! m/s, perturbation to the average vertical velocity calculated by quic
		local_w,											& ! m/s, local vertical velocity = average + perturbation (wprime)
		hwindmag,										& ! m/s, horizontal wind magnitude (sqrt(u**2+v**2))
		windmag											  ! m/s, total wind magnitude (sqrt(u**2+v**2+w**2))

	real ::												&
		alpha												  ! N/A, Raupach factor to calculate w

	if(matching_grids /= 1) then
		! calc local_w
		alpha = 4.7*(GRAVITY*curr_cellvol/(p_t0*c_cp*RHO_AIR*zb))**ONE_THIRD	 ! w/(q^(1/3)) from a point Raupach
		wprime = p_w * alpha*(ww_q**ONE_THIRD - avg_q)    !perturbation  in w
		local_w = (wprime + avg_w)   !avg_comes from QUIC
	else
		local_w = avg_w
		wprime = 0.
	endif

	hwindmag = sqrt(avg_u**2 + avg_v**2)					! m/s, horizontal wind magnitude
	windmag = sqrt(avg_u**2 + avg_v**2 + local_w**2)	! m/s, total wind magnitude

	end
!======================================================================================
!======================================================================================
	subroutine ComputeTurbulence(ew,ns,ud,windmag,wprime,lengthscale,delta, &
		weather_sigma,k,mixing,k2sqrt)

	use grid_module
	use winds_module
	use fireca_constants_module

	implicit none


	real,parameter ::					&
		p_mixing_rate = 0.1			  ! N/A, controls turbulent mixing in mass loss rate calculation		

	integer,intent(IN) :: ew,ns,ud ! N/A, cell indexes

	real,intent(IN) ::				&
		windmag,							& ! m/s, total wind magnitude (sqrt(u**2+v**2+w**2))
		wprime,							& ! m/s, perturbation to the average vertical velocity calculated by quic
		delta,							& ! m, 4*(cell_volume)**(1/3)
		lengthscale						  ! m, flame length

	real,intent(OUT) ::				&
		k2sqrt,							&
		k,									&
		mixing,							& ! m2/s, turbulent viscosity
		weather_sigma

	real ::								&
		shear,							& ! (m/s)**2, wind shear
		calc_shear						  ! function to calculate the shear

	! -- wind shear in (m/s)^2
	shear = calc_shear(fcawinds%u,fcawinds%v,fcawinds%w,ew,ns,ud)
	! -- squared turbulent viscosity factor, in (m2/s)^2  =  (m^2) * (m/s)^2
	weather_sigma = min( max( (SMAG_COEFF*delta)**2*sqrt(shear), 0.001), 0.4*windmag)

	! Turbulent kinetic energy ((m2/s)^2 )
	k = max(.25,1.5 * (0.0001+weather_sigma*max(1., abs(wprime)/(windmag+.0001)))) !experimental slg
    !write (*,*) 'k=',k
	k2sqrt = sqrt(2.*k)

	! rrl burnout_time is hardcoded above and should be put in with the parameters in the QUIC_fire.inp.
	! This value should be a nominal time for a fine fuel element
	! Calc mixing term => NEED TO FIX THIS
	! m2/s, turbulent viscosity
	mixing = max(p_mixing_rate * min(2.,lengthscale) * sqrt(k) , 0.001)

	end
!======================================================================================
!======================================================================================
	subroutine UpdateBurnCenter(nbtime,fire_timedecay, &
		w_rr,w_rr_o2_turb,sumdecay,w_depl,					&
		sum_dcw_dcxy,sum_dcw_dcxy2,							&
		fire_burncenter_info,fire_burncenter,				&
		fire_spat_var,fuels_depl_var,fuels_depl_center)

	use grid_module
	use fireca_constants_module

	implicit none

	integer,intent(IN) :: nbtime

	real,intent(IN) ::		&
		w_rr,						& ! kg/m3/s, current reaction rate
		sumdecay,				& ! N/A, cumulative sum over fire%timedecay
		w_rr_o2_turb,			& ! ???
		w_depl

	real,dimension(nbtime),intent(IN) :: fire_timedecay

	real,dimension(nbtime,3),intent(INOUT) :: fire_burncenter_info
	
	real,dimension(2),intent(INOUT) ::		&
		fire_burncenter,							&
		fuels_depl_var,							&
		fuels_depl_center,						&
		fire_spat_var,								&
		sum_dcw_dcxy,sum_dcw_dcxy2

	real ::											&
		rr_sum,										& ! kg/m3/s, sum of the reaction rate of the previous times
		dcx, dcy, dcw,								& ! N/A, normalized centroid of ignitions, weighted by their reaction rate
		numerator,denominator,					&
		sum_diff,sum_wgts

	integer :: icount,ii

	real,dimension(2) :: dc

	! this flags the instances where cells are artificially ignited (such as initial ignition)
	! or the first time an ignition lands in a cell
	if (fire_burncenter_info(1,1) .le. 0. .and. &
		fire_burncenter_info(1,2) .le. 0.) then

		! a) Fill burncener_Info array  sum of all elements should be w_rr/w_rr_o2_turb
		! b) set initial burn center and variance in this cell
		! (different from subsequent determinations of these variables because there is no "past" centers and var

		! normalized centroid of ignitions, weighted by their reaction rate
		dcx = 0.0
		dcy = 0.0
		dcw = 0.0

		do icount = 1,nbtime
			fire_burncenter_info(icount,1) = fire_burncenter(1)
			fire_burncenter_info(icount,2) = fire_burncenter(2)
			fire_burncenter_info(icount,3) = w_rr*fire_timedecay(icount)/(w_rr_o2_turb*sumdecay)
		enddo
	   
		fuels_depl_center = fire_burncenter
		sum_dcw_dcxy = w_rr*fire_burncenter

		sum_dcw_dcxy2(1) = w_rr*fire_burncenter(1)**2
		sum_dcw_dcxy2(2) = w_rr*fire_burncenter(2)**2
	else

		! a) sum rr for elements 2 through n using updated decay and O2/turb
		rr_sum = sum(fire_burncenter_info(2:nbtime,3)) * w_rr_o2_turb
		if (rr_sum > w_rr) then
			dcw = 0.0
			fire_burncenter_info(:,3) = w_rr/ rr_sum * fire_burncenter_info(:,3)
			fire_burncenter_info(1,3) = 1e-6
		else
			! b) update reaction rate for first element of array
			fire_burncenter_info(1,3) = (1e-6 + w_rr - rr_sum) / w_rr_o2_turb
		endif

		! c) calculate burn center and variance using weighted average of elements 1 through n
		! rrl p_rem down to **1
		do ii = 1,2   !x and y directions
			numerator = sum(fire_burncenter_info(:,ii) * fire_burncenter_info(:,3))
			denominator = sum(fire_burncenter_info(:,3))
			fire_burncenter(ii) = numerator/denominator     !cell wide burn center

			sum_diff = 0.
			sum_wgts = 0.
			do icount = 1,nbtime
				sum_diff = sum_diff + (fire_burncenter(ii)-fire_burncenter_info(icount,ii))**2 * &
					fire_burncenter_info(icount,3)
				sum_wgts = sum_wgts + fire_burncenter_info(icount,3)
			enddo

			if(CELL_SIZE_MODE == CELL_SIZE_MODE_SMALL) then
				fire_spat_var(ii) = min(0.5*ONE_OVER_ROOTTHREE, &
					sqrt(sum_diff/(sum_wgts + 1e-6)))
			else
				fire_spat_var(ii) = min(0.5*ONE_OVER_ROOTTHREE, &
					max(0.5*ONE_OVER_ROOTTHREE*firegrid%reciprocal(ii), &
					sqrt(sum_diff/(sum_wgts + 1e-6))))
				
			endif
			! shifting burn center so that it is inside the cell
			fire_burncenter(ii) = max(fire_spat_var(ii)*rootthree, &
				fire_burncenter(ii))

			fire_burncenter(ii) = min(1.-fire_spat_var(ii)*rootthree, &
				fire_burncenter(ii))

		enddo

		! f) use discarded array element to update depleted center and variance
		if (w_depl .gt. 0.) then
			dc = 0.0
			dcw = 0.0

			do icount = 1,nbtime
				dc(1) = dc(1) + fire_burncenter_info(icount,1)*fire_burncenter_info(icount,3)
				dc(2) = dc(2) + fire_burncenter_info(icount,2)*fire_burncenter_info(icount,3)
				dcw = dcw + fire_burncenter_info(icount,3)
			enddo
			!dcw can be set to w_rr (in therory)
			dc = dc/dcw			
			dcw = dcw*w_rr_o2_turb

			sum_wgts = w_depl + dcw
			
			do ii = 1,2
				!working on the depletion center and var for x and y
				fuels_depl_center(ii) = ( w_depl*fuels_depl_center(ii) + &
					dc(ii)*dcw) / sum_wgts

				sum_diff = sqrt(max(w_depl*fuels_depl_center(ii)**2- &
					2.*sum_dcw_dcxy(ii)*fuels_depl_center(ii)+  &
					sum_dcw_dcxy2(ii),0.)+ &
					dcw*(dc(ii)-fuels_depl_center(ii))**2)

				fuels_depl_var(ii) = sqrt(sum_diff/sum_wgts)
				fuels_depl_var(ii) = min(0.5*ONE_OVER_ROOTTHREE,fuels_depl_var(ii))

				! shifting burn center so that it is inside the cell
				fuels_depl_center(ii) = max(fuels_depl_var(ii)*rootthree,fuels_depl_center(ii))

				fuels_depl_center(ii) =	min(1.-fuels_depl_var(ii)*rootthree, fuels_depl_center(ii))

				! preparing for next time steps calculation of depl variables in
				! this cell (excep that w_depl is updated at the end of the time step
				sum_dcw_dcxy(ii) = sum_dcw_dcxy(ii) + dcw*dc(ii)

				sum_dcw_dcxy2(ii) = sum_dcw_dcxy2(ii) + dcw*dc(ii)**2
			enddo
		endif
	endif

	

	end
!======================================================================================
!======================================================================================
	subroutine GenerateNewIgnitions(TID,nbtime,fire_ignitions,						&
		curr_cellvol,energy_to_fuels,w_rr_o2_turb,w_rr,fire_burncenter_info,		&
		w_ignitions,new_ignitions,new_smold_ignitions,supported_ignitions)

	use time_module
	use rnd_num_vars
	use fireca_constants_module

	implicit none

	integer,intent(IN) ::	&
		fire_ignitions,		& ! N/A, number of ignition per cell
		TID,						& ! N/A, threads counter
		nbtime

	real,dimension(nbtime,3),intent(IN) :: fire_burncenter_info

	real,intent(IN) ::									&
		energy_to_fuels,									& ! kW/m3, energy released to the fuel
		w_rr_o2_turb,										&
		w_rr,													& ! kg/m3/s, reaction rate in the cell
		curr_cellvol										  ! m3, cell volume

	integer,intent(INOUT) :: w_ignitions			  ! N/A, total number of ignitions in a cell (old + new)

	integer,intent(OUT) ::								&
		new_ignitions,										& ! N/A, how many flaming ignitions are generated with the available energy
		new_smold_ignitions,								& ! N/A, how many flaming ignitions are generated with the available energy
		supported_ignitions								  ! N/A, how many ignitions could be generated with this much

	real ::													&
		extra_energy,										& ! kW, energy remaining after the energy to fuel is partitioned over the ignitions
		new_ignitions_temp,								& ! N/A, how many flaming combustions ignitions are actually generated, assuming that only
																  !      the first 1/3 most recent ignited fuel is in flaming combustion
		new_smold_ignitions_temp						  ! N/A, how many creeping combustions ignitions are actually generated, assuming that only
															     !      the final 2/3 of burn time is contributing to this type of ignition

	integer :: ii

	! determine ignitions per cell and new ignitions
	call random_number(harvest,stat(TID))

	! -- energy remaining after the energy to fuel is partitioned over the ignitions, kW
	extra_energy = max(0.,mod(energy_to_fuels*curr_cellvol,ENERGY_PER_IGNITION))

	! -- how many ignitions could be generated with this much energy to the fuel with extra_energy assigned stocastically
	supported_ignitions = (int( energy_to_fuels*curr_cellvol/ ENERGY_PER_IGNITION )  + &
		int(extra_energy/ENERGY_PER_IGNITION + harvest(1)))

	! -- how many ignitions are actually generated, assuming that only the first 1/3 most recent ignited fuel is in flaming combustion
	new_ignitions_temp = 0.
	new_smold_ignitions_temp = 0.
	if (supported_ignitions > 0) then
		do ii = 1,int(BURNOUT_TIME/fcatime%dt_real*ONE_THIRD) 		
     			new_ignitions_temp = new_ignitions_temp +							&
				(real(supported_ignitions)*fire_burncenter_info(ii,3)*		&
				w_rr_o2_turb/(0.999*w_rr))			
		enddo
	endif
	new_ignitions = int(new_ignitions_temp)
	new_smold_ignitions = supported_ignitions - new_ignitions
	w_ignitions = w_ignitions + max(supported_ignitions , fire_ignitions)

	end
!======================================================================================
!======================================================================================
	subroutine FirePropagationCreeping(TID,it,ew_temp,ns_temp,ud_temp,	&
		partical_burn_vel_dt,start_coord,fire_spat_var,u,v,					&
		w_ew,w_ns,w_ud,landlocation,launchlocation,flame_components)
	
	use rnd_num_vars
	use fireca_constants_module
	use constants
	use grid_module

	implicit none

	integer,intent(IN) ::					&
		TID,									& ! N/A, thread counter
		it,										& ! N/A, time step counter
		ew_temp,ns_temp,ud_temp			  ! N/A, current rabbit position

	real,intent(IN) ::						&
		u,v,										& ! m/s, wind speed
		partical_burn_vel_dt

	real,intent(IN),dimension(2) ::		&
		fire_spat_var
	
	real,intent(IN),dimension(3) ::		&
		start_coord
		
	real,dimension(3),intent(OUT) ::		&
		landlocation,							& ! N/A, where the ignition lands [0,1]
		launchlocation,						& ! N/A, where the ignition is launched from [0,1] in the cell
		flame_components                   ! m, flame components in the x,y,z directions

	integer,intent(OUT) ::					&
		w_ew,w_ns,w_ud							  ! N/A, new cell where the ignitions land

	integer ::	ii								  ! N/A, counter

	real ::										&
		umag,										& ! m/s, perturbed local winds to compute the fire spread with the mean wind
		ax,										& ! rad, angles in the horizontal and vertical
		d,											& ! m, current flame length
		magnitude_scale


	! Fire propagation with local winds.
	! This is intended to capture the effects of the bimodal fireline dynamics seen in grass fires.
	call random_number(harvest3,stat(TID))

	! Perturbed local winds
	umag = sqrt(u**2 + v**2)

	ax = 2. * GREEK_PI * harvest3(1)

   magnitude_scale = 0.5 * (1. - (u*cos(ax) + v*sin(ax)) / (umag+1e-6))
   !magnitude_scale = .75
	
   !d = partical_burn_vel_dt *2.4*exp(-2.*fuels_density)* (1. - sqrt(harvest3(2))) * magnitude_scale
   d = partical_burn_vel_dt *1.* (1. - sqrt(harvest3(2))) * magnitude_scale
	flame_components(1) = d*cos(ax)*firegrid%dxi
	flame_components(2) = d*sin(ax)*firegrid%dyi
	flame_components(3) = 0.0

	if(it == 1) then
		call random_number(harvest2,stat(TID))
		do ii = 1,2
			launchlocation(ii) = start_coord(ii) + &
				(harvest2(ii)*2. - 1.)*fire_spat_var(ii)*ROOTTHREE
		enddo
		launchlocation(3) = start_coord(3)
	else
		launchlocation = start_coord
	endif

	landlocation(1) = mod(launchlocation(1) + flame_components(1),1.)
	if (landlocation(1).lt.0) landlocation(1) = 1. + landlocation(1)

	landlocation(2) = mod(launchlocation(2) + flame_components(2),1.)
	if (landlocation(2).lt.0) landlocation(2) = 1. + landlocation(2)

	landlocation(3) = launchlocation(3)

	w_ew = int(float(ew_temp) + launchlocation(1) + flame_components(1))
	w_ns = int(float(ns_temp) + launchlocation(2) + flame_components(2))

   !w_ud = 1+(float(ud_temp)-1)*nint(harvest3(3))
   !write (*,*) 'ud_temp',ud_temp,w_ud,harvest3(3)
	w_ud = int(float(ud_temp) + launchlocation(3) + flame_components(3))

	end 
!======================================================================================
!======================================================================================
	subroutine ComputeMoistureDepletion(fuels_moisture,fuels_moisture_initial, &
		fuels_moist_depl_var,fuels_moist_depl_center)

	use constants
	use fireca_constants_module
	use coord_module

	implicit none

	real,intent(IN) ::							&
		fuels_moisture,							& ! kg/m3
		fuels_moisture_initial					  ! kg/m3

	real,dimension(2),intent(INOUT) ::		&
		fuels_moist_depl_var,					&
		fuels_moist_depl_center

	real ::											&
		moist_depl_local,							&
		expansion_factor

	integer :: i

	moist_depl_local = (fuels_moisture_initial - fuels_moisture) / &
		(12.*fuels_moist_depl_var(1)*fuels_moist_depl_var(2))

	if (moist_depl_local > fuels_moisture_initial) then
		expansion_factor = sqrt(moist_depl_local/fuels_moisture_initial)

		do i = 1,2
			fuels_moist_depl_var(i) = fuels_moist_depl_var(i)*expansion_factor
			fuels_moist_depl_var(i) = min(0.5*ONE_OVER_ROOTTHREE , fuels_moist_depl_var(i))
		enddo
	endif

	! shifting moisture depletion center so that it is inside the cell
	do i = 1,2
		fuels_moist_depl_center(i) = max(fuels_moist_depl_var(i)*ROOTTHREE, &
			fuels_moist_depl_center(i))
		fuels_moist_depl_center(i) = min(1.-fuels_moist_depl_var(i)*ROOTTHREE, &
			fuels_moist_depl_center(i))
	enddo

	end
!======================================================================================
!======================================================================================
	subroutine ComputeNewCentroid(nbtime,sum_pos,sum_pos_squared,n_new_starts, &
		fire_burncenter_info,fire_burncenter,fire_spat_var)

	use constants
	use grid_module
	use fireca_constants_module

	implicit none

	integer,intent(IN) :: nbtime

	real,intent(IN) :: n_new_starts

	real,dimension(2),intent(IN) ::					&
		sum_pos,												&
		sum_pos_squared

	real,dimension(nbtime,3),intent(OUT) :: fire_burncenter_info

	real,dimension(2),intent(INOUT) :: fire_burncenter, fire_spat_var	

	integer :: ivar

	real ::													&
		w_new_spat_var,									&
		w_burncenter

	! e) calculate centroid of new ignitions and insert as first element of array
	do ivar = 1,2
		w_burncenter = sum_pos(ivar)/n_new_starts
		w_new_spat_var = sum_pos_squared(ivar)/n_new_starts - w_burncenter**2
		
		if(CELL_SIZE_MODE == CELL_SIZE_MODE_SMALL) then
			w_new_spat_var = max(w_new_spat_var,1e-6)
			w_new_spat_var = max(0.1*ONE_OVER_ROOTTHREE*firegrid%dxi,  &
			sqrt(w_new_spat_var)/n_new_starts)
		else
			w_new_spat_var = max(w_new_spat_var,0.)
			w_new_spat_var = max(0.5*ONE_OVER_ROOTTHREE*firegrid%dxi,  &
            sqrt(w_new_spat_var)/n_new_starts)
		endif
		
		w_new_spat_var = min(0.5*ONE_OVER_ROOTTHREE,w_new_spat_var)

		! shifting burn center so that it is inside the cell
		w_burncenter = max(w_new_spat_var*ROOTTHREE , w_burncenter)
		w_burncenter = min(1. - w_new_spat_var*ROOTTHREE , w_burncenter)

		fire_burncenter_info(1,ivar) = w_burncenter
		if(fire_burncenter(ivar) <= BURNCENTER_INIT) then
			fire_burncenter(ivar) = w_burncenter
			fire_spat_var(ivar) = w_new_spat_var
		endif

	enddo

	end
!======================================================================================
!======================================================================================
	subroutine SpreadFireSpotting(curr_cellvol,												&
		fuels_density,fuels_moisture,fuels_moisture_initial,fuels_moist_depl_var,	&
		fuels_moist_depl_center,ww_depl,ww_ignitions,sum_pos,								&
		sum_pos_squared,fuels_density_initial,TID,nn_new_starts,num_spotfires)

	use constants
	use rnd_num_vars
	use fireca_constants_module

	implicit none

	real,parameter :: fuels_moist_threshold = 0.25
	integer,intent(IN) ::		&
		num_spotfires,				&
		TID							  ! N/A, thread counter
	real,intent(IN) ::			&
		fuels_density,				& ! kg/m3, fuel density in the cell
		fuels_density_initial,	& ! kg/m3, initial fuel density in the cell (before ignition)
		fuels_moisture_initial,	& ! kg water/kg air, initial moisture in the cell (before ignition)
		curr_cellvol				  ! m3, cell volume

	real ::							&
		fuels_moisture,			&  ! kg water/kg air, moisture in the cell
		ww_depl,						&
		nn_new_starts
	
	integer ::						&
		ww_ignitions,				& ! N/A, number of ignitions in a cell
		new_ignitions_spotting, &
		ii								  ! N/A, counter
	real,dimension(2) ::			&
		fuels_moist_depl_var,	&
		fuels_moist_depl_center,&
		sum_pos,						&
		sum_pos_squared
	
	real ::							&
		ignition_success_prob,	&
		fuels_moist_old,			&
		lengthscale_area,			&
		moist_depl_local,			&
		current_moisture_loss,	&
		past_moisture_loss,		&
		past_moist_depl_center


	ignition_success_prob = 1.

	fuels_moist_old = fuels_moisture

	if (fuels_moisture .gt. MIN_MOISTURE) then
		ignition_success_prob = 0.
		lengthscale_area = 12.*fuels_moist_depl_var(1)*fuels_moist_depl_var(2)

		if (lengthscale_area .gt. 1) then
			write (msgoutfile,*) 'moisture area bigger than a cell', &
				lengthscale_area,fuels_moist_depl_var(1),fuels_moist_depl_var(2)
			call TerminateProgram()
		endif

		moist_depl_local = (fuels_moisture_initial - fuels_moisture) / lengthscale_area
		if (moist_depl_local .ge. fuels_moisture_initial) then
			ignition_success_prob = lengthscale_area
			moist_depl_local = fuels_moisture_initial
		else
			ignition_success_prob = lengthscale_area*max(0.0,  &
				(1.-(fuels_moisture_initial-moist_depl_local)/fuels_moist_threshold))
		endif

		fuels_moisture = max(0., &
			fuels_moisture - num_spotfires* (1.-ignition_success_prob)*&
			ENERGY_PER_IGNITION/(LATENT_HEAT_WATER_KG* &
			fuels_density*curr_cellvol))
	endif

	current_moisture_loss = fuels_moist_old - fuels_moisture
	past_moisture_loss = fuels_moisture_initial - fuels_moist_old


	do ii = 1,2
		past_moist_depl_center = fuels_moist_depl_center(ii)

		fuels_moist_depl_center(ii)= &
			(fuels_moist_depl_center(ii)* (past_moisture_loss+1.e-6) + &
			.5*current_moisture_loss) / &
			(current_moisture_loss+past_moisture_loss+1.e-6)

		fuels_moist_depl_var(ii)= &
			sqrt((current_moisture_loss* (.5*ONE_OVER_ROOTTHREE)**2+&
			(past_moisture_loss+1.e-6)*(fuels_moist_depl_var(ii))**2)/ &
			(current_moisture_loss+past_moisture_loss+1.e-6))

		fuels_moist_depl_var(ii) = min(fuels_moist_depl_var(ii), &
			0.5 *ONE_OVER_ROOTTHREE)
	enddo

	if (ww_depl .gt. -1.) then
		ignition_success_prob = ignition_success_prob * &
			(fuels_density_initial-fuels_density)/fuels_density_initial
	endif

	call random_number(harvest,stat(TID))

	new_ignitions_spotting = int(ignition_success_prob*float(num_spotfires)+harvest(1))
	ww_ignitions = ww_ignitions + new_ignitions_spotting
	nn_new_starts = nn_new_starts + real(new_ignitions_spotting)

	do ii = 1,2
		sum_pos(ii) = sum_pos(ii) + 0.5 * new_ignitions_spotting
		sum_pos_squared(ii) = sum_pos_squared(ii) + 0.5 * new_ignitions_spotting**2
	enddo

	end
!======================================================================================
!======================================================================================

