	subroutine SpreadFire(w_ew,w_ns,w_ud, TID,landlocation, launchlocation,	&
		fuel_height_launched,fuel_height_landed,										&					
		flame_components, l_gap, l_patch, curr_cellvol,								&
		fuels_density_initial,fuels_density,fuels_moisture_initial,				&
		fuels_moisture,fuels_depl_var,fuels_depl_center,n_new_starts,			&
		w_ignitions,nosubseqent_jump,ew_temp,ns_temp,ud_temp,sum_pos,			&
		sum_pos_squared,fuels_moist_depl_var,fuels_moist_depl_center,			&
		start_coord,num_ignitions_not_absorbed,										&
		energy_used_to_evaporate_water,fire_burncenter, fire_spat_var,			&
		function_switch,supported_ignitions, vert_wind, windmag, w_depl)
	
	use constants
	use rnd_num_vars
	use fireca_constants_module
	use time_module
	use grid_module

	implicit none

	real,parameter ::											&
		FUEL_MOIST_THRESHOLD = 0.0001/50.,				&
		FUELDENSITY_IGNITION_SUCCESS_CONTACT = 6.,	& 
		SAPV = 2/0.0005,										& !this needs to be put in the input conditions either as SAPV or sizescale which is what the firetec file is.  TBD so it is hardcoded for the moment.
		SAPV_NORM_INV = 1./SAPV,                     &
		windmag_norm_sqrd = 6.

	integer,intent(IN) ::								&
		w_ew,w_ns,w_ud,									&
		TID,													& ! N/A, thread counter
		function_switch, 									&
		supported_ignitions

	real,intent(IN),dimension(2) ::					&
		fuels_depl_var,									&
		fuels_depl_center,								&
		fire_spat_var
		
	real,intent(IN),dimension(3) ::					&
		landlocation,										& ! N/A, where the ignition lands [0,1] in the cell
		launchlocation,									&
		flame_components,									&		
		fire_burncenter
		
	real,intent(IN) ::									&		
		w_depl,												&
		fuels_density,										& ! kg/m3, fuel density in the cell
		fuels_density_initial,							& ! kg/m3, initial fuel density in the cell (before ignition)
		fuels_moisture_initial,							& ! kg water/kg air, initial moisture in the cell (before ignition)
		curr_cellvol,										& ! m3, cell volume
		windmag,												& ! m/s, horizontal wind 
		vert_wind,											& ! m/s, vertical wind component in this cell
		l_gap,												&
		l_patch,                                  &
		fuel_height_launched,fuel_height_landed	  ! m, height of the fuel in cells where energy originated and deposited
		
	real,intent(INOUT) ::								&
		fuels_moisture,									&  ! kg water/kg air, moisture in the cell
		n_new_starts

	integer,intent(INOUT) ::							&
		num_ignitions_not_absorbed,					&
		w_ignitions,										& ! N/A, number of ignitions in a cell
		ew_temp,ns_temp,ud_temp,						&
		nosubseqent_jump

	real,dimension(2),intent(INOUT) ::				&
		fuels_moist_depl_var,							&
		fuels_moist_depl_center,						&
		sum_pos,												&
		sum_pos_squared

	real,dimension(3),intent(INOUT) ::				&		
		start_coord

	real,dimension(2)::									&
		depl_vicinity,										&
		fire_vicinity,										&
		fuel_depl_vicinity,xa

	real ::													&
		ignition_success_prob,							&
		fuels_moist_old,									&
		creeping_adjustment,                      &
		lengthscale_area,									&
		previous_moist_depl_local,						&
		fuel_density_local,                       &
		current_moisture_loss,							&
		past_moisture_loss,                       &
		energy_used_to_evaporate_water,				&
		energy_evap_per_moisture_fraction,			&
		previous_num_ignitions_local_real,			&
		ignition_satuation_level,						&
		overshoot,											&
		perc_ignition,										&
		d,														&
		prob_hitting_patch,								&
		in_out_switch,										&
		temp,													&
		negposswitch,fuel_density_local_funct

	integer :: ii, computing_option,idir			 ! N/A, counter

	ignition_success_prob=1.   !assuming dense fuel and already dried

	if (function_switch == 2) then
		creeping_adjustment = .05
	else
		creeping_adjustment = 1.
	endif

	!compare fuel heights from launched to landed cells if in first vertical cell
	If (w_ud.eq.1) then
		if (fuel_height_landed.lt.fuel_height_launched)					&
			ignition_success_prob = ignition_success_prob*				&
			   fuel_height_landed / fuel_height_launched    !here we have assumed that if the energy packet is launced from first cell and is landing in first cell then it is moving parllel to the ground and it chance of hitting fuels depends on the random launch height within the launching fuel bed
	endif
	!check to see if in region where fuel has been depleted
	do ii = 1,2
		fuel_depl_vicinity(ii) = fuels_depl_var(ii)*ROOTTHREE -		&
			abs(landlocation(ii) - fuels_depl_center(ii))
	enddo

	lengthscale_area = 12.*fuels_depl_var(1)*fuels_depl_var(2)     !lengthscale_area is normalized by the area of the cell since var is already normalized by dx or dy  (it is a ratio of fire area to cell area)

	if (fuel_depl_vicinity(1) > 0.0 .and. fuel_depl_vicinity(2) > 0.0)  then  !inside region where fuels have been depleted

		fuel_density_local = max(0.0,(fuels_density_initial -						&
			(fuels_density_initial-fuels_density)/ lengthscale_area+1e-6))		

	else  !outside region where fuels have been depleted
		!the term that is added (but its value will be negative) here compensates for any fuel depletion that does not fit within the fuel depletion zone.  In theory if should be small to zero

		fuel_density_local = fuels_density_initial+(min(0.0,(fuels_density_initial-				&
			(fuels_density_initial-fuels_density)/(lengthscale_area+1.e-6))))*lengthscale_area	&
			/(1.-lengthscale_area + 1e-6)        !the term that is added (but its value will be negative) here compensates for any fuel depletion that does not fit within the fuel depletion zone.  In theory if should be small to zero
	endif

	if (fuel_density_local.le.0) then
		ignition_success_prob = 0.0
	else
       fuel_density_local_funct=.7*(fuel_density_local/.7)**.3333
		if (w_ud == 1) then
		!if (w_ud .le. 20) then
		   ignition_success_prob = ignition_success_prob *							&
			   (1. - exp(-fuel_density_local /						               &
			   FUELDENSITY_IGNITION_SUCCESS_CONTACT/creeping_adjustment*		&			   
				SAPV*SAPV_NORM_inv))*														&			   
			   fuel_density_local/ fuels_density_initial

      else
			 
			ignition_success_prob = ignition_success_prob *							&
			(1. - exp(-fuel_density_local_funct /						               &
			FUELDENSITY_IGNITION_SUCCESS_CONTACT* & !/creeping_adjustment*		&			   
			!1.*firegrid%z(w_ud)*    &
			SAPV*SAPV_NORM_inv*														&		! needs )) if remove the next line
			(windmag_norm_sqrd+max(0.0,vert_wind)**2)/windmag_norm_sqrd)) *	&	   
			fuel_density_local/ fuels_density_initial
			
			
			!ignition_success_prob = ignition_success_prob*                   &
          !(1.-exp(-fuel_density_local_funct / creeping_adjustment/         &
          !FUELDENSITY_IGNITION_SUCCESS_CONTACT*1.*firegrid%z(w_ud)*			&   !try replacing 1-exp with 1
          !SAPV*SAPV_NORM_inv*																&
          !(windmag_norm_sqrd+max(0.0,vert_wind)**2)/windmag_norm_sqrd)) *	&    !rrl_4_2 (changed inv to sqrd)
          !fuel_density_local/ fuels_density_initial
		endif

		call random_number(harvest3,stat(TID))

		do ii = 1,2
			fire_vicinity(ii) = fire_spat_var(ii)*ROOTTHREE -						&
				abs(landlocation(ii) - fire_burncenter(ii))
		enddo

		if (fire_vicinity(1) > 0.0 .and. fire_vicinity(2) > 0.0)  then  !inside active fire zone

			lengthscale_area = 12.*fire_spat_var(1)*fire_spat_var(2)

			!lengthscale_area is normalized by the area of the cell since var is already normalized by dx or dy  (it is a ratio of fire area to cell area)
			previous_num_ignitions_local_real = 										&
				real(supported_ignitions) / (lengthscale_area+1e-6)

			ignition_satuation_level=fuel_density_local*								&
				curr_cellvol*300.*															&
				CP_WOOD/(fcatime%dt_real*0.5*ENERGY_PER_IGNITION)

			ignition_success_prob = ignition_success_prob *							&
				exp(-previous_num_ignitions_local_real/ignition_satuation_level)
		endif

		if (l_gap.gt.0.0) then
			prob_hitting_patch=1
			computing_option=2

			if (computing_option.eq.1) then
				do idir = 1,2
					if (flame_components(idir).gt.0.0) then
						if (launchlocation(idir) 									&
							+flame_components(idir).gt.1) then 
							xa(idir)=mod(launchlocation(idir)+flame_components(idir),1.)
						else 
							xa(idir)=flame_components(idir)
						endif
					else
						if (1.-launchlocation(idir)-flame_components(idir).gt.1) then 
							xa(idir)=mod(1-launchlocation(idir)-flame_components(idir),1.)
						else 
							xa(idir)=-flame_components(idir)
						endif
					endif
				enddo
				d=sqrt((xa(1)*firegrid%dx)**2+(xa(2)*firegrid%dy)**2)
                
				if (mod(d-harvest3(1)*l_patch+l_patch+l_gap,l_patch+l_gap).lt.l_gap)  &
				    prob_hitting_patch=0.0
			elseif (computing_option.eq.2) then	
				do idir=1,2
					temp=launchlocation(idir)+flame_components(idir)
					negposswitch=(1-.5*(abs(temp)+temp)/temp)
					in_out_switch=float(int     &
						((float(int(abs(temp) +  negposswitch  )))    &
						/(float(int(abs(temp) +negposswitch ))-.000001)))
						
					negposswitch=.5*(flame_components(idir)+abs(flame_components(idir))) &
										/flame_components(idir)
					xa(idir)=negposswitch* &
						 mod(launchlocation(idir)*in_out_switch+flame_components(idir),1.)+ &
						 (1.-negposswitch)* &
						 mod(1-(launchlocation(idir)*in_out_switch+flame_components(idir)),1.)
				enddo
				d=sqrt((xa(1)*firegrid%dx)**2+(xa(2)*firegrid%dy)**2)
				temp=mod(d-harvest3(1)*l_patch+l_patch+l_gap,l_patch+l_gap)-l_gap
				prob_hitting_patch=.5*(abs(temp)+temp)/temp
			endif
			ignition_success_prob=ignition_success_prob*prob_hitting_patch
		endif

	endif
	call random_number(harvest3,stat(TID))

	if  (harvest3(1) + ignition_success_prob .lt. 1) then
		num_ignitions_not_absorbed = num_ignitions_not_absorbed + 1
		ignition_success_prob = 0.
	elseif (harvest3(1) + ignition_success_prob .ge. 1) then  ! ignition landed where we still have fuel and fire not saturated
		ignition_success_prob = 1.0

		if (fuels_moisture > MIN_MOISTURE) then

			fuels_moist_old = fuels_moisture
			past_moisture_loss = fuels_moisture_initial - fuels_moist_old

			! check to see if in area where moisture has already been driven off
			do ii = 1,2
				depl_vicinity(ii) = fuels_moist_depl_var(ii)*ROOTTHREE -			&
					abs(landlocation(ii) - fuels_moist_depl_center(ii))
			enddo
			
			lengthscale_area = 12.*fuels_moist_depl_var(1)*fuels_moist_depl_var(2)
			if (depl_vicinity(1) > 0.0 .and. depl_vicinity(2) > 0.0)  then  !inside moisture depletion zone
				!lengthscale_area is normalized by the area of the cell since var is already normalized by dx or dy  (it is a ratio of depeted area to cell area)

				previous_moist_depl_local = min(fuels_moisture_initial,			&
					(fuels_moisture_initial - fuels_moisture) / lengthscale_area)

			else  !if outside region where moisture has been partially evaporated
				!the term that is added (but its value will be negative) here compensates for any fuel depletion that does not fit within the fuel depletion zone.  In theory if should be small to zero

				overshoot=max(fuels_moisture_initial,			&
					(fuels_moisture_initial - fuels_moisture) / (lengthscale_area+1e-6))-    &
					fuels_moisture_initial

				previous_moist_depl_local = min(fuels_moisture_initial,			&
					(overshoot*lengthscale_area)/(1.-lengthscale_area+1e-6))
			endif

			if (previous_moist_depl_local+MIN_MOISTURE.lt.fuels_moisture_initial) then
            ignition_success_prob= &
               0.0*(1.-exp(-.005*previous_moist_depl_local/ &
                       (fuels_moisture_initial-previous_moist_depl_local+1e-6) * &
                       SAPV*SAPV_NORM_inv *windmag_norm_sqrd/(windmag_norm_sqrd+windmag)))   !rrl_4_2 changed inv to sqrd
				ignition_success_prob = 0.
			else
				ignition_success_prob = 1.
			endif


			if (previous_moist_depl_local+MIN_Moisture.lt.fuels_moisture_initial) then   
				!if there is still moisture where energy is deposited
				energy_evap_per_moisture_fraction =	LATENT_HEAT_WATER_KG * fuels_density * curr_cellvol  ! kW

				fuels_moisture = fuels_moisture - 									&
					(1.-ignition_success_prob)*fcatime%dt_real*					&
					ENERGY_PER_IGNITION/(energy_evap_per_moisture_fraction)

				If (fuels_moisture < 0.) then
					perc_ignition = -fuels_moisture*(energy_evap_per_moisture_fraction)/ &
						(fcatime%dt_real*ENERGY_PER_IGNITION)
					if (harvest3(3) + perc_ignition > 1) then
						ignition_success_prob = 1.

						!else consider if this needs to be sending energy to atm
					endif
					fuels_moisture = 0.0
				endif

				current_moisture_loss = fuels_moist_old - fuels_moisture

				energy_used_to_evaporate_water =											&
					energy_used_to_evaporate_water+current_moisture_loss *		&
					energy_evap_per_moisture_fraction
			else
				ignition_success_prob = 1
				current_moisture_loss = 0.0
			endif

			do ii = 1,2
				fuels_moist_depl_center(ii) =												&
					(fuels_moist_depl_center(ii) * (past_moisture_loss+1.e-6) + &
					landlocation(ii)*current_moisture_loss) /                   &
					(current_moisture_loss+past_moisture_loss+1.e-6)

				fuels_moist_depl_var(ii) =                                     &
					sqrt((current_moisture_loss*                                &
					(max(0.0001,landlocation(ii)-fuels_moist_depl_center(ii)))**2 +    &
					(past_moisture_loss+1.e-6)*(fuels_moist_depl_var(ii))**2)/	&
					(current_moisture_loss+past_moisture_loss+1.e-6))

				fuels_moist_depl_var(ii) = min(fuels_moist_depl_var(ii),       &
					0.5 * ONE_OVER_ROOTTHREE)
			enddo

		endif
	endif

	if (harvest3(2) + ignition_success_prob > 1) then
		w_ignitions = w_ignitions + 1
		n_new_starts = n_new_starts + 1.

		do ii = 1,2
			sum_pos(ii) = sum_pos(ii) + landlocation(ii)
			sum_pos_squared(ii) = sum_pos_squared(ii) + landlocation(ii)**2
		enddo

		start_coord = landlocation

		ew_temp = w_ew
		ns_temp = w_ns
		ud_temp = w_ud
	else
		nosubseqent_jump = 1
	endif

		end subroutine SpreadFire
!======================================================================================
!======================================================================================
	subroutine UpdateMoisture()

	use fireca_module
	use grid_module

	implicit none

	integer :: index

	!$OMP parallel do private(index)
	do index = 1,firegrid%num_fuel_cells
		! d) shift array elements down
		if (w_rr(index) > 0.0 .and. w_depl(index) <= 0.0) then
			fuels%depl_center(index)%val(1) = fire%burncenter(index)%val(1)
			fuels%depl_center(index)%val(2) = fire%burncenter(index)%val(2)
		endif

		if (n_new_starts(index) > 0.0) then
			if (fuels%moisture(index) > MIN_MOISTURE) then
				call ComputeMoistureDepletion(fuels%moisture(index),								&
					fuels%moisture_initial(index),fuels%moist_depl_var(index)%val,				&
					fuels%moist_depl_center(index)%val)
			endif
			call ComputeNewCentroid(nbtime,firecell(index)%sum_pos,firecell(index)%sum_pos_squared,	&
				n_new_starts(index),fire%burncenter_info(index)%val2,								&
				fire%burncenter(index)%val,fire%spat_var(index)%val)
		endif
	enddo
	!$OMP end parallel do

	end subroutine UpdateMoisture
!======================================================================================
!======================================================================================
	subroutine ComputeReactionVariables(firegrid_cellarea,lsa,cellvol,windmag,	&
		avg_rr,local_w,mixing,fuel_dens,fuel_dens0,ignitions,							&
		wstar,burnfraction,w_fraction_to_atmos,net_energy_rate,						&
		w_rr,w_rr_o2_turb,w_fueldensity,w_mass,energy_to_fuels,						&
		lengthscale,o2density,mass_burned,fuelheight,firegrid_dz,TID)

	use time_module
	use fireca_constants_module	
	use constants
	
	implicit none

	real,parameter ::									&
		O2_FLOOR = 0.1,								& ! N/A, fraction of available O2 that can be consumed
		RHO_HYDRO_THRESHOLD = 0.4,					& ! ??, the cf values were taken from FIRETEC
		CF_HYDRO = 32.,								&
		CF_CHAR = 10.,									&
		CHARRING_FRACT_TO_FUEL = 0.75,			& ! N/A, fraction of fire output that goes to fuels
		FLAMING_FRACT_TO_FUEL = 0.25,				& ! N/A, fraction of fire output that goes to atmosphere
		P_O2 = 2.0                               ! s*m3/kg, has to be the inverse of the reaction rate to calculate the oxygen density,
	! o2density = (X_O2MAX-O2_FLOOR)*exp(-P_O2*avg_rr) + O2_FLOOR
	
	integer,intent(IN) :: TID

	real,intent(IN) ::								&
		windmag,											& ! m/s, wind speed
		cellvol,											& ! m3, cell volume
		firegrid_cellarea,							& ! m2, horizontal cell area
		local_w,											& ! m/s, local vertical velocity = average + perturbation (wprime)
		mixing,											& ! m2/s, turbulent viscosity
		fuel_dens,										& ! kg/m3, fuel density
		fuel_dens0,										& ! kg/m3, fuel density before ignition
		avg_rr,											& ! kg/m3/s, average reaction rate when the quic and
		fuelheight,										& ! m, fuel height
		firegrid_dz,									& ! m, height of firegrid cell
		lsa												  ! length scale area

	integer,intent(IN) :: ignitions				  ! N/A, number of ignitions

	real,intent(INOUT) :: o2density				  ! kg/m3, density of oxygen close to the combustion zone

	real,intent(OUT) ::					&
		w_mass,								& ! kg, mass in the cell
		wstar,								& ! m/s
		burnfraction,						& ! N/A, fraction of mass in a cell burning because of the ignition energy
		w_fraction_to_atmos,				& ! N/A, fraction of energy to the atmosphere
		net_energy_rate,					& ! KW/m3, energy released by the combustion
		energy_to_fuels,					& ! KW/m3, energy released by the combustion that goes to the fuel
		lengthscale,						& ! m
		w_rr,									& ! kg/m3/s, current reaction rate
		w_rr_o2_turb,						& ! added to track o2 and turbulence contribution to reaction rate
		w_fueldensity,						& ! kg/m3, current fuel density
		mass_burned							  ! g, mass burned

	real ::									&
		wind2d_sq,							& ! (m/s)^2 , u**2 + v**2
		wind3d,								& ! m/s^2 , wind speed including all 3 components
		lambda,								& ! N/A, stoichiometry ratio for 1 unit of product
		n_e_per_ign,						& ! kW, total energy of ignitions in a cell
		perchydroremaining,				& ! N/A, percentage of hydrocarbons remaining
		cf,									& ! N/A
		w_fraction_to_fuels				  ! N/A, fraction of energy to the fuel


	! -- density of oxygen close to the combustion zone, kg/m3
	! In theory this still needs a rho_air term in front of it to make the units correct.
	! If density of gas is close to 1 then minimal impact
	o2density = RHO_AIR * ((X_O2MAX-O2_FLOOR)*exp(-P_O2*firegrid_cellarea*avg_rr*0.5448/(mixing+1.e-8)) + O2_FLOOR)

	! -- From Drysdale, "Fire Dynamics", stoichiometry ratio to use in the calculation of the reaction rate
	! 0.4552 (1/=2.1968366) and 0.5448 (1/=1.835536) are stoichiometry coefficients from Drysdale ("Fire Dynamics")
	! These coefficients are the stoichiometric coefficients for wood and oxygen respectively for 1 unit of product
	lambda = fuel_dens * o2density / (fuel_dens*2.1968366 + o2density*1.835536)**2

	! -- Total energy of ignitions in a cell, kW
	n_e_per_ign = real(ignitions)*ENERGY_PER_IGNITION

	! -- Fraction of mass in a cell burning because of the ignition energy
	burnfraction = min(1.0 , max(0.0 , n_e_per_ign)/ &
		(fuel_dens*cellvol*NET_ENERGY_WOOD_KG/BURNOUT_TIME))  !rrltoday

	! -- m/s
	wstar = max(local_w , max(0.,local_w/(burnfraction**ONE_THIRD+0.0001)))
	wstar = max(0.,local_w)
	! write (*,*) 'in computereactionvar', wstar, local_w, burnfraction**ONE_THIRD, &
	!     n_e_per_ign, ignitions,fuel_dens*cellvol*NET_ENERGY_WOOD_KG/BURNOUT_TIME

	! Calc updated reaction rate and updated fuel density
	perchydroremaining = max(0.,(fuel_dens-RHO_HYDRO_THRESHOLD*fuel_dens0) / &
		(fuel_dens0*(1.-RHO_HYDRO_THRESHOLD)))
	!!!rrl  this might be unneeded or even double counting now with the burndout paradigm
	cf = CF_HYDRO*perchydroremaining + CF_CHAR*(1.-perchydroremaining)
	! -- fraction of energy to fuel
	w_fraction_to_fuels = CHARRING_FRACT_TO_FUEL*(1.-perchydroremaining) + &
		FLAMING_FRACT_TO_FUEL*perchydroremaining
	! -- fraction of energy to the atmosphere
	w_fraction_to_atmos = 1. - w_fraction_to_fuels

	! -- current reaction rate, kg/m3/s
	! rrl try limiting condition here to avoid needing to below
	w_rr = max(0., min(fuel_dens/fcatime%dt_real, &
		cf * fuel_dens * o2density * burnfraction * mixing * lambda ))
	o2density = RHO_AIR * ((X_O2MAX-O2_FLOOR)*exp(-P_O2*firegrid_cellarea*w_rr*0.5448/(mixing+1.e-8))+ O2_FLOOR)
	w_rr = max(0., min(fuel_dens/fcatime%dt_real, &
		cf * fuel_dens * o2density * burnfraction * mixing * lambda ))
	w_rr_o2_turb = o2density * mixing * lambda

	! -- current fuel density, kg/m3
	w_fueldensity = fuel_dens - fcatime%dt_real * w_rr

	! -- energy released by the combustion, KJ/m^3/s
	!!!rrl consider using net_energy_rate instead of w_rr in burncenter_info array
	net_energy_rate = w_rr * NET_ENERGY_WOOD_KG

	! -- energy to the fuel, kJ/m3/s
	! This radiation loss should probably become an exponential function at some point
	energy_to_fuels = net_energy_rate * (1. - .5* RADIATON_LOSS_FRACT_COEFF)

	! -- recalculate flame length with current reaction rate
	if(w_fueldensity > 0) then
		wind2d_sq = windmag**2
		wind3d = sqrt(windmag**2 + local_w**2)
		call fire_length_scale_calc(wind2d_sq, wind3d,								&
				lsa,ignitions,o2density,fuelheight,							&
				firegrid_cellarea,firegrid_dz, TID,lengthscale)		
	else
		lengthscale = 0.
	endif
	! -- update mass in cell
	w_mass = w_fueldensity * cellvol

	!get the mass burned in GRAMS this time step and cell
	mass_burned = w_rr * cellvol * fcatime%dt_real * 1000.

	end subroutine ComputeReactionVariables
!======================================================================================
!======================================================================================
	subroutine FirePropagationMainWind(TID,it,ew_temp,ns_temp,ud_temp,	&
		lengthscale,k2sqrt,wstar,start_coord,fire_spat_var,					&
		w_ew,w_ns,w_ud,d,landlocation,launchlocation,flame_components,		&
		fuel_height_launched)

	use rnd_num_vars
	use constants
	use grid_module
	use winds_module
	use fireca_constants_module
	use time_module

	implicit none

	integer,intent(IN) ::						&
		TID,											& ! N/A, thread counter
		it,											& ! N/A, time step counter
		ew_temp,ns_temp,ud_temp

	real,intent(IN) ::							&
		fuel_height_launched,					& ! m, height of the fuel in the cell
		lengthscale,								& ! m, flame length
		k2sqrt,										&
		wstar											  ! m/s

	real,dimension(2),intent(IN) ::			&
		fire_spat_var
	
	real,dimension(3),intent(IN) ::			&		
		start_coord		

	real, intent(OUT) :: d
	real,dimension(3),intent(OUT) ::			&
		landlocation,								& ! N/A, where the ignition lands [0,1]
		launchlocation,							& ! N/A, where the ignition is launched from [0,1] in the cell
		flame_components

	integer,intent(OUT) ::						&
		w_ew,w_ns,w_ud 							  ! N/A, new cell where the ignitions land

	integer :: ii									  ! N/A, counter

	real ::											&
		u_local,v_local,							& ! m/s, perturbed local winds to compute the fire spread with the mean wind
		ah,ax,										& ! rad, angles in the horizontal and vertical
		sinah,sinax,cosah,cosax,				& ! N/A, trigonometry
		ahwindcomp,									&
		lengthscale_area_inv,					&
		minimum_area
	
	real :: erfinv

	! Fire propagation with local winds.
	! This is intended to capture the effects of the bimodal fireline dynamics seen in grass fires.
	call random_number(harvest3,stat(TID))

	minimum_area = 2./firegrid%cellarea

	lengthscale_area_inv = 1./max(minimum_area,12.*fire_spat_var(1)*fire_spat_var(2))
	
	! Perturbed local winds
	u_local = fcawinds%u(ew_temp,ns_temp,ud_temp) + k2sqrt * erfinv(2.*harvest3(1)-1.)*0.5
	v_local = fcawinds%v(ew_temp,ns_temp,ud_temp) - k2sqrt * erfinv(2.*harvest3(2)-1.)*0.5

	if (ud_temp.eq.1) then
		u_local = u_local*fuel_height_launched/firegrid%dz_array(1)
		v_local = v_local*fuel_height_launched/firegrid%dz_array(1)
  	endif

	ah = HORIZONTAL_VARIATION*aint(max(0.0,   &
		1. + lengthscale_area_inv*wstar/          &
		(lengthscale_area_inv*wstar+sqrt(u_local**2+v_local**2)) - harvest3(3)))
	
	ahwindcomp = (ah*(lengthscale_area_inv*wstar)+  &
		(HORIZONTAL_VARIATION-ah)*sqrt(u_local**2+v_local**2)) &
		/(HORIZONTAL_VARIATION* &
		sqrt((lengthscale_area_inv*wstar)**2  &
		+u_local**2+v_local**2))

		!write (*,*) wstar*lengthscale_area_inv, sqrt(u_local**2+v_local**2),ah, &
		 !  lengthscale_area_inv,fuel_height_launched,firegrid%dz_array(1)
	ax = atan2(v_local, u_local)

	sinah = sin(ah)
	cosah = cos(ah)
	sinax = sin(ax)
	cosax = cos(ax)

	call random_number(harvest3,stat(TID))

	d = lengthscale * ahwindcomp * (1. - sqrt(1-harvest3(3)))

	! determine x', y', z' and determine impacted grid cell
	flame_components(1) = d*cosah*cosax*firegrid%dxi
	flame_components(2) = d*cosah*sinax*firegrid%dyi
	flame_components(3) = d*sinah/firegrid%dz_array(ud_temp)

	if(it == 1) then
		do ii = 1,2
			launchlocation(ii) = start_coord(ii) + &
				(harvest3(ii)*2. - 1.)*fire_spat_var(ii)*ROOTTHREE
		enddo
		launchlocation(3) = start_coord(3)
	else
		launchlocation = start_coord
	endif

	landlocation(1) = mod(launchlocation(1) + flame_components(1), 1.)
	if (landlocation(1) < 0.) landlocation(1) = 1. + landlocation(1)

	landlocation(2) = mod(launchlocation(2) + flame_components(2),1.)
	if (landlocation(2) < 0.) landlocation(2) = 1. + landlocation(2)

	landlocation(3) = mod(launchlocation(3) + flame_components(3),1.)	

	w_ew = int(float(ew_temp) + launchlocation(1) + flame_components(1))
	w_ns = int(float(ns_temp) + launchlocation(2) + flame_components(2))
	w_ud = int(float(ud_temp) + launchlocation(3) + flame_components(3))
	
	end
!======================================================================================
!======================================================================================
	subroutine fire_length_scale_calc(windmag2d_sq,windmag3d,lsa,n,o2density,fuelheight, & 
		firegrid_cellarea,dz,TID,fire_length_scale)

	! This function contains three options for calculating the final result (fire_length_scale).
	! - The first option is from Scott Goodric, and is experimental for forest type fuels,
	!   but inappropriate for large urban fires.
	! - The second option has two portions: 1) the mass flow rate of volatiles in the system
	!	 per unit cell volume multiplied by the ratio of the bulk volume of a burning collection
	!   of fines (a bush, branch, office chair, etc.) divided by the surface are of the fines.
	!   Concetually, for a give flow rate of volatiles/cell_volume, the larger the bulk volume, the
	!   the bulk volume, the more flame merging effects we have to lengthen the flame, and the
	!   smaller the surface are of the fines the more rapid the devolatilization (for a give
	!	 devolatilization rate). The second part of the second equation is a balance between
	!   momentum and turbulance, normalized by ambient air density.
	! - The last equation is adapte from Stephen Turn's Itro to Combustion (2nd editing page 495-500).
	!   This equation onceptually balances impacts of buoyancy, initial momentum, turbulence,
	!   stoichiometry, and volatile density, and applies empirical corellations and scaling
	!   arguments to incorporate flow and combustion characteristics.

	use constants
	use fireca_constants_module
	use rnd_num_vars

	implicit none

	real,parameter ::							&
		RHO_GAS_RATIO = 0.813,				& ! N/A, ratio of heated volatile density to ambient air density, assuming
													  !      T_ambient=300 K, T_pyrolysis=474 K, and mw of volatiles is 37.26 gm/mol
		TEMP_GRAD_RATIO = 6.51,				& ! N/A, The ratio of the temperature gradient to ambient temperature
													  ! assuming ambient, stoichiometric flame is 2253 K (2253-300)/300
		CONST_A = 0.04648,					& ! moles/g, moles of O2 needed to consume 1 gram of fuel assuming douglas fir
		MIN_HEIGHT = 0.05,					& ! m, minimum fuel scale
		FUEL_HEIGHT_NORM = 30.				  ! N/A, normalization factor in the exponential decay of fuel height

	integer,intent(IN) :: TID, n            ! n = number of ignitions

	real, intent(IN) ::						&
		lsa,										& ! lengthscale area
		o2density,								& ! kg/m^3,density of O2
		fuelheight,								& ! m, height of the fuel
		firegrid_cellarea,					& ! m2, cell area
		dz,										& ! m, height of a firegrid cell
		windmag3d, 								& ! m/s, wind speed
		windmag2d_sq							  ! ! (m/s)^2, horizontal wind speed squared

	real,intent(OUT) ::						&
		fire_length_scale						  ! m
	
	real ::										&
		froud,									& ! N/A, froude number
		Lnd,										& ! N/A
		stoich,									& ! N/A, stoichiometric equivilance ratio needed for the flame length calculation
		b,											& ! moles/g, moles of N2 per mole of O2 in depleted air times CONST_A
		c,											& ! moles/g, moles of CO2 per mole of O2 in depleted air times CONST_A
		depletion,								& ! N/A, fraction of O2 remaining in ambient air
		ma,										& ! g/g, mass of air in grams per gram of fuel
		x_fuel,									& ! N/A, random variable
		fuel_length_scale,					& ! m
		fuel_height_var,						& ! m, fuel height used in the calculations
		I, X,										& ! kW m^-1, fireline intensity		
		wi

	integer,parameter ::  FuelScaleChoice=4    ! 1=constant fraction of fuel height
	! 2=random fraction of fuel height
	! 3=specified distribution of fuel heights used for Atlanta test case
	! 4=flame height based on fire intentisty following Nelson, Butler & Weise 2012 IJWF


	! Integrate in the impacts of depletion and calculate the stoichiometric equivalence ratio
	depletion = min(o2density / INITIAL_O2DENSITY , 1.0) ! The fraction of O2 remaining in ambient air
	b = (0.79/(0.21*depletion))*CONST_A  ! The moles of N2 per mole of O2 in depleted air times CONST_A
	c = CONST_A*0.21*(1. - depletion) ! The moles of CO2 per mole of O2 in depleted air times CONST_A
	ma = 32. * CONST_A + 28. * b + 44. * c ! The mass of air in grams per gram of fuel
	stoich = 1./(ma+1.) ! The stoichiometric equivilance ratio needed for the flame length calculation

	! Calculate the length scale for this flame
	! (this covers individual flames on various fuel sizes and clumps)
	fuel_height_var = min(fuelheight,dz)

	if ( FuelScaleChoice < 4 ) then
		call random_number(harvest,stat(TID))
		x_fuel = harvest(1)

		if ( FuelScaleChoice==1 ) then
			fuel_length_scale = 0.5*fuel_height_var
		elseif ( FuelScaleChoice==2 ) then
			fuel_length_scale = 0.25 * fuel_height_var * (1.+2.*x_fuel)
		elseif ( FuelScaleChoice==3 ) then
			if (NINT(x_fuel) == 0 .or. fuelheight <= 1.0) then
				fuel_height_var = fuelheight
			else
				fuel_height_var = 0.1
			endif
			fuel_length_scale = MIN_HEIGHT - log(x_fuel)*min(fuel_height_var,dz) / FUEL_HEIGHT_NORM ! meters
		endif

		froud = windmag3d*stoich**1.5/(RHO_GAS_RATIO**0.25*sqrt(TEMP_GRAD_RATIO*GRAVITY*fuel_length_scale)) ! N/A
		Lnd = 13.5 * froud**0.4 / (1. + 0.07*froud**2)**0.2

		fire_length_scale = max(0.01 , Lnd*fuel_length_scale*sqrt(RHO_GAS_RATIO)/stoich)
	else
		if(FuelScaleChoice == 4) then			
			I = ENERGY_PER_IGNITION * real(n) / max(0.2, sqrt(lsa*firegrid_cellarea))       ! Estimated fireline intensity in kW m^-1, the length in the denominator must be greater than 5 cm
			X = 1. + 0.0155 * I**0.667         ! flame height from Nelson, Butler & Weise eq in fig 5b			
			wi = min(1.,dz*0.1)*0.377*I**ONE_THIRD					 ! equation 18 in Nelson et al. provides a charcteristic vertical velocity for the whole fire, factor in front is to adjust for velocity close to ground
			fire_length_scale =  X * max(1., ROOTTHREE * (windmag2d_sq/wi**2)**0.25 )
		elseif(FuelScaleChoice == 5) then
			I = ENERGY_PER_IGNITION * real(n) / (RHO_AIR*1.004*AMBIENT_TEMPERATURE*sqrt(GRAVITY))
			X = windmag2d_sq/GRAVITY	
			fire_length_scale = sqrt((X+1.)*I**0.4)
		endif
	endif
	
	end subroutine fire_length_scale_calc
!======================================================================================
!======================================================================================
	subroutine lengthscale_area_calc(spat_var, cellarea, lengthscale_area)
		
		implicit none
		
		real, dimension(2),intent(IN) :: spat_var
		real,intent(IN) :: cellarea
		real,intent(OUT) :: lengthscale_area
				
		lengthscale_area = max(0.01, min(1., 12.*spat_var(1)*spat_var(2)))		
		
	end subroutine lengthscale_area_calc
!======================================================================================
!======================================================================================
	SUBROUTINE SetParamsFireSim()
	
	use plume_const_module
	use fireca_constants_module
	use grid_module
	
	implicit none
	
	plume_const%do_w_fill = 1
	moist_depl_var_init = 0.0001
	init_var_length = min(2., min(firegrid%dx, firegrid%dy))
	CELL_SIZE_MODE = CELL_SIZE_MODE_SMALL
		
	END 
!======================================================================================
!======================================================================================
	SUBROUTINE CheckVer(invar)
	
	implicit none
	
	integer, intent(IN) :: invar
	
	if(invar /= 1)then
		print*, 'Incorrect compilation mode'
		call TerminateProgram()
	endif
	
	END 
!======================================================================================
!======================================================================================

	
