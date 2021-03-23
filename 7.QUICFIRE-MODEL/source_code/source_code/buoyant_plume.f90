	subroutine buoyant_plume(fire_energy_rate)

	use time_module
	use file_handling_module
	use sor_module
	use flags_module
	use plume_module
	use plume_const_module
	use grid_module
	use winds_module
	use time_module
	use constants

	implicit none

	! 3D array that stores energy that fire emits to the atmosphere the array is stored on the FIRE grid
	real, dimension(firegrid%nx, firegrid%ny, firegrid%nz_en2atmos),intent(IN) :: fire_energy_rate	
	integer :: i,j,k
	
	! ---- Initialize variables
	to_delete = 0   ! no plumes to delete initially
		
	! Clear fire from previous time
   quwinds%wplume = 0.
	
	! --- Initialize time 	
	plume_const%total_time = 0.	
	
   ! ---- Generate the new plumes
	if(plume_const%time_step_flag == 1) then
		plume_const%time_step = fcatime%dt_real
	endif
	call InitializeNewPlumes(fire_energy_rate)
	
	!do it = 1, n_plume_steps
	do while(plume_const%total_time < fcatime%dt_real)
		plume_const%total_time = plume_const%total_time + plume_const%time_step
		write(*, '(a, f8.2, a, f8.2)') 'Plume time step ', plume_const%time_step, '; Plume tot time ', plume_const%total_time
		
		! For stable conditions, check if max height has been reached
		if(plume_const%BRUNT_VAISALA_FREQ2 /= 0) then
			call CheckMaxHeight()
			if(sum(to_delete(1:total_fires)) > 0) call DeletePlumes()
		endif
		
		! Move the plumes with the wind and pull them together
		call MovePlumes() ! after this step, plumes are at k+1
		if(sum(to_delete(1:total_fires)) > 0) call DeletePlumes()
		
      if(total_fires > 0)then
			call UpdatePlumeProperties()
         if(sum(to_delete(1:total_fires)) > 0) call DeletePlumes()
		endif
								
		! Merge plumes
		if(total_fires > 1)then
			call MergePlumes()
			if(sum(to_delete(1:total_fires)) > 0) call DeletePlumes()
		endif
		
		! Compute W 
		call WWIndexCalc()
		call WWCalc()
		
		! Determine plume time, kfire if the plume has stopped rising
      if(flag%plume_traj_file > 0) then			
			do i = 1,total_fires				
				call WriteOutPlume(plume(i), 0)
         enddo
		endif
		
		call PlumeFinalization()
		if(sum(to_delete(1:total_fires)) > 0) call DeletePlumes()
		
		
		do i = 1,total_fires
			plume_diagnostics%max_w = max(plume_diagnostics%max_w, plume(i)%traj_vel(3))
		enddo
		
		!write(ID_FILE_TIMELOG,*)'t =  ', fcatime%current_int - fcatime%dt_int+it*plume_const%time_step, &
		!	'  # plumes: ', total_fires
	enddo
	
	sor%alpha2_fire = sor%alpha2
	sor%denom = sor%denom0
   !$OMP parallel do private(i,j,k)
	do k = 1,qugrid%nz
      do j = 1,qugrid%ny-1
         do i = 1,qugrid%nx-1
				if(quwinds%wplume(i,j,k) /= 0) then
					quwinds%wplume(i,j,k) = quwinds%wplume(i,j,k)**plume_const%INV_WC_EXPONENT
					sor%alpha2_fire(i,j,k) = sor%ALPHA2_FIRE_VAL
					sor%denom(i,j,k) = sor%omegarelax/			&
						(sor%bc%e(i,j,k) + sor%bc%f(i,j,k) + sor%bc%g(i,j,k) + sor%bc%h(i,j,k) +			&						
						(sor%bc%m(i,j,k) + sor%bc%n(i,j,k)) * sor%alpha2sq / sor%alpha2_fire(i,j,k)**2 )
				endif
         enddo
      enddo
	enddo
	!$OMP end parallel do 
		
	if(total_fires == 0)then
		call ResetSORMCVars(0)		
	else
		sor%mc_is = qugrid%nx+10
		sor%mc_ie = -1
		sor%mc_js = qugrid%ny+10
		sor%mc_je = -1
		sor%mc_ks = qugrid%nz+10
		sor%mc_ke = -1
		do i = 1,total_fires
			sor%mc_is = min(sor%mc_is,plume(i)%istart)
			sor%mc_ie = max(sor%mc_ie,plume(i)%iend)
			sor%mc_js = min(sor%mc_js,plume(i)%jstart)
			sor%mc_je = max(sor%mc_je,plume(i)%jend)
			sor%mc_ks = min(sor%mc_ks,plume(i)%kstart)
			sor%mc_ke = max(sor%mc_ke,plume(i)%kend)
		enddo
		sor%mc_is = max(sor%mc_is-sor%mc_buffer,1)
		sor%mc_ie = min(sor%mc_ie+sor%mc_buffer,qugrid%nx)
		sor%mc_js = max(sor%mc_js-sor%mc_buffer,1)
		sor%mc_je = min(sor%mc_je+sor%mc_buffer,qugrid%ny)
		sor%mc_ks = max(sor%mc_ks-sor%mc_buffer,1)
		sor%mc_ke = min(qugrid%nz+2,sor%mc_ke+sor%mc_buffer)
	endif
	
	return
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE ResetSORMCVars(reset_alpha2)
	
	use sor_module
	use grid_module
	
	implicit none
	
	integer, intent(IN) :: reset_alpha2
	sor%mc_is = 1
	sor%mc_ie = qugrid%nx
	sor%mc_js = 1
	sor%mc_je = qugrid%ny
	sor%mc_ks = 1
	sor%mc_ke = qugrid%nz+2

	if(reset_alpha2 == 1) sor%alpha2_fire = sor%alpha2
	
	end	
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE ComputePerimeter(radius_h, gamma, perimeter, perimeter_over_rph)
	
	use constants
	
	implicit none
	
	real, intent(IN) :: radius_h, gamma
	real, intent(OUT) :: perimeter, perimeter_over_rph
	
	perimeter_over_rph = PI  * (3. * (1. + gamma) - sqrt((3.*gamma + 1.) * (3. + gamma)) )
	perimeter = radius_h * perimeter_over_rph
	
	END
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE InitializeNewPlumes(fire_energy_rate)
 
	use plume_module
	use constants
	use mersenne_twister
	use rnd_num_vars
	use time_module
	use file_handling_module
	use flags_module
	use plume_const_module
	use interpolation_module	
	use grid_module
	use fireca_module
	
	implicit none
	
	real, parameter ::				&
		cp = 1.004,						& ! kJ/kg/K, air heat capacity
		air_temp = 300					  ! K, air temperature		
		
	real, dimension(firegrid%nx, firegrid%ny, firegrid%nz_en2atmos),intent(IN) :: fire_energy_rate
	
	real ::								&
		uc, vc, wc,						& ! m/s, wind components at the center of the QU cell
		xc, yc, zc,						& ! m, plume centerline coordinates
		ztot,								& ! m, cumulative z
		ws,								& ! m/s, wind speed at the center of the QU cell
		sum_fire_energy_rate,		& ! energy provided by the fire		
		mult,								&		
		one_over_z						  ! 1/m, 1/z
	
	integer ::							&
		is_replace,						&
		curr_plume,						& ! N/A, plume to replace in the arrays
		niter,							& ! N/A, max number of iterations in the do while
		iii,jjj,kkk,					& ! N/A, counter
		i,j,k,							& ! N/A, counter
		old_fire_num,					& ! N/A, total number of plumes before eliminating some
		i1,j1,k1,						& ! N/A, loop start index
		init_val							  ! N/A, index used to decide where to start the loops
 
		
	mult = GRAVITY / (air_temp * cp * RHO_AIR)  ! Fb = g/T *(H/rho/cp)  => Fb = mult * H * volume   (dvol used to convert units)		
	
	i1 = loopidx%list(loopidx%list_count,1)
	j1 = loopidx%list(loopidx%list_count,2)
	k1 = loopidx%list(loopidx%list_count,3)
	
	do k = loopidx%kstart(k1),loopidx%kend(k1),loopidx%increm(k1)
      old_fire_num = total_fires
		one_over_z = 1./ qugrid%dz_array(k)
		do j = loopidx%jstart(j1),loopidx%jend(j1),loopidx%increm(j1)
			do i = loopidx%istart(i1),loopidx%iend(i1),loopidx%increm(i1)
				sum_fire_energy_rate = 0.0
 
				! Compute the energy emitted to the atmosphere by the fire	(on the quic grid) and count number of new fires
				do kkk = fca_2_qu_kstart(k),fca_2_qu_kend(k)
					do jjj = 1 + (j-1) * yratio_int, j * yratio_int
						do iii = 1 + (i-1) * xratio_int , i * xratio_int
							sum_fire_energy_rate = sum_fire_energy_rate + & 
								fire_energy_rate(iii,jjj,kkk) * fca2quic(k)%volRatio(kkk-fca_2_qu_kstart(k)+1) * firegrid%cellvol_en2atmos(kkk)
						enddo
					enddo
				enddo
 
				if (mult * sum_fire_energy_rate > 0.0) then  ! we have a fire in this QU cell and consequently a plume is born
					sum_fire_energy_rate = mult * sum_fire_energy_rate
 
					! - Winds
					xc = (real(i)-0.5)*qugrid%dx
					yc = (real(j)-0.5)*qugrid%dy
					zc = qugrid%z(k)
					ztot = qugrid%z(k) + plume_const%zvirt0 - qugrid%z(k-1)
					call ComputeWinds(i, j, k,	xc, yc, zc, uc, vc, ws)
 
					! - Updraft
					call ComputeUpdraft(sum_fire_energy_rate, ws, ztot, wc)
 
					if(wc >= plume_const%SPEEDS_RATIO * ws .and. wc >= plume_const%WC_MIN) then
						! New plume is valid
 
						! Check if the arrays can host all the fires
						if(total_fires + 1 > plume_const%MAX_NUM_PLUMES_TIMESTEP) then
							! Find the plume to replace
							
							call FindPlumeToReplace(wc, is_replace, curr_plume)
							
							if(is_replace == 1) then
								! Replace existing plume, otherwise the new plume is disregarded beacause too weak
								call SpecPropertiesNewPlume(curr_plume, xc, yc, zc, ztot, ws,  &
									uc, vc, wc, i, j, k, sum_fire_energy_rate)
							endif
						else
							! Create new plume
							total_fires = total_fires + 1
							call SpecPropertiesNewPlume(total_fires, xc, yc, zc, ztot, ws,  &
									uc, vc, wc, i, j, k, sum_fire_energy_rate)
							
							if(minw_plume > wc) then 
								minw_plume = wc
								idx_minw_plume = total_fires
							endif
						endif
											
					endif				
				endif
			enddo
		enddo
   enddo
	
   init_val = loopidx%list_count
	niter = 0
	do while(init_val	== loopidx%list_count .and. niter < 3)
		niter = niter + 1
		call random_index(8,iharvest)
		loopidx%list_count = iharvest(1)
	enddo
 
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE FindPlumeToReplace(wc, is_replace, idx)
	
	use plume_module
	
	implicit none

	real, intent(IN) :: wc			 ! m/s, new plume updraft
	
	real ::								&
		wmin								  ! m/s, min value of wc

	integer :: 							&
		is_replace,						& ! N/A, 0 = new plume is the weakest, 1 = a plume should be replaced with the new one
		ip,								& ! N/A, counter
		idx								  ! N/A, index of the plume with the lowest w
	
	if(idx_minw_plume < 0) then
		wmin = 1e8
		idx = -1
		do ip = 1,total_fires-1
			if(plume(ip)%traj_vel(3) < wmin)then
				wmin = plume(ip)%traj_vel(3)
				idx = ip
			endif
		enddo
		if(wmin < wc) then
			is_replace = 1			
			minw_plume = wc
		else
			is_replace = 0			
			minw_plume = wmin
		endif
		idx_minw_plume = ip
	else
		! We already know which plume has the minw
		if(minw_plume < wc)then
			is_replace = 1
			minw_plume = wc
		else
			is_replace = 0
		endif
		idx = idx_minw_plume
	endif
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE SpecPropertiesNewPlume(curr_plume, xc, yc, zc, ztot, ws, uc, vc, wc, &
		i, j, k,	sum_fire_energy_rate)
			
	use constants
	use grid_module
	use time_module
	use plume_const_module
	use plume_module
	use flags_module
	
	implicit none
	
	integer, intent(IN) :: i,j,k,curr_plume
	real, intent(IN) :: xc, yc, zc, ztot, ws, uc, vc, wc, sum_fire_energy_rate
	
	! Initialize plume properties
	! - Buoyancy factor
	plume(curr_plume)%heat_of_fire = sum_fire_energy_rate	! convert from the kW/m3 into kW and then into m4/s3 = FB in Briggs

	! - Location = center of the cell in the horizontal, top of the cell in the vertical
	plume(curr_plume)%traj_coord = (/xc, yc, zc/)	
					
	plume(curr_plume)%traj_coord_prev = (/xc, yc, qugrid%z(k-1)/)
	
	call ComputeTrajectoryVector(plume(curr_plume))
	call ComputeTrajectoryLength(plume(curr_plume))
					
	plume(curr_plume)%i = i
	plume(curr_plume)%j = j
	plume(curr_plume)%k = k
					
	plume(curr_plume)%i0 = i
	plume(curr_plume)%j0 = j
	plume(curr_plume)%k0 = k
					
	plume(curr_plume)%id = max_id_plume + 1
	max_id_plume = max_id_plume + 1
					
	! - Virtual source
	plume(curr_plume)%zfire_orig = qugrid%z(k-1)  ! where the fire is considered to originate. So for a ground fire, we can assign w at z=dz
	plume(curr_plume)%zvirt = plume_const%zvirt0
	plume(curr_plume)%total_z = ztot
	
	! - Winds
	plume(curr_plume)%traj_vel = (/uc, vc, wc/)	
	plume(curr_plume)%ws = ws
	
	! - Updraft
	plume(curr_plume)%wc_prev = wc
					
	! Base vectors
	call ComputeBaseVectors(plume(curr_plume))
	
	! Plume size from cell horizontal size
	plume(curr_plume)%gamma = 1.	
	plume(curr_plume)%radius_h = sqrt(qugrid%dx*qugrid%dy / PI)
	plume(curr_plume)%radius_n = plume(curr_plume)%radius_h
	call ComputePerimeter(plume(curr_plume)%radius_h, plume(curr_plume)%gamma, &
		plume(curr_plume)%perimeter, plume(curr_plume)%perimeter_over_rph)
	
	! - Maximum height reached in stable conditions
	if(plume_const%BRUNT_VAISALA_FREQ2 /= 0) then
		plume(curr_plume)%zmax = (6. * plume(curr_plume)%heat_of_fire / &
			(plume_const%BETA_ENTR2 * plume_const%BRUNT_VAISALA_FREQ2 * plume(curr_plume)%ws))**(1./3.)
	endif
					
	plume(curr_plume)%start_time	= real(fcatime%current_int - fcatime%dt_int)
	plume(curr_plume)%total_time = 0.
	
	! - Loop index to compute interactions among plumes
	plume(curr_plume)%istart = plume(curr_plume)%i - ceiling(plume(curr_plume)%radius_h * qugrid%dxi)
	plume(curr_plume)%istart = max(plume(curr_plume)%istart, 1)
						
	plume(curr_plume)%iend = plume(curr_plume)%i + ceiling(plume(curr_plume)%radius_h * qugrid%dxi)
	plume(curr_plume)%iend = min(plume(curr_plume)%iend, qugrid%nx-1)

	plume(curr_plume)%jstart = plume(curr_plume)%j - ceiling(plume(curr_plume)%radius_h * qugrid%dyi)
	plume(curr_plume)%jstart = max(plume(curr_plume)%jstart, 1)
						
	plume(curr_plume)%jend = plume(curr_plume)%j + ceiling(plume(curr_plume)%radius_h * qugrid%dyi)
	plume(curr_plume)%jend = min(plume(curr_plume)%jend, qugrid%ny-1)

	plume(curr_plume)%kstart = plume(curr_plume)%k
	plume(curr_plume)%kend = min(plume(curr_plume)%k + 1, qugrid%nz-1)
	
   if(flag%plume_traj_file > 0) then
		call WriteOutPlume(plume(curr_plume), 0)
	endif	
	
	if(plume_const%time_step_flag == 1) then
		call PlumeTimeStepCalc(plume(curr_plume))		
	endif
	
	end
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE PlumeTimeStepCalc(p)
	
	use grid_module
	use time_module
	use plume_const_module
	use plume_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: p		
	real :: umin, vmin, wmin
		
	umin = qugrid%dx / (abs(p%traj_vel(1)) + 1e-6)
	vmin = qugrid%dy / (abs(p%traj_vel(2)) + 1e-6)
	wmin = qugrid%dz_array(p%k) / (abs(p%traj_vel(3)) + 1e-6)
	plume_const%time_step = min(plume_const%time_step, 0.7 * min(umin, vmin, wmin))
	plume_const%time_step = max(plume_const%time_step, plume_const%min_time_step)
	
	if(plume_const%total_time + plume_const%time_step > fcatime%dt_real) then
		plume_const%time_step = fcatime%dt_real - plume_const%total_time
	endif
		
	
	END
	!=====================================================================================================================
	!=====================================================================================================================						
	SUBROUTINE ComputeTrajectoryVector(currplume)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType) :: currplume
	
	currplume%s = currplume%traj_coord - currplume%traj_coord_prev
	
	end	
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE ComputeTrajectoryLength(currplume)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType) :: currplume
	
	currplume%traj_len = sqrt(sum(currplume%s**2))
	
	end	
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE ComputeBaseVectors(currplume)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType) :: currplume
	real:: den, phi
		
	if(currplume%s(1) == 0 .and. currplume%s(2) == 0) then
		! Plume is straight up
		currplume%eh(1) = currplume%traj_vel(2) / currplume%ws
		currplume%eh(2) = -currplume%traj_vel(1) / currplume%ws
		
		den = 1./ (currplume%traj_len * currplume%ws)
		currplume%en(1) = currplume%traj_vel(1) * currplume%s(3) * den
		currplume%en(2) = currplume%traj_vel(2) * currplume%s(3) * den
		currplume%en(3) = 0.
		
	else
		! Angle between sz and s
		phi = acos(max(min(currplume%s(3) / &
		                   currplume%traj_len, 0.999999), -0.999999))
		den = 1. / (sin(phi) * currplume%traj_len)
	
		currplume%eh(1) = currplume%s(2) * den
		currplume%eh(2) = -currplume%s(1) * den
		
		den = 1./ (currplume%traj_len**2 * sin(phi))
		currplume%en(1) = currplume%s(1) * currplume%s(3)
		currplume%en(2) = currplume%s(2) * currplume%s(3)
		currplume%en(3) = -(currplume%s(1)**2 + currplume%s(2)**2)
		currplume%en = currplume%en * den
	endif
	
	currplume%eh(3) = 0.
	
	end	
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE CalcZVirt(gamma, plume_radius_h, plume_perimeter_over_rph, z, zorig, zVirt)
	
	use plume_const_module
	use constants
	
	implicit none

	real,intent(IN) :: plume_radius_h, z, zorig, gamma, plume_perimeter_over_rph
	real :: zVirt
	
	! this is the quantity that added to the plume z gives you the z to use in the Briggs eqn	
	zVirt = 2 * PI * gamma * plume_radius_h / (plume_perimeter_over_rph * plume_const%BETA_ENTR) + zorig - z

	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE DeletePlumes()

	use plume_module
	use file_handling_module
	use time_module
	use flags_module
	
	implicit none
	
	integer :: i,j
	
	i = 0
	j = total_fires
	do while(i < j .and. i < total_fires)
		i = i + 1
		if(to_delete(i) == 1) then
			if(flag%plume_traj_file > 0) call WriteOutPlume(plume(i), 0)
			
			do while(j > i .and. to_delete(j) == 1)
				! Track maximum height reached by the plumes
				plume_diagnostics%max_height = max(plume_diagnostics%max_height, plume(j)%traj_coord(3))
				if(flag%plume_traj_file > 0) call WriteOutPlume(plume(j), 0)
				j = j - 1
				total_fires = total_fires - 1				
			enddo
			if(to_delete(j) == 0) then
				! Swap plume to delete with valid plumes
				plume(i) = plume(j)
				to_delete(i) = to_delete(j)
				j = j - 1
				total_fires = total_fires - 1
				if(idx_minw_plume == j) then
					idx_minw_plume = i
				endif
			else
				total_fires	= i-1				
				if(flag%plume_traj_file > 0) call WriteOutPlume(plume(j), 0)
			endif			
		endif
	enddo
	to_delete = 0
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE ComputeWinds(i, j , k,xplume,yplume,zplume,uc,vc,ws)

	use grid_module
	
	use winds_module
	
	implicit none
		
	integer,intent(IN) :: i,j,k  ! of the plume
	real,intent(IN) :: xplume,yplume,zplume	
	real,intent(OUT) :: uc, vc, ws
	real :: u1, u2
	
	! Interpolate u
	call linear_interpolation(qugrid%zm(k),qugrid%zm(k+1),	&
		quwinds%u(i,j,k), quwinds%u(i,j,k+1), zplume, u1)  ! vertical interpolation at the cell top

	call linear_interpolation(qugrid%zm(k),qugrid%zm(k+1),	&
		quwinds%u(i+1,j,k), quwinds%u(i+1,j,k+1), zplume, u2)  ! vertical interpolation at the cell top
	uc = u1 + (u2-u1) * (xplume - real(i-1)*qugrid%dx) * qugrid%dxi
	
	! Interpolate v
	call linear_interpolation(qugrid%zm(k),qugrid%zm(k+1),	&
		quwinds%v(i,j,k), quwinds%v(i,j,k+1), zplume, u1)  ! vertical interpolation at the cell top
	call linear_interpolation(qugrid%zm(k),qugrid%zm(k+1),	&	
		quwinds%v(i,j+1,k), quwinds%v(i,j+1,k+1), zplume, u2)  ! vertical interpolation at the cell top
	vc = u1 + (u2-u1) * (yplume - real(j-1)*qugrid%dy) * qugrid%dyi

	ws = max(sqrt(uc**2 + vc**2) , 0.1)
	

	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE ComputeUpdraft(heat_of_fire,ws,ztot,updraft)

	use plume_const_module
	use fireca_constants_module
	
	implicit none
	
	real :: sqrt_fac
	real,intent(IN) :: heat_of_fire,ws,ztot
	real,intent(out) :: updraft

	! Updraft for ws << wc
	!updraft = 5./(6. * plume_const%ALPHA_ENTR) * (0.9 * plume_const%ALPHA_ENTR * heat_of_fire / ztot)**ONE_THIRD 
	!if(ws > updraft) then
	!	sqrt_fac = 2. - plume_const%BRUNT_VAISALA_FREQ2 * ws *&
	!		plume_const%BETA_ENTR2 * ztot**3 / 3. / heat_of_fire
	!	if(plume_const%BRUNT_VAISALA_FREQ2 <= 0 .or. sqrt_fac < 0.) then
	!		updraft = sqrt(2. * heat_of_fire / (3. * plume_const%BETA_ENTR2 * ws * ztot))
	!	else
	!		updraft = sqrt(heat_of_fire / (3. * plume_const%BETA_ENTR2 * ws * ztot) * sqrt_fac)
	!	endif
	!endif
	
	sqrt_fac = 2. - plume_const%BRUNT_VAISALA_FREQ2 * ws *&
		plume_const%BETA_ENTR2 * ztot**3 / 3. / heat_of_fire
	if(plume_const%BRUNT_VAISALA_FREQ2 == 0 .or. sqrt_fac < 0.) then
		updraft = sqrt(2. * heat_of_fire / (3. * plume_const%BETA_ENTR2 * ws * ztot))
	else
		updraft = sqrt(heat_of_fire / (3. * plume_const%BETA_ENTR2 * ws * ztot) * sqrt_fac)
	endif

	end
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE WWIndexCalc()
 
	use plume_module
	use grid_module

	implicit none
		
	real :: arg, theta, phi, dx, dy, dz
	integer :: arg_int, kkk, ip
 	
	!$OMP parallel do private(ip, theta, phi, dx, dy, dz, arg, arg_int, kkk)
   do ip = 1, total_fires
		
		theta = atan2(plume(ip)%s(3), sqrt(plume(ip)%s(1)**2 + plume(ip)%s(2)**2))
		phi = atan2(plume(ip)%s(2), plume(ip)%s(1))
		
		! - x
		dx = abs(plume(ip)%radius_h * sin(phi)) + abs(plume(ip)%radius_n * sin(theta) * cos(phi))
		
		arg = min(plume(ip)%traj_coord_prev(1), plume(ip)%traj_coord(1)) - dx
		arg_int = ceiling( arg * qugrid%dxi )
		plume(ip)%istart = max(arg_int, 1)
				 
		arg = max(plume(ip)%traj_coord_prev(1), plume(ip)%traj_coord(1)) + dx
		arg_int = ceiling( arg * qugrid%dxi )
		plume(ip)%iend = min( arg_int , qugrid%nx-1)
   
		! - y
		dy = abs(plume(ip)%radius_h * cos(phi)) + abs(plume(ip)%radius_n * sin(theta) * sin(phi))
		
		arg = min(plume(ip)%traj_coord_prev(2), plume(ip)%traj_coord(2)) - dy
		arg_int = ceiling( arg * qugrid%dyi )
		plume(ip)%jstart = max(arg_int, 1)
				 
		arg = max(plume(ip)%traj_coord_prev(2), plume(ip)%traj_coord(2)) + dy
		arg_int = ceiling( arg * qugrid%dyi )
		plume(ip)%jend = min( arg_int , qugrid%ny-1)
			
		! - z
		dz = abs(plume(ip)%radius_n * cos(theta))
		
		arg = min(plume(ip)%traj_coord_prev(3), plume(ip)%traj_coord(3)) - dz		
		kkk = 3
		do while(qugrid%zm(kkk) < arg .and. kkk < qugrid%nz-2)
			kkk = kkk + 1
		enddo
		if(qugrid%zm(kkk) >  arg) kkk = max(kkk - 1,3)
		plume(ip)%kstart = kkk
		
		arg = max(plume(ip)%traj_coord_prev(3), plume(ip)%traj_coord(3)) + dz
		do while(qugrid%zm(kkk) < arg .and. kkk < qugrid%nz-2)
			kkk = kkk + 1
		enddo
		plume(ip)%kend = min(max(kkk, plume(ip)%kstart+1), qugrid%nz-1)
	enddo
	!$OMP end parallel do
 
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE WWCalc()

	! This routine returns quwinds%wplume**WC_EXPONENT
	
	! wc is the top hat value
	
	!https://www.youtube.com/watch?v=9wznbg_aKOo
		
	use time_module
	use constants
	use rnd_num_vars
	use omp_lib
	use constants
	use plume_module
	use plume_const_module
	use grid_module
	use winds_module
	
	implicit none
	
	integer :: iii, jjj, kkk, ip
	real, dimension(3) :: dtest  ! m, segment (xprev, yprev, zprev)-(point in the QU grid)
	real ::							&
		dotres,						&
		qu_point_dist,				& ! N/A, distance of point in the QU grid from plume centerline, normalized
		dtest_proj,					& ! m, projection of the segment (xprev, yprev, zprev)-(point in the QU grid) on the trajectory
		dh,							& ! m, segment (xprev, yprev, zprev)-(point in the QU grid) DOT eh
		dn,							& ! m, segment (xprev, yprev, zprev)-(point in the QU grid) DOT eh								
		w								  ! m/s, linear interpolation of (wprev, w) on the trajectory based on dtest
		
	!!$OMP parallel do private(ip, iii, jjj, kkk, dtest, dtest_proj, dh, dn, w, qu_point_dist, dotres)
   do ip = 1,total_fires
		do kkk = plume(ip)%kstart,plume(ip)%kend
			do jjj = plume(ip)%jstart,plume(ip)%jend
				do iii = plume(ip)%istart,plume(ip)%iend
					dtest = (/																	&
						qugrid%xcenters(iii) - plume(ip)%traj_coord_prev(1),		&
						qugrid%ycenters(jjj) - plume(ip)%traj_coord_prev(2),		&
						(qugrid%z(kkk-1)		- plume(ip)%traj_coord_prev(3))/)
										
					call dot(dtest, plume(ip)%s, dotres)
					dtest_proj = dotres / plume(ip)%traj_len
					
					if(dtest_proj >= 0 .and. dtest_proj <= plume(ip)%traj_len) then
					
						call dot(dtest, plume(ip)%eh, dh)
						call dot(dtest, plume(ip)%en, dn)
												
						!dist_AC = sqrt(sum(dtest)**2)
						!dist_BC = sqrt(												&
						!	(qugrid%xcenters(iii) - plume(ip)%x)**2 +			&
						!	(qugrid%ycenters(jjj) - plume(ip)%y)**2 +			&
						!	(qugrid%z(kkk-1)		 - plume(ip)%z)**2)
				
						qu_point_dist = (dn / plume(ip)%radius_n)**2 + (dh / plume(ip)%radius_h)**2
												
						!if(dist_AC <= 1e-2) then
						!	! To avoid divisions by zero
						!	call omp_set_lock(plume_lock(iii, jjj, kkk))
						!	call ComputeW(plume(ip)%wc_prev, quwinds%wplume(iii,jjj,kkk), qu_point_dist)
						!	call omp_unset_lock(plume_lock(iii, jjj, kkk))
						!elseif(dist_BC <= 1e-2) then
						!	! To avoid divisions by zero
						!	call omp_set_lock(plume_lock(iii, jjj, kkk))
						!	call ComputeW(plume(ip)%wc, quwinds%wplume(iii,jjj,kkk), qu_point_dist)
						!	call omp_unset_lock(plume_lock(iii, jjj, kkk))
						!else
						
						if(qu_point_dist <= 1.) then
							! Check if point is within the perimeter of the ellipse
							w = plume(ip)%wc_prev + (plume(ip)%traj_vel(3) - plume(ip)%wc_prev) * dtest_proj / plume(ip)%traj_len
							!call omp_set_lock(plume_lock(iii, jjj, kkk))
							call ComputeW(w, quwinds%wplume(iii,jjj,kkk), qu_point_dist)
							!call omp_unset_lock(plume_lock(iii, jjj, kkk))
						endif
							
					endif
				enddo
			enddo
		enddo 
	enddo
	!!$OMP end parallel do 
	
	!if(plume_const%do_w_fill == 1) call FillWW()
	
	end
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE ComputeW(w, wplume, qu_point_dist)
	
	use plume_const_module
	
	implicit none
	
	real,parameter :: TOP_HAT_TO_GAUSS = 1.6
	real, intent(IN) :: w, qu_point_dist
	real :: wplume
	
	wplume = wplume +	(TOP_HAT_TO_GAUSS * w)**plume_const%WC_EXPONENT *		&
		exp(- plume_const%WC_EXPONENT * qu_point_dist)
	
	end
	!=====================================================================================================================
	!=====================================================================================================================	
	SUBROUTINE ComputePlumeRadius(currplume) 
	
	use plume_const_module
	use plume_module
	use constants
	
	implicit none
	
	TYPE(PlumeType) :: currplume	
		
	currplume%radius_h = plume_const%BETA_ENTR * currplume%total_z / (2. * PI * currplume%gamma) * &
		currplume%perimeter_over_rph
	
	currplume%radius_n = currplume%radius_h * currplume%gamma
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE PlumeFinalization()
	
	use plume_module
	use plume_const_module
	
	implicit none
		
   integer :: i

	to_delete = 0
		   
	if(plume_const%time_step_flag == 1) then
		plume_const%time_step = 1e8
	endif
	!$OMP parallel do private(i)
	do i = 1,total_fires ! after the initialization, it will do only the new fires		
		! Plume time
		plume(i)%total_time = plume(i)%total_time + plume_const%time_step
				
		if(plume_const%time_step_flag == 1) then
		plume_const%time_step = 1e8
		endif
		
		if(plume_const%BRUNT_VAISALA_FREQ2 > 0 .and. plume(i)%total_time > plume_const%max_time)then
			to_delete(i) = 1		
		endif		
	enddo
	!$OMP end parallel do 
	
	if(plume_const%time_step_flag == 1) then
		do i = 1,total_fires
			if(to_delete(i) == 0)call PlumeTimeStepCalc(plume(i))			
		enddo
	endif
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE MovePlumes()

	use plume_module
	use plume_const_module
	use time_module
	use grid_module
	
	implicit none
	
	integer :: i
		
	!$OMP parallel do private(i)
	do i = 1,total_fires
		! Store previous position of the trajectory
		plume(i)%traj_coord_prev = plume(i)%traj_coord		
	enddo
	!$OMP end parallel do
	
	! Pull plumes toward each other
	call PullPlumesTogether()
	
	! Move plumes with the wind
	!$OMP parallel do private(i)
	do i = 1,total_fires
					
		plume(i)%traj_coord = plume(i)%traj_coord + plume(i)%traj_vel * plume_const%time_step
		plume(i)%total_z = plume(i)%total_z + plume(i)%traj_vel(3) * plume_const%time_step
				
		! Cell index
		call GetPlumeCellIndex(plume(i))
			
		! Check if plume is in the grid
		if(plume(i)%traj_coord(1) < 0. .or. plume(i)%traj_coord(1) > qugrid%Lx .or. &
			plume(i)%traj_coord(2) < 0. .or. plume(i)%traj_coord(2) > qugrid%Ly .or. &
			plume(i)%traj_coord(3) < 0. .or. plume(i)%traj_coord(3) > qugrid%Lz) then
			to_delete(i) = 1
		endif
		
	enddo
	!$OMP end parallel do 
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE GetPlumeCellIndex(p)
	
	use plume_module
	use grid_module
	
	TYPE(PlumeType) :: p
	
	p%i = ceiling(p%traj_coord(1) * qugrid%dxi)
	p%j = ceiling(p%traj_coord(2) * qugrid%dyi)
	call GetPlumeK(p%traj_coord(3), p%k)		
	
	END
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE GetPlumeK(plume_z, k)
	
	use grid_module
	
	implicit none
	
	integer :: k
	real,intent(IN) :: plume_z
	
	k = 1
	do while(qugrid%z(k) < plume_z .and. k < qugrid%nz - 1)
		k = k + 1
	enddo		
	
	END
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE UpdatePlumeProperties()
	
	use plume_const_module
	use plume_module
	use file_handling_module
	
	implicit none
	
	integer :: i
	
	!$OMP parallel do private(i)
	do i = 1,total_fires
		plume(i)%wc_prev = plume(i)%traj_vel(3)
		! Winds - we have moved the trajectory from the point kk to the point kk+1			
		call ComputeWinds(plume(i)%i,plume(i)%j,plume(i)%k,								&
			plume(i)%traj_coord(1), plume(i)%traj_coord(2), plume(i)%traj_coord(3),	&
			plume(i)%traj_vel(1), plume(i)%traj_vel(2), plume(i)%ws)
			
		! Updraft	
		call ComputeUpdraft(plume(i)%heat_of_fire,plume(i)%ws,plume(i)%total_z,plume(i)%traj_vel(3))
			
		if(plume(i)%traj_vel(3) < plume_const%SPEEDS_RATIO * plume(i)%ws .or. plume(i)%traj_vel(3) < plume_const%WC_MIN) then
			to_delete(i) = 1
		else
			! Plume trajectory properties
			call ComputeTrajectoryVector(plume(i))
			call ComputeTrajectoryLength(plume(i))
			
			! Plume base vectors
			call ComputeBaseVectors(plume(i))
			! Plume radius
			call ComputePlumeRadius(plume(i))
			
			! Plume perimeter
			call ComputePerimeter(plume(i)%radius_h, plume(i)%gamma, &
				plume(i)%perimeter, plume(i)%perimeter_over_rph)
		endif	
	enddo
	!$OMP end parallel do 
	
	END
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE PullPlumesTogether()
	
	use constants
	use plume_module
	use plume_const_module
	use rnd_num_vars
	use grid_module
	
	implicit none
		
	real,dimension(total_fires, 3, numthreads) ::		&
		mov					 ! m, total movement of a point
	real,dimension(total_fires) ::		&		
		dist_moved
	
	real, dimension(3) ::					&
		delta,									& ! m
		mincoord,								& ! m, min location of the plumes
		maxcoord									  ! m, max location of the plumes
	real ::										&
		const,									& 
		plume_dist,								& ! m, distance between two plumes centerlines at time t
		plume_dist2,							& ! m^2, distance between two plumes centerlines at time t, squared
		dist_1,									& ! m, total distance travelled by plume ip1 because of the updraft at point ip2
		dist_2,									& ! m, total distance travelled by plume ip2 because of the updraft at point ip1
		sum_dist 								  ! m, sum of the distance travelled by plume ip1 and plume ip2
	
	integer ::									&	
		TID,										&
		i,											&
		ip1,ip2							  	     ! N/A, plumes 
	
	! Initilization
   mov = 0.	
	
	const = 0.5 * plume_const%BETA_ENTR * plume_const%time_step / PI
	do ip1 = 1,total_fires
		dist_moved(ip1) =  plume(ip1)%traj_vel(3) * plume(ip1)%traj_len * plume(ip1)%perimeter
	enddo
	dist_moved = dist_moved * const

	mincoord = (/qugrid%Lx, qugrid%Ly, qugrid%Lz/)
	maxcoord = (/0., 0., 0./)
	
	do ip1 = 1,total_fires		
		do i = 1,3
			mincoord(i) = min(plume(ip1)%traj_coord(i), mincoord(i))
			maxcoord(i) = max(plume(ip1)%traj_coord(i), maxcoord(i))
		enddo		
	enddo
		
	! Calc pull
	!$OMP parallel do private(ip1, ip2, plume_dist, plume_dist2, dist_1, dist_2, sum_dist, &
	!$OMP delta, TID, i)
	do ip1 = 1,total_fires-1
		TID = 1
		!$ TID = TID + OMP_GET_THREAD_NUM()
		do ip2 = ip1+1,total_fires
			
			! Only if the two plumes do not originate from the same cell
			if(plume(ip1)%i0 /= plume(ip2)%i0 .or. plume(ip1)%j0 /= plume(ip2)%j0 .or. &
				plume(ip1)%k0 /= plume(ip2)%k0) then
			
				! 1e-2 added so they are always at some distance
				plume_dist2 = sum((plume(ip1)%traj_coord - plume(ip2)%traj_coord)**2) + 1e-2
				plume_dist = sqrt(plume_dist2)
					
				! Distance traveled by ip1 because of ip2
				dist_1 = dist_moved(ip2) / plume_dist2
				! Distance traveled by ip2 because of ip1
				dist_2 = dist_moved(ip1) / plume_dist2
                     					
				sum_dist = dist_1 + dist_2
				! Plumes cannot go past each other
				if(sum_dist > plume_dist) then
					dist_1 = dist_1 * plume_dist / sum_dist
					dist_2 = dist_2 * plume_dist / sum_dist
				endif
		
				! Project in x, y, z directions
				delta = (plume(ip2)%traj_coord - plume(ip1)%traj_coord) / plume_dist
				
				do i = 1,3
					! Move plume 1
					mov(ip1,i,TID) = mov(ip1,i,TID) + dist_1 * delta(i)
					
					! Move plume 2
					mov(ip2,i,TID) = mov(ip2,i,TID) - dist_2 * delta(i)
				enddo
			endif
		enddo		
	enddo
	!$OMP end parallel do 
				
	!$OMP parallel do private(ip1, i)
	do ip1 = 1,total_fires
		
		! Move plumes
		do i = 1,3
			plume(ip1)%traj_coord(i) = plume(ip1)%traj_coord(i) + sum(mov(ip1,i,:))		
		
			! Keep plumes within inital box
			plume(ip1)%traj_coord(i) = min(max(plume(ip1)%traj_coord(i), mincoord(i)), maxcoord(i))
		enddo
		
		! Find cells for plumes
		plume(ip1)%i = ceiling(plume(ip1)%traj_coord(1) * qugrid%dxi)
		plume(ip1)%j = ceiling(plume(ip1)%traj_coord(2) * qugrid%dyi)
		
		call GetPlumeK(plume(ip1)%traj_coord(3), plume(ip1)%k)
	enddo	
	!$OMP end parallel do 
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE dot(v1, v2, outvar)
	
	implicit none
	
	real,dimension(3),intent(IN) :: v1, v2 ! input vectors
	real,intent(OUT) :: outvar  ! dot product result
	integer :: i
	
	outvar = 0.
	do i = 1,3
		outvar = outvar + v1(i) * v2(i)
	enddo
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE CheckLen(v, check)
	
	implicit none

	real, dimension(3) :: v
	real, parameter :: TOL = 1e-4			  ! m^2, minimum distance to say that two points are coincident squared (1 cm)	
	real :: vlen
	integer, intent(OUT) :: check
	
	vlen = sum(v**2)
	if(vlen <= TOL) then
		check = 1
	else
		check = 0
	endif
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE CheckPointCoincidence(p1, p2, p3, p4, check)
	
	implicit none
	
	real, dimension(3), intent(IN) ::	&
		p1, p2,									& ! m, coordinates of beginning and end of a plume trajectory
		p3, p4									  ! m, coordinates of beginning and end of anoter plume trajectory
	integer, intent(OUT) :: check	
	real, dimension(3) :: v

	check = 0
	
	v = p1 - p3
	call CheckLen(v, check)
	if(check == 1) return
		
	v = p1 - p4
	call CheckLen(v, check)
	if(check == 1) return
	
	v = p2 - p3
	call CheckLen(v, check)
	if(check == 1) return
		
	v = p2 - p4
	call CheckLen(v, check)
	if(check == 1) return
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE DistBetween2Segment(p1, p2, p3, p4, ip1, ip2, distance, rpw, rqw)

	! http://geomalgorithms.com/a07-_distance.html
	
	use plume_const_module
	use plume_module
	
	implicit none
	
	integer, intent(IN) ::					&
		ip1, ip2									  ! N/A, indexes of the plumes that could merge
	real, dimension(3), intent(IN) ::	&
		p1, p2,									& ! m, coordinates of beginning and end of a plume trajectory
		p3, p4									  ! m, coordinates of beginning and end of anoter plume trajectory
	real, intent(out) ::						&
		rpw, rqw,								& ! m, projection of the ellipse radius on wc
		distance									  ! m, minimum distance between two line segments	
	real, dimension(3) ::					&
		wc,										& ! m, vector connecting closest points
		u, v, w,									& ! m, vectors built from trajectory points
		dp,										& ! m, vector difference of the two closest points
		pc_contact, qc_contact				  ! m, points at the extremes of a vector connecting the input segments		
		!cross_val								  ! N/A, cross product result
	real ::										&		
		dot_wcehp, dot_wcenp,				& ! m, dot product between eh and en of Pc and wc
		dot_wcehq, dot_wcenq,				& ! m, dot product between eh and en of Qc and wc
		a,b,c,d,e,								& ! m^2, dot product result
		den,										& ! m^4, denominator in the expressions for the distance segment		
		sN,sD,tN,tD,sc,tc						  ! N/A, temporary parameters to determine the distance segment
	real,parameter ::							&
		SMALL_NUM = 1e-6,						& ! m^4, tolerance on denominator == 0 (parallel segments)					
		TOL = 1e-4								  ! N/A, tolerance on sc, tc being 0 or 1
	integer :: check
	
	! Check if points coincide
	!call CheckPointCoincidence(p1, p2, p3, p4, check)
	!if(check == 1) then
	!	distance = 0
	!	return
	!endif
	
	! Compute distance
	u = p1 - p2
	v = p3 - p4
	w = p2 - p4
    
	call dot(u,u,a)
	call dot(u,v,b)
	call dot(v,v,c)
	call dot(u,w,d)
	call dot(v,w,e)
	den = a*c - b**2
	sD = den
	tD = den
      
	check = 0
	! compute the line parameters of the two closest points
	if (den < SMALL_NUM)then  ! the lines are almost parallel
		sN = 0.0			! force using point P0 on segment S1
      sD = 1.0			! to prevent possible division by 0.0 later
      tN = e
      tD = c
	else                ! get the closest points on the infinite lines
		sN = (b*e - c*d)
		tN = (a*e - b*d)
		if (sN < 0.0) then  ! sc < 0 => the s=0 edge is visible       
			sN = 0.0
			tN = e
			tD = c
		elseif (sN > sD) then! sc > 1 => the s=1 edge is visible
			sN = sD
			tN = e + b
			tD = c
		endif
	endif
    
	if (tN < 0.0)then            ! tc < 0 => the t=0 edge is visible
		tN = 0.0
		! recompute sc for this edge
		if (-d < 0.0)then
			sN = 0.0
		elseif (-d > a)then
			sN = sD
		else
			sN = -d
			sD = a
		endif
	elseif (tN > tD) then      ! tc > 1 => the t=1 edge is visible
		tN = tD
		! recompute sc for this edge
		if ((-d + b) < 0.0)then
			sN = 0
		elseif ((-d + b) > a)then
			sN = sD
		else 
			sN = (-d + b)
			sD = a
		endif
	endif
    
	! finally do the division to get sc and tc
	if(abs(sN) < SMALL_NUM)then
		sc = 0.0
	else
		sc = sN / sD
	endif
    
	if(abs(tN) < SMALL_NUM)then
		tc = 0.0
	else
		tc = tN / tD
	endif
    
	! get the difference of the two closest points
	dP = w + (sc * u) - (tc * v)  ! = S1(sc) - S2(tc)

	distance = sqrt(dP(1)**2 + dp(2)**2 + dp(3)**2)
	
	pc_contact = p2 + sc*u
	qc_contact = p4 + tc*v
	wc = pc_contact - qc_contact
	
	call dot(wc, plume(ip1)%eh, dot_wcehp)
	call dot(wc, plume(ip1)%en, dot_wcenp)
	call dot(wc, plume(ip2)%eh, dot_wcehq)
	call dot(wc, plume(ip2)%en, dot_wcenq)
	rpw = sqrt( (dot_wcehp * plume(ip1)%radius_h)**2 + (dot_wcenp * plume(ip1)%radius_n)**2) / distance
	rqw = sqrt( (dot_wcehq * plume(ip2)%radius_h)**2 + (dot_wcenq * plume(ip2)%radius_n)**2) / distance	

	end
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE CheckOverlap(p1, p2, p3, p4, is_valid)

	use plume_const_module
	use constants
	
	implicit none 
	
	real, dimension(3), intent(IN) ::	&
		p1, p2,									& ! m, coordinates of beginning and end of a plume trajectory
		p3, p4									  ! m, coordinates of beginning and end of anoter plume trajectory
	integer, intent(OUT) :: is_valid
	real, dimension(3) :: u, v	
	real ::										&
		dotres,									& ! m, result of the dot product
		u_len,									& ! m, length of vector
		max_angle								  ! rad, max angle between segments
	real, dimension(4) ::	angles		  ! rad, angles between segments
	integer, dimension(1) :: max_idx

	is_valid	= 1
	angles = 0.
	
	! Check angles to keep only plumes overlapping a lot
	! case #1
	u = p2 - p1
	u_len = sum(u**2)	
	v = p4 - p1
	call calc_angle(u, v, u_len, angles(1))
	
	! case #2
	v = p3 - p1
	call calc_angle(u, v, u_len, angles(2))	
	
	! case #3
	u = -u
	v = p3 - p2
	call calc_angle(u, v, u_len, angles(3))
	
	! case #4
	v = p4 - p2
	call calc_angle(u, v, u_len, angles(4))
    
	max_angle = maxval(angles)
	max_idx = maxloc(angles)
	
	if(max_idx(1) == 1) then
		u = p2 - p1
		v = p4 - p1
	elseif(max_idx(1) == 2) then
		u = p2 - p1
		v = p3 - p1
	elseif(max_idx(1) == 3) then
		u = p1 - p2
		v = p3 - p2
	else
		u = p1 - p2
		v = p4 - p2
	endif
	call dot(u, v, dotres)
	u_len = sqrt(u_len)
	if(max_angle > PI * 0.5 .and. abs(dotres) / u_len > plume_const%max_fraction_traj * u_len)then
		is_valid = 0
	endif

	end
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE calc_angle(u, v, u_len, angles)
	
	use constants
	
	implicit none
	
	real, dimension(3), intent(IN) :: u, v
	real, intent(IN) :: u_len
	real, intent(OUT) :: angles	
	real :: outdot, v_len
	real, parameter :: TOL = 1.2e-7
	
	v_len = sum(v**2)
   call dot(u, v, outdot)
	! To avoid floating point errors
	outdot = outdot / sqrt(u_len*v_len)
	outdot = max(outdot, -1.+TOL)
	outdot = min(1.-TOL, outdot)	
   angles = acos( max(min(outdot, 0.999999), -0.999999))
	
	end
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE CheckParallelSegments(p1, p2, p3, p4, check)
	
	implicit none
	
	real,dimension(3), intent(IN) :: p1, p2, p3, p4  ! p1,p2 = segment 1; p3,p4 = segment 2
	integer,intent(OUT) :: check
	real :: alpha, sina, cosa, beta, sinb, cosb, gamma, sing, cosg, lh, p11, p22, p33, p44
	real,dimension(3) :: rot
	real, dimension(3) :: direction
	
	integer :: i
	
	
	! Rotate to align x-axis with the segment p1-p2, origin in p1
	direction = p2 - p1
	! Yaw
	alpha = atan2(direction(2), direction(1))
	sina = sin(alpha)
	cosa = cos(alpha)
	! pitch
	Lh = sqrt( direction(1)**2 + direction(2)**2)
	beta = atan2( direction(3), Lh)	
	sinb = sin(beta)
	cosb = cos(beta)
	! roll
	gamma = 0
	sing = sin(gamma)
	cosg = cos(gamma)

	rot(1) = cosa*cosb
	rot(2) = cosa*sinb*sing - sina*cosg
	rot(3) = cosa*sinb*cosg - sina*sing
	!rot(1,1) = cosa*cosb
	!rot(1,2) = cosa*sinb*sing - sina*cosg
	!rot(1,3) = cosa*sinb*cosg - sina*sing
	!rot(2,1) = sina*cosb
	!rot(2,2) = sina*sinb*sing + cosa*cosg
	!rot(2,3) = sina*sinb*cosg - cosa*sing
	!rot(3,1) = -sinb
	!rot(3,2) = cosb*sing
	!rot(3,3) = cosb*cosg

	! Rotate points into new coordinate system
	p11 = 0	
	p22 = 0
	p33 = 0
	p44 = 0
	do i = 1,3
		p22 = p22 + rot(i)*(p2(i)-p1(i))
		p33 = p33 + rot(i)*(p3(i)-p1(i))
		p44 = p44 + rot(i)*(p4(i)-p1(i))
	enddo
	!Check if there could be overlap
	! c1: p1 <= p3 <= p2
	! c2: p1 <= p4 <= p2
	! c3: p3 <= p1 & p4 >= p2
	check	= 0
	if(p33 > p11 .and. p33 < p22) then
		check	= 1
		return
	elseif(p44 > p11 .and. p44 < p22) then
		check	= 1
		return
	elseif(p33 <= p11 .and. p44 >= p22)then
		check	= 1
		return
	endif
	
	end
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE WriteOutPlume(plume_in, it)
	
	use plume_module
	use file_handling_module
	use time_module
	use constants
	use flags_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: plume_in
	integer, intent(IN) :: it  ! 0 = print current time step, 1 = print prev variables
	
	
	if(it == 0) then
		if(flag%plume_traj_file == 1) then
			write(ID_FILE_TRAJ,'(2(i8,a,2x),15(g22.7,a,2x))')	&
				plume_in%id,',',											&
				fcatime%current_int - fcatime%dt_int,',',			&
				plume_in%total_time,',',								&
				plume_in%traj_coord(1),',',							&
				plume_in%traj_coord(2),',',							&
				plume_in%traj_coord(3),',',							&
				plume_in%traj_vel(1),',',								&
				plume_in%traj_vel(2),',',								&
				plume_in%traj_vel(3),',',								&
				plume_in%radius_h,',',									&
				plume_in%eh(1),',',										&
				plume_in%eh(2),',',										&
				plume_in%eh(3),',',										&
				plume_in%radius_n,',',									&
				plume_in%en(1),',',										&
				plume_in%en(2),',',										&
				plume_in%en(3),','
		else
		write(ID_FILE_TRAJ)										&
			plume_in%id,											&
			fcatime%current_int - fcatime%dt_int,			&
			plume_in%total_time,									&
			plume_in%traj_coord(1),								&
			plume_in%traj_coord(2),								&
			plume_in%traj_coord(3),								&
			plume_in%traj_vel(1),								&
			plume_in%traj_vel(2),								&
			plume_in%traj_vel(3),								&
			plume_in%radius_h,									&
			plume_in%eh(1),										&
			plume_in%eh(2),										&
			plume_in%eh(3),										&
			plume_in%radius_n,									&
			plume_in%en(1),										&
			plume_in%en(2),										&
			plume_in%en(3)
		endif
	else
		if(flag%plume_traj_file == 1) then
			write(ID_FILE_TRAJ,'(2(i8,a,2x),15(g22.7,a,2x))')	&
				plume_in%id,',',											&
				fcatime%current_int - fcatime%dt_int,',',			&
				plume_in%total_time,',',								&
				plume_in%traj_coord_prev(1),',',						&
				plume_in%traj_coord_prev(2),',',						&
				plume_in%traj_coord_prev(3),',',						&
				plume_in%traj_vel(1),',',								&
				plume_in%traj_vel(2),',',								&
				plume_in%traj_vel(3),',',								&
				plume_in%radius_h,',',									&
				plume_in%eh(1),',',										&
				plume_in%eh(2),',',										&
				plume_in%eh(3),',',										&
				plume_in%radius_n,',',									&
				plume_in%en(1),',',										&
				plume_in%en(2),',',										&
				plume_in%en(3),','
	else
		write(ID_FILE_TRAJ)										&	
			plume_in%id,											&
			fcatime%current_int - fcatime%dt_int,			&
			plume_in%total_time,									&
			plume_in%traj_coord_prev(1),						&
			plume_in%traj_coord_prev(2),						&
			plume_in%traj_coord_prev(3),						&
			plume_in%traj_vel(1),								&
			plume_in%traj_vel(2),								&
			plume_in%traj_vel(3),								&
			plume_in%radius_h,									&
			plume_in%eh(1),										&
			plume_in%eh(2),										&
			plume_in%eh(3),										&
			plume_in%radius_n,									&
			plume_in%en(1),										&
			plume_in%en(2),										&
			plume_in%en(3)
	endif
	endif

	end
!=====================================================================================================================
!===================================================================================================================== 
	SUBROUTINE MergePlumes()
	
	use plume_module
	use rnd_num_vars
	use time_module
	
	implicit none
			
	real, parameter :: epsilon = 1.
	integer,dimension(total_fires) :: is_used	
	integer,dimension(:),allocatable :: to_merge
	integer :: nplumes_new,ip1,ip2,count,j,k,temp,is_valid,jp
	real ::										&
		rpw, rqw,								& ! m, projection of the ellipse radius on wc
		plume_dist								  ! m, minimum distance between the two plumes
	integer,dimension(total_fires) :: plume_index
	real, dimension(1) :: myrand
	real,dimension(3) :: pp1,pp2,pp3,pp4
	
	! https://tekpool.wordpress.com/2006/10/06/shuffling-shuffle-a-deck-of-cards-knuth-shuffle/
	plume_index = (/(ip1,ip1=1,total_fires)/)
	do j = total_fires,1,-1
		call random_number(myrand,stat(1))
		k = int(real(j)*myrand(1) + 1)
		temp = plume_index(k)
		plume_index(k) = plume_index(j)
		plume_index(j) = temp
	enddo
	
	is_used = 0
	nplumes_new = 0	
   
	allocate(to_merge(total_fires))
	
	to_merge = 0
	
	do j = 1,total_fires
		ip1 = plume_index(j)		
		count = 0
		if(is_used(ip1) == 0) then		
			is_used(ip1) = 1
			pp1 = plume(ip1)%traj_coord_prev
			pp2 = plume(ip1)%traj_coord
			
			do ip2 = 1, total_fires	
				
				if(is_used(ip2) == 0 .and. (plume(ip1)%i0 /= plume(ip2)%i0 .or. &
					plume(ip1)%j0 /= plume(ip2)%j0 .or. plume(ip1)%k0 /= plume(ip2)%k0)) then
					pp3 = plume(ip2)%traj_coord_prev
					pp4 = plume(ip2)%traj_coord
					
					!call AngleWithWind(pp2, pp4, plume(ip1)%traj_vel, is_valid)
					!if(is_valid == 1) then
						call AngleBetweenSegments(pp1, pp2, pp3, pp4, is_valid)

						if(is_valid == 1) then
							call DistBetween2Segment(pp1, pp2, pp3, pp4, ip1, ip2, plume_dist, rpw, rqw)
						
							! Either they are closer than the merging distance
							if(plume_dist < epsilon * min(rpw, rqw) + max(rpw, rqw)) then
							
								call CheckOverlap(pp1, pp2, pp3, pp4, is_valid)
								if(is_valid == 1) then
									
									! Check condition on initial location on all plumes already in the list
									is_valid = 1
									jp = 0
									do while(jp < count .and. is_valid == 1)
										jp = jp + 1
										
										if(plume(to_merge(jp))%i0 == plume(ip2)%i0 .and. &
											plume(to_merge(jp))%j0 == plume(ip2)%j0 .and. &
											plume(to_merge(jp))%k0 == plume(ip2)%k0)then
											is_valid = 0
										endif
									enddo
									if(is_valid == 1) then
										is_used(ip2) = 1
										count = count + 1
										to_merge(count) = ip2
										to_delete(ip2) = 1
									endif
								endif
							endif
						endif
					!endif
				endif
			enddo
			if(count > 0)then				
				call DoMerging(count,to_merge(1:count),ip1)
				to_merge = 0
			endif
		endif
	enddo
   deallocate(to_merge)
	
	end
!=====================================================================================================================
!===================================================================================================================== 
	SUBROUTINE AngleWithWind(pp2, pp4, traj_vel, is_valid)
	
	use constants
	
	implicit none
	
	real,dimension(3), intent(IN) :: pp2, pp4, traj_vel
	real,dimension(3) :: p1, p2
	integer, intent(OUT) :: is_valid
	real :: dotres, theta, p1_len, p2_len

	p1 = traj_vel - pp2
	p2 = pp4 - pp2
	
	dotres = p1(1)*p2(1) + p1(2)*p2(2)
	
	p1_len = p1(1)**2 + p1(2)**2
	p2_len = p2(1)**2 + p2(2)**2
	! Acos does not work with 1
	theta = acos(max(min(dotres / sqrt(p1_len*p2_len), 0.99999),-0.99999))
	
	if(abs(theta) < 30. * PI / 180.) then
		is_valid = 1
	else
		is_valid = 0
	endif
	
	END
!=====================================================================================================================
!===================================================================================================================== 
	SUBROUTINE AngleBetweenSegments(pp1, pp2, pp3, pp4, is_valid)

	use plume_const_module
	use constants
	
	implicit none
	
	real,dimension(3), intent(IN) :: pp1,pp2,pp3,pp4
	real,dimension(3) :: pp12, pp34
	integer, intent(OUT) :: is_valid
	real ::				&
		arg,				&
		dot_res,			& ! Result of the dot product
		len_1, len_2,	& ! m, length of the two segments (squared)
		theta				  ! rad, angle between the two segments
	integer :: i
	
	pp12 = pp2-pp1
	pp34 = pp4-pp3
	call dot(pp12, pp34, dot_res)
	len_1 = 0
	len_2 = 0
	do i = 1,3
		len_1 = len_1 + (pp2(i) - pp1(i))**2
		len_2 = len_2 + (pp4(i) - pp3(i))**2
	enddo
	arg = max(min(dot_res / sqrt(len_1 * len_2), 1.), -1.)
	theta = acos(min(max(-0.999999,arg),0.999999))
	
	! See if this angle is below max amount (30 deg)
	is_valid	= 0	
	if(abs(theta) < plume_const%max_angle) is_valid = 1 
	
	end
!=====================================================================================================================
!===================================================================================================================== 
	SUBROUTINE Cross(ab, ac, outres)
	
	implicit none
	
	real,intent(IN),dimension(3) :: ab, ac
	real,intent(OUT),dimension(3) :: outres
	
	outres(1) = ac(2) * ab(3) - ac(3) * ab(2)
	outres(2) = ac(3) * ab(1) - ac(1) * ab(3)
	outres(3) = ac(1) * ab(2) - ac(2) * ab(1)
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE ComputeRpComponents(plumeToMergeInto, ip1, to_merge, nelem, Qtot)
	
	use plume_module
	use constants
	
	implicit none
		
	TYPE(PlumeType) :: plumeToMergeInto
	integer, intent(IN) :: nelem, ip1
	integer, dimension(nelem), intent(IN) :: to_merge	
	real, intent(IN) :: Qtot
	integer :: i, ip2
	real :: max_rph, max_rpn
	
	! Compute rpi for the plume we merge into
	max_rph = 0.
	max_rpn = 0.
	plume(ip1)%radius_h = 0.
	plume(ip1)%radius_n = 0.
	
	call MaxValCalc(plumeToMergeInto, plume(ip1), max_rph, max_rpn)
	call CalcOtherRpTerm(plumeToMergeInto, plume(ip1))
	
	! Compute rpi for other plumes
	do i = 1, nelem
		ip2 = to_merge(i)
		call MaxValCalc(plume(ip2), plume(ip1), max_rph, max_rpn)
		call CalcOtherRpTerm(plume(ip2), plume(ip1))
	enddo
	
	max_rph = max_rph / real(nelem + 1)
	max_rpn = max_rpn / real(nelem + 1)
	
	! New primed values for the merged plume
	plume(ip1)%radius_h = max_rph + sqrt(plume(ip1)%radius_h / plume(ip1)%heat_of_fire)
	plume(ip1)%radius_n = max_rpn + sqrt(plume(ip1)%radius_n / plume(ip1)%heat_of_fire)
	
	! Recompte gamma
	plume(ip1)%gamma = abs(plume(ip1)%radius_n / plume(ip1)%radius_h)
	
	! Recompute plume radii
	plume(ip1)%radius_h = sqrt(Qtot / (plume(ip1)%gamma * plume(ip1)%traj_vel(3)))  ! Without PI because Qtot does not include it either
	plume(ip1)%radius_n = plume(ip1)%radius_h * plume(ip1)%gamma
		
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE CalcOtherRpTerm(p, pmerge)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: p
	TYPE(PlumeType) :: pmerge
	integer :: i
	real :: temp_h, temp_n, dff
	
	temp_h = 0.	
	temp_n = 0.
	do i = 1, 3
		dff = p%traj_coord(i) - pmerge%traj_coord(i)
		temp_h = temp_h + dff * pmerge%eh(i)
		temp_n = temp_n + dff * pmerge%en(i)
	enddo

	pmerge%radius_h = pmerge%radius_h + p%heat_of_fire * temp_h**2
	pmerge%radius_n = pmerge%radius_n + p%heat_of_fire * temp_n**2
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE MaxValCalc(p, pmerge, max_rph, max_rpn)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: p, pmerge
	real, intent(INOUT) :: max_rph, max_rpn
	real :: dot1, dot2
		
	call dot(pmerge%eh, p%eh, dot1)
	call dot(pmerge%eh, p%en, dot2)
	max_rph = max_rph	+							&
		sqrt(											&
		(p%radius_h * abs(dot1) )**2 +		&
		(p%radius_n * abs(dot2) )**2 )
	
	call dot(pmerge%en, p%eh, dot1)
	call dot(pmerge%en, p%en, dot2)

	max_rpn = max_rpn +							&
		sqrt(											&
		(p%radius_h * abs(dot1) )**2 +		&
		(p%radius_n * abs(dot2) )**2 )
	
	end
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE DoMerging(nelem,to_merge,ip1)
	
	use time_module
	use plume_module
	use file_handling_module
	use flags_module
	use plume_const_module
	use constants
	use grid_module
	
	implicit none
	
	TYPE(PlumeType) :: plumeToMergeInto
	integer,intent(IN) :: ip1, nelem
	integer,intent(IN),dimension(nelem) :: to_merge	
	integer :: ip2, count
	real :: Qtot	

	if(flag%plume_traj_file > 0) call WriteOutPlume(plume(ip1), 0)
	
	! Save copy of this plume for the computation of rpx, rpy, rpz
	plumeToMergeInto = plume(ip1)
	
	Qtot = plume(ip1)%radius_n * plume(ip1)%radius_h * plume(ip1)%traj_vel(3)

	plume(ip1)%traj_vel(3) = plume(ip1)%traj_vel(3)**plume_const%WC_EXPONENT
	plume(ip1)%wc_prev = plume(ip1)%wc_prev**plume_const%WC_EXPONENT
	plume(ip1)%traj_coord = plume(ip1)%traj_coord * plume(ip1)%heat_of_fire	
	plume(ip1)%traj_coord_prev = plume(ip1)%traj_coord_prev * plume(ip1)%heat_of_fire
		
	if(flag%plume_traj_file == 1) then
		write(ID_FILE_TRAJ_MERGE,'(3(i10, a))') plume(ip1)%id,',',max_id_plume + 1,',',fcatime%current_int - fcatime%dt_int, ','
	elseif(flag%plume_traj_file == 2) then
		write(ID_FILE_TRAJ_MERGE) plume(ip1)%id, max_id_plume + 1,fcatime%current_int - fcatime%dt_int
	endif
	
	max_id_plume = max_id_plume + 1
	plume(ip1)%id = max_id_plume
		
	do count = 1,nelem
		ip2 = to_merge(count)
		
		if(flag%plume_traj_file > 0) then			
			call WriteOutPlume(plume(ip2), 1)
			call WriteOutPlume(plume(ip2), 0)
		endif
		
		Qtot									= Qtot + plume(ip2)%radius_n * plume(ip2)%radius_h * plume(ip2)%traj_vel(3)
		
		plume(ip1)%heat_of_fire			= plume(ip1)%heat_of_fire + plume(ip2)%heat_of_fire			
		plume(ip1)%traj_coord			= plume(ip1)%traj_coord + plume(ip2)%traj_coord * plume(ip2)%heat_of_fire
			
		plume(ip1)%traj_coord_prev		= plume(ip1)%traj_coord_prev + plume(ip2)%traj_coord_prev * plume(ip2)%heat_of_fire
			
		plume(ip1)%traj_vel(3)			= plume(ip1)%traj_vel(3) + plume(ip2)%traj_vel(3)**plume_const%WC_EXPONENT
		plume(ip1)%wc_prev				= plume(ip1)%wc_prev + plume(ip2)%wc_prev**plume_const%WC_EXPONENT
		
		plume(ip1)%total_time			= max(plume(ip1)%total_time,plume(ip2)%total_time)
				   		
		if(flag%plume_traj_file == 1) then
			write(ID_FILE_TRAJ_MERGE,'(3(i10, a))')plume(ip2)%id, ',', plume(ip1)%id, ',', fcatime%current_int - fcatime%dt_int, ','
		elseif(flag%plume_traj_file == 2) then
			write(ID_FILE_TRAJ_MERGE)plume(ip2)%id, plume(ip1)%id, fcatime%current_int - fcatime%dt_int
		endif
	enddo
	
	
	plume(ip1)%traj_coord			= plume(ip1)%traj_coord / plume(ip1)%heat_of_fire
	
	plume(ip1)%traj_coord_prev		= plume(ip1)%traj_coord_prev / plume(ip1)%heat_of_fire	
		
	plume(ip1)%i						= ceiling(plume(ip1)%traj_coord(1) * qugrid%dxi)
	plume(ip1)%j						= ceiling(plume(ip1)%traj_coord(2) * qugrid%dyi)
	
	call GetPlumeK(plume(ip1)%traj_coord(3), plume(ip1)%k)
	
	if(plume_const%BRUNT_VAISALA_FREQ2 == 0) then		
		plume(ip1)%total_time = 0.			! plume is created like new		
	endif
	
	plume(ip1)%i0						= plume(ip1)%i
	plume(ip1)%j0						= plume(ip1)%j
	plume(ip1)%k0						= plume(ip1)%k
	
	plume(ip1)%traj_vel(3)			= plume(ip1)%traj_vel(3)**plume_const%INV_WC_EXPONENT
	plume(ip1)%wc_prev				= plume(ip1)%wc_prev**plume_const%INV_WC_EXPONENT
	
	call ComputeWinds(plume(ip1)%i, plume(ip1)%j, plume(ip1)%k,									&
		plume(ip1)%traj_coord(1), plume(ip1)%traj_coord(2), plume(ip1)%traj_coord(3),		&
      plume(ip1)%traj_vel(1), plume(ip1)%traj_vel(2), plume(ip1)%ws)
	
	! Plume trajectory properties
	call ComputeTrajectoryVector(plume(ip1))
	call ComputeTrajectoryLength(plume(ip1))
	! Plume base vectors
	call ComputeBaseVectors(plume(ip1))
				
	call ComputeRpComponents(plumeToMergeInto, ip1, to_merge, nelem, Qtot)
	
	! Recompute zvirt	and zfire_orig: set zvirt to zero and use zfire_orig at the point of origin based on the radius
	!plume(ip1)%zfire_orig = 2. *	 PI * plume(ip1)%gamma * plume(ip1)%radius_h / &
	!	(plume(ip1)%perimeter_over_rph * plume_const%BETA_ENTR)
	!plume(ip1)%total_z = plume(ip1)%zfire_orig
	!plume(ip1)%zvirt = 0.
	
	! Diffuse along the trajectory
	call Diffuse(plumeToMergeInto, ip1, nelem, to_merge)
		
	! Print both origin and destination of the merged plume
	if(flag%plume_traj_file > 0) then
		call WriteOutPlume(plume(ip1), 1)
		call WriteOutPlume(plume(ip1), 0)
	endif	
	
	end	
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE Diffuse(plumeToMergeInto, ip1, nelem, to_merge)
	
	use plume_module
	use grid_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: plumeToMergeInto	
	integer, intent(IN) :: ip1, nelem
	integer, intent(IN), dimension(nelem) :: to_merge	
	integer :: ip2, count, has_changed
	real, dimension(3) :: temp_coord, temp_coord_prev
	real, dimension(2) :: temp
	real, dimension(2, nelem+1) :: projval
	real :: minv, maxv	
	
	call ComputeSegmProj(plumeToMergeInto, plume(ip1), temp)
	projval(1, 1) = temp(1)
	projval(2, 1) = temp(2)	
	do count = 1,nelem
		ip2 = to_merge(count)
	
		call ComputeSegmProj(plume(ip2), plume(ip1), temp)
		projval(:, count + 1) = temp
	enddo
	minv = minval(projval(1,:))
	maxv = maxval(projval(2,:))
	has_changed	= 0
	
	temp_coord = plume(ip1)%traj_coord
	temp_coord_prev = plume(ip1)%traj_coord_prev
	if(minv < 0) then
		has_changed = 1
		plume(ip1)%traj_coord_prev = temp_coord - plume(ip1)%s * (1. - minv / plume(ip1)%traj_len)
	endif
	if(maxv > plume(ip1)%traj_len) then
		has_changed = 1
		plume(ip1)%traj_coord = temp_coord_prev + plume(ip1)%s * (maxv / plume(ip1)%traj_len)
	endif
	
	if(has_changed == 1) then
		! Prevent changed plume from going outside the domain
		call BoundPlumeCoord(plume(ip1)%traj_coord)
		call BoundPlumeCoord(plume(ip1)%traj_coord_prev)
			
		call GetPlumeCellIndex(plume(ip1))
		
		call ComputeWinds(plume(ip1)%i, plume(ip1)%j, plume(ip1)%k,									&
			plume(ip1)%traj_coord(1), plume(ip1)%traj_coord(2), plume(ip1)%traj_coord(3),		&
			plume(ip1)%traj_vel(1), plume(ip1)%traj_vel(2), plume(ip1)%ws)
	
		call ComputeTrajectoryVector(plume(ip1))
		call ComputeTrajectoryLength(plume(ip1))
		! Plume base vectors
		call ComputeBaseVectors(plume(ip1))
		
	endif
	
	
	end	
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE BoundPlumeCoord(traj_coord)
	
	use grid_module
	
	implicit none
	
	real, dimension(3) :: traj_coord
	real, parameter :: delta = 0.1
	
	traj_coord(1) = min(traj_coord(1), qugrid%Lx-delta*qugrid%dx)
	traj_coord(1) = max(traj_coord(1), delta*qugrid%dx)
	traj_coord(2) = min(traj_coord(2), qugrid%Ly-delta*qugrid%dy)
	traj_coord(2) = max(traj_coord(2), delta*qugrid%dy)
	traj_coord(3) = min(traj_coord(3), qugrid%Lz-delta*qugrid%dz_array(qugrid%nz))
	traj_coord(3) = max(traj_coord(3), delta*delta*qugrid%dz_array(1))		
		
	end	
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE ComputeSegmProj(p, pmerged, projval)
	
	use plume_module
	
	implicit none
	
	TYPE(PlumeType), intent(IN) :: p, pmerged
	real :: proj
	real, dimension(3) :: ph
	real, dimension(2) :: projval
	
	projval(1) = +1e8  ! min
	projval(2) = -1e8  ! max
	
	! ------------------- Direction of rh
	
	! (1) ph: end of the merged plume traj, along rh, one side
	ph = p%traj_coord + p%eh * p%radius_h - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (2) ph: end of the merged plume traj, along rh, other side
	ph = p%traj_coord - p%eh * p%radius_h - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (3) ph: beginning of the merged plume traj, along rh, one side
	ph = p%traj_coord_prev + p%eh * p%radius_h - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (4) ph: beginning of the merged plume traj, along rh, othr side
	ph = p%traj_coord_prev - p%eh * p%radius_h - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! ------------------- Direction of rn
	
	! (1) ph: end of the merged plume traj, along rh, one side
	ph = p%traj_coord + p%en * p%radius_n - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (2) ph: end of the merged plume traj, along rh, other side
	ph = p%traj_coord - p%en * p%radius_n - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (3) ph: beginning of the merged plume traj, along rh, one side
	ph = p%traj_coord_prev + p%en * p%radius_n - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	! (4) ph: beginning of the merged plume traj, along rh, othr side
	ph = p%traj_coord_prev - p%en * p%radius_n - pmerged%traj_coord_prev
	call dot(pmerged%s, ph, proj)
	proj = proj / pmerged%traj_len
		
	projval(1) = min(projval(1), proj)
	projval(2) = max(projval(2), proj)
	
	END
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE CheckMaxHeight()
	
	use plume_module
	
	implicit none
		
	integer :: ip
	
	do ip = 1,total_fires
		if(plume(ip)%zmax < plume(ip)%traj_coord(3)) then
			to_delete(ip) = 1
		endif
	enddo
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	SUBROUTINE FillWW()
	
	use grid_module
	use plume_module
	use winds_module
	
	implicit none
	
	integer ::				&
		imin, imax,			&
		jmin, jmax,			&
		kmin, kmax,			&
		ip, i, j, k
	
	imin = qugrid%nx + 1
	imax = 0
	jmin = qugrid%ny + 1
	jmax = 0
	kmin = qugrid%nz + 2
	kmax = 0
	
	do ip = 1,total_fires
		imin = min(imin, plume(ip)%istart)
		jmin = min(jmin, plume(ip)%jstart)
		kmin = min(kmin, plume(ip)%kstart)
		
		imax = max(imax, plume(ip)%iend)
		jmax = max(jmax, plume(ip)%jend)
		kmax = max(kmax, plume(ip)%kend)
	enddo
	
	imin = max(imin, 2)
	imax = min(imax, qugrid%nx-2)
	jmin = max(jmin, 2)
	jmax = min(jmax, qugrid%ny-2)
	kmin = max(kmin, 3)
	kmax = min(kmax, qugrid%nz-2)
	
	!$OMP parallel do private(i,j,k,ip)
	do ip = 1,2 ! Fill 2-cell gaps
		do k = kmin, kmax
			do j = jmin, jmax
				do i = imin, imax
					if(quwinds%wplume(i,j,k) == 0) then
						if(quwinds%wplume(i-1,j,k) /= 0. .and. quwinds%wplume(i+1,j,k) /= 0.) then
							quwinds%wplume(i,j,k) = 0.5 * (quwinds%wplume(i-1,j,k) + quwinds%wplume(i+1,j,k))
						elseif(quwinds%wplume(i,j-1,k) /= 0. .and. quwinds%wplume(i,j+1,k) /= 0.) then
							quwinds%wplume(i,j,k) = 0.5 * (quwinds%wplume(i,j-1,k) + quwinds%wplume(i,j+1,k))
						elseif(quwinds%wplume(i,j,k-1) /= 0. .and. quwinds%wplume(i,j,k+1) /= 0.) then
							! w(k) is in z(k-1)
							quwinds%wplume(i,j,k) = quwinds%wplume(i,j,k-1) + (quwinds%wplume(i,j,k+1)-quwinds%wplume(i,j,k-1))/ 	&
								(qugrid%z(k)-qugrid%z(k-2))*(qugrid%z(k-1)-qugrid%z(k-2))
						endif
					endif
				enddo
			enddo
		enddo
	enddo	
	!$OMP end parallel do
	
	end
	!=====================================================================================================================
	!=====================================================================================================================
	