	subroutine update_canopy()

	
	use fireca_module
	use fireca_constants_module
	use omp_lib
	use rnd_num_vars
	use canopy_module
	use constants
	use plume_module
	use bld_module
	use sor_module
	use grid_module
	use interpolation_module
	use winds_module

	implicit none	

	real :: vegvelfrac,vegvelfrac0,avg_atten,num_atten,ncell,ucan,ulog,density0,ave_density
	integer ::index,i,j,k,TID,isrun
	integer,dimension(numThreads) :: iimin,iimax,jjmin,jjmax,kkmin,kkmax
		
	iimin = qugrid%nx + 1
	iimax = 0
	jjmin = qugrid%ny + 1
	jjmax = 0
	kkmin = qugrid%nz+3
	kkmax = 0
	isrun = 0
	!$OMP parallel do private(i,j,k,index, &
	!$OMP ave_density,density0,ncell, &
	!$OMP avg_atten,num_atten,vegvelfrac,ucan,ulog,vegvelfrac0,TID)
	do index = 1,qugrid%num_fuel_cells
		TID = 1
		!$ TID = TID + OMP_GET_THREAD_NUM()
		
		i = qugrid%ijk_cell_index(index,1)
		j = qugrid%ijk_cell_index(index,2)
		k = qugrid%ijk_cell_index(index,3)
		
		ave_density = sum(fuels%density(grid_conv(i,j,k)%index) * grid_conv(i,j,k)%en_ratio)
		density0 = init_dens_canopy(index) !sum(fuels%density_initial(grid_conv(i,j,k)%index) * grid_conv(i,j,k)%en_ratio)
			
		ncell = 1. / real(grid_conv(i,j,k)%nelem)
		ave_density = ave_density * ncell
		density0 = density0 * ncell
			
		if(ave_density < density0) then
			if(ave_density >= 5e-6)then
				
				avg_atten = canopy%atten(i,j,k)
				if(canopy%atten(i,j,k+1) .ne. canopy%atten(i,j,k) &
						.or. canopy%atten(i,j,k-1) .ne. canopy%atten(i,j,k))then
					num_atten=1.
					if(canopy%atten(i,j,k+1) .gt. 0.)then
						avg_atten = avg_atten + canopy%atten(i,j,k+1)
						num_atten=num_atten+1.
					endif
					if(canopy%atten(i,j,k-1) .gt. 0.)then
						avg_atten = avg_atten + canopy%atten(i,j,k-1)
						num_atten=num_atten+1.
					endif
					avg_atten=avg_atten/num_atten
				endif
					
				if(avg_atten > 0) then
					if(bld%icellflag(i,j,k) == 0)then
						! ---- Solid building
						ulog = canopy%ustar(i,j)/vk*log(qugrid%zm(k)/canopy%zo(i,j))
						vegvelfrac = exp(canopy%atten(i,j,k)*(qugrid%zm(k)/canopy%top(i,j)-1.))
						ucan = quwinds%u_mc(i,j,canopy%ktop(i,j)-1) * vegvelfrac
						quwinds%uo(i,j,k) = ulog + (ucan - ulog) * ave_density / density0
						ucan = quwinds%v_mc(i,j,canopy%ktop(i,j)-1) * vegvelfrac
						quwinds%vo(i,j,k) = ulog + (ucan - ulog) * ave_density / density0						
					else
						! ---- Vegetation or building canopy
						! reduce speed of the log law
						if(grid_conv(i,j,k)%nelem > 0) then
							vegvelfrac0 = log((canopy%top(i,j)-canopy%d(i,j))/canopy%zo(i,j))*&
											exp(avg_atten*((qugrid%zm(k)/canopy%top(i,j))-1.))/&
											log(qugrid%zm(k)/canopy%zo(i,j))
							
							vegvelfrac = vegvelfrac0 + (vegvelfrac0-1.)*(ave_density / density0-1.)
						
							!		current deviation from log prof    vel after bld params       deviation due to full vegetation
							! uo = (u0_bef_params * vegvelfrac)    +     u_bld_aware        -   (u0_bef_params * vegvelfrac0)
							!quwinds%uo(i,j,k) = quwinds%u_mc(i,j,k) - quwinds%uo_before_params(index)*(vegvelfrac0 - vegvelfrac)
							!quwinds%vo(i,j,k) = quwinds%v_mc(i,j,k) - quwinds%vo_before_params(index)*(vegvelfrac0 - vegvelfrac)
														
							quwinds%uo(i,j,k) = quwinds%u_mc(i,j,k) / vegvelfrac0 * vegvelfrac
							quwinds%vo(i,j,k) = quwinds%v_mc(i,j,k) / vegvelfrac0 * vegvelfrac
						else
							canopy%atten(i,j,k) = 0.
							quwinds%uo(i,j,k) = quwinds%u_mc(i,j,k) !quwinds%uo_before_params(index)
							quwinds%vo(i,j,k) = quwinds%v_mc(i,j,k) !quwinds%vo_before_params(index)
						endif						
					endif
					isrun = 1
					iimin(TID) = min(iimin(TID), i)
					iimax(TID) = max(iimax(TID), i)
					jjmin(TID) = min(jjmin(TID), j)
					jjmax(TID) = max(jjmax(TID), j)
					kkmin(TID) = min(kkmin(TID), k)
					kkmax(TID) = max(kkmax(TID), k)
				endif
			else
				! remove canopy
				canopy%atten(i,j,k) = 0.
				quwinds%uo(i,j,k) = quwinds%u_mc(i,j,k) !quwinds%uo_before_params(index)
				quwinds%vo(i,j,k) = quwinds%v_mc(i,j,k) !quwinds%vo_before_params(index)
			endif
		endif
	enddo
	!$OMP END PARALLEL DO
	
	if(total_fires > 0 .and. isrun == 1) then
		sor%mc_is = max(minval(iimin)-sor%mc_buffer,1)
		sor%mc_ie = min(maxval(iimax)+sor%mc_buffer,qugrid%nx)
		sor%mc_js = max(minval(jjmin)-sor%mc_buffer,1)
		sor%mc_je = min(maxval(jjmax)+sor%mc_buffer,qugrid%ny)
		sor%mc_ks = max(minval(kkmin)-sor%mc_buffer,1)
		sor%mc_ke = min(maxval(kkmax)+sor%mc_buffer,qugrid%nz+2)		
	endif
	
	return

	end
!============================================================================
!============================================================================
	subroutine ConvertBldgToCanopy()

	
	use constants
	use canopy_module
	use bld_module
	use grid_module
	use winds_module
	use interpolation_module
	
	implicit none

	integer :: k_above,i,j,k
	real :: theta_site,mag_site,bisect,avg_atten,vegvelfrac,num_atten

	if(allocated(canopy%atten) .eqv. .false.)then
		allocate(canopy%atten(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
		canopy%atten = 0.
	endif
	if(allocated(canopy%zo) .eqv. .false.)then
		allocate(canopy%zo(qugrid%nx-1,qugrid%ny-1))
		canopy%zo = 0.
	endif		
	if(allocated(canopy%d) .eqv. .false.)then
		allocate(canopy%d(qugrid%nx-1,qugrid%ny-1))
		canopy%d = 0
	endif	
	if(allocated(canopy%ustar) .eqv. .false.)then
		allocate(canopy%ustar(qugrid%nx-1,qugrid%ny-1))
		canopy%ustar = 0
	endif
	if(allocated(canopy%ktop) .eqv. .false.)then
		allocate(canopy%ktop(qugrid%nx-1,qugrid%ny-1))
		canopy%ktop = 0
	endif
	if(allocated(canopy%top) .eqv. .false.)then
		allocate(canopy%top(qugrid%nx-1,qugrid%ny-1))
		canopy%top = 0
	endif
	
	!$OMP parallel do private(i,j,k_above,mag_site,theta_site, &
	!$OMP avg_atten,num_atten,vegvelfrac)
	do j = 1,qugrid%ny-1
		do i = 1,qugrid%nx-1
			if(any(bld%icellflag(i,j,2:qugrid%nz-1) == 0)) then   ! at least one cell is in a building
				
				! - Define canopy%top
				k_above = 2
				do while(k_above < qugrid%nz-1  .and. bld%icellflag(i,j,k_above) == 0)
					k_above = k_above + 1
				enddo
				canopy%ktop(i,j) = k_above-1
				canopy%top(i,j) = qugrid%z(k_above)
				
				! - specify z0
				canopy%zo(i,j) = canopy%BLDG_canopy_Z0
				
				! - specify attenuation coefficient
				canopy%atten(i,j,2:k_above-1) = canopy%BLDG_ATTEN_COEFF
				
				! - use the first cell that is not a building to define the reference velocity
				mag_site = sqrt(quwinds%u(i,j,k_above)**2 + quwinds%v(i,j,k_above)**2)
				do while(mag_site == 0 .and. k_above < qugrid%nz+2)
					k_above = k_above + 1
					mag_site = sqrt(quwinds%u(i,j,k_above)**2 + quwinds%v(i,j,k_above)**2)
				enddo
				if(mag_site == 0) then
					write(msgoutfile,*) 'No winds above building at i=',i,', j=',j
					call TerminateProgram()
				endif
				
				theta_site = atan2(quwinds%v(i,j,k_above) , quwinds%u(i,j,k_above))
				! - the canopy ustar is calculated from the log-profile
				canopy%ustar(i,j) = mag_site * vk / (log(canopy%top(i,j)/canopy%zo(i,j)))
				
				! - define canopy displacement height for neutral stability
				canopy%d(i,j) = bisect(canopy%ustar(i,j),canopy%zo(i,j),	&
					canopy%top(i,j),canopy%BLDG_ATTEN_COEFF,0.)
				if(canopy%d(i,j) .gt. 0.99*canopy%top(i,j))then
					canopy%d(i,j) = 0.7*canopy%top(i,j)
				   canopy%zo(i,j) = min(0.1*canopy%top(i,j),0.5*qugrid%dz_array(1))
				endif
				
				! - profile
            do k = 2,canopy%ktop(i,j)
					if(grid_conv(i,j,k)%nelem > 0) then
						avg_atten = canopy%atten(i,j,k)
						if(canopy%atten(i,j,k+1) .ne. canopy%atten(i,j,k) &
								.or. canopy%atten(i,j,k-1) .ne. canopy%atten(i,j,k))then
							num_atten=1.
							if(canopy%atten(i,j,k+1) .gt. 0.)then
								avg_atten = avg_atten + canopy%atten(i,j,k+1)
								num_atten=num_atten+1.
							endif
							if(canopy%atten(i,j,k-1) .gt. 0.)then
								avg_atten = avg_atten + canopy%atten(i,j,k-1)
								num_atten=num_atten+1
							endif
							avg_atten=avg_atten/num_atten
						endif

						vegvelfrac = min(exp(avg_atten*((qugrid%zm(k)/canopy%top(i,j))-1.)),1.)
						quwinds%u(i,j,k) = quwinds%u(i,j,k_above)*vegvelfrac
						quwinds%v(i,j,k) = quwinds%v(i,j,k_above)*vegvelfrac
					endif
				enddo
			endif
		enddo
	enddo
	!$OMP end parallel do 
	
	end
!============================================================================
!============================================================================
	subroutine ConvertFuelToCanopy()

	use fireca_module
	
	use constants
	use canopy_module
	use grid_module
	use winds_module
	use interpolation_module

	implicit none

	integer :: k_above,index,i,j,k
	real :: theta_site,mag_site,bisect,avg_atten,vegvelfrac,vegvelfrac0,num_atten

	if(allocated(canopy%atten) .eqv. .false.)then
		allocate(canopy%atten(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
		canopy%atten = 0.
	endif
	if(allocated(canopy%zo) .eqv. .false.)then
		allocate(canopy%zo(qugrid%nx-1,qugrid%ny-1))
		canopy%zo = 0.
	endif		
	if(allocated(canopy%d) .eqv. .false.)then
		allocate(canopy%d(qugrid%nx-1,qugrid%ny-1))
		canopy%d = 0
	endif	
	if(allocated(canopy%ustar) .eqv. .false.)then
		allocate(canopy%ustar(qugrid%nx-1,qugrid%ny-1))
		canopy%ustar = 0
	endif
	if(allocated(canopy%ktop) .eqv. .false.)then
		allocate(canopy%ktop(qugrid%nx-1,qugrid%ny-1))
		canopy%ktop = 0
	endif
	if(allocated(canopy%top) .eqv. .false.)then
		allocate(canopy%top(qugrid%nx-1,qugrid%ny-1))
		canopy%top = 0
	endif
	
	!$OMP parallel do private(i,j,k,k_above,mag_site,theta_site, &
	!$OMP avg_atten,num_atten,vegvelfrac)
	do j = 1,qugrid%ny-1
		do i = 1,qugrid%nx-1
			!if(canopy%top(i,j) == 0) then				
				! is there fuel there?
				do k = 2,qugrid%kfire_top
					if(grid_conv(i,j,k)%nelem > 0) canopy%ktop(i,j) = k
				enddo
						
				if(canopy%ktop(i,j) > 0) then
					k_above = min(canopy%ktop(i,j)+1,qugrid%nz)
					canopy%top(i,j) = qugrid%z(canopy%ktop(i,j))
					! - specify z0
					canopy%zo(i,j) = canopy%FUEL_CANOPY_Z0
					! - specify attenuation coefficient (only if not already specified via veg canopy)
					do k = 2,canopy%ktop(i,j)
						if(canopy%atten(i,j,k) <= 0) canopy%atten(i,j,k) = canopy%FUEL_ATTEN_COEFF
					enddo
					
					! - use the first cell that is not a building to define the reference velocity
					mag_site = sqrt(quwinds%u(i,j,k_above)**2 + quwinds%v(i,j,k_above)**2)
					do while(mag_site == 0 .and. k_above < qugrid%nz+2)
						k_above = k_above + 1
						mag_site = sqrt(quwinds%u(i,j,k_above)**2 + quwinds%v(i,j,k_above)**2)
					enddo
					if(mag_site == 0) then
						write(msgoutfile,*) 'No winds above building at i=',i,', j=',j
						call TerminateProgram()
					endif
				
					theta_site = atan2(quwinds%v(i,j,k_above) , quwinds%u(i,j,k_above))
					! - the canopy ustar is calculated from the log-profile
					canopy%ustar(i,j) = mag_site * vk / (log(canopy%top(i,j)/canopy%zo(i,j)))
				
					! - define canopy displacement height for neutral stability
					canopy%d(i,j) = bisect(canopy%ustar(i,j),canopy%zo(i,j),canopy%top(i,j),canopy%FUEL_ATTEN_COEFF,0.)
					if(canopy%d(i,j) .gt. 0.99*canopy%top(i,j))then
						canopy%d(i,j) = 0.7*canopy%top(i,j)
						canopy%zo(i,j) = min(0.1*canopy%top(i,j),0.5*qugrid%dz_array(1))
					endif
				
					! - profile
					do k = 2,canopy%ktop(i,j)
						if(grid_conv(i,j,k)%nelem > 0) then
							avg_atten = canopy%atten(i,j,k)
							if(canopy%atten(i,j,k+1) .ne. canopy%atten(i,j,k) &
									.or. canopy%atten(i,j,k-1) .ne. canopy%atten(i,j,k))then
								num_atten=1.
								if(canopy%atten(i,j,k+1) .gt. 0.)then
									avg_atten = avg_atten + canopy%atten(i,j,k+1)
									num_atten=num_atten+1.
								endif
								if(canopy%atten(i,j,k-1) .gt. 0.)then
									avg_atten = avg_atten + canopy%atten(i,j,k-1)
									num_atten=num_atten+1.
								endif
								avg_atten=avg_atten/num_atten
							endif

							vegvelfrac = min(exp(avg_atten*((qugrid%zm(k)/canopy%top(i,j))-1.)),1.)
						
							quwinds%u(i,j,k) = quwinds%u(i,j,k_above)*vegvelfrac
							quwinds%v(i,j,k) = quwinds%v(i,j,k_above)*vegvelfrac
						endif
					enddo
				endif
			!endif
		enddo
	enddo
	!$OMP end parallel do	
	
	do index = 1,qugrid%num_fuel_cells
		i = qugrid%ijk_cell_index(index,1)
		j = qugrid%ijk_cell_index(index,2)
		k = qugrid%ijk_cell_index(index,3)
			
		call AssignInitVel(i,j,k,vegvelfrac0,canopy%atten,canopy%top(i,j), &
			canopy%d(i,j),canopy%zo(i,j),qugrid%zm(k),qugrid%nx,qugrid%ny,qugrid%nz)
				
		quwinds%uo_before_params(index) = quwinds%u(i,j,k) / vegvelfrac0
		quwinds%vo_before_params(index) = quwinds%v(i,j,k) / vegvelfrac0
	enddo
	
	end
!============================================================================
!============================================================================
	subroutine AssignInitVel(i,j,k,vegvelfrac0,canopy_atten,canopy_top, &
		canopy_d,canopy_zo,zm,nx,ny,nz)	
	
	implicit none
	integer,intent(IN) :: i,j,k,nx,ny,nz
	real,intent(IN) :: canopy_top,canopy_d,canopy_zo,zm
	real :: vegvelfrac0
	real :: avg_atten,num_atten
	real,intent(IN),dimension(nx-1,ny-1,nz-1) :: canopy_atten
	
	avg_atten = canopy_atten(i,j,k)
	if(canopy_atten(i,j,k+1) .ne. canopy_atten(i,j,k) &
			.or. canopy_atten(i,j,k-1) .ne. canopy_atten(i,j,k))then
		num_atten=1.
		if(canopy_atten(i,j,k+1) .gt. 0.)then
			avg_atten = avg_atten + canopy_atten(i,j,k+1)
			num_atten=num_atten+1.
		endif
		if(canopy_atten(i,j,k-1) .gt. 0.)then
			avg_atten = avg_atten + canopy_atten(i,j,k-1)
			num_atten=num_atten+1.
		endif
		avg_atten=avg_atten/num_atten
	endif
		
	vegvelfrac0 = log((canopy_top-canopy_d)/canopy_zo)*	&
		exp(avg_atten*((zm/canopy_top)-1.))/					&
		log(zm/canopy_zo)
	
	end
!============================================================================
!============================================================================
