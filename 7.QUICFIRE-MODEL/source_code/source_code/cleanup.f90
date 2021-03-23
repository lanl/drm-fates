	SUBROUTINE CLEANUP()

	use fireca_module
	use rnd_num_vars
	use file_handling_module
	use interpolation_module
	use sor_module
	use canopy_module
	use ship_module
	use flags_module
	use wind_profile_module
	use bld_module
	use landuse_module
	use diffusion_module
	use ignitions_module
	use grid_module
	use plume_module
	use flags_module
	use winds_module
	use time_module

	implicit none
	
	integer :: ew, ns, ud, i

	IF(ALLOCATED(lock))then
		!$OMP parallel do private(ew, ns, ud)
		do ud = 1, firegrid%nz_en2atmos
			do ns = 1, firegrid%ny
				do ew	= 1, firegrid%nx
					CALL omp_destroy_lock (lock(ew, ns, ud))
				enddo
			enddo
		enddo
		!$OMP end parallel do
		DEALLOCATE(lock)
	endif
	if(allocated(fire_spread_lock))then
		!$OMP parallel do private(ud)
		do ud = 1, firegrid%num_fuel_cells
			call omp_destroy_lock(fire_spread_lock(ud))
		enddo
		!$OMP end parallel do 
	endif 
	
	IF(ALLOCATED(bld%icellflag))DEALLOCATE(bld%icellflag)
	IF(ALLOCATED(bld%flag%isbld))DEALLOCATE(bld%flag%isbld)
	IF(ALLOCATED(bld%istart))DEALLOCATE(bld%istart)
	IF(ALLOCATED(bld%iend))DEALLOCATE(bld%iend)
	IF(ALLOCATED(bld%jstart))DEALLOCATE(bld%jstart)
	IF(ALLOCATED(bld%jend))DEALLOCATE(bld%jend)
	IF(ALLOCATED(bld%kstart))DEALLOCATE(bld%kstart)
	IF(ALLOCATED(bld%kend))DEALLOCATE(bld%kend)
	IF(ALLOCATED(bld%Lf))DEALLOCATE(bld%Lf)
	IF(ALLOCATED(bld%Lr))DEALLOCATE(bld%Lr)
	IF(ALLOCATED(bld%Weff))DEALLOCATE(bld%Weff)
	IF(ALLOCATED(bld%Leff))DEALLOCATE(bld%Leff)
	IF(ALLOCATED(bld%Lt))DEALLOCATE(bld%Lt)
	IF(ALLOCATED(bld%Rscale))DEALLOCATE(bld%Rscale)
	IF(ALLOCATED(bld%Rcx))DEALLOCATE(bld%Rcx)
	IF(ALLOCATED(bld%Wt))DEALLOCATE(bld%Wt)
	IF(ALLOCATED(bld%num))DEALLOCATE(bld%num)
	IF(ALLOCATED(bld%btype))DEALLOCATE(bld%btype)
	IF(ALLOCATED(bld%geometry))DEALLOCATE(bld%geometry)
	IF(ALLOCATED(bld%invnum))DEALLOCATE(bld%invnum)
	IF(ALLOCATED(bld%startidx))DEALLOCATE(bld%startidx)
	IF(ALLOCATED(bld%stopidx))DEALLOCATE(bld%stopidx)
	IF(ALLOCATED(bld%numpolygons))DEALLOCATE(bld%numpolygons)
	IF(ALLOCATED(bld%roof))DEALLOCATE(bld%roof)
	IF(ALLOCATED(bld%x))DEALLOCATE(bld%x)
	IF(ALLOCATED(bld%y))DEALLOCATE(bld%y)
	IF(ALLOCATED(bld%cx))DEALLOCATE(bld%cx)
	IF(ALLOCATED(bld%cy))DEALLOCATE(bld%cy)
	IF(ALLOCATED(bld%wall))DEALLOCATE(bld%wall)
	IF(ALLOCATED(quwinds%uo))DEALLOCATE(quwinds%uo)
	IF(ALLOCATED(quwinds%vo))DEALLOCATE(quwinds%vo)
	IF(ALLOCATED(quwinds%wo))DEALLOCATE(quwinds%wo)
	IF(ALLOCATED(bld%ufarwake))DEALLOCATE(bld%ufarwake)
	IF(ALLOCATED(bld%vfarwake))DEALLOCATE(bld%vfarwake)
	IF(ALLOCATED(sor%p1))DEALLOCATE(sor%p1)
	IF(ALLOCATED(sor%p2))DEALLOCATE(sor%p2)
	IF(ALLOCATED(sor%r))DEALLOCATE(sor%r)
	IF(ALLOCATED(sor%bc%e))DEALLOCATE(sor%bc%e)
	IF(ALLOCATED(sor%bc%f))DEALLOCATE(sor%bc%f)
	IF(ALLOCATED(sor%bc%g))DEALLOCATE(sor%bc%g)
	IF(ALLOCATED(sor%bc%h))DEALLOCATE(sor%bc%h)
	IF(ALLOCATED(sor%bc%m))DEALLOCATE(sor%bc%m)
	IF(ALLOCATED(sor%bc%n))DEALLOCATE(sor%bc%n)
	IF(ALLOCATED(sor%denom))DEALLOCATE(sor%denom)
	IF(ALLOCATED(sor%denom0))DEALLOCATE(sor%denom0)
	IF(ALLOCATED(quwinds%u))DEALLOCATE(quwinds%u)
	IF(ALLOCATED(quwinds%v))DEALLOCATE(quwinds%v)
	IF(ALLOCATED(quwinds%w))DEALLOCATE(quwinds%w)
	IF(ALLOCATED(bld%Ht))DEALLOCATE(bld%Ht)
	IF(ALLOCATED(bld%Wti))DEALLOCATE(bld%Wti)
	IF(ALLOCATED(bld%Lti))DEALLOCATE(bld%Lti)
	IF(ALLOCATED(bld%aa))DEALLOCATE(bld%aa)
	IF(ALLOCATED(bld%bb))DEALLOCATE(bld%bb)
	IF(ALLOCATED(bld%xfo))DEALLOCATE(bld%xfo)
	IF(ALLOCATED(bld%yfo))DEALLOCATE(bld%yfo)
	IF(ALLOCATED(bld%zfo))DEALLOCATE(bld%zfo)
	IF(ALLOCATED(bld%gamma))DEALLOCATE(bld%gamma)
	IF(ALLOCATED(windProfile%ws_data))DEALLOCATE(windProfile%ws_data)
	IF(ALLOCATED(windProfile%wd_data))DEALLOCATE(windProfile%wd_data)
	IF(ALLOCATED(windProfile%z_data))DEALLOCATE(windProfile%z_data)
	IF(ALLOCATED(qutime%unix))DEALLOCATE(qutime%unix)
	IF(ALLOCATED(windProfile%site%xcoord))DEALLOCATE(windProfile%site%xcoord)
	IF(ALLOCATED(windProfile%site%ycoord))DEALLOCATE(windProfile%site%ycoord)
	IF(ALLOCATED(windProfile%site%u_prof))DEALLOCATE(windProfile%site%u_prof)
	IF(ALLOCATED(windProfile%site%v_prof))DEALLOCATE(windProfile%site%v_prof)
	IF(ALLOCATED(windProfile%site%uoint))DEALLOCATE(windProfile%site%uoint)
	IF(ALLOCATED(windProfile%site%voint))DEALLOCATE(windProfile%site%voint)
	IF(ALLOCATED(windProfile%site%wm))DEALLOCATE(windProfile%site%wm)
	IF(ALLOCATED(windProfile%site%wms))DEALLOCATE(windProfile%site%wms)
	IF(ALLOCATED(windProfile%site%z_data))DEALLOCATE(windProfile%site%z_data)
	IF(ALLOCATED(windProfile%site%ws_data))DEALLOCATE(windProfile%site%ws_data)
	IF(ALLOCATED(windProfile%site%wd_data))DEALLOCATE(windProfile%site%wd_data)
	IF(ALLOCATED(windProfile%site%u_data))DEALLOCATE(windProfile%site%u_data)
	IF(ALLOCATED(windProfile%site%v_data))DEALLOCATE(windProfile%site%v_data)
	IF(ALLOCATED(windProfile%site%nz_data))DEALLOCATE(windProfile%site%nz_data)
	IF(ALLOCATED(windProfile%site%blayer_flag))DEALLOCATE(windProfile%site%blayer_flag)
	IF(ALLOCATED(windProfile%site%pp))DEALLOCATE(windProfile%site%pp)
	IF(ALLOCATED(windProfile%site%H))DEALLOCATE(windProfile%site%H)
	IF(ALLOCATED(windProfile%site%ac))DEALLOCATE(windProfile%site%ac)
	IF(ALLOCATED(windProfile%site%rL))DEALLOCATE(windProfile%site%rL)
	IF(ALLOCATED(canopy%ktop))DEALLOCATE(canopy%ktop)
	IF(ALLOCATED(canopy%top))DEALLOCATE(canopy%top)
	IF(ALLOCATED(canopy%atten))DEALLOCATE(canopy%atten)
	IF(ALLOCATED(canopy%zo))DEALLOCATE(canopy%zo)
	IF(ALLOCATED(canopy%ustar))DEALLOCATE(canopy%ustar)
	IF(ALLOCATED(canopy%d))DEALLOCATE(canopy%d)
	IF(ALLOCATED(bld%atten))DEALLOCATE(bld%atten)
	IF(ALLOCATED(bld%group_id))DEALLOCATE(bld%group_id)
	IF(ALLOCATED(bld%zfo_actual))DEALLOCATE(bld%zfo_actual)
	IF(ALLOCATED(diff%visc))DEALLOCATE(diff%visc)
	IF(ALLOCATED(diff%Fxd))DEALLOCATE(diff%Fxd)
	IF(ALLOCATED(diff%Fyd))DEALLOCATE(diff%Fyd)
	IF(ALLOCATED(diff%Fzd))DEALLOCATE(diff%Fzd)
	IF(ALLOCATED(bld%rooftop_flag))DEALLOCATE(bld%rooftop_flag)
	IF(ALLOCATED(quwinds%uo_roof))DEALLOCATE(quwinds%uo_roof)
	IF(ALLOCATED(quwinds%vo_roof))DEALLOCATE(quwinds%vo_roof)
	IF(ALLOCATED(qugrid%z))DEALLOCATE(qugrid%z)
	IF(ALLOCATED(qugrid%zm))DEALLOCATE(qugrid%zm)
	IF(ALLOCATED(qugrid%dz_array))DEALLOCATE(qugrid%dz_array)
	IF(ALLOCATED(bld%damage))DEALLOCATE(bld%damage)
	IF(ALLOCATED(landuse%val))DEALLOCATE(landuse%val)
	IF(ALLOCATED(landuse%height))DEALLOCATE(landuse%height)
	IF(ALLOCATED(landuse%atten))DEALLOCATE(landuse%atten)
	IF(ALLOCATED(bld%LrNode))DEALLOCATE(bld%LrNode)
	IF(ALLOCATED(bld%LrFace))DEALLOCATE(bld%LrFace)
	IF(ALLOCATED(bld%FaceRelWindDir))DEALLOCATE(bld%FaceRelWindDir)
	IF(ALLOCATED(quwinds%uint))DEALLOCATE(quwinds%uint)
	IF(ALLOCATED(quwinds%vint))DEALLOCATE(quwinds%vint)
	IF(ALLOCATED(quwinds%undisturbed))DEALLOCATE(quwinds%undisturbed)
	IF(ALLOCATED(ship%speed))DEALLOCATE(ship%speed)
	IF(ALLOCATED(ship%bearing))DEALLOCATE(ship%bearing)
	IF(ALLOCATED(ship%currentSpeed))DEALLOCATE(ship%currentSpeed)
	IF(ALLOCATED(ship%currentDirection))DEALLOCATE(ship%currentDirection)
	IF(ALLOCATED(bld%icellwake))DEALLOCATE(bld%icellwake)
	IF(ALLOCATED(bld%uProfile))DEALLOCATE(bld%uProfile)
	IF(ALLOCATED(bld%vProfile))DEALLOCATE(bld%vProfile)
	IF(ALLOCATED(bld%FaceWakeDir))DEALLOCATE(bld%FaceWakeDir)
	IF(ALLOCATED(bld%speedProfile))DEALLOCATE(bld%speedProfile)

	close(ID_FILE_IGNITE_SEL)
	close(ID_FILE_TRAJ)
	close(ID_FILE_TRAJ_MERGE)

	if(flag%output_initialized .gt. 0)then
		if(flag%uofield.eq.1) close(ID_FILE_QU_UOFIELD)
		if(flag%frm.eq.1 .or. flag%frm.eq.3)then
			close(ID_FILE_QU_CELLTYPE_BIN)
			close(ID_FILE_QU_VELOCITY)
		endif
		if(flag%frm.eq.2 .or. flag%frm.eq.3)then
			close(ID_FILE_QU_VELOCITY_BIN)
			close(ID_FILE_QU_CELLTYPE_BIN)
		endif

		close(ID_FILE_QU_BUILDOUT)
		if(flag%staggered .eq. 1)close(ID_FILE_QU_STAGGERED)
		if(canopy%flag .gt. 0)close(ID_FILE_QU_VEG_BIN)
		close(ID_FILE_QU_UNDIST_BIN)
		if(flag%errorWrite .gt. 0)close(ID_FILE_QU_ERROR)
	endif

	if(flag%isfire == 1) then


		if(fb%flag == 1) close(ID_FILE_FB_OUT)

		! winds
		IF(ALLOCATED(fcawinds%u)) deallocate(fcawinds%u)
		IF(ALLOCATED(fcawinds%v)) deallocate(fcawinds%v)
		IF(ALLOCATED(fcawinds%w)) deallocate(fcawinds%w)
		IF(ALLOCATED(quwinds%uo_before_params)) deallocate(quwinds%uo_before_params)
		IF(ALLOCATED(quwinds%vo_before_params)) deallocate(quwinds%vo_before_params)
		IF(ALLOCATED(quwinds%u_ave)) deallocate(quwinds%u_ave)
		IF(ALLOCATED(quwinds%v_ave)) deallocate(quwinds%v_ave)
		IF(ALLOCATED(quwinds%w_ave)) deallocate(quwinds%w_ave)

		IF(ALLOCATED(mass_int_0)) deallocate(mass_int_0)

		! Winds
		IF(ALLOCATED(fcawinds%sigma)) deallocate(fcawinds%sigma)

		! Fuels
		IF(ALLOCATED(fuels%density)) deallocate(fuels%density)
		IF(ALLOCATED(fuels%density_initial)) deallocate(fuels%density_initial)
		IF(ALLOCATED(fuels%moisture)) deallocate(fuels%moisture)
		IF(ALLOCATED(fuels%depl_center)) deallocate(fuels%depl_center)
		IF(ALLOCATED(fuels%depl_var)) deallocate(fuels%depl_var)		
		IF(ALLOCATED(fuels%moisture_initial)) deallocate(fuels%moisture_initial)
		IF(ALLOCATED(fuels%moist_depl_center)) deallocate(fuels%moist_depl_center)
		IF(ALLOCATED(fuels%moist_depl_var)) deallocate(fuels%moist_depl_var)

		IF(ALLOCATED(fuels%mu_soot)) deallocate(fuels%mu_soot)
		IF(ALLOCATED(fuels%sigma_soot)) deallocate(fuels%sigma_soot)
		IF(ALLOCATED(fuels%conv_human)) deallocate(fuels%conv_human)
		IF(ALLOCATED(fuels%therm_dose)) deallocate(fuels%therm_dose)
		IF(ALLOCATED(fuels%windmag)) deallocate(fuels%windmag)
		IF(ALLOCATED(fuels%o2density)) deallocate(fuels%o2density)

		! Fire
		IF(ALLOCATED(fire%ignitions)) deallocate(fire%ignitions)
		IF(ALLOCATED(fire%reaction_rate)) deallocate(fire%reaction_rate)
		IF(ALLOCATED(fire%energy_to_atmos)) deallocate(fire%energy_to_atmos)
		IF(ALLOCATED(fire%burncenter)) deallocate(fire%burncenter)
		IF(ALLOCATED(fire%spat_var)) deallocate(fire%spat_var)
		IF(ALLOCATED(fire%burncenter_info)) deallocate(fire%burncenter_info)
		IF(ALLOCATED(fire%timedecay)) deallocate(fire%timedecay)
		if(ALLOCATED(fire%energy_used_to_evaporate_water)) deallocate(fire%energy_used_to_evaporate_water)

		! Working arrays
		IF(ALLOCATED(w_rr_o2_turb)) deallocate(w_rr_o2_turb)
		IF(ALLOCATED(w_rr)) deallocate(w_rr)
		IF(ALLOCATED(w_mass)) deallocate(w_mass)
		IF(ALLOCATED(w_moisture)) deallocate(w_moisture)
		IF(ALLOCATED(w_fueldensity)) deallocate(w_fueldensity)
		IF(ALLOCATED(w_q)) deallocate(w_q)
		IF(ALLOCATED(n_new_starts)) deallocate(n_new_starts)
		IF(ALLOCATED(w_depl)) deallocate(w_depl)
		IF(ALLOCATED(w_ignitions)) deallocate(w_ignitions)

		! FB
		IF(ALLOCATED(fb%time_delay)) deallocate(fb%time_delay)
		IF(ALLOCATED(fb%num_ignitions)) deallocate(fb%num_ignitions)


		IF(ALLOCATED(fuels%height_initial)) deallocate(fuels%height_initial)
		IF(ALLOCATED(fuels%actual_height)) deallocate(fuels%actual_height)

		! Other FireCA
		IF(ALLOCATED(firegrid%cell_index)) deallocate(firegrid%cell_index)
		IF(ALLOCATED(firegrid%ijk_cell_index)) deallocate(firegrid%ijk_cell_index)
		IF(ALLOCATED(firegrid%idx)) deallocate(firegrid%idx)
		IF(ALLOCATED(firegrid%z_en2atmos)) deallocate(firegrid%z_en2atmos)
		IF(ALLOCATED(firegrid%zm_en2atmos)) deallocate(firegrid%zm_en2atmos)
		IF(ALLOCATED(firegrid%dz_array_en2atmos)) deallocate(firegrid%dz_array_en2atmos)
		IF(ALLOCATED(firegrid%cellvol_en2atmos)) deallocate(firegrid%cellvol_en2atmos)


		IF(ALLOCATED(grid_conv)) deallocate(grid_conv)
		IF(ALLOCATED(grid_where)) deallocate(grid_where)
		IF(ALLOCATED(init_dens_canopy)) deallocate(init_dens_canopy)

		! Mapping
		IF(ALLOCATED(qugrid%ijk_cell_index)) deallocate(qugrid%ijk_cell_index)

		! Plumes
		IF(ALLOCATED(quwinds%wplume)) deallocate(quwinds%wplume)
		IF(ALLOCATED(quwinds%u_mc)) deallocate(quwinds%u_mc)
		IF(ALLOCATED(quwinds%v_mc)) deallocate(quwinds%v_mc)
		IF(ALLOCATED(quwinds%w_mc)) deallocate(quwinds%w_mc)
		IF(ALLOCATED(plume)) deallocate(plume)
		IF(ALLOCATED(to_delete)) deallocate(to_delete)
		IF(ALLOCATED(qugrid%xcenters)) deallocate(qugrid%xcenters)
		IF(ALLOCATED(qugrid%ycenters)) deallocate(qugrid%ycenters)
		IF(ALLOCATED(qugrid%xedge)) deallocate(qugrid%xedge)
		IF(ALLOCATED(qugrid%yedge)) deallocate(qugrid%yedge)
		IF(ALLOCATED(fca_2_qu_kstart)) deallocate(fca_2_qu_kstart)
		IF(ALLOCATED(fca_2_qu_kend)) deallocate(fca_2_qu_kend)
		IF(ALLOCATED(en2atm_2_qu_kstart)) deallocate(en2atm_2_qu_kstart)
		IF(ALLOCATED(en2atm_2_qu_kend)) deallocate(en2atm_2_qu_kend)

		! Other datamodule
		IF(ALLOCATED(qu_2fca_kuv)) deallocate(qu_2fca_kuv)
		IF(ALLOCATED(qu_2fca_w)) deallocate(qu_2fca_w)
		IF(ALLOCATED(qugrid%dzmi)) deallocate(qugrid%dzmi)
		IF(ALLOCATED(kmap_start)) deallocate(kmap_start)
		IF(ALLOCATED(kmap_end)) deallocate(kmap_end)

		! Fire grid
		IF(ALLOCATED(firegrid%dz_array)) deallocate(firegrid%dz_array)
		IF(ALLOCATED(firegrid%z)) deallocate(firegrid%z)
		IF(ALLOCATED(firegrid%zm)) deallocate(firegrid%zm)
		IF(ALLOCATED(firegrid%cellvol)) deallocate(firegrid%cellvol)
		IF(ALLOCATED(firegrid%dzmi)) deallocate(firegrid%dzmi)
		IF(ALLOCATED(firegrid%xcenters)) deallocate(firegrid%xcenters)
		IF(ALLOCATED(firegrid%ycenters)) deallocate(firegrid%ycenters)

		! Ignitions
		IF(ALLOCATED(fire_ignition%x)) deallocate(fire_ignition%x)
		IF(ALLOCATED(fire_ignition%y))	deallocate(fire_ignition%y)
		IF(ALLOCATED(fire_ignition%z))	deallocate(fire_ignition%z)
		IF(ALLOCATED(fire_ignition%time))	deallocate(fire_ignition%time)
		IF(ALLOCATED(fire_ignition%radius))	deallocate(fire_ignition%radius)
		IF(ALLOCATED(fire_ignition%new_num))	deallocate(fire_ignition%new_num)

		! Conversion between grids
		if(allocated(fca2quic))deallocate(fca2quic)
		if(allocated(ft_2_fca_conv_x))then
			do i = 1, firegrid%nx
				if(ft_2_fca_conv_x(i)%nelem > 0) then
					deallocate(ft_2_fca_conv_x(i)%index)
					deallocate(ft_2_fca_conv_x(i)%length)
				endif
			enddo
			deallocate(ft_2_fca_conv_x)
	
			do i = 1, firegrid%ny
				if(ft_2_fca_conv_y(i)%nelem > 0) then
					deallocate(ft_2_fca_conv_y(i)%index)
					deallocate(ft_2_fca_conv_y(i)%length)
				endif
			enddo
			deallocate(ft_2_fca_conv_y)
	
			do i = 1, firegrid%nz
				if(ft_2_fca_conv_z(i)%nelem > 0) then
					deallocate(ft_2_fca_conv_z(i)%index)
					deallocate(ft_2_fca_conv_z(i)%length)
				endif
			enddo
			deallocate(ft_2_fca_conv_z)
		endif
		! Random numbers
		IF(ALLOCATED(stat)) deallocate(stat)

	endif

	END SUBROUTINE CLEANUP