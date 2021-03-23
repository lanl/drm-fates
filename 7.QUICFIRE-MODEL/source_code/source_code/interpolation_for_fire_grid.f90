subroutine interpolate_for_fire_grid !(fire_time)
    
	! Trilinear interpolation of winds from the QU grid to the QFire grid (same or higher resolution than the QU grid)
	! The quic fire quantities are cell centered
	
	use constants
   
	use fireca_module
	use interpolation_module
	use grid_module
	use winds_module
    
   implicit none
	
	real,dimension(qugrid%nx,qugrid%ny) :: 				&
		vert_interpolation_u,									&
		vert_interpolation_v,									&
		vert_interpolation_w
	real,dimension(firegrid%nx,firegrid%ny) :: temp_u,temp_v,temp_w
	integer :: i,j,k

	do k = 1,firegrid%nz+1
		
		call VerticalInterpolation(qu_2fca_kuv(k),qu_2fca_w(k),k,			&
			vert_interpolation_u,vert_interpolation_v,vert_interpolation_w)
		
		call XYInterp(vert_interpolation_u,vert_interpolation_v,vert_interpolation_w,		&
			temp_u,temp_v,temp_w,jstart_interp,jend_interp,istart_interp,iend_interp)
		
		!$OMP parallel do private(i,j)
		do j = 1,firegrid%ny
			do i = 1,firegrid%nx
				fcawinds%u(i,j,k) = temp_u(i,j)
				fcawinds%v(i,j,k) = temp_v(i,j)
				fcawinds%w(i,j,k) = temp_w(i,j)
			enddo
		enddo
		!$OMP end parallel do 
		
	enddo

	end subroutine interpolate_for_fire_grid                             

!=====================================================================================================================
!=====================================================================================================================		
	subroutine VerticalInterpolation(kuv,kw,kfire, &
		vert_interpolation_u,vert_interpolation_v,vert_interpolation_w)

	use grid_module
	use winds_module
	
	implicit none
	
	integer,intent(IN) :: kuv,kw,kfire	
	real,dimension(qugrid%nx,qugrid%ny) :: vert_interpolation_u,vert_interpolation_v,vert_interpolation_w
	integer :: i,j

	
	! Initialization 
	!$OMP parallel do private(i,j)
	do j = 1,qugrid%ny
		do i = 1,qugrid%nx
			vert_interpolation_u(i,j) = 0.
			vert_interpolation_v(i,j) = 0.
			vert_interpolation_w(i,j) = 0.
		enddo
	enddo
	!$OMP end parallel do
	
	! Interpolate	
	if(kuv > 2)then
		!$OMP parallel do private(i,j)
		do j = 1,qugrid%ny
			do i = 1,qugrid%nx
				
				call linear_interpolation(qugrid%zm(kuv-1),qugrid%zm(kuv), &
					quwinds%u(i,j,kuv-1),quwinds%u(i,j,kuv),firegrid%zm(kfire), vert_interpolation_u(i,j))
				call linear_interpolation(qugrid%zm(kuv-1),qugrid%zm(kuv), &
					quwinds%v(i,j,kuv-1),quwinds%v(i,j,kuv),firegrid%zm(kfire), vert_interpolation_v(i,j))
			enddo
		enddo
		!$OMP end parallel do
	else
		!$OMP parallel do private(i,j)
		do j = 1,qugrid%ny
			do i = 1,qugrid%nx
				call linear_interpolation(0.,qugrid%zm(kuv),0., &
					quwinds%u(i,j,kuv),firegrid%zm(kfire),vert_interpolation_u(i,j))
				call linear_interpolation(0.,qugrid%zm(kuv),0., &
					quwinds%v(i,j,kuv),firegrid%zm(kfire),vert_interpolation_v(i,j))
			enddo
		enddo
		!$OMP end parallel do
	endif
		
	if(kw > 3)then
		!$OMP parallel do private(i,j)
		do j = 1,qugrid%ny
			do i = 1,qugrid%nx
				call linear_interpolation(qugrid%z(kw-1),qugrid%z(kw),	&
					quwinds%w(i,j,kw-1),quwinds%w(i,j,kw),firegrid%zm(kfire),vert_interpolation_w(i,j))
			enddo
		enddo
		!$OMP end parallel do
	else
		!$OMP parallel do private(i,j)
		do j = 1,qugrid%ny
			do i = 1,qugrid%nx
				call linear_interpolation(0.,qugrid%z(kw),0.,quwinds%w(i,j,kw),	&
				firegrid%zm(kfire),vert_interpolation_w(i,j))
			enddo
		enddo
		!$OMP end parallel do
	endif	
	
	return
	
	end	
!=====================================================================================================================
!=====================================================================================================================		
	subroutine XYInterp(vert_interpolation_u,vert_interpolation_v,vert_interpolation_w,	&
		weather_u,weather_v,weather_w,jstart,jend,istart,iend)
	
	use grid_module
	
	implicit none
	
	integer,intent(IN) :: jstart,jend,istart,iend
	real,intent(IN),dimension(qugrid%nx,qugrid%ny) :: vert_interpolation_u,vert_interpolation_v,vert_interpolation_w
	real,dimension(firegrid%nx,firegrid%ny) :: weather_u,weather_v,weather_w
	integer :: i,j,k,i_qu,j_qu,i_qu_0,j_qu_0,i1,i2,j1,j2
	real :: x1,x2,x,y1,y2,y,x_qu,y_qu,temp_1,temp_2
	integer,dimension(2) :: low_extr,up_extr,j_extr
		
	! if ratio = 4, then ceil(ratio*0.5) = 2 but this cell center is below the center of the qu cell so add 1
	! if ratio = 5, then ceil(ratio*0.5) = 3 this cell center is at the center of the qu cell so add 0
	
	!!$OMP parallel do private(y,y_qu,j_qu,j_qu_0,x,x_qu,i_qu,i_qu_0,i,j,i1,i2,x1,x2,j1,j2,y1,y2,temp_1,temp_2)
	do j = jstart,jend
 
		y = firegrid%ycenters(j)			! center of the fire grid cell in the y direction
		j_qu = ceiling(y * qugrid%dyi)	! QU cell in the y direction
		y_qu = qugrid%ycenters(j_qu)			! center of the QU grid cell in the y direction
		
		j_qu_0 = j_qu
		if(y_qu > y) then
			j_qu = j_qu-1  ! Interpolate with the cell below
		endif
 
		do i = istart,iend
			x = firegrid%xcenters(i)			! center of the fire grid cell in the x direction
			i_qu = ceiling(x * qugrid%dxi)	! QU cell in the x direction
			x_qu = qugrid%xcenters(i_qu)				! center of the QU grid cell in the x direction
		
			i_qu_0 = i_qu
			if(x_qu > x) then
				i_qu = i_qu-1  ! Interpolate with the cell below
			endif
 
			! -- u
			! interpolation in the x-direction at two different y-levels
			i1 = i_qu_0					! side of the QU cell before the center of the QF cell
			x1 = qugrid%xedge(i1)
			i2 = i_qu_0 + 1			! side of the QU cell after  the center of the QF cell			
			x2 = qugrid%xedge(i2)!x1 + dx
			j1 = j_qu					! center of the QU cell before the center of the QF cell
			y1 = qugrid%ycenters(j1)
			j2 = j_qu + 1				! center of the QU cell after  the center of the QF cell			
			y2 = qugrid%ycenters(j2) !y1 + dy
			call linear_interpolation( x1 , x2 , vert_interpolation_u(i1,j1),vert_interpolation_u(i2,j1),x, temp_1)
			call linear_interpolation( x1 , x2 , vert_interpolation_u(i1,j2),vert_interpolation_u(i2,j2),x, temp_2)
			! interpolation in the y-direction
			call linear_interpolation( y1 , y2 , temp_1 , temp_2 , y, weather_u(i,j))
 
			! -- v
			! interpolation in the x-direction at two different y-levels
			i1 = i_qu					! center of the QU cell before the center of the QF cell
			x1 = qugrid%xcenters(i1)
			i2 = i_qu + 1				! center of the QU cell after  the center of the QF cell
			x2 = qugrid%xcenters(i2) !x1 + dx
			j1 = j_qu_0					! side of the QU cell before the center of the QF cell
			y1 = qugrid%yedge(j1)
			j2 = j_qu_0	+ 1			! side of the QU cell after  the center of the QF cell
			y2 = qugrid%yedge(j2) !y1 + dy
			call linear_interpolation( x1 , x2 , vert_interpolation_v(i1,j1),vert_interpolation_v(i2,j1),x, temp_1)
			call linear_interpolation( x1 , x2 , vert_interpolation_v(i1,j2),vert_interpolation_v(i2,j2),x, temp_2)
			! interpolation in the y-direction
			call linear_interpolation( y1 , y2 , temp_1, temp_2 , y, weather_v(i,j))
 
			! -- w
			! interpolation in the x-direction at two different y-levels
			i1 = i_qu				! center of the QU cell before the center of the QF cell
			x1 = qugrid%xcenters(i1)
			i2 = i_qu + 1			! center of the QU cell after  the center of the QF cell			
			x2 = qugrid%xcenters(i2)!x1 + dx
			j1 = j_qu				! center of the QU cell before the center of the QF cell
			y1 = qugrid%ycenters(j1)
			j2 = j_qu + 1		! center of the QU cell after  the center of the QF cell
			y2 = qugrid%ycenters(j2) !y1 + dy
			call linear_interpolation( x1 , x2 , vert_interpolation_w(i1,j1),vert_interpolation_w(i2,j1),x,temp_1)
			call linear_interpolation( x1 , x2 , vert_interpolation_w(i1,j2),vert_interpolation_w(i2,j2),x,temp_2)
			! interpolation in the y-direction
			call linear_interpolation( y1 , y2 , temp_1, temp_2 , y,weather_w(i,j))
 
		enddo	
	enddo
	!!$OMP end parallel do
	
	! boundaries (top & bottom)
	low_extr = (/1,jend+1/)
	up_extr = (/jstart-1,firegrid%ny/)
	j_extr = (/1,qugrid%ny-1/)
	!!$OMP parallel do private(i,x,i_qu,x_qu,i_qu_0,k,j_qu,j_qu_0,j,y,i1,x1,i2,x2,j1,y1,y2,temp_1,temp_2)
	do i = 1,firegrid%nx
		x = firegrid%xcenters(i)		! center of the fire grid cell in the x direction
		i_qu = ceiling(x * qugrid%dxi)	! QU cell in they direction
		x_qu = qugrid%xcenters(i_qu)			! center of the QU grid cell in the x direction
		
		i_qu_0 = i_qu
		if(x_qu > x) then
			i_qu = max(i_qu-1,1)  ! Interpolate with the cell below
		endif
 
		do k = 1,2
			j_qu = j_extr(k)
			j_qu_0 = j_qu
			do j = low_extr(k),up_extr(k)
				y = firegrid%ycenters(j)	! center of the fire grid cell in the x direction
				
				! -- u
				i1 = i_qu_0			! side of the QU cell before the center of the QF cell
				x1 = qugrid%xedge(i1)
				i2 = i_qu_0 + 1	! side of the QU cell after  the center of the QF cell
				x2 = qugrid%xedge(i2) !x1 + dx
				j1 = j_qu			! center of the QU cell before the center of the QF cell
				call linear_interpolation( x1 , x2 , vert_interpolation_u(i1,j1),vert_interpolation_u(i2,j1),x,weather_u(i,j) )
			
				! -- v
				! interpolation in the x-direction at two different y-levels
				i1 = i_qu			! center of the QU cell before the center of the QF cell
				x1 = qugrid%xcenters(i1)
				i2 = i_qu + 1		! center of the QU cell after  the center of the QF cell
				x2 = qugrid%xcenters(i2) !x1 + dx
				j1 = j_qu_0			! side of the QU cell before the center of the QF cell
				y1 = qugrid%yedge(j1)
				j2 = j_qu_0 + 1	! side of the QU cell after  the center of the QF cell
				y2 = qugrid%yedge(j2) !y1 + dy
				call linear_interpolation( x1 , x2 , vert_interpolation_v(i1,j1),vert_interpolation_v(i2,j1),x,temp_1)
				call linear_interpolation( x1 , x2 , vert_interpolation_v(i1,j2),vert_interpolation_v(i2,j2),x,temp_2)
				! interpolation in the y-direction
				call linear_interpolation( y1 , y2 , temp_1, temp_2 , y, weather_v(i,j))
			
				! -- w
				i1 = i_qu			! center of the QU cell before the center of the QF cell
				x1 = qugrid%xcenters(i1)
				i2 = i_qu + 1		! center of the QU cell after  the center of the QF cell
				x2 = qugrid%xcenters(i2) !x1 + dx
				j1 = j_qu			! center of the QU cell before the center of the QF cell
				call linear_interpolation( x1 , x2 , vert_interpolation_w(i1,j1),vert_interpolation_w(i2,j1),x,weather_w(i,j))
			enddo
		enddo
	enddo
	!!$OMP end parallel do 
	
	! boundaries (right & left)
	low_extr = (/1,iend+1/)
	up_extr = (/istart-1,firegrid%nx/)
	j_extr = (/1,qugrid%nx-1/)	
	!!$OMP parallel do private(j,y,j_qu,y_qu,j_qu_0,k,i_qu,i_qu_0,i,x,i1,i2,x1,x2,y1,y2,temp_1,temp_2)
	do j = 2,firegrid%ny-1  ! 1 and firegrid%ny already done in the previous pass
		y = firegrid%ycenters(j)	! center of the fire grid cell in the y direction
		j_qu = ceiling(y * qugrid%dyi)					! QU cell in the y direction
		y_qu = qugrid%ycenters(j_qu)			! center of the QU grid cell in the y direction
		
		j_qu_0 = j_qu
		if(y_qu > y) then
			j_qu = max(j_qu-1,1)  ! Interpolate with the cell below
		endif
	
		do k = 1,2
			i_qu = j_extr(k)
			i_qu_0 = i_qu
			do i = low_extr(k),up_extr(k)
				x = firegrid%xcenters(i)	! center of the fire grid cell in the x direction
				
				! -- u
				! interpolation in the x-direction at two different y-levels
				i1 = i_qu_0					! side of the QU cell before the center of the QF cell
				x1 = qugrid%xedge(i1)
				i2 = i_qu_0	+ 1			! side of the QU cell after  the center of the QF cell				
				x2 = qugrid%xedge(i2) !x1 + dx
				j1 = j_qu					! center of the QU cell before the center of the QF cell
				y1 = qugrid%ycenters(j1)
				j2 = j_qu + 1				! center of the QU cell after  the center of the QF cell
				y2 = qugrid%ycenters(j2) !y1 + dy
				call linear_interpolation( x1 , x2 , vert_interpolation_u(i1,j1),vert_interpolation_u(i2,j1),x,temp_1)
				call linear_interpolation( x1 , x2 , vert_interpolation_u(i1,j2),vert_interpolation_u(i2,j2),x,temp_2)
				! interpolation in the y-direction
				call linear_interpolation( y1 , y2 , temp_1 , temp_2 , y, weather_u(i,j))
			
				! -- v
				i1 = i_qu			! center of the QU cell before the center of the QF cell
				j1 = j_qu_0			! side of the QU cell before the center of the QF cell
				y1 = qugrid%yedge(j1)
				j2 = j_qu_0	+ 1	! side of the QU cell after  the center of the QF cell				
				y2 = qugrid%yedge(j2) !y1 + dy
				call linear_interpolation( y1 , y2 , vert_interpolation_v(i1,j1),vert_interpolation_v(i1,j2),y,weather_v(i,j))
							
				! -- w
				i1 = i_qu					! center of the QU cell before the center of the QF cell
				j1 = j_qu	! center of the QU cell before the center of the QF cell
				y1 = qugrid%ycenters(j1)
				j2 = j_qu + 1	! center of the QU cell after  the center of the QF cell
				y2 = qugrid%ycenters(j2) !y1 + dy
				call linear_interpolation( y1 , y2 , vert_interpolation_w(i1,j1),vert_interpolation_w(i1,j2),y,weather_w(i,j))
			enddo
		enddo
	enddo
	!!$OMP end parallel do 
	

	end
!=====================================================================================================================
!=====================================================================================================================		
	subroutine linear_interpolation(x1,x2,y1,y2,xinterpolated,outvar) 
	! used for interpolation of wind data from the QUIC grid to the fire grid 
	! in the interpolation_for-fire_frid.f90 routine which is called in main.f90
    
   implicit none   
         
   real, intent(in) :: x1,x2,y1,y2,xinterpolated
	real, intent(out) :: outvar

	outvar = (xinterpolated - x1)*(y2 - y1)/(x2 - x1) + y1	
	
	end 
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE AverageFiretechFuelDepth(field_fca,	field_ft)

	! fca = FireCA
	! ft = FireTech
		
	use grid_module
	use interpolation_module

	implicit none
	
	real,dimension(ft%nx,ft%ny), intent(IN) :: field_ft
	real,dimension(firegrid%nx,firegrid%ny) :: field_fca
	real :: ntot
	integer :: i,j,ii,jj

	! Initialization
	field_fca = 0.
		
	do j = 2,firegrid%ny-1		
		do i = 2,firegrid%nx-1
			ntot = 0
			do jj = ft_2_fca_conv_y(j)%index(1), ft_2_fca_conv_y(j)%index(ft_2_fca_conv_y(j)%nelem)
				do ii = ft_2_fca_conv_x(i)%index(1), ft_2_fca_conv_x(i)%index(ft_2_fca_conv_x(i)%nelem)			
					ntot = ntot + 1				
					field_fca(i,j) = field_fca(i,j) + field_ft(ii, jj)
				enddo
			enddo
			field_fca(i,j) = min(field_fca(i,j) / ntot, firegrid%dz_array(1))
		enddo
	enddo
	
	end 
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE InterpolateFiretechFuelDepth(field_fca,	field_ft)

	! fca = FireCA
	! ft = FireTech
		
	use grid_module

	implicit none

	real,dimension(ft%nx,ft%ny),intent(IN) :: field_ft
	real,dimension(firegrid%nx,firegrid%ny) :: field_fca
		
	call XYInterpFiretech(field_ft,field_fca)			
	
	end	
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE DefineInterpolationArrays()
	
	use interpolation_module
	use fireca_module
	use grid_module
	
	implicit none
	
	integer ::														&
		i,j,k,														& ! N/A, loop indexes
		k_start,k_end												  ! N/A, indexes of the first and last cell in the Firetec domain inside a FireCA domain
		
	real ::															&
		ftz,															& ! m, z for firetec, either the cell top or the fuel top for the first layer
		z1,z2,														& ! m, part of a Firetec cell inside a FireCA cell
		dutmx,dutmy,												& ! m, Fireca south-west corner - Firetec south-west corner 
		xs,xe,														& ! m, start and end of location in the Firetec domain (x-dir)
		xsfca,xefca,												& ! m, left border of the cell in the FireCA domain, relative coordinates (from SW corner)
		ys,ye,														& ! m, start and end of location in the Firetec domain (y-dir)
		ysfca,yefca												    ! m, bottom border of the cell in the FireCA domain, relative coordinates (from SW corner)
		
	allocate(														&
		ft_2_fca_conv_x(firegrid%nx),							&
		ft_2_fca_conv_y(firegrid%ny),							&
		ft_2_fca_conv_z(firegrid%nz))

	! Initialization	
	dutmx = firegrid%utmx - ft%utmx
	dutmy = firegrid%utmy - ft%utmy
	
	!!! X_DIRECTION
	do i = 2,firegrid%nx-1
		! Left border of the cell in the FireCA domain, relative coordinates (from SW corner)
		xsfca = real(i-1)*firegrid%dx
		! Location in the Firetec domain, relative coordinates (from SW corner)
		xs = xsfca + dutmx
		! Cell in the Firetec domain
		k_start = max(floor(xs * ft%dxi), 1)
		if(mod(xs, ft%dx) == 0) k_start = k_start + 1
		
		! Right border of the cell in the FireCA domain, relative coordinates (from SW corner)
		xefca = real(i)*firegrid%dx
		! Location in the Firetec domain, relative coordinates (from SW corner)
		xe = xefca + dutmx
		! Cell in the Firetec domain
		k_end = min(ceiling(xe * ft%dxi),ft%nx)
		
		! How many Firetec cells per FireCA cells
		ft_2_fca_conv_x(i)%nelem = k_end - k_start + 1
		allocate(																&
			ft_2_fca_conv_x(i)%index(ft_2_fca_conv_x(i)%nelem),	&
			ft_2_fca_conv_x(i)%length(ft_2_fca_conv_x(i)%nelem))
		
		do k = 1,ft_2_fca_conv_x(i)%nelem
			! Firetec cell index
			j = k_start + k - 1
			
			ft_2_fca_conv_x(i)%index(k) = j
			
			! Coordinates
			z1 = max(xs, real(j-1)*ft%dx)
			z2 = min(xe, real(j)*ft%dx)
			
			! Part of the Firetec cell in the FireCA cell
			ft_2_fca_conv_x(i)%length(k) = (z2-z1)
		enddo
	enddo
	
	!!! Y-DIRECTION
	do i = 2,firegrid%ny-1
		! Left border of the cell in the FireCA domain, relative coordinates (from SW corner)
		ysfca = real(i-1)*firegrid%dy
		! Location in the Firetec domain, relative coordinates (from SW corner)
		ys = ysfca + dutmy		
		! Cell in the Firetec domain
		k_start = max(floor(ys * ft%dyi), 1)
		if(mod(ys, ft%dy) == 0) k_start = k_start + 1
		
		! Right border of the cell in the FireCA domain, relative coordinates (from SW corner)
		yefca = real(i)*firegrid%dy
		! Location in the Firetec domain, relative coordinates (from SW corner)
		ye = yefca + dutmy
		! Cell in the Firetec domain
		k_end = min(ceiling(ye * ft%dyi),ft%ny)
		
		! How many Firetec cells per FireCA cells
		ft_2_fca_conv_y(i)%nelem = k_end - k_start + 1
		allocate(																&
			ft_2_fca_conv_y(i)%index(ft_2_fca_conv_y(i)%nelem),	&
			ft_2_fca_conv_y(i)%length(ft_2_fca_conv_y(i)%nelem))
		
		do k = 1,ft_2_fca_conv_y(i)%nelem
			! Firetec cell index
			j = k_start + k - 1
			
			ft_2_fca_conv_y(i)%index(k) = j
			
			! Coordinates
			z1 = max(ys, real(j-1)*ft%dy)
			z2 = min(ye, real(j)*ft%dy)
			
			! Part of the Firetec cell in the FireCA cell
			ft_2_fca_conv_y(i)%length(k) = (z2-z1)
		enddo
	enddo	
	
	!!! Z-DIRECTION
	do k = 1,firegrid%nz
		! cells in the firetech domain
		k_start = 1
		do while(k_start < ft%nz .and. ft%z(k_start) <= firegrid%z(k))
			k_start = k_start + 1
		enddo
		k_start = max(1,k_start-1)
		
		k_end = k_start
		do while(k_end < ft%nz .and. ft%z(k_end) < firegrid%z(k+1))
			k_end = k_end + 1
		enddo
		k_end = max(k_end-1,1)
		
		ft_2_fca_conv_z(k)%nelem = k_end - k_start + 1
		allocate(																&
			ft_2_fca_conv_z(k)%index(ft_2_fca_conv_z(k)%nelem),	&
			ft_2_fca_conv_z(k)%length(ft_2_fca_conv_z(k)%nelem))
		
		do i = 1,ft_2_fca_conv_z(k)%nelem
			j = k_start + i - 1
			ft_2_fca_conv_z(k)%index(i) = j
			
			if(j == 1)then
				z1 = max(0. , firegrid%z(k))
			else
				z1 = max(ft%z(j) , firegrid%z(k))
			endif
			z2 = min(ft%z(j+1) , firegrid%z(k+1))
			ft_2_fca_conv_z(k)%length(i) = z2 - z1
		enddo
	enddo
	
	END
!=====================================================================================================================
!=====================================================================================================================	
	SUBROUTINE InterpolateFiretechFile(field_fca, field_ft)

	! fca = FireCA
	! ft = FireTech
		
	use fireca_module
	use grid_module
	use interpolation_module
	
	implicit none
	
	real,dimension(ft%nx,ft%ny,ft%nz),intent(IN) :: field_ft	! field to interpolate from Firetec grid to FireCA grid		
	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz),intent(OUT) :: field_fca	! FireCA interpolated field
	integer ::														&
		i,j,k,														& ! N/A, loop indexes
		ii,jj,kk														  ! N/A, loop indexes
	
	real :: mass,mass1,													& ! depends on the field interpolated
		frz															  ! m, fraction of cell with mass	
	
	! Initialization
	field_fca = 0.
	
	!$OMP parallel do private(i,j,k,mass,kk,ii,jj,frz)
	do k = 1,firegrid%nz
		do j = 2,firegrid%ny-1
			do i = 2,firegrid%nx-1
				mass = 0.				
				do kk = ft_2_fca_conv_z(k)%index(1), ft_2_fca_conv_z(k)%index(ft_2_fca_conv_z(k)%nelem)
					do jj = ft_2_fca_conv_y(j)%index(1), ft_2_fca_conv_y(j)%index(ft_2_fca_conv_y(j)%nelem)
						do ii = ft_2_fca_conv_x(i)%index(1), ft_2_fca_conv_x(i)%index(ft_2_fca_conv_x(i)%nelem)
							if(kk == 1 .and. fuels%actual_height(i, j) < ft%dz_array(1)) then
								if(k == 1) then
									frz = firegrid%z(k+1) / max(fuels%actual_height(i, j),firegrid%z(k+1)) * ft%dz_array(kk)
								else
									frz = max((min(fuels%actual_height(i, j),firegrid%z(k+1)) - firegrid%z(k)), 0.) / &
										fuels%actual_height(i, j) * ft%dz_array(kk)
								endif
							else
								frz = ft_2_fca_conv_z(k)%length(kk-ft_2_fca_conv_z(k)%index(1)+1)
							endif
							mass = mass + field_ft(ii,jj,kk) *											&
								ft_2_fca_conv_x(i)%length(ii-ft_2_fca_conv_x(i)%index(1)+1)*	&
								ft_2_fca_conv_y(j)%length(jj-ft_2_fca_conv_y(j)%index(1)+1)*	&
								frz
								!ft_2_fca_conv_z(k)%length(kk-ft_2_fca_conv_z(k)%index(1)+1)
						enddo
					enddo
				enddo
				field_fca(i,j,k) = mass / firegrid%cellvol(k)
			enddo
		enddo
	enddo
	!$OMP end parallel do
	
	!open(111, file='temp.bin', form='unformatted')
	!mass = 0.
	!do k = 1,firegrid%nz
	!	mass = mass + sum(field_fca(:,:,k)) * firegrid%cellvol(k)
	!enddo
	!print*,'FCA mass=',mass
	!write(111) firegrid%nx,firegrid%ny,firegrid%nz
	!write(111) firegrid%cellvol(1:firegrid%nz)
	!write(111) field_fca
	!write(111) ft%nx,ft%ny,ft%nz
	!write(111) ft%dz_array * ft%dx * ft%dy
	!write(111) field_ft
	!close(111)
	!mass1 = 0.
	!do k = 1, ft%nz
	!	mass1 = mass1 + sum(field_ft(:,:,k)) * ft%dz_array(k)
	!enddo	
	!mass1 = mass1 * ft%dx * ft%dy
	!print*,'FT mass=',mass1
	!print*,'diff=',mass1-mass
	!close(111)
	!read(*,*)
	!	
	!xs = dutmx
	!if(mod(xs,dx_ft) == 0.)then
	!	k_start = int(xs/dx_ft)+1
	!else
	!	k_start = int(xs/dx_ft)
	!endif
	!xefca = real(firegrid%nx)*firegrid%dx
	!xe = xefca+dutmx
	!if(mod(xe,dx_ft) == 0.)then
	!	k_end = int(xe/dx_ft)
	!else
	!	k_end = int(xe/dx_ft)+1
	!endif
 !
	!xs = dutmy
	!if(mod(xs,dy_ft) == 0.)then
	!	ii = int(xs/dy_ft)+1
	!else
	!	ii = int(xs/dy_ft)
	!endif
	!xefca = real(firegrid%ny)*firegrid%dy
	!xe = xefca+dutmy
	!if(mod(xe,dy_ft) == 0.)then
	!	jj = int(xe/dy_ft)
	!else
	!	jj = int(xe/dy_ft)+1
	!endif
	!	
	!z1 = 0.
	!do k = 1,nz_ft
	!	z1 = z1 + sum(field_ft(k_start:k_end,ii:jj,k)) * dx_ft*dy_ft*(z_ft_top(k+1)-z_ft_top(k))
	!enddo
	!print*,'FT  mass=',z1
	!print*,'%=',(z1-mass)/z1*100.
	!print*,k_start,k_end,ii,jj
	!read(*,*)
	!do k = 1,firegrid%nz
	!	do j = 1,firegrid%ny
	!		do i = 1,firegrid%nx
	!			if(field_fca(i,j,k)/=field_ft(i,j,k))then
	!				print*,field_fca(i,j,k),field_ft(i,j,k),field_fca(i,j,k)-field_ft(i,j,k)
	!			endif
	!		enddo
	!	enddo
	!enddo
	
	!call NearestNeighbourInterpolation(field_fca,x_sw_fca,y_sw_fca, &
	!		field_ft,nx_ft,ny_ft,nz_ft,dx_ft,dy_ft,z_ft,x_sw_ft,y_sw_ft)
	
	!do k = 1,firegrid%nz
	!		
	!	if(z_ft_top(nz_ft+1) > firegrid%z(k))then
 !
	!		if(nz_ft > 1) then
	!			call VerticalInterpolationFiretech(field_ft,nx_ft,ny_ft,nz_ft,z_ft, &
	!				vert_interpolation_ft,k)
	!		else
	!			vert_interpolation_ft = field_ft(:,:,1)
	!		endif
	!	
	!		call XYInterpFiretech(x_sw_fca,y_sw_fca, &
	!			nx_ft,ny_ft,dx_ft,dy_ft,x_sw_ft,y_sw_ft,vert_interpolation_ft,outfield)
 !
	!		!$OMP parallel do private(i,j)
	!		do j = 1,firegrid%ny
	!			do i = 1,firegrid%nx
	!				field_fca(i,j,k) = outfield(i,j)
	!			enddo
	!		enddo
	!		!$OMP end parallel do
	!	endif
	!enddo
	
	end	
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE VerticalInterpolationFiretech(field_ft,nx_ft,ny_ft,nz_ft,z_ft, &
			vert_interpolation_ft,kfire)
	
	use grid_module
	
	implicit none
	
	integer,intent(IN) :: nx_ft,ny_ft,nz_ft,kfire	
	real,dimension(nz_ft),intent(IN) :: z_ft	
	real,dimension(nx_ft,ny_ft) :: vert_interpolation_ft
	real,dimension(nx_ft,ny_ft,nz_ft) :: field_ft
	integer :: found,k,i,j
	
	! Initialization 
	!$OMP parallel do private(i,j)
	do j = 1,ny_ft
		do i = 1,nx_ft
			vert_interpolation_ft(i,j) = 0.
		enddo
	enddo
	!$OMP end parallel do
	
	! Linear interpolation between layers	
	found = 0
	k = 0	
	do while(found == 0 .and. k < nz_ft)
		k = k + 1
		if(z_ft(k) >= firegrid%zm(kfire)) then
			found = 1
		endif
	enddo
	if(found == 0)then
		goto 10
	endif
		
	! Interpolate
	if(k >= 2)then
		!$OMP parallel do private(i,j)
		do j = 1,ny_ft
			do i = 1,nx_ft
				call linear_interpolation(z_ft(k-1),z_ft(k),field_ft(i,j,k-1), &
					field_ft(i,j,k),firegrid%zm(kfire),vert_interpolation_ft(i,j))				
			enddo
		enddo
		!$OMP end parallel do
	else
		!$OMP parallel do private(i,j)
		do j = 1,ny_ft
			do i = 1,nx_ft
				vert_interpolation_ft(i,j) = field_ft(i,j,k)
			enddo
		enddo
		!$OMP end parallel do
	endif
		
	return
	
10 print*,'Error in the interpolation of the wind field on the fire grid.'
	print*,'Press any key to terminate.'
	read(*,*)
	stop
	
	END
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE XYInterpFiretech(vert_interpolation_ft,outfield)
	
   use grid_module
	
	implicit none
	
	real,intent(IN),dimension(ft%nx,ft%ny) :: vert_interpolation_ft
	real,dimension(firegrid%nx,firegrid%ny) :: outfield
	real :: x,x1,x2,y,y1,y2,temp_1,temp_2,dutmx,dutmy
	integer :: i,j,jft1,jft2,ift1,ift2,borderx,bordery
	
	outfield = 0.
	
	dutmx = firegrid%utmx - ft%utmx
	dutmy = firegrid%utmy - ft%utmy
	
	do j = 1,firegrid%ny
		! fireca coordinate
		y = (real(j)-0.5)*firegrid%dy
		! firetech cell
		jft1 = ceiling((y+dutmy)/ft%dy)
		! firetech coordinate
		y1 = (real(jft1)-0.5)*ft%dy
		
		bordery = 0
		if(y+dutmy >= y1)then
			if(jft1 == ft%ny)then
				bordery = 1
				jft2 = jft1
			else
				y2 = (real(jft1)+0.5)*ft%dy
				jft2 = jft1 + 1
			endif
		else
			if(jft1 == 1)then
				bordery = 1
				jft2 = jft1
			else
				y2 = (real(jft1)-1.5)*ft%dy
				jft2 = jft1 - 1
			endif
		endif
		
		do i = 1,firegrid%nx
			x = (real(i)-0.5)*firegrid%dx
			ift1 = ceiling((x+dutmx)/ft%dx)
			x1 = (real(ift1)-0.5)*ft%dx
			borderx = 0
			if(x+dutmx >= x1)then
				if(ift1 == ft%nx)then
					borderx = 1
					ift2 = ift1
				else
					x2 = (real(ift1)+0.5)*ft%dx
					ift2 = ift1 + 1
				endif
			else
				if(ift1 == 1)then
					borderx = 1
					ift2 = ift1
				else
					x2 = (real(ift1)-1.5)*ft%dx
					ift2 = ift1 - 1
				endif
			endif
			! Interpolation in the x-direction at two different y levels
			if(borderx == 1) then
				temp_1 = vert_interpolation_ft(ift1,jft1)
				temp_2 = vert_interpolation_ft(ift1,jft2)
			else
				call linear_interpolation(x1,x2,vert_interpolation_ft(ift1,jft1),vert_interpolation_ft(ift2,jft1),x+dutmx,temp_1)
				call linear_interpolation(x1,x2,vert_interpolation_ft(ift1,jft2),vert_interpolation_ft(ift2,jft2),x+dutmx,temp_2)
			endif
			if(bordery == 1) then
				if(jft1 == 1)then
					outfield(i,j) = temp_1
				else
					outfield(i,j) = temp_2
				endif
			else
				call linear_interpolation(y1,y2,temp_1,temp_2,y+dutmy,outfield(i,j))
			endif
		enddo
	enddo
	
	END
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE NearestNeighbourInterpolation(field_fca,x_sw_fca,y_sw_fca, &
			field_ft,nx_ft,ny_ft,nz_ft,dx_ft,dy_ft,z_ft,x_sw_ft,y_sw_ft)

	! fca = FireCA
	! ft = FireTech
		
	use grid_module

	implicit none
	
	integer,intent(IN) :: nx_ft,ny_ft,nz_ft
	real,intent(IN) :: dx_ft,dy_ft,x_sw_ft,y_sw_ft,x_sw_fca,y_sw_fca	
	real,dimension(nz_ft),intent(IN) :: z_ft
	real,dimension(nx_ft,ny_ft,nz_ft),intent(IN) :: field_ft
	real,dimension(firegrid%nx,firegrid%ny,firegrid%nz) :: field_fca
	real :: dutmx,dutmy
	integer :: i,j,k			 ! N/A, cell in the fuel domain
	integer,dimension(firegrid%nx) :: ii
	integer,dimension(firegrid%ny) :: jj			 ! N/A
	integer,dimension(1) :: kk
	
	dutmx = x_sw_fca - x_sw_ft
	dutmy = y_sw_fca - y_sw_ft
	
	do j = 1, firegrid%ny
		jj(j) = ceiling(((real(j)-0.5)*firegrid%dy + dutmy) / dy_ft)
	enddo
	do i = 1, firegrid%nx
		ii(i) = ceiling(((real(i)-0.5)*firegrid%dx + dutmx) / dx_ft)
	enddo
	
	!$OMP parallel do private(i,j,k,ii,jj,kk)
	do k = 1, firegrid%nz
		kk = minloc(abs(z_ft - firegrid%zm(k)))
		do j = 1, firegrid%ny
			do i = 1, firegrid%nx
				field_fca(i,j,k) = field_ft(ii(i),jj(j),kk(1))
			enddo
		enddo
	enddo
	!$OMP end parallel do 
	
	END
!	
!=====================================================================================================================
!=====================================================================================================================
	SUBROUTINE GetWinds(coord_fb,winds,shear,is_simple,w_mult)
	
	
	use grid_module
	use winds_module
	
	implicit none
	
	real,intent(IN) :: w_mult
 	real,intent(IN),dimension(3) :: coord_fb	! m, coordinates of the firebrands
	real,dimension(3) :: winds		! m/s, wind components at the firebrand location
	integer,intent(IN) :: is_simple
	integer :: iw,jw,kw
	real :: shear
	real :: dudx,dvdy,dwdz,dudy,dvdx,dudz,dwdx,dvdz,dwdy
	
	! Get u,v
	iw = nint(coord_fb(1)*qugrid%dxi)
	iw = min(max(iw,1),qugrid%nx)
	jw = nint(coord_fb(2)*qugrid%dyi)
	jw = min(max(jw,1),qugrid%ny)
	kw = 2
	do while(qugrid%z(kw) < coord_fb(3))
		kw = kw + 1
	enddo
	kw = min(max(kw,2),qugrid%nz)
	
	winds(1) = quwinds%u(iw,jw,kw)
	winds(2) = quwinds%v(iw,jw,kw)
		
	! Get w
	iw = ceiling(coord_fb(1)*qugrid%dxi)
	jw = ceiling(coord_fb(2)*qugrid%dyi)
	kw = kw + 1
	kw = min(max(kw,3),qugrid%nz)
	winds(3) = quwinds%w(iw,jw,kw)
	
	winds(3) = winds(3) * w_mult
	
	! --------- shear	
	if(is_simple == 1) then
		shear = 0.
	else
		iw = ceiling(coord_fb(1)*qugrid%dxi)
		jw = ceiling(coord_fb(2)*qugrid%dyi)
	
		dudx = (quwinds%u(iw+1,jw,kw) - quwinds%u(iw,jw,kw) ) * qugrid%dxi
		dvdy = (quwinds%v(iw,jw+1,kw) - quwinds%v(iw,jw,kw) ) * qugrid%dyi
		dwdz = (quwinds%w(iw,jw,kw)   - quwinds%w(iw,jw,kw-1)) / qugrid%dz_array(kw)
	
		dudy = 0.5 * (														&
			(quwinds%u(iw  ,jw+1,kw) - quwinds%u(iw  ,jw,kw)) * qugrid%dyi +	&
			(quwinds%u(iw+1,jw+1,kw) - quwinds%u(iw+1,jw,kw)) * qugrid%dyi)
	
		dvdx = 0.5 * (														&
			(quwinds%v(iw  ,jw+1,kw) - quwinds%v(iw  ,jw,kw)) * qugrid%dxi +	&
			(quwinds%v(iw+1,jw+1,kw) - quwinds%v(iw+1,jw,kw)) * qugrid%dxi)
		
		dudz = 0.5 * (																	&
			(quwinds%u(iw  ,jw,kw) - quwinds%u(iw  ,jw,kw-1)) / qugrid%dz_array(kw) +	&
			(quwinds%u(iw+1,jw,kw) - quwinds%u(iw+1,jw,kw-1)) / qugrid%dz_array(kw))
		
		dwdx = (quwinds%w(iw+1,jw,kw) - quwinds%w(iw,jw,kw)) * qugrid%dxi
	
		dvdz = 0.5 * (																	&
			(quwinds%v(iw,jw  ,kw) - quwinds%v(iw,jw  ,kw-1)) / qugrid%dz_array(kw) +	&
			(quwinds%v(iw,jw+1,kw) - quwinds%v(iw,jw+1,kw-1)) / qugrid%dz_array(kw))
	
		dwdy = (quwinds%w(iw,jw+1,kw) - quwinds%w(iw,jw,kw)) * qugrid%dyi
	
		shear = 2.0*(dudx**2+dvdy**2+ dwdz**2) + (dudy+dvdx) + (dudz+dwdx) + (dvdz+dwdy)
	endif
	
	END
!=====================================================================================================================
!=====================================================================================================================