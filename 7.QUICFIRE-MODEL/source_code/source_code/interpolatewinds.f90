	SUBROUTINE INTERPOLATEWINDS()

	
	use constants
	use ship_module
	use wind_profile_module
	use bld_module
	use grid_module
	use winds_module
	use time_module
	
	implicit none

	integer :: i,j,k,ii,kk
	integer iwork, jwork
	real ustar
	real site_zref, theta_site, mag_site, umult_site, vmult_site

	real yc, xc, rc, rcsum, rcval, d
	real sumw, sumwv, sumwu
	real sgamma, lamda, deln
	real dxx, dyy, u12, u34, v12, v34
	! MAN 04/05/2007 Time varying profile parameters
	real bisect
	real xtemp,psi_m !erp 9/18/2006 stability variables

	! MAN 10/23/2013 advanced data point interpolation
	integer iter,logflag, i_time
	real zonew,zolow,zohigh,unew,unewl,unewh
	real dwinddir
	real a1,a2,a3, uH
	integer InterpAverage
	real Usum,Vsum,ZlowerBound,ZupperBound


	i_time = qutime%current_iter
	! ---------------------------------------------------------------------------------
	! This do loop defines each vertical velocity profile at each sensor site so the
	! Barnes Mapping interpolation scheme has a velocity at each site for each height
	lp003:   do kk=1,windProfile%num_sites
		if(windProfile%site%blayer_flag(kk,i_time) .lt. 4)then
			theta_site=(270.-windProfile%site%wd_data(kk,i_time,1))*pi/180.
			mag_site = windProfile%site%ws_data(kk,i_time,1)
			site_zref = windProfile%site%z_data(kk,i_time,1)
			umult_site=cos(theta_site)
			vmult_site=sin(theta_site)
			! power law profile
			if(windProfile%site%blayer_flag(kk,i_time) .eq. 2)then
				! MAN 07/25/2008 stretched vertical grid
				do k=2,qugrid%nz      ! loop through vertical cell indices
					windProfile%site%u_prof(kk,k)=umult_site*mag_site*	&
						((qugrid%zm(k)/site_zref)**windProfile%site%pp(kk,i_time))
					windProfile%site%v_prof(kk,k)=vmult_site*mag_site*	&
						((qugrid%zm(k)/site_zref)**windProfile%site%pp(kk,i_time))
				enddo
			endif !erp 2/6/2003 end power law profile
			! logrithmic velocity profile
			if(windProfile%site%blayer_flag(kk,i_time) .eq. 1)then
				! MAN 05/15/2007 adjust for stability
				do k=2,qugrid%nz      ! loop through vertical cell indices
					if(k .eq. 2)then
						if(site_zref*windProfile%site%rL(kk,i_time) .ge. 0)then
							psi_m=4.7*(site_zref+windProfile%site%pp(kk,i_time))*windProfile%site%rL(kk,i_time)
						else
							xtemp=(1.-15*(site_zref+windProfile%site%pp(kk,i_time))*windProfile%site%rL(kk,i_time))**(0.25)
							psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
						endif
						ustar=mag_site*vk/(log((site_zref+windProfile%site%pp(kk,i_time))/windProfile%site%pp(kk,i_time))+psi_m)
					endif
					! end MAN 05/15/2007
					! MAN 07/25/2008 stretched vertical grid
					if(windProfile%site%rL(kk,i_time) .ge. 0)then
						psi_m=4.7*(qugrid%zm(k)+windProfile%site%pp(kk,i_time))*windProfile%site%rL(kk,i_time)
					else
						xtemp=(1.-15*(qugrid%zm(k)+windProfile%site%pp(kk,i_time))*windProfile%site%rL(kk,i_time))**(0.25)
						psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
					endif
					windProfile%site%u_prof(kk,k)=(umult_site*ustar/vk)*(log((qugrid%zm(k)+windProfile%site%pp(kk,i_time)) &
						/windProfile%site%pp(kk,i_time))+psi_m)
					windProfile%site%v_prof(kk,k)=(vmult_site*ustar/vk)*(log((qugrid%zm(k)+windProfile%site%pp(kk,i_time)) &
						/windProfile%site%pp(kk,i_time))+psi_m)
				enddo
			endif ! erp 2/6/2003 end log law profile
			! Canopy profile
			if(windProfile%site%blayer_flag(kk,i_time) .eq. 3)then
				do k=2,qugrid%nz      ! loop through vertical cell indices
					if(k .eq. 2)then ! only calculate d once
						! MAN 05/15/2007 adjust for stability
						if(site_zref*windProfile%site%rL(kk,i_time) .ge. 0)then
							psi_m=4.7*site_zref*windProfile%site%rL(kk,i_time)
						else
							xtemp=(1.-15*site_zref*windProfile%site%rL(kk,i_time))**(0.25)
							psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
						endif
						ustar = mag_site*vk/(log(site_zref/windProfile%site%pp(kk,i_time))+psi_m)
						d = bisect(ustar,windProfile%site%pp(kk,i_time),windProfile%site%H(kk,i_time), &
							windProfile%site%ac(kk,i_time),psi_m)
						! end MAN 05/15/2007
						if(windProfile%site%H(kk,i_time)*windProfile%site%rL(kk,i_time) .ge. 0)then
							psi_m=4.7*(windProfile%site%H(kk,i_time)-d)*windProfile%site%rL(kk,i_time)
						else
							xtemp=(1.-15*(windProfile%site%H(kk,i_time)-d)*windProfile%site%rL(kk,i_time))**(0.25)
							psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
						endif
						uH = (ustar/vk)*(log((windProfile%site%H(kk,i_time)-d)/windProfile%site%pp(kk,i_time))+psi_m);
						if(site_zref .le. windProfile%site%H(kk,i_time))then
							mag_site=mag_site/(uH*exp(windProfile%site%ac(kk,i_time)*(site_zref/windProfile%site%H(kk,i_time) -1.)))
						else
							if(site_zref*windProfile%site%rL(kk,i_time) .ge. 0)then
								psi_m=4.7*(site_zref-d)*windProfile%site%rL(kk,i_time)
							else
								xtemp=(1.-15*(site_zref-d)*windProfile%site%rL(kk,i_time))**(0.25)
								psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
							endif
							mag_site=mag_site/((ustar/vk)*(log((site_zref-d)/bld%zo)+psi_m))
						endif
						ustar=mag_site*ustar
						uH=mag_site*uH
					endif
					! MAN 07/25/2008 stretched vertical grid
					if(qugrid%zm(k) .le. windProfile%site%H(kk,i_time))then ! lower canopy profile
						windProfile%site%u_prof(kk,k) = umult_site * uH*exp(windProfile%site%ac(kk,i_time)&
							*(qugrid%zm(k)/windProfile%site%H(kk,i_time) -1))
						windProfile%site%v_prof(kk,k) = vmult_site * uH*exp(windProfile%site%ac(kk,i_time)&
							*(qugrid%zm(k)/windProfile%site%H(kk,i_time) -1))
					endif
					if(qugrid%zm(k) .gt. windProfile%site%H(kk,i_time))then ! upper canopy profile
						if(qugrid%zm(k)*windProfile%site%rL(kk,i_time) .ge. 0)then
							psi_m=4.7*(qugrid%zm(k)-d)*windProfile%site%rL(kk,i_time)
						else
							xtemp=(1.-15*(qugrid%zm(k)-d)*windProfile%site%rL(kk,i_time))**(0.25)
							psi_m=-2*log(0.5*(1+xtemp))-log(0.5*(1+xtemp**2.))+2*atan(xtemp)-0.5*pi
						endif
						windProfile%site%u_prof(kk,k)=(umult_site*ustar/vk)*&
							(log((qugrid%zm(k)-d)/windProfile%site%pp(kk,i_time))+psi_m)

						windProfile%site%v_prof(kk,k)=(vmult_site*ustar/vk)*&
							(log((qugrid%zm(k)-d)/windProfile%site%pp(kk,i_time))+psi_m)
					endif !end urban canopy TMB 6/16/03
				enddo
			endif
		endif
		! new 2/7/2005 velocity profile entry
		if(windProfile%site%blayer_flag(kk,i_time) .eq. 4)then  ! data entry profile
			ii=0
			theta_site=(270.-windProfile%site%wd_data(kk,i_time,1))*pi/180.
			do k=2,qugrid%nz ! loop through vertical cell indices
				if(qugrid%zm(k) .lt. windProfile%site%z_data(kk,i_time,1) .or. windProfile%site%nz_data(kk,i_time) .eq. 1)then
					windProfile%site%u_prof(kk,k)= (windProfile%site%ws_data(kk,i_time,1)*cos(theta_site)/&
						(log((windProfile%site%z_data(kk,i_time,1)+windProfile%site%pp(kk,i_time)) &
						/windProfile%site%pp(kk,i_time))))*log((qugrid%zm(k)+windProfile%site%pp(kk,i_time))/ &
						windProfile%site%pp(kk,i_time))
					windProfile%site%v_prof(kk,k)= (windProfile%site%ws_data(kk,i_time,1)*sin(theta_site)/ &
					(log((windProfile%site%z_data(kk,i_time,1)+windProfile%site%pp(kk,i_time)) &
						/windProfile%site%pp(kk,i_time))))*log((qugrid%zm(k)+windProfile%site%pp(kk,i_time))/ &
						windProfile%site%pp(kk,i_time))
				elseif( qugrid%zm(k) .ge. windProfile%site%z_data(kk,i_time,windProfile%site%nz_data(kk,i_time)))then
					ii=windProfile%site%nz_data(kk,i_time)
					mag_site=windProfile%site%ws_data(kk,i_time,ii)
					theta_site=(270.-(windProfile%site%wd_data(kk,i_time,ii)))*pi/180.
					windProfile%site%u_prof(kk,k)= mag_site*cos(theta_site)
					windProfile%site%v_prof(kk,k)= mag_site*sin(theta_site)
				else
					InterpAverage=0
					if(ii .lt. windProfile%site%nz_data(kk,i_time)-1 .and. &
						qugrid%zm(k) .ge. windProfile%site%z_data(kk,i_time,ii+1))then
						ii=ii+1
						if(qugrid%z(k) .gt. windProfile%site%z_data(kk,i_time,ii+1))then
							InterpAverage=1
							! grid cells are larger than the distance between data points: average data in grid cell
							ZlowerBound=qugrid%z(k-1) !start at the bottom of the wind cell
							ZupperBound=0.5*(windProfile%site%z_data(kk,i_time,ii)+	&
								windProfile%site%z_data(kk,i_time,ii+1)) !influence of ii extends to halfway between points
							theta_site=(270.-windProfile%site%wd_data(kk,i_time,ii))*pi/180.
							Usum=windProfile%site%ws_data(kk,i_time,ii)*cos(theta_site)*(ZupperBound-ZlowerBound)
							Vsum=windProfile%site%ws_data(kk,i_time,ii)*sin(theta_site)*(ZupperBound-ZlowerBound)
							do while(ii .lt. windProfile%site%nz_data(kk,i_time)-1 .and.	&
								 qugrid%z(k) .ge. windProfile%site%z_data(kk,i_time,ii+1))
								ii=ii+1
								ZlowerBound=ZupperBound
								ZupperBound=0.5*(windProfile%site%z_data(kk,i_time,ii)+windProfile%site%z_data(kk,i_time,ii+1))
								if(qugrid%z(k) .lt. ZupperBound)then
									ZupperBound = qugrid%z(k)
								endif
								theta_site=(270.-windProfile%site%wd_data(kk,i_time,ii))*pi/180.
								Usum=windProfile%site%ws_data(kk,i_time,ii)*cos(theta_site)*(ZupperBound-ZlowerBound)+Usum
								Vsum=windProfile%site%ws_data(kk,i_time,ii)*sin(theta_site)*(ZupperBound-ZlowerBound)+Vsum
							enddo
							if(ii .eq. windProfile%site%nz_data(kk,i_time)-1 .and. &
								qugrid%z(k) .ge. windProfile%site%z_data(kk,i_time,ii+1))then
								ii=ii+1
								ZlowerBound=ZupperBound
								ZupperBound = qugrid%z(k)
								theta_site=(270.-windProfile%site%wd_data(kk,i_time,ii))*pi/180.
								Usum=windProfile%site%ws_data(kk,i_time,ii)*cos(theta_site)*(ZupperBound-ZlowerBound)+Usum
								Vsum=windProfile%site%ws_data(kk,i_time,ii)*sin(theta_site)*(ZupperBound-ZlowerBound)+Vsum
							endif
							windProfile%site%u_prof(kk,k)= Usum/qugrid%dz_array(k)
							windProfile%site%v_prof(kk,k)= Vsum/qugrid%dz_array(k)

						else
							! grid cells are smaller than data points: interpolate
							InterpAverage=0
						endif
						if(abs(windProfile%site%wd_data(kk,i_time,ii+1)-windProfile%site%wd_data(kk,i_time,ii)) .gt. 180.)then
							if(windProfile%site%wd_data(kk,i_time,ii+1) .gt. windProfile%site%wd_data(kk,i_time,ii))then
								dwinddir=(windProfile%site%wd_data(kk,i_time,ii+1)-360.-windProfile%site%wd_data(kk,i_time,ii)) &
									/(windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii))
							else
								dwinddir=(windProfile%site%wd_data(kk,i_time,ii+1)-windProfile%site%wd_data(kk,i_time,ii)+360.) &
									/(windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii))
							endif
						else
							dwinddir=(windProfile%site%wd_data(kk,i_time,ii+1)-windProfile%site%wd_data(kk,i_time,ii)) &
								/(windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii))
						endif
						! a1=(windProfile%site%ws_data(kk,i_time,ii+1)-windProfile%site%ws_data(kk,i_time,ii))/(windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii))

						zohigh=20.
						ustar=vk*windProfile%site%ws_data(kk,i_time,ii)/log((windProfile%site%z_data(kk,i_time,ii)+zohigh)/zohigh)
						unewh=(ustar/vk)*log((windProfile%site%z_data(kk,i_time,ii)+zohigh)/zohigh)
						zolow=1.e-9
						ustar=vk*windProfile%site%ws_data(kk,i_time,ii)/log((windProfile%site%z_data(kk,i_time,ii)+zolow)/zolow)
						unewl=(ustar/vk)*log((windProfile%site%z_data(kk,i_time,ii+1)+zolow)/zolow)
						if(windProfile%site%ws_data(kk,i_time,ii+1) .gt. unewl .and. windProfile%site%ws_data(kk,i_time,ii+1) .lt. unewh)then
							logflag=1
							iter=0
							zonew=windProfile%site%pp(kk,i_time)
							ustar=vk*windProfile%site%ws_data(kk,i_time,ii)/log((windProfile%site%z_data(kk,i_time,ii)+zonew)/zonew)
							unew=(ustar/vk)*log((windProfile%site%z_data(kk,i_time,ii+1)+zonew)/zonew)
							do while(iter .lt. 200 .and. abs(unew-windProfile%site%ws_data(kk,i_time,ii+1)) &
								.gt. 0.0001*windProfile%site%ws_data(kk,i_time,ii+1))
								iter=iter+1
								zonew=0.5*(zolow+zohigh)
								ustar=vk*windProfile%site%ws_data(kk,i_time,ii)/log((windProfile%site%z_data(kk,i_time,ii)+zonew)/zonew)
								unew=(ustar/vk)*log((windProfile%site%z_data(kk,i_time,ii+1)+zonew)/zonew)
								if(unew .gt. windProfile%site%ws_data(kk,i_time,ii+1))then
									zohigh=zonew
								else
									zolow=zonew
								endif
							enddo
						else
							logflag=0
							if(ii .lt. windProfile%site%nz_data(kk,i_time)-1)then
								a1=((windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii)) &
									*(windProfile%site%ws_data(kk,i_time,ii+2)-windProfile%site%ws_data(kk,i_time,ii)) &
									+(windProfile%site%z_data(kk,i_time,ii)-windProfile%site%z_data(kk,i_time,ii+2)) &
									*(windProfile%site%ws_data(kk,i_time,ii+1)-windProfile%site%ws_data(kk,i_time,ii))) &
									/((windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii)) &
									*((windProfile%site%z_data(kk,i_time,ii+2)**2)-(windProfile%site%z_data(kk,i_time,ii)**2)) &
									+((windProfile%site%z_data(kk,i_time,ii+1)**2)-(windProfile%site%z_data(kk,i_time,ii)**2)) &
									*(windProfile%site%z_data(kk,i_time,ii)-windProfile%site%z_data(kk,i_time,ii+2)))
							else
								a1=0
							endif
							a2=((windProfile%site%ws_data(kk,i_time,ii+1)-windProfile%site%ws_data(kk,i_time,ii)) &
								-a1*((windProfile%site%z_data(kk,i_time,ii+1)**2)-(windProfile%site%z_data(kk,i_time,ii)**2))) &
								/(windProfile%site%z_data(kk,i_time,ii+1)-windProfile%site%z_data(kk,i_time,ii))
							a3=windProfile%site%ws_data(kk,i_time,ii)-a1*(windProfile%site%z_data(kk,i_time,ii)**2) &
								-a2*windProfile%site%z_data(kk,i_time,ii);
						endif
					endif
					if(InterpAverage .eq. 0)then
						if(logflag .eq. 1)then
							mag_site=(ustar/vk)*log((qugrid%zm(k)+zonew)/zonew);
						else
							mag_site=a1*(qugrid%zm(k)**2)+a2*qugrid%zm(k)+a3
						endif
						!mag_site=a1*(zm(k)-windProfile%site%z_data(kk,i_time,ii))+windProfile%site%ws_data(kk,i_time,ii)
						theta_site=(270.-(windProfile%site%wd_data(kk,i_time,ii)+ &
							dwinddir*(qugrid%zm(k)-windProfile%site%z_data(kk,i_time,ii))))*pi/180.
						windProfile%site%u_prof(kk,k)= mag_site*cos(theta_site)
						windProfile%site%v_prof(kk,k)= mag_site*sin(theta_site)
					endif
				endif
			enddo			
		endif ! erp 2/6/2003 end data entry
	enddo   lp003       ! num_sites kk
	if(ship%movingCoordsFlag .gt. 0)then
		qugrid%domain_rotation=ship%bearing(i_time)-ship%relativeBearing
		! Adjust the effective velocity field from the moving coordinate system
		print*,'Effective Wind Dir',windProfile%site%wd_data(1,i_time,1)
		print*,'Effective Wind Vector Bottom',windProfile%site%u_prof(kk,1),windProfile%site%v_prof(kk,1)
		print*,'Effective Wind Vector Top',windProfile%site%u_prof(kk,qugrid%nz-1),windProfile%site%v_prof(kk,qugrid%nz-1)
		print*,'Domain rotation',qugrid%domain_rotation
		print*,'Ship',ship%speed(i_time),ship%bearing(i_time)
		print*,'Current',ship%currentSpeed(i_time),ship%currentDirection(i_time)
		theta_site = (90.0 - ship%relativeBearing)*pi/180.
		ship%shipU = ship%speed(i_time)*cos(theta_site)
		ship%shipV = ship%speed(i_time)*sin(theta_site)
		theta_site = (90.0 - ship%currentDirection(i_time) + qugrid%domain_rotation)*pi/180.
		!print*,90.0 - currentDirection(i_time) + domain_rotation
		ship%currentU = ship%currentSpeed(i_time)*cos(theta_site)
		ship%currentV = ship%currentSpeed(i_time)*sin(theta_site)
		windProfile%site%u_prof(:,:) = windProfile%site%u_prof(:,:) - ship%shipU - ship%currentU
		windProfile%site%v_prof(:,:) = windProfile%site%v_prof(:,:) - ship%shipV - ship%currentV
		print*,'Ship',ship%shipU,ship%shipV
		print*,'Current',ship%currentU,ship%currentV
	endif
	! end MAN 04/05/2007
	! ---------------------------------------------------------------------------------
	! find average distance of each measuring station relative to all other
	! measuring stations

	if(windProfile%num_sites .eq. 1)then
		do k=2,qugrid%nz
			quwinds%uo(:,:,k)=windProfile%site%u_prof(1,k)
			quwinds%vo(:,:,k)=windProfile%site%v_prof(1,k)
		enddo
	else ! Barnes Mapping Scheme
		rcsum=0.
		do kk=1,windProfile%num_sites
			rcval=1000000. ! ignore anything over 1000 kilometers
			do k=1,windProfile%num_sites
				xc=windProfile%site%xcoord(k)-windProfile%site%xcoord(kk)
				yc=windProfile%site%ycoord(k)-windProfile%site%ycoord(kk)
				rc=sqrt(xc**2+yc**2)
				if(rc .lt. rcval .and. k .ne. kk) rcval=rc ! shortest distance
			enddo
			rcsum=rcval+rcsum ! sum of shortest distances
		enddo
		deln=rcsum/real(windProfile%num_sites)  ! average Radius (i.e. computed data spacing)
		lamda=5.052*(2.*deln/pi)**2 ! weight parameter
		! numerical convergence parameter
		sgamma = 0.2         ! gamma=.2 max detail, gamma=1 min detail
		! ------------------------------------------------------------------------------
		! first and second barnes pass done for each cell level in z direction

		! compute weight of each site on point (i,j)
		!$omp parallel do default(shared) private(i,j,kk)
		do j=1,qugrid%ny
			do i=1,qugrid%nx
				do kk=1,windProfile%num_sites
					windProfile%site%wm(kk,i,j)=exp(-1/lamda*(windProfile%site%xcoord(kk)-qugrid%xcenters(i))**2   &
						-1/lamda*(windProfile%site%ycoord(kk)-qugrid%ycenters(j))**2)
					windProfile%site%wms(kk,i,j)=exp(-1/(sgamma*lamda)*(windProfile%site%xcoord(kk)-qugrid%xcenters(i))**2   &
						-1/(sgamma*lamda)*(windProfile%site%ycoord(kk)-qugrid%ycenters(j))**2)
				enddo
				if(sum(windProfile%site%wm(:,i,j)) .eq. 0.)then
					windProfile%site%wm(:,i,j)=1e-20
				endif
			enddo
		enddo
		!$omp end parallel do
		lp004:      do k=2,qugrid%nz
			! interpolate onto the grid
			! do first Barnes pass
			!$omp parallel do private(i,sumwu,sumwv,sumw)
			do j=1,qugrid%ny
				do i=1,qugrid%nx
					sumwu=sum(windProfile%site%wm(:,i,j)*windProfile%site%u_prof(:,k))
					sumwv=sum(windProfile%site%wm(:,i,j)*windProfile%site%v_prof(:,k))
					sumw=sum(windProfile%site%wm(:,i,j))
					quwinds%uo(i,j,k)=sumwu/sumw
					quwinds%vo(i,j,k)=sumwv/sumw
				enddo ! i=nx
			enddo ! j=ny
			!$omp end parallel do
			! before doing the 2nd pass for the Barnes Method
			! use a 4-point bilinear interpolation
			! scheme to get estimated values at measured point (+)
			! using the 1st pass calculated data at grid points (*)
			!
			!     *       *            1        2
			!
			!          +                    +        !definition of points
			!
			!     *       *            3        4
			!     |    |  |
			!     | dxx|  |
			!     |       |  !definition of measurements
			!     |  ddx  |
			!
			do kk=1,windProfile%num_sites
				if(windProfile%site%xcoord(kk) .gt. 0. .and. &
					windProfile%site%xcoord(kk) .lt. real(qugrid%nx-1)*qugrid%dx .and. &
					windProfile%site%ycoord(kk) .gt. 0. .and. &
					windProfile%site%ycoord(kk) .lt. real(qugrid%ny-1)*qugrid%dy)then
					do j=1,qugrid%ny
						!find closest grid location on lower side of site
						if(qugrid%ycenters(j).lt.windProfile%site%ycoord(kk)) jwork=j
					enddo
					do i=1,qugrid%nx
						!find closest grid location on left side of site
						if(qugrid%xcenters(i).lt.windProfile%site%xcoord(kk)) iwork=i
					enddo !i=nx
					! distance to site point from lower and left sides
					dxx=windProfile%site%xcoord(kk)-qugrid%xcenters(iwork)
					dyy=windProfile%site%ycoord(kk)-qugrid%ycenters(jwork)
					! MAN 7/7/2005 fixed interpolation of velocities and var dz conversion
					! upper u interpolated velocity
					u12 = (1-dxx/qugrid%dx)*quwinds%uo(iwork,jwork+1,k)+(dxx/qugrid%dx)*quwinds%uo(iwork+1,jwork+1,k)
					! lower u interplotaed velocity
					u34 = (1-dxx/qugrid%dx)*quwinds%uo(iwork,jwork,k)+(dxx/qugrid%dx)*quwinds%uo(iwork+1,jwork,k)
					! total interpolated u velocity
					windProfile%site%uoint(kk)=(dyy/qugrid%dy)*u12+(1-dyy/qugrid%dy)*u34

					! upper v interpolated velocity
					v12 = (1-dxx/qugrid%dx)*quwinds%vo(iwork,jwork+1,k)+(dxx/qugrid%dx)*quwinds%vo(iwork+1,jwork+1,k)
					! lower v interplotaed velocity
					v34 = (1-dxx/qugrid%dx)*quwinds%vo(iwork,jwork,k)+(dxx/qugrid%dx)*quwinds%vo(iwork+1,jwork,k)
					! total interpolated u velocity
					windProfile%site%voint(kk)=(dyy/qugrid%dy)*v12+(1-dyy/qugrid%dy)*v34
				else
					windProfile%site%uoint(kk)=windProfile%site%u_prof(kk,k)
					windProfile%site%voint(kk)=windProfile%site%v_prof(kk,k)
				endif
				! end MAN 7/7/2005
			enddo ! kk=num_sites
			! end bilinear interpolation section
			! Begin 2nd Barnes pass
			!$omp parallel do private(i,sumwu,sumwv,sumw)
			do j=1,qugrid%ny
				do i=1,qugrid%nx
					sumwu=sum(windProfile%site%wms(:,i,j)*(windProfile%site%u_prof(:,k)-windProfile%site%uoint(:)))
					sumwv=sum(windProfile%site%wms(:,i,j)*(windProfile%site%v_prof(:,k)-windProfile%site%voint(:)))
					sumw=sum(windProfile%site%wms(:,i,j))
					if(sumw .ne. 0.)then
						quwinds%uo(i,j,k)=quwinds%uo(i,j,k) + sumwu/sumw
						quwinds%vo(i,j,k)=quwinds%vo(i,j,k) + sumwv/sumw
					endif
				enddo
			enddo
			!$omp end parallel do
		enddo   lp004       ! k=1:nz
	endif
	!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
	do k=1,qugrid%nz-1
		do j=1,qugrid%ny-1
			do i=1,qugrid%nx-1
				quwinds%uint(i,j,k)=(0.5*(quwinds%uo(i,j,k)+quwinds%uo(i+1,j,k)))
				quwinds%vint(i,j,k)=(0.5*(quwinds%vo(i,j,k)+quwinds%vo(i,j+1,k)))
				if(quwinds%uint(i,j,k) .ne. quwinds%uint(i,j,k) .or. quwinds%vint(i,j,k) .ne. quwinds%vint(i,j,k))then
					print*,'Interpolation scheme NaN at ',i,j,k
				endif
			enddo
		enddo
	enddo
	!$OMP END PARALLEL DO

	quwinds%u = 0.
	quwinds%v = 0.
	quwinds%w = 0.
	! erp initialize upwind array erp 6/18/04
	print*,'Finished interpolating winds'
	
	END SUBROUTINE INTERPOLATEWINDS