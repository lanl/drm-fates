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
      subroutine read_quic_met

         
			use constants
			use file_handling_module
         use wind_profile_module
         use grid_module
			use time_module
! operations  done a priori to speed up the code AAG & IS  07/03/06
         implicit none
         
			integer ii,kk,tt,coordFlag, utm_zone
         character*128 f_name    ! filename for a data profile at one of the sites.
         real*8 convergence, lon, lat, utm_x, utm_y, theta
         
         ! blayer_flag has already been read in
         read(228,*)windProfile%num_sites			! total number of sites
         read(228,*)windProfile%num_vert_points
         ! allocate arrays
         allocate(windProfile%site%blayer_flag(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%pp(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%H(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%rL(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%ac(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%xcoord(windProfile%num_sites))
         allocate(windProfile%site%ycoord(windProfile%num_sites))
         allocate(windProfile%site%u_prof(windProfile%num_sites,qugrid%nz),   &
            windProfile%site%v_prof(windProfile%num_sites,qugrid%nz))
         allocate(windProfile%site%uoint(windProfile%num_sites),windProfile%site%voint(windProfile%num_sites))
         allocate(windProfile%site%wm(windProfile%num_sites,qugrid%nx,qugrid%ny),   &
            windProfile%site%wms(windProfile%num_sites,qugrid%nx,qugrid%ny))
         allocate(windProfile%site%nz_data(windProfile%num_sites,qutime%nsteps))
         ! MAN 4/5/2007 number of data points to allocate in data points profiles
         allocate(windProfile%site%z_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%ws_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%wd_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))

         ! initialize velocity profiles at each site
         windProfile%site%u_prof = 0
         windProfile%site%v_prof = 0
         
         theta=real(qugrid%domain_rotation*pi/180.,8)

         ! initialize the vectors
         windProfile%site%xcoord(1:windProfile%num_sites) = 0
         windProfile%site%ycoord(1:windProfile%num_sites) = 0
         windProfile%site%nz_data(:,:)=1

lp001:   do kk=1,windProfile%num_sites
            read(228,*)						! The site name
            read(228,*)						! Description line
            read(228,*)f_name				! File name of the individual file
            ! Reading each profile/sensor file
            open(unit=ID_FILE_SENSOR,file=TRIM(workingDirectory)//TRIM(fileSeparator)//TRIM(f_name),status='old')
            read(ID_FILE_SENSOR,*)						! The site name
				read(ID_FILE_SENSOR,*)
				read(ID_FILE_SENSOR,*)
            read(ID_FILE_SENSOR,*)coordFlag
            read(ID_FILE_SENSOR,*)windProfile%site%xcoord(kk)		! x coordinate of site location (meters)
            read(ID_FILE_SENSOR,*)windProfile%site%ycoord(kk)		! y coordinate of site location (meters)
            select case(coordFlag)
               case(2) ! UTM coordinates
                  read(ID_FILE_SENSOR,*) utm_x      ! UTMX position of site location (meters)
                  read(ID_FILE_SENSOR,*) utm_y      ! UTMY position of site location (meters)
                  read(ID_FILE_SENSOR,*) utm_zone   ! UTM zone number of site location
                  read(ID_FILE_SENSOR,*) ! UTM zone letter of site location (1=A, 2=B, etc.)
               case(3)
                  read(ID_FILE_SENSOR,*) lat ! latitude of site location (dd)
                  read(ID_FILE_SENSOR,*) lon ! longitude of site location (dd)
            end select
            convergence = 0.0
            if(qugrid%utmx .ne. 0. .and. qugrid%utmy .ne. 0.) then
            	select case(coordFlag)
            		case(1) ! QUIC coordinates
                     utm_x = real(windProfile%site%xcoord(kk),8)*dcos(theta) + &
                        real(windProfile%site%ycoord(kk),8)*dsin(theta) + real(qugrid%utmx,8)
                     utm_y = -real(windProfile%site%xcoord(kk),8)*dsin(theta) + &
                        real(windProfile%site%ycoord(kk),8)*dcos(theta) + real(qugrid%utmy,8)
            			call utm_geo(lon, lat, utm_x, utm_y, qugrid%utmzone, 1)
                  case(2) ! UTM coordinates
                     call utm_geo(lon, lat, utm_x, utm_y, utm_zone, 1)
               end select
               call getConvergence(lon, lat, qugrid%utmzone, convergence)
            endif
            do tt=1,qutime%nsteps
               read(ID_FILE_SENSOR,*) !time stamp
               read(ID_FILE_SENSOR,*)windProfile%site%blayer_flag(kk,tt) ! boundary layer flag for each site (1=log,2=exp,3=canopy,4=data)e
	            read(ID_FILE_SENSOR,*)windProfile%site%pp(kk,tt)			! if blayer = 2 windProfile%site%pp = exp else windProfile%site%pp = zo
               select case(windProfile%site%blayer_flag(kk,tt))
                  case(1)! logarithmic profile
                     read(ID_FILE_SENSOR,*)windProfile%site%rL(kk,tt)			! reciprocal Monin-Obukhov length
                  case(3)! urban canopy
                     read(ID_FILE_SENSOR,*)windProfile%site%rL(kk,tt)			! reciprocal Monin-Obukhov length
                     read(ID_FILE_SENSOR,*)windProfile%site%H(kk,tt)			   ! canopy height
                     read(ID_FILE_SENSOR,*)windProfile%site%ac(kk,tt)			! atenuation coefficient
                  case(4)! data points
                     read(ID_FILE_SENSOR,*)windProfile%site%nz_data(kk,tt)		! number of data points in vertical wind profile
               end select
               read(ID_FILE_SENSOR,*)! skip line			!"height  direction   magnitude"  Label
               do ii=1,windProfile%site%nz_data(kk,tt)
                  read(ID_FILE_SENSOR,*)windProfile%site%z_data(kk,tt,ii),&
                     windProfile%site%ws_data(kk,tt,ii),windProfile%site%wd_data(kk,tt,ii)
! MAN 02/05/2007 Domain Rotation
                  ! Scot - Modify here... 
                  windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii) - &
							qugrid%domain_rotation - real(convergence)
                  if(windProfile%site%wd_data(kk,tt,ii).lt. 0.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)+360.
                  elseif(windProfile%site%wd_data(kk,tt,ii).ge. 360.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)-360.
                  endif
               enddo
            enddo
            close(ID_FILE_SENSOR)	! closing the profile/sensor data file
         enddo   lp001       !kk = windProfile%num_sites
! MAN 10/10/2007 if there is only one measurement make sure it is in the domain
         if(windProfile%num_sites .eq. 1)then
            windProfile%site%xcoord(1)=0.5*real(qugrid%nx-1)*qugrid%dx
            windProfile%site%ycoord(1)=0.5*real(qugrid%ny-1)*qugrid%dy
         endif
! end MAN 10/10/2007
         return
      end


subroutine getConvergence(lon, lat, zone, convergence)
    implicit none
    integer zone
    double precision, parameter :: PI = 3.141592653589793d0
    double precision, parameter :: degrad=PI/180.d0, raddeg=180.d0/PI
    double precision lat, lon, tempLon, convergence
    
    tempLon = (zone * 6.0) - 183.0 - lon
    
    convergence = datan(dtan(tempLon * degrad) * dsin(lat * degrad)) * raddeg
end subroutine