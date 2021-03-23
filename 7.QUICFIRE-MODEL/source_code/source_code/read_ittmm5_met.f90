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
      subroutine read_ittmm5_met

          
			use constants
         use file_handling_module
         use grid_module
			use wind_profile_module
         use time_module
         
! operations  done a priori to speed up the code AAG & IS  07/03/06
         implicit none
         
			integer ii,jj,kk,tt
         character*128 f_name,dum    !filename for a data profile at one of the sites.
         integer nittx,nitty,nittz
         real ittutmx,ittutmy,ittutmzone,rad2deg
         real dumx,dumy,dumz,dumu,dumv,costheta,sintheta
         real dist1,dist2,dcx,dcy,meso_search_distance
         real, allocatable :: ittx(:,:,:),itty(:,:,:),ittz(:,:,:)
         real, allocatable :: ittzo(:,:,:),ittws(:,:,:,:),ittwd(:,:,:,:)
         integer, allocatable :: x_idx(:),y_idx(:)
         
         rad2deg=180./pi
         read(228,*)meso_search_distance
         read(228,*) ! Description line
			allocate(ittx(nittx,nitty,nittz),&
                  itty(nittx,nitty,nittz),&
                  ittz(nittx,nitty,nittz))
         allocate(ittzo(nittx,nitty,qutime%nsteps),&
                  ittws(nittx,nitty,nittz,qutime%nsteps),&
                  ittwd(nittx,nitty,nittz,qutime%nsteps))
         ittzo(:,:,:)=0.
         ittws(:,:,:,:)=0.
         ittwd(:,:,:,:)=0.
         do tt=1,qutime%nsteps
            read(228,*)f_name				!File name of the individual file
            !Reading each profile/sensor file
            open(unit=ID_FILE_SENSOR,file=TRIM(workingDirectory)//TRIM(fileSeparator)//TRIM(f_name),status='old')
            read(ID_FILE_SENSOR,*) !description line
            read(ID_FILE_SENSOR,*)dum !lat long of sw grid cell center
            read(ID_FILE_SENSOR,*)ittutmx,ittutmy,ittutmzone !utmx of sw grid cell center
            read(ID_FILE_SENSOR,*)dum !time stamp
            read(ID_FILE_SENSOR,*)nittx !num cells in east-west direction for MM5 data
            read(ID_FILE_SENSOR,*)nitty !num cells in north-south direction for MM5 data
            read(ID_FILE_SENSOR,*)nittz !num cells in vertical direction for MM5 data
            read(ID_FILE_SENSOR,*) !description line
            read(ID_FILE_SENSOR,*) !description line
            if(tt .eq. 1)then
               if(ittutmzone .ne. qugrid%utmzone)then
                  print*,'Project UTM Zone and Met data UTM Zone are not equal'
               endif
               
            endif
            ittx(:,:,:)=0.
            itty(:,:,:)=0.
            ittz(:,:,:)=0.
            do ii=1,nittx
               do jj=1,nitty
                  read(ID_FILE_SENSOR,*)dumx,dumy,dumz,ittzo(ii,jj,tt)
               enddo
            enddo
            read(ID_FILE_SENSOR,*) ! description line
            read(ID_FILE_SENSOR,*) ! description line
            read(ID_FILE_SENSOR,*) ! description line
            do kk=1,nittz
               do jj=1,nitty
                  do ii=1,nittx
                     read(ID_FILE_SENSOR,*)ittx(ii,jj,kk),itty(ii,jj,kk),ittz(ii,jj,kk),&
                               dumu,dumv
                     ittws(ii,jj,kk,tt)=sqrt((dumu**2.)+(dumv**2.))
                     ittwd(ii,jj,kk,tt)=270.-rad2deg*atan2(dumv,dumu)
                  enddo
               enddo
            enddo
            close(ID_FILE_SENSOR)	! closing the profile/sensor data file
         enddo
         costheta=cos(-qugrid%domain_rotation*pi/180.)
         sintheta=sin(-qugrid%domain_rotation*pi/180.)
         do kk=1,nittz
            do jj=1,nitty
               do ii=1,nittx
                  ittx(ii,jj,kk)=ittx(ii,jj,kk)-qugrid%utmx
                  itty(ii,jj,kk)=itty(ii,jj,kk)-qugrid%utmy
                  dumx=costheta*ittx(ii,jj,kk)+sintheta*itty(ii,jj,kk)
                  dumy=-sintheta*ittx(ii,jj,kk)+costheta*itty(ii,jj,kk)
                  ittx(ii,jj,kk)=dumx
                  itty(ii,jj,kk)=dumy
               enddo
            enddo
         enddo
         windProfile%num_sites=0
         allocate(x_idx(nittx*nitty),y_idx(nittx*nitty))
         do jj=1,nitty
            do ii=1,nittx
               if(ittx(ii,jj,1) .ge. -meso_search_distance .and. &
                     ittx(ii,jj,1) .le. (qugrid%nx-1)*qugrid%dx+meso_search_distance .and. &
                     itty(ii,jj,1) .ge. -meso_search_distance .and. &
                     itty(ii,jj,1) .le. (qugrid%ny-1)*qugrid%dy+meso_search_distance)then
                  windProfile%num_sites=windProfile%num_sites+1
                  x_idx(windProfile%num_sites)=ii
                  y_idx(windProfile%num_sites)=jj
               endif
            enddo
         enddo
         if(windProfile%num_sites .eq. 0)then
            windProfile%num_sites=1
            dist1=0.5*((ittx(nittx,1,1)-ittx(1,1,1))+(ittx(1,nitty,1)-ittx(1,1,1)))
            dcx=0.5*qugrid%dx*(qugrid%nx-1)
            dcy=0.5*qugrid%dy*(qugrid%ny-1)
            do jj=1,nitty
               do ii=1,nittx
                  dist2=sqrt(((ittx(ii,jj,1)-dcx)**2.)+((itty(ii,jj,1)-dcy)**2.))
                  if(dist2 .le. dist1)then
                     dist1=dist2
                     x_idx(windProfile%num_sites)=ii
                     y_idx(windProfile%num_sites)=jj
                  endif
               enddo
            enddo
            ittx(x_idx(1),y_idx(1),1)=dcx
            itty(x_idx(1),y_idx(1),1)=dcy
         endif
         windProfile%num_vert_points=nittz
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
         allocate(windProfile%site%z_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%ws_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%wd_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         windProfile%site%blayer_flag = 4
         windProfile%site%u_prof = 0
         windProfile%site%v_prof = 0
         windProfile%site%xcoord = 0
         windProfile%site%ycoord = 0
         windProfile%site%nz_data=nittz

lp001:   do kk=1,windProfile%num_sites
            windProfile%site%xcoord(kk)=ittx(x_idx(kk),y_idx(kk),1) ! x coordinate of site location (meters)
            windProfile%site%ycoord(kk)=itty(x_idx(kk),y_idx(kk),1) ! y coordinate of site location (meters)
            do tt=1,qutime%nsteps
               windProfile%site%pp(kk,tt)=ittzo(x_idx(kk),y_idx(kk),tt)		! if blayer = 2 windProfile%site%pp = exp else windProfile%site%pp = zo
               do ii=1,windProfile%site%nz_data(kk,tt)
                  windProfile%site%z_data(kk,tt,ii)=ittz(x_idx(kk),y_idx(kk),ii)
                  windProfile%site%ws_data(kk,tt,ii)=ittws(x_idx(kk),y_idx(kk),ii,tt)
! MAN 02/05/2007 Domain Rotation
                  windProfile%site%wd_data(kk,tt,ii)=ittwd(x_idx(kk),y_idx(kk),ii,tt)-qugrid%domain_rotation
                  if(windProfile%site%wd_data(kk,tt,ii).lt. 0.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)+360.
                  elseif(windProfile%site%wd_data(kk,tt,ii).ge. 360.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)-360.
                  endif
! end MAN 02/05/2007
               enddo
            enddo
         enddo   lp001       ! kk = windProfile%num_sites
         deallocate(ittx,itty,ittz)
         deallocate(ittzo,ittws,ittwd)
         deallocate(x_idx,y_idx)
			
         return
      end
