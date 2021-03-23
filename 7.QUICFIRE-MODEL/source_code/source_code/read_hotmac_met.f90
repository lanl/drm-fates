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
      subroutine read_hotmac_met

         
			use constants
			use file_handling_module
         use wind_profile_module
         use grid_module
         use time_module
			
         implicit none
			
         integer ii,jj,kk,tt,t_idx_hotmac,skip_idx_hotmac,i,j,k
         character*128 f_name    !filename for a data profile at one of the sites.
         integer nx_hotmac,ny_hotmac,nz_hotmac
         integer nx_hotmacp1,ny_hotmacp1,nz_hotmacp1,nzs_hotmac
         real dx_hotmac,dy_hotmac,meso_search_distance
         real utmx_hotmac,utmy_hotmac,rad2deg
         real dumx,dumy,dumu,dumv,costheta,sintheta
         real dist1,dist2,dcx,dcy,gmt,delgmt,gmtday,hotmac_time
         real, allocatable :: x_hotmac(:,:),y_hotmac(:,:),z_hotmac(:,:,:)
         real, allocatable :: u_hotmac(:,:,:),v_hotmac(:,:,:)
         real, allocatable :: zo_hotmac(:,:,:),ws_hotmac(:,:,:,:),wd_hotmac(:,:,:,:)
         integer, allocatable :: x_idx(:),y_idx(:)
         character*1 adum
         real, allocatable :: zgrnd(:,:),ztop(:,:),ustar_hotmac(:,:)
         real, allocatable :: hmz(:),hmzm(:)
         
         rad2deg=180./pi
         
         read(ID_FILE_QU_METPARAMS,*)meso_search_distance
         read(ID_FILE_QU_METPARAMS,*) !Description line
         read(ID_FILE_QU_METPARAMS,*)f_name				!File name of the individual file
         read(ID_FILE_QU_METPARAMS,*)skip_idx_hotmac  !Number of HOTMAC time steps to skip
         !Reading each profile/sensor file
         open(unit=ID_FILE_SENSOR,file=TRIM(workingDirectory)//TRIM(fileSeparator)//TRIM(f_name),form='unformatted',status='old')
         rewind(ID_FILE_SENSOR)
         do i=1,10
            read(ID_FILE_SENSOR)adum
         enddo
         read(ID_FILE_SENSOR) gmt,delgmt,gmtday,hotmac_time,utmx_hotmac,utmy_hotmac,&
            nx_hotmac,ny_hotmac,nz_hotmac,nzs_hotmac,dx_hotmac,dy_hotmac
         nx_hotmacp1=nx_hotmac+1
         ny_hotmacp1=ny_hotmac+1
         nz_hotmacp1=nz_hotmac+1
         allocate(hmz(nz_hotmacp1),hmzm(nz_hotmacp1))
         allocate(zgrnd(nx_hotmacp1,ny_hotmacp1),ztop(nx_hotmacp1,ny_hotmacp1))
         allocate(ustar_hotmac(nx_hotmacp1,ny_hotmacp1))
         allocate(u_hotmac(nx_hotmacp1,ny_hotmacp1,nz_hotmacp1),&
                  v_hotmac(nx_hotmacp1,ny_hotmacp1,nz_hotmacp1))
         allocate(x_hotmac(nx_hotmac,ny_hotmac),&
                  y_hotmac(nx_hotmac,ny_hotmac),&
                  z_hotmac(nx_hotmac,ny_hotmac,nz_hotmac))
         allocate(zo_hotmac(nx_hotmac,ny_hotmac,qutime%nsteps),&
                  ws_hotmac(nx_hotmac,ny_hotmac,nz_hotmac,qutime%nsteps),&
                  wd_hotmac(nx_hotmac,ny_hotmac,nz_hotmac,qutime%nsteps))
         zo_hotmac(:,:,:)=0.
         u_hotmac(:,:,:)=0.
         v_hotmac(:,:,:)=0.
         ws_hotmac(:,:,:,:)=0.
         wd_hotmac(:,:,:,:)=0.
         x_hotmac(:,:)=0.
         y_hotmac(:,:)=0.
         z_hotmac(:,:,:)=0.
         rewind(ID_FILE_SENSOR)
         read(ID_FILE_SENSOR) (hmz(k),k=1,nz_hotmacp1)
         read(ID_FILE_SENSOR) (hmzm(k),k=1,nz_hotmacp1)
         read(ID_FILE_SENSOR) adum !(zsoil(k),k=1,nzs_hotmac) !
         read(ID_FILE_SENSOR) adum !(zmsoil(k),k=1,nzs_hotmac)
         read(ID_FILE_SENSOR) adum !((iwater(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
         read(ID_FILE_SENSOR) adum !((ftree(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
         read(ID_FILE_SENSOR) ((zgrnd(i,j),i=1,nx_hotmacp1),j=1,ny_hotmacp1)
         read(ID_FILE_SENSOR) adum !((dzgdx(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
         read(ID_FILE_SENSOR) adum !((dzgdy(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
         read(ID_FILE_SENSOR) ((ztop(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
         do k=2,nz_hotmacp1
            do j=1,ny_hotmac
               do i=1,nx_hotmac
                  z_hotmac(i,j,k-1)=(ztop(i,j)-zgrnd(i,j))*hmzm(k)/hmz(nz_hotmacp1)
               enddo
            enddo
         enddo
         costheta=cos(-qugrid%domain_rotation*pi/180.)
         sintheta=sin(-qugrid%domain_rotation*pi/180.)
         do j=1,ny_hotmac
            do i=1,nx_hotmac
               x_hotmac(i,j)=(real(i-1)+0.5)*dx_hotmac+utmx_hotmac*1000.-qugrid%utmx
               y_hotmac(i,j)=(real(j-1)+0.5)*dy_hotmac+utmy_hotmac*1000.-qugrid%utmy
               dumx=costheta*x_hotmac(i,j)+sintheta*y_hotmac(i,j)
               dumy=-sintheta*x_hotmac(i,j)+costheta*y_hotmac(i,j)
               x_hotmac(i,j)=dumx
               y_hotmac(i,j)=dumy
            enddo
         enddo
         t_idx_hotmac=0
         tt=0
         do while(tt .lt. qutime%nsteps)
            t_idx_hotmac=t_idx_hotmac+1
            if(t_idx_hotmac .gt. skip_idx_hotmac)then
               read(ID_FILE_SENSOR,end=41111) adum !gmt,delgmt,gmtday
               read(ID_FILE_SENSOR) ((ustar_hotmac(i,j),i=1,nx_hotmacp1),j=1,ny_hotmacp1)
               read(ID_FILE_SENSOR) adum !((theta_star(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((qstar(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((bowen(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((solar(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((shortw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((uplongw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((dlongw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((sensib(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((latent(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((soilfl(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((dfxtree(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((ufxtree(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((wopt(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum !((preg(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) (((u_hotmac(i,j,k),i=1,nx_hotmacp1),j=1,ny_hotmacp1),k=1,nz_hotmacp1)
               read(ID_FILE_SENSOR) (((v_hotmac(i,j,k),i=1,nx_hotmacp1),j=1,ny_hotmacp1),k=1,nz_hotmacp1)
               tt=tt+1
               do k=2,nz_hotmacp1
                  do j=1,ny_hotmac
                     do i=1,nx_hotmac
                        dumu=(0.5*(u_hotmac(i,j,k)+u_hotmac(i+1,j,k)))
                        dumv=(0.5*(v_hotmac(i,j,k)+v_hotmac(i,j+1,k)))
                        ws_hotmac(i,j,k-1,tt)=sqrt((dumu**2.)+(dumv**2.))
                        wd_hotmac(i,j,k-1,tt)=270.-rad2deg*atan2(dumv,dumu)
                     enddo
                  enddo
               enddo
               do j=1,ny_hotmac
                  do i=1,nx_hotmac
                     zo_hotmac(i,j,tt)=z_hotmac(i,j,1)*exp(-0.4*ws_hotmac(i,j,1,tt)/ustar_hotmac(i,j))
                     if(zo_hotmac(i,j,tt) .lt. 1e-5)then
                        zo_hotmac(i,j,tt)=1e-5
                     elseif(zo_hotmac(i,j,tt) .gt. 0.45*qugrid%dz )then
                        zo_hotmac(i,j,tt)=0.45*qugrid%dz
                     endif
                  enddo
               enddo
            else
               read(ID_FILE_SENSOR,end=41111) adum ! gmt,delgmt,gmtday
               read(ID_FILE_SENSOR) adum ! ((ustar_hotmac(i,j),i=1,nx_hotmacp1),j=1,ny_hotmacp1)
               read(ID_FILE_SENSOR) adum ! ((theta_star(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((qstar(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((bowen(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((solar(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((shortw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((uplongw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((dlongw(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((sensib(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((latent(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((soilfl(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((dfxtree(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((ufxtree(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((wopt(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! ((preg(i,j),i=1,nx_hotmac),j=1,ny_hotmac)
               read(ID_FILE_SENSOR) adum ! (((u_hotmac(i,j,k),i=1,nx_hotmacp1),j=1,ny_hotmacp1),k=1,nz_hotmacp1)
               read(ID_FILE_SENSOR) adum ! (((v_hotmac(i,j,k),i=1,nx_hotmacp1),j=1,ny_hotmacp1),k=1,nz_hotmacp1)
            endif
            read(ID_FILE_SENSOR) adum ! (((w_hotmac(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((thet(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((wvap(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((qsq(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((q2l(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((tsoil(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,ksmax)
            read(ID_FILE_SENSOR) adum ! (((avthet(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((qliq(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((eddy(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((black(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((edkh(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((hedyx(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((hedyy(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((ratio(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
            read(ID_FILE_SENSOR) adum ! (((radcl(i,j,k),i=1,nx_hotmac),j=1,ny_hotmac),k=1,nz_hotmac)
         enddo
41111    continue
         close(ID_FILE_SENSOR)	!closing the profile/sensor data file
         if(tt .eq. 0)then
            print*,'Error: too many HOTMAC time steps were skipped.  The HOTMAC simulation only had',&
               t_idx_hotmac,'time steps.'
         elseif(tt .lt. qutime%nsteps)then
            print*,'Error: there were an insufficient number of HOTMAC time steps to yield all',qutime%nsteps,&
               'QUIC-URB time steps.  The last HOTMAC time step will be used for all subsequent QUIC-URB time steps.'
            t_idx_hotmac=tt
            do tt=t_idx_hotmac+1,qutime%nsteps
               ws_hotmac(:,:,:,tt)=ws_hotmac(:,:,:,t_idx_hotmac)
               wd_hotmac(:,:,:,tt)=wd_hotmac(:,:,:,t_idx_hotmac)
               zo_hotmac(:,:,tt)=zo_hotmac(:,:,t_idx_hotmac)
            enddo
         endif
         windProfile%num_sites=0
         allocate(x_idx(nx_hotmac*ny_hotmac),y_idx(nx_hotmac*ny_hotmac))
         do jj=1,ny_hotmac
            do ii=1,nx_hotmac
               if(x_hotmac(ii,jj) .ge. -meso_search_distance .and. &
                     x_hotmac(ii,jj) .le. (qugrid%nx-1)*qugrid%dx+meso_search_distance .and. &
                     y_hotmac(ii,jj) .ge. -meso_search_distance .and. &
                     y_hotmac(ii,jj) .le. (qugrid%ny-1)*qugrid%dy+meso_search_distance)then
                  windProfile%num_sites=windProfile%num_sites+1
                  x_idx(windProfile%num_sites)=ii
                  y_idx(windProfile%num_sites)=jj
               endif
            enddo
         enddo
         if(windProfile%num_sites .eq. 0)then
            windProfile%num_sites=1
            dist1=0.5*((x_hotmac(nx_hotmac,1)-x_hotmac(1,1))+(x_hotmac(1,ny_hotmac)-x_hotmac(1,1)))
            dcx=0.5*qugrid%dx*(qugrid%nx-1)
            dcy=0.5*qugrid%dy*(qugrid%ny-1)
            do jj=1,ny_hotmac
               do ii=1,nx_hotmac
                  dist2=sqrt(((x_hotmac(ii,jj)-dcx)**2.)+((y_hotmac(ii,jj)-dcy)**2.))
                  if(dist2 .le. dist1)then
                     dist1=dist2
                     x_idx(windProfile%num_sites)=ii
                     y_idx(windProfile%num_sites)=jj
                  endif
               enddo
            enddo
            x_hotmac(x_idx(1),y_idx(1))=dcx
            y_hotmac(x_idx(1),y_idx(1))=dcy
         endif
         windProfile%num_vert_points=nz_hotmac
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
         allocate(windProfile%site%wm(windProfile%num_sites,qugrid%nx,qugrid%ny),  &
            windProfile%site%wms(windProfile%num_sites,qugrid%nx,qugrid%ny))
         allocate(windProfile%site%nz_data(windProfile%num_sites,qutime%nsteps))
         allocate(windProfile%site%z_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%ws_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         allocate(windProfile%site%wd_data(windProfile%num_sites,qutime%nsteps,windProfile%num_vert_points))
         windProfile%site%blayer_flag=4
         windProfile%site%u_prof = 0
         windProfile%site%v_prof = 0
         windProfile%site%xcoord = 0
         windProfile%site%ycoord = 0
         windProfile%site%nz_data=nz_hotmac

         do kk=1,windProfile%num_sites
            windProfile%site%xcoord(kk)=x_hotmac(x_idx(kk),y_idx(kk)) !x coordinate of site location (meters)
            windProfile%site%ycoord(kk)=y_hotmac(x_idx(kk),y_idx(kk)) !y coordinate of site location (meters)
            do tt=1,qutime%nsteps
               windProfile%site%pp(kk,tt)=min(zo_hotmac(x_idx(kk),y_idx(kk),tt),0.45*qugrid%dz)		!if blayer = 2 windProfile%site%pp = exp else windProfile%site%pp = zo
               do ii=1,windProfile%site%nz_data(kk,tt)
                  windProfile%site%z_data(kk,tt,ii)=z_hotmac(x_idx(kk),y_idx(kk),ii)
                  windProfile%site%ws_data(kk,tt,ii)=ws_hotmac(x_idx(kk),y_idx(kk),ii,tt)
! MAN 02/05/2007 Domain Rotation
                  windProfile%site%wd_data(kk,tt,ii)=wd_hotmac(x_idx(kk),y_idx(kk),ii,tt)-qugrid%domain_rotation
                  if(windProfile%site%wd_data(kk,tt,ii).lt. 0.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)+360.
                  elseif(windProfile%site%wd_data(kk,tt,ii).ge. 360.)then
                       windProfile%site%wd_data(kk,tt,ii)=windProfile%site%wd_data(kk,tt,ii)-360.
                  endif
! end MAN 02/05/2007
               enddo
            enddo
         enddo    !kk = num_sites
         deallocate(hmz,hmzm,ustar_hotmac,zgrnd,ztop)
         deallocate(x_hotmac,y_hotmac,z_hotmac)
         deallocate(zo_hotmac,ws_hotmac,wd_hotmac)
         deallocate(x_idx,y_idx)
         return
      end
