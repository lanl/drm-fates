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
      subroutine upwind(ibuild)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
         use constants
         use bld_module
         use flags_module
         use grid_module
         use winds_module
         
         implicit none
			
			integer,intent(IN) :: ibuild
         integer perpendicular_flag,ns_flag,i,j,k
         integer upIstart,upIstop,upJstart,upJstop
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real xf1,yf1,xf2,yf2,tol,ynorm,lfcoeff
         real zf,x_u,y_u,x_v,y_v,x_w,y_w
         real xs_u,xs_v,xs_w,xv_u,xv_v,xv_w,xrz_u,xrz_v
         real urot,vrot,uhrot,vhrot,vel_mag
         real vortex_height,build_width,retarding_factor
         real length_factor,height_factor,rz_end,retarding_height,eff_height
         real totalLength,perpendicularDir,bld%gamma_eff
         integer ktop,kbottom,iface,ivert
         real, allocatable :: LfFace(:),LengthFace(:)
         
         if(bldgeometry(ibuild) .eq. 4)then
            eff_height=0.8*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild)
         else
            eff_height=bld%Ht(ibuild)
         endif
         if(bldgeometry(ibuild) .eq. 6)then
            xco = bld%cx(ibuild)
            yco = bld%cy(ibuild)
         else
            xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
            yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         endif
         ! find upwind direction and determine the type of flow regime
         uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         upwind_rel=upwind_dir-bld%gamma(ibuild)
         uhrot=uo_h*cos(bld%gamma(ibuild))+vo_h*sin(bld%gamma(ibuild))
         vhrot=-uo_h*sin(bld%gamma(ibuild))+vo_h*cos(bld%gamma(ibuild))
         vel_mag=sqrt((uo_h**2.)+(vo_h**2.))
         tol=bld%params%upwindAngleRange*pi/180.
         
         ! retarding_factor=0.4
         retarding_factor=0.4
         length_factor=0.4
         height_factor=0.6
         if(bld%flag%upwind .eq. 1)then
            lfcoeff=2
         else
            lfcoeff=1.5
         endif
         if(upwind_rel .gt. pi)upwind_rel=upwind_rel-2*pi
         if(upwind_rel .le. -pi)upwind_rel=upwind_rel+2*pi
         if(bldgeometry(ibuild) .eq. 6)then
            allocate(LfFace(bld%stopidx(ibuild)-bld%startidx(ibuild)),LengthFace(bld%stopidx(ibuild)-bld%startidx(ibuild)))
            iface=0
            do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
               x1=0.5*(bld%x(ivert)+bld%x(ivert+1))
               y1=0.5*(bld%y(ivert)+bld%y(ivert+1))
               xf1=(bld%x(ivert)-x1)*cos(upwind_dir)+(bld%y(ivert)-y1)*sin(upwind_dir)
               yf1=-(bld%x(ivert)-x1)*sin(upwind_dir)+(bld%y(ivert)-y1)*cos(upwind_dir)
               xf2=(bld%x(ivert+1)-x1)*cos(upwind_dir)+(bld%y(ivert+1)-y1)*sin(upwind_dir)
               yf2=-(bld%x(ivert+1)-x1)*sin(upwind_dir)+(bld%y(ivert+1)-y1)*cos(upwind_dir)
               upwind_rel=atan2(yf2-yf1,xf2-xf1)+0.5*pi
               if(upwind_rel .gt. pi)upwind_rel=upwind_rel-2*pi
               ! print*,upwind_rel,pi-tol
               if(abs(upwind_rel) .gt. pi-tol)then
                  perpendicularDir=atan2(bld%y(ivert+1)-bld%y(ivert),bld%x(ivert+1)-bld%x(ivert))+0.5*pi
                  if(perpendicularDir.le.-pi)perpendicularDir=perpendicularDir+2*pi
                  if(abs(perpendicularDir) .ge. 0.25*pi .and. abs(perpendicularDir) .le. 0.75*pi)then
                     ns_flag=1
                  else
                     ns_flag=0
                  endif
                  bld%gamma_eff=perpendicularDir
                  if(bld%gamma_eff .ge. 0.75*pi)then
                     bld%gamma_eff=bld%gamma_eff-pi
                  elseif(bld%gamma_eff .ge. 0.25*pi)then
                     bld%gamma_eff=bld%gamma_eff-0.5*pi
                  elseif(bld%gamma_eff .lt. -0.75*pi)then
                     bld%gamma_eff=bld%gamma_eff+pi
                  elseif(bld%gamma_eff .lt. -0.25*pi)then
                     bld%gamma_eff=bld%gamma_eff+0.5*pi
                  endif
                  uhrot=uo_h*cos(bld%gamma_eff)+vo_h*sin(bld%gamma_eff)
                  vhrot=-uo_h*sin(bld%gamma_eff)+vo_h*cos(bld%gamma_eff)
                  iface=iface+1
                  LengthFace(iface)=sqrt(((xf2-xf1)**2.)+((yf2-yf1)**2.))
                  LfFace(iface)=abs(lfcoeff*LengthFace(iface)*cos(upwind_rel)/(1+0.8*LengthFace(iface)/eff_height))
                  if(bld%flag%upwind .eq. 3)then
                     vortex_height=min(LengthFace(iface),eff_height)
                     retarding_height=min(LengthFace(iface),eff_height)
                     ! retarding_height=eff_height
                  else
                     vortex_height=eff_height
                     retarding_height=eff_height
                  endif
                  ! MAN 07/25/2008 stretched vertical grid
                  do k=2,bld%kstart(ibuild)
                     kbottom=k
                     if(bld%zfo(ibuild) .le. zm(k))exit
                  enddo
                  do k=bld%kstart(ibuild),qugrid%nz-1
                     ktop=k
                     if(height_factor*retarding_height+bld%zfo_actual(ibuild) .le. z(k))exit
                  enddo
                  upIstart=max(nint(min(bld%x(ivert),bld%x(ivert+1))/qugrid%dx)-nint(1.5*LfFace(iface)/qugrid%dx),2)
                  upIstop=min(nint(max(bld%x(ivert),bld%x(ivert+1))/qugrid%dx)+nint(1.5*LfFace(iface)/qugrid%dx),qugrid%nx-1)
                  upJstart=max(nint(min(bld%y(ivert),bld%y(ivert+1))/qugrid%dy)-nint(1.5*LfFace(iface)/qugrid%dy),2)
                  upJstop=min(nint(max(bld%y(ivert),bld%y(ivert+1))/qugrid%dy)+nint(1.5*LfFace(iface)/qugrid%dy),qugrid%ny-1)
                  ynorm=abs(yf2)
                  do k=kbottom,ktop
                     zf=zm(k)-bld%zfo(ibuild)
                     do j=upJstart,upJstop
                        do i=upIstart,upIstop
                           x_u=((real(i)-1)*qugrid%dx-x1)*cos(upwind_dir)+ &
                                        ((real(j)-0.5)*qugrid%dy-y1)*sin(upwind_dir)
                           y_u=-((real(i)-1)*qugrid%dx-x1)*sin(upwind_dir)+ &
                                        ((real(j)-0.5)*qugrid%dy-y1)*cos(upwind_dir)
                           x_v=((real(i)-0.5)*qugrid%dx-x1)*cos(upwind_dir)+ &
                                        ((real(j)-1)*qugrid%dy-y1)*sin(upwind_dir)
                           y_v=-((real(i)-0.5)*qugrid%dx-x1)*sin(upwind_dir)+	&
                                        ((real(j)-1)*qugrid%dy-y1)*cos(upwind_dir)
                           x_w=((real(i)-0.5)*qugrid%dx-x1)*cos(upwind_dir)+ &
                                        ((real(j)-0.5)*qugrid%dy-y1)*sin(upwind_dir)
                           y_w=-((real(i)-0.5)*qugrid%dx-x1)*sin(upwind_dir)+	&
                                        ((real(j)-0.5)*qugrid%dy-y1)*cos(upwind_dir)
!u values
                           if(abs(y_u) .le. ynorm .and. height_factor*vortex_height .gt. zf)then
                              xs_u=((xf2-xf1)/(yf2-yf1))*(y_u-yf1)+xf1
                              xv_u=-LfFace(iface)*sqrt((1-((y_u/ynorm)**2.))*(1-((zf/(height_factor*vortex_height))**2.)))
                              xrz_u=-LfFace(iface)*sqrt((1-((y_u/ynorm)**2.))*(1-((zf/(height_factor*retarding_height))**2.)))
                              if(zf .gt. height_factor*vortex_height)then
                                 rz_end=0.
                              else
                                 rz_end=length_factor*xv_u
                              endif
                              if(bld%flag%upwind .eq. 1)then
                                 if(x_u-xs_u .ge. xv_u .and. x_u-xs_u .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                    quwinds%uo(i,j,k)=0.
                                 endif
                              else
                                 if(x_u-xs_u .ge. xrz_u .and. x_u-xs_u .lt. rz_end &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    if(bld%flag%upwind .eq. 3)then
                                       quwinds%uo(i,j,k)=((x_u-xs_u-xrz_u)*(retarding_factor-1.)/(rz_end-xrz_u)+1.)*quwinds%uo(i,j,k)
                                    else
                                       quwinds%uo(i,j,k)=retarding_factor*quwinds%uo(i,j,k)
                                    endif
                                    if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized U exceeds max in upwind',&
                                          quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                 endif
                                 if(x_u-xs_u .ge. length_factor*xv_u .and. x_u-xs_u .le. 0.1*qugrid%dxy &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    urot=quwinds%uo(i,j,k)*cos(bld%gamma_eff)
                                    vrot=-quwinds%uo(i,j,k)*sin(bld%gamma_eff)
                                    if(ns_flag .eq. 1)then
                                       vrot=-vhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                            *(-height_factor*sin(((pi*abs(x_u-xs_u))/(length_factor*LfFace(iface)))+0))
                                    else
                                       urot=-uhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                            *(-height_factor*sin(((pi*abs(x_u-xs_u))/(length_factor*LfFace(iface)))+0))
                                    endif
                                    quwinds%uo(i,j,k)=urot*cos(-bld%gamma_eff)+vrot*sin(-bld%gamma_eff)
                                    if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized U exceeds max in upwind',&
                                          quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                 endif
                              endif
                           endif
!v values
                           if(abs(y_v) .le. ynorm .and. height_factor*vortex_height .gt. zf)then
                              xs_v=((xf2-xf1)/(yf2-yf1))*(y_v-yf1)+xf1
                              xv_v=-LfFace(iface)*sqrt((1-((y_v/ynorm)**2.))*(1-((zf/(height_factor*vortex_height))**2.)))
                              xrz_v=-LfFace(iface)*sqrt((1-((y_v/ynorm)**2.))*(1-((zf/(height_factor*retarding_height))**2.)))
                              if(zf .ge. height_factor*vortex_height)then
                                 rz_end=0.
                              else
                                 rz_end=length_factor*xv_v
                              endif
                              if(bld%flag%upwind .eq. 1)then
                                 if(x_v-xs_v .ge. xv_v .and. x_v-xs_v .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                    quwinds%vo(i,j,k)=0.
                                 endif
                              else
                                 if(x_v-xs_v .ge. xrz_v .and. x_v-xs_v .lt. rz_end &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    if(bld%flag%upwind .eq. 3)then
                                       quwinds%vo(i,j,k)=((x_v-xs_v-xrz_v)*(retarding_factor-1.)/(rz_end-xrz_v)+1.)*quwinds%vo(i,j,k)
                                    else
                                       quwinds%vo(i,j,k)=retarding_factor*quwinds%vo(i,j,k)
                                    endif
                                    if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized V exceeds max in upwind',&
                                          quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                 endif
                                 if(x_v-xs_v .ge. length_factor*xv_v .and. x_v-xs_v .le. 0.1*qugrid%dxy &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    urot=quwinds%vo(i,j,k)*sin(bld%gamma_eff)
                                    vrot=quwinds%vo(i,j,k)*cos(bld%gamma_eff)
                                    if(ns_flag .eq. 1)then
                                       vrot=-vhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                            *(-height_factor*sin(((pi*abs(x_v-xs_v))/(length_factor*LfFace(iface)))+0))
                                    else
                                       urot=-uhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                            *(-height_factor*sin(((pi*abs(x_v-xs_v))/(length_factor*LfFace(iface)))+0))
                                    endif
                                    quwinds%vo(i,j,k)=-urot*sin(-bld%gamma_eff)+vrot*cos(-bld%gamma_eff)
                                    if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized V exceeds max in upwind',&
                                          quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                 endif
                              endif
                           endif
!w values
                           if(abs(y_w) .le. ynorm .and. height_factor*vortex_height .gt. zf)then
                              xs_w=((xf2-xf1)/(yf2-yf1))*(y_w-yf1)+xf1
                              xv_w=-LfFace(iface)*sqrt((1-((y_w/ynorm)**2.))*(1-((zf/(height_factor*vortex_height))**2.)))
                              if(bld%flag%upwind .eq. 1)then
                                 if(x_w-xs_w .ge. xv_w .and. x_w-xs_w .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                    quwinds%wo(i,j,k)=0.
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       bld%icellflag(i,j,k)=21 ! 2 is the old value
                                    endif
                                 endif
                              else
                                 if(x_w-xs_w .ge. xv_w .and. x_w-xs_w .lt. length_factor*xv_w &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    quwinds%wo(i,j,k)=retarding_factor*quwinds%wo(i,j,k)
                                    if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized W exceeds max in upwind',&
                                          quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       bld%icellflag(i,j,k)=21 ! 2 is the old value
                                    endif
                                 endif
                                 if(x_w-xs_w .ge. length_factor*xv_w .and. x_w-xs_w .le. 0.1*qugrid%dxy &
                                       .and. bld%icellflag(i,j,k) .ne. 0)then
                                    quwinds%wo(i,j,k)=-vel_mag*(0.1*cos(((pi*abs(x_w-xs_w))/(length_factor*LfFace(iface))))-0.05)
                                    if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized W exceeds max in upwind',&
                                          quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                    endif
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       bld%icellflag(i,j,k)=20 ! 2 is the old value
                                    endif
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               endif
               if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                     .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))exit
            enddo
            if(iface .gt. 0)then
               totalLength=0.
               bld%lf(ibuild)=0.
               do ivert=1,iface
                  bld%lf(ibuild)=bld%lf(ibuild)+LfFace(ivert)*LengthFace(ivert)
                  totalLength=totalLength+LengthFace(ivert)
               enddo
               bld%lf(ibuild)=bld%lf(ibuild)/totalLength
            else
               bld%lf(ibuild)=-999.0
            endif
            deallocate(LfFace,LengthFace)
         else
            !Location of corners relative to the center of the building
            x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
            y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
            x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
            y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
            x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            perpendicular_flag=0
            
            ! if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
            !    xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xw2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yw2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yf2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    perpendicular_flag=0
            ! elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
            !    xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    
            !    xf1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yf1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xf3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yf3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    perpendicular_flag=1
            ! elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
            !    xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yf2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    perpendicular_flag=0
            ! elseif(abs(upwind_rel) .le. tol)then
            !    xf1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yf1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xf3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yf3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    
            !    
            !    xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    perpendicular_flag=1
            ! elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
            !    xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xf2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yf2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    perpendicular_flag=0
            ! elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
            !    xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    
            !    
            !    xf1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yf1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xf3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yf3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    perpendicular_flag=1
            ! elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
            !    xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yf2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    perpendicular_flag=0
            ! else
            !    xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            !    yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            !    xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            !    yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            !    xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    
            !    xf1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            !    yf1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            !    xf3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            !    yf3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            !    perpendicular_flag=1
            ! endif
            
            if(upwind_rel .gt. 0.5*pi-tol .and. upwind_rel .lt. 0.5*pi+tol)then
               xf1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
               yf1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
               xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
               yf2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
               perpendicular_flag=1
               ns_flag=1
               bld%lf(ibuild)=abs(lfcoeff*bld%Lti(ibuild)*sin(upwind_rel)/(1+0.8*bld%Lti(ibuild)/eff_height))
               build_width=bld%Lti(ibuild)
            elseif(upwind_rel .gt. -tol .and. upwind_rel .lt. tol)then
               xf1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
               yf1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
               xf2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
               yf2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
               perpendicular_flag=1
               ns_flag=0
               bld%lf(ibuild)=abs(lfcoeff*bld%Wti(ibuild)*cos(upwind_rel)/(1+0.8*bld%Wti(ibuild)/eff_height))
               build_width=bld%Wti(ibuild)
            elseif(upwind_rel .gt. -0.5*pi-tol .and. upwind_rel .lt. -0.5*pi+tol)then
               xf1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
               yf1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
               xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
               yf2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
               perpendicular_flag=1
               ns_flag=1
               bld%lf(ibuild)=abs(lfcoeff*bld%Lti(ibuild)*sin(upwind_rel)/(1+0.8*bld%Lti(ibuild)/eff_height))
               build_width=bld%Lti(ibuild)
            elseif(upwind_rel .gt. pi-tol .or. upwind_rel .lt. -pi+tol)then
               xf1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
               yf1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
               xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
               yf2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
               perpendicular_flag=1
               ns_flag=0
               bld%lf(ibuild)=abs(lfcoeff*bld%Wti(ibuild)*cos(upwind_rel)/(1+0.8*bld%Wti(ibuild)/eff_height))
               build_width=bld%Wti(ibuild)
            endif
            if(perpendicular_flag .eq. 1)then
               ynorm=max(abs(yf2),abs(yf1))
               if(bld%flag%upwind .eq. 3)then
                  vortex_height=min(build_width,eff_height)
                  retarding_height=min(build_width,eff_height)
                  ! retarding_height=eff_height
               else
                  vortex_height=eff_height
                  retarding_height=eff_height
               endif
               ! MAN 07/25/2008 stretched vertical grid
               do k=2,bld%kstart(ibuild)
                  kbottom=k
                  if(bld%zfo(ibuild) .le. zm(k))exit
               enddo
               do k=bld%kstart(ibuild),qugrid%nz-1
                  ktop=k
                  if(height_factor*retarding_height+bld%zfo_actual(ibuild) .le. z(k))exit
               enddo
               ! print*,ibuild,kbottom,ktop
               ! print*,vortex_height,retarding_height
               upIstart=max(bld%istart(ibuild)-nint(1.5*bld%lf(ibuild)/qugrid%dx),2)
               upIstop=min(bld%iend(ibuild)+nint(1.5*bld%lf(ibuild)/qugrid%dx),qugrid%nx-1)
               upJstart=max(bld%jstart(ibuild)-nint(1.5*bld%lf(ibuild)/qugrid%dy),2)
               upJstop=min(bld%jend(ibuild)+nint(1.5*bld%lf(ibuild)/qugrid%dy),qugrid%ny-1)
lp003:         do k=kbottom,ktop
                  zf=zm(k)-bld%zfo(ibuild)
lp002:            do j=upJstart,upJstop
lp001:               do i=upIstart,upIstop
                        x_u=((real(i)-1)*qugrid%dx-xco)*cos(upwind_dir)+ &
                                     ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                        y_u=-((real(i)-1)*qugrid%dx-xco)*sin(upwind_dir)+ &
                                     ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)
                        x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+ &
                                     ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)
                        y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&
                                     ((real(j)-1)*qugrid%dy-yco)*cos(upwind_dir)
                        x_w=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+ &
                                     ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                        y_w=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&
                                     ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)
! u values
                        if(y_u .ge. yf1 .and. y_u .le. yf2)then
                           xs_u=((xf2-xf1)/(yf2-yf1))*(y_u-yf1)+xf1
                           
                           if(zf .gt. height_factor*vortex_height)then
                              rz_end=0.
                              xv_u=0.
                              xrz_u=0.
                           else
                           	xv_u=-bld%lf(ibuild)*sqrt((1-((y_u/ynorm)**2.))*(1-((zf/(height_factor*vortex_height))**2.)))
                              xrz_u=-bld%lf(ibuild)*sqrt((1-((y_u/ynorm)**2.))*(1-((zf/(height_factor*retarding_height))**2.)))
                              rz_end=length_factor*xv_u
                           endif
                           if(bld%flag%upwind .eq. 1)then
                              if(x_u-xs_u .ge. xv_u .and. x_u-xs_u .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                 quwinds%uo(i,j,k)=0.
                              endif
                           else
                              if(x_u-xs_u .ge. xrz_u .and. x_u-xs_u .lt. rz_end &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 if(bld%flag%upwind .eq. 3)then
                                    quwinds%uo(i,j,k)=((x_u-xs_u-xrz_u)*(retarding_factor-1.)/(rz_end-xrz_u)+1.)*quwinds%uo(i,j,k)
                                 else
                                    quwinds%uo(i,j,k)=retarding_factor*quwinds%uo(i,j,k)
                                 endif
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in upwind',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                              endif
                              if(x_u-xs_u .ge. length_factor*xv_u .and. x_u-xs_u .le. 0.1*qugrid%dxy &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 urot=quwinds%uo(i,j,k)*cos(bld%gamma(ibuild))
                                 vrot=-quwinds%uo(i,j,k)*sin(bld%gamma(ibuild))
                                 if(ns_flag .eq. 1)then
                                    vrot=-vhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                         *(-height_factor*sin(((pi*abs(x_u-xs_u))/(length_factor*bld%lf(ibuild)))+0))
                                 else
                                    urot=-uhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                         *(-height_factor*sin(((pi*abs(x_u-xs_u))/(length_factor*bld%lf(ibuild)))+0))
                                 endif
                                 quwinds%uo(i,j,k)=urot*cos(-bld%gamma(ibuild))+vrot*sin(-bld%gamma(ibuild))
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in upwind',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                              endif
                           endif
                        endif
!v values
                        if(y_v .ge. yf1 .and. y_v .le. yf2)then
                           xs_v=((xf2-xf1)/(yf2-yf1))*(y_v-yf1)+xf1
                           
                           if(zf .ge. height_factor*vortex_height)then
                              rz_end=0.
                              xv_v=0.
                              xrz_v=0.
                           else
                              xv_v=-bld%lf(ibuild)*sqrt((1.-((y_v/ynorm)**2.))*(1.-((zf/(height_factor*vortex_height))**2.)))

                              xrz_v=-bld%lf(ibuild)*sqrt((1.-((y_v/ynorm)**2.))*(1.-((zf/(height_factor*retarding_height))**2.)))
                              rz_end=length_factor*xv_v
                           endif
                           if(bld%flag%upwind .eq. 1)then
                              if(x_v-xs_v .ge. xv_v .and. x_v-xs_v .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                 quwinds%vo(i,j,k)=0.
                              endif
                           else
                              if(x_v-xs_v .ge. xrz_v .and. x_v-xs_v .lt. rz_end &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 if(bld%flag%upwind .eq. 3)then
                                    quwinds%vo(i,j,k)=((x_v-xs_v-xrz_v)*(retarding_factor-1.)/(rz_end-xrz_v)+1.)*quwinds%vo(i,j,k)
                                 else
                                    quwinds%vo(i,j,k)=retarding_factor*quwinds%vo(i,j,k)
                                 endif
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in upwind',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                              endif
                              if(x_v-xs_v .ge. length_factor*xv_v .and. x_v-xs_v .le. 0.1*qugrid%dxy &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 urot=quwinds%vo(i,j,k)*sin(bld%gamma(ibuild))
                                 vrot=quwinds%vo(i,j,k)*cos(bld%gamma(ibuild))
                                 if(ns_flag .eq. 1)then
                                    vrot=-vhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                         *(-height_factor*sin(((pi*abs(x_v-xs_v))/(length_factor*bld%lf(ibuild)))+0))
                                 else
                                    urot=-uhrot*(-height_factor*cos(((pi*zf)/(0.5*vortex_height)))+0.05)   &
                                         *(-height_factor*sin(((pi*abs(x_v-xs_v))/(length_factor*bld%lf(ibuild)))+0))
                                 endif
                                 quwinds%vo(i,j,k)=-urot*sin(-bld%gamma(ibuild))+vrot*cos(-bld%gamma(ibuild))
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in upwind',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                              endif
                           endif
                        endif
!w values
                        if(y_w .ge. yf1 .and. y_w .le. yf2)then
                           xs_w=((xf2-xf1)/(yf2-yf1))*(y_w-yf1)+xf1
                           
                           if(zf .ge. height_factor*vortex_height)then
							         xv_w=0.
                           else
                              xv_w=-bld%lf(ibuild)*sqrt((1-((y_w/ynorm)**2.))*(1-((zf/(height_factor*vortex_height))**2.)))

                           endif
                           if(bld%flag%upwind .eq. 1)then
                              if(x_w-xs_w .ge. xv_w .and. x_w-xs_w .le. 0.1*qugrid%dxy .and. bld%icellflag(i,j,k) .ne. 0)then
                                 quwinds%wo(i,j,k)=0.
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    bld%icellflag(i,j,k)=21 ! 2 is the old value
                                 endif
                              endif
                           else
                              if(x_w-xs_w .ge. xv_w .and. x_w-xs_w .lt. length_factor*xv_w &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 quwinds%wo(i,j,k)=retarding_factor*quwinds%wo(i,j,k)
                                 if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized W exceeds max in upwind',&
                                       quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    bld%icellflag(i,j,k)=21 ! 2 is the old value
                                 endif
                              endif
                              if(x_w-xs_w .ge. length_factor*xv_w .and. x_w-xs_w .le. 0. &
                                    .and. bld%icellflag(i,j,k) .ne. 0)then
                                 quwinds%wo(i,j,k)=-vel_mag*(0.1*cos(((pi*abs(x_w-xs_w))/(length_factor*bld%lf(ibuild))))-0.05)
                                 if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized W exceeds max in upwind',&
                                       quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    bld%icellflag(i,j,k)=20 ! 2 is the old value
                                 endif
                              endif
                           endif
                        endif
                     enddo   lp001      
                  enddo   lp002      
               enddo   lp003
            else
               bld%lf(ibuild)=-999.0
            endif
         endif
         return
      end
