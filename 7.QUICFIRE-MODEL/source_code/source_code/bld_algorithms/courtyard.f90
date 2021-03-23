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

      subroutine courtyard(ibuild)
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
			integer circle_flag,perpendicular_flag,k_top
         real uo_h,vo_h,upwind_dir,xco,yco,tol
         real zb,x_u,y_u,x_v,y_v,ymin,ymax
         real dNu,dNv,farwake_exponent,farwake_factor,farwake_velocity
         real thetamax,thetamin,thetai
         real xnorm,ynorm,radius,xnorm_bisect,cav_fac,wake_fac
         real x1,y1,x2,y2,x3,y3,x4,y4
         real xw1,yw1,xw2,yw2,xw3,yw3,xw4,yw4
         real beta,upwind_rel,upwind_rel_norm
         real theta_u,theta_v,theta_corner,velmag,usign,xmax
         real roof_ratio,bld%zfo_roof,court_frac,ro_u,ro_v,ri_u,ri_v,rr_u,rr_v,r_u,r_v
         real indoor_vel_frac,Lt_roof,Wt_roof
         real xu_rel,xv_rel,recirc_frac,courtyard_Lr,eff_height
         integer ivert,ipoly,iu,ju,iv,jv,uwakeflag,vwakeflag,x_idx,y_idx,x_idx_min,ktop,kbottom
         real xwall,xwallu,xwallv,xu,yu,xv,yv,xp,yp,ynormm,ynormp,xc,yc
         real epsilon
			integer i,j,k
         
         epsilon = 10e-10
         
         recirc_frac=0.5
         farwake_exponent=1.5
         farwake_factor=3
         if(bld%geometry(ibuild) .eq. 6)then
            xco=bld%cx(ibuild)
            yco=bld%cy(ibuild)
         else
            xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
            yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         endif
         uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         indoor_vel_frac=0.75
! find upwind direction and determine the type of flow regime
         upwind_dir=atan2(vo_h,uo_h)
         select case(bld%geometry(ibuild))
            case(4) !rectangular stadium
               upwind_rel=upwind_dir-bld%gamma(ibuild)
               if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
               if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
               upwind_rel_norm=upwind_rel+0.5*pi
               if(upwind_rel_norm .gt. pi)upwind_rel_norm=upwind_rel_norm-2*pi
               beta=abs(atan2(bld%Lti(ibuild),bld%Wti(ibuild)))
               if(bld%flag%wake .eq. 2)then
                  cav_fac=1.1
                  wake_fac=0.1
               else
                  cav_fac=1.
                  wake_fac=0.
               endif
               tol=0.01*pi/180.
               farwake_exponent=1.5
               farwake_factor=3
               !Location of corners relative to the center of the building
               x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
               y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
               x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
               y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
               x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
               y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
               x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
               y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
               if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
                  xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  xw4=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw4=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  perpendicular_flag=0
                  usign=1
               elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
                  xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  xmax=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  perpendicular_flag=1
               elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
                  xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  xw4=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw4=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  perpendicular_flag=0
                  usign=-1
               elseif(abs(upwind_rel) .le. tol)then
                  xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  xmax=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  perpendicular_flag=1
               elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
                  xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  xw2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  xw4=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw4=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  perpendicular_flag=0
                  usign=1
               elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
                  xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  xmax=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  perpendicular_flag=1
               elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
                  xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                  xw4=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                  yw4=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                  perpendicular_flag=0
                  usign=-1
               else
                  xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                  yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                  xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                  yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                  xmax=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                  perpendicular_flag=1
               endif
               if(perpendicular_flag .eq. 0)then
                  theta_corner=atan2(yw2,xw2)
                  velmag=sqrt((uo_h**2)+(vo_h**2))
               endif
               courtyard_Lr=7.6*bld%Ht(ibuild)
               if(bld%roof(ibuild) .eq. 0)then
                  do k=bld%kstart(ibuild),bld%kend(ibuild)  !erp 7/23/03
                     zb=zm(k)-bld%zfo(ibuild)
                     do j=bld%jstart(ibuild),bld%jend(ibuild)
                        do i=bld%istart(ibuild),bld%iend(ibuild)
                           x_u=((real(i)-1)*qugrid%dx-xco)*cos(upwind_dir)+		&
                                        ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                           y_u=-((real(i)-1)*qugrid%dx-xco)*sin(upwind_dir)+		&
                                        ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)
                           x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+	&
                                        ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)
                           y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&
                                        ((real(j)-1)*qugrid%dy-yco)*cos(upwind_dir)
                           theta_u=atan2(y_u,x_u)
                           theta_v=atan2(y_v,x_v)
! u values           
                           if(y_u.ge.0.)then
                              ynorm=yw1
                           else
                              ynorm=yw3
                           endif         
                           if(perpendicular_flag .gt. 0)then
                              xnorm=xw1
                           else
                              if(y_u.ge.yw2)then
                                 xnorm=((xw2-xw1)/(yw2-yw1))*(y_u-yw1)+xw1
                              else
                                 xnorm=((xw3-xw2)/(yw3-yw2))*(y_u-yw2)+xw2
                              endif
                              if(y_u.ge.yw4)then
                                 xmax=((xw4-xw1)/(yw4-yw1))*(y_u-yw1)+xw1
                              else
                                 xmax=((xw3-xw4)/(yw3-yw4))*(y_u-yw4)+xw4
                              endif
                           endif
! Far wake           
                           if(y_u .ge. yw3 .and. y_u .le. yw1 .and. k .le. bld%kend(ibuild) &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i-1,j,k) .ne. 0
                              if(abs(ynorm) > epsilon .and. bld%Ht(ibuild) > epsilon)then
                                 dNu=sqrt(exp(-10.*(y_u/ynorm)**2.)*   &
                                          (1.-((zb)/bld%Ht(ibuild))**2.)*courtyard_Lr**2)
                              else
                              	dNu=0.
                              endif           
                              if(x_u .gt. xnorm+dNu .and. x_u-xnorm .le. &
                                    farwake_factor*dNu .and. x_u .lt. xmax)then
                                 farwake_velocity=quwinds%uo(i,j,k)*&
                                    (1.-(dNu/(x_u-xnorm+wake_fac*dNu))**(farwake_exponent))
                                 if(abs(farwake_velocity).lt.abs(quwinds%uo(i,j,k)))then
                                    quwinds%uo(i,j,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       if(bld%icellflag(i,j,k) .ne. 0)then
                                          bld%icellflag(i,j,k)=5
                                       endif
                                    endif
                                 endif
                              endif
                           endif
! Cavity                      
                           if(y_u .ge. yw3 .and. y_u .le. yw1 .and. k .le. bld%kend(ibuild) &
                                 .and. x_u .gt. xnorm .and. x_u .le. dNu+xnorm &
                                 .and. x_u .lt. xmax &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i-1,j,k) .ne. 0
                              quwinds%uo(i,j,k)=-uo_h*min(sqrt(1.-abs(y_u/ynorm)),1.)*(1.-(x_u-xnorm)/(cav_fac*dNu))**2.
                              if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 write(325,*)'Parameterized U exceeds max in courtyard',&
                                    quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                              endif
                              quwinds%wo(i,j,k)=0.
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                 if(bld%icellflag(i,j,k) .ne. 0)then
                                    bld%icellflag(i,j,k)=4
                                 endif
                              endif
                           endif
! v values           
                           if(y_v.ge.0.)then
                              ynorm=yw1
                           else
                              ynorm=yw3
                           endif
                           if(perpendicular_flag .gt. 0)then
                              xnorm=xw1
                           else
                              if(y_v.ge.yw2)then
                                 xnorm=((xw2-xw1)/(yw2-yw1))*(y_v-yw1)+xw1
                              else
                                 xnorm=((xw3-xw2)/(yw3-yw2))*(y_v-yw2)+xw2
                              endif
                              if(y_v.ge.yw4)then
                                 xmax=((xw4-xw1)/(yw4-yw1))*(y_v-yw1)+xw1
                              else
                                 xmax=((xw3-xw4)/(yw3-yw4))*(y_v-yw4)+xw4
                              endif
                           endif
! Far wake           
                           if(y_v .ge. yw3 .and. y_v .le. yw1 .and. k .le. bld%kend(ibuild) &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i,j-1,k) .ne. 0
                              if(abs(ynorm) > epsilon .and. bld%Ht(ibuild) > epsilon)then
                              	dNv=sqrt(exp(-10.*(y_v/ynorm)**2.)*   &
                                          (1.-((zb)/bld%Ht(ibuild))**2.)*courtyard_Lr**2)
                              else
                              	dNv=0.
                              endif
                              if(x_v .gt. xnorm+dNv .and. x_v-xnorm .le. &
                                    farwake_factor*dNv .and. x_v .lt. xmax)then
                                 farwake_velocity=quwinds%vo(i,j,k)*&
                                    (1.-(dNv/(x_v-xnorm+wake_fac*dNv))**(farwake_exponent))
                                 if(abs(farwake_velocity).lt.abs(quwinds%vo(i,j,k)))then
                                    quwinds%vo(i,j,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       if(bld%icellflag(i,j,k) .ne. 0)then
                                          bld%icellflag(i,j,k)=5
                                       endif
                                    endif
                                 endif
                              endif
                           endif
! Cavity             
                           if(y_v .ge. yw3 .and. y_v .le. yw1 .and. k .le. bld%kend(ibuild) &
                                 .and. x_v .gt. xnorm .and. x_v .le. dNv+xnorm &
                                 .and. x_v .lt. xmax &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i,j-1,k) .ne. 0
                              quwinds%vo(i,j,k)=-vo_h*min(sqrt(1.-abs(y_v/ynorm)),1.)*(1.-(x_v-xnorm)/(cav_fac*dNv))**2.
                              if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 write(325,*)'Parameterized V exceeds max in courtyard',&
                                    quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                              endif
                              quwinds%wo(i,j,k)=0.
                              if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                 if(bld%icellflag(i,j,k) .ne. 0)then
                                    bld%icellflag(i,j,k)=4
                                 endif
                              endif
                           endif
                        enddo         
                     enddo         
                  enddo
               else
                  roof_ratio=0.8
                  bld%zfo_roof=(bld%Ht(ibuild)-bld%zfo_actual(ibuild))*roof_ratio
                  court_frac=bld%wall(ibuild)
                  do j=bld%jstart(ibuild),bld%jend(ibuild)
                     do i=bld%istart(ibuild),bld%iend(ibuild)
                        x_u=((real(i)-1)*qugrid%dx-xco)*cos(bld%gamma(ibuild))+		&
                                     ((real(j)-0.5)*qugrid%dy-yco)*sin(bld%gamma(ibuild))
                        y_u=-((real(i)-1)*qugrid%dx-xco)*sin(bld%gamma(ibuild))+		&
                                     ((real(j)-0.5)*qugrid%dy-yco)*cos(bld%gamma(ibuild))
                        x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(bld%gamma(ibuild))+	&
                                     ((real(j)-1)*qugrid%dy-yco)*sin(bld%gamma(ibuild))
                        y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(bld%gamma(ibuild))+	&
                                     ((real(j)-1)*qugrid%dy-yco)*cos(bld%gamma(ibuild))
                        xu_rel=((real(i)-1)*qugrid%dx-xco)*cos(upwind_rel)+		&
                                     ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                        xv_rel=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_rel)+	&
                                     ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)
                        do k=bld%kstart(ibuild),bld%kend(ibuild)
                           zb=zm(k)-bld%zfo(ibuild)
                           if(zb .le. bld%zfo_roof)then
                              if(abs(x_u) .lt. bld%Lt(ibuild) &
                                    .and. abs(y_u) .lt. bld%Wt(ibuild) &
                                    .and. (abs(x_u) .gt. recirc_frac*(bld%Lt(ibuild)-court_frac) &
                                    .or. abs(y_u) .gt. recirc_frac*(bld%Wt(ibuild)-court_frac)) &
                                    .and. bld%icellflag(i,j,k) .ne. 0  &
                                    .and. xu_rel .lt. 0)then !.and. bld%icellflag(i-1,j,k) .ne. 0
                                 quwinds%uo(i,j,k)=-indoor_vel_frac*uo_h*(1-zb/bld%zfo_roof)
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in partial roof courtyard',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              if(abs(x_v) .lt. bld%Lt(ibuild) &
                                    .and. abs(y_v) .lt. bld%Wt(ibuild) &
                                    .and. (abs(x_v) .gt. recirc_frac*(bld%Lt(ibuild)-court_frac) &
                                    .or. abs(y_v) .gt. recirc_frac*(bld%Wt(ibuild)-court_frac)) &
                                    .and. bld%icellflag(i,j,k) .ne. 0 &
                                    .and. xv_rel .lt. 0)then !.and. bld%icellflag(i,j-1,k) .ne. 0
                                 quwinds%vo(i,j,k)=-indoor_vel_frac*vo_h*(1-zb/bld%zfo_roof)
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in partial roof courtyard',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           else
                              Lt_roof=bld%Lt(ibuild)-court_frac*((zb-bld%zfo_roof)/(0.25*bld%zfo_roof))**2.
                              Wt_roof=bld%Wt(ibuild)-court_frac*((zb-bld%zfo_roof)/(0.25*bld%zfo_roof))**2.
                              if(abs(x_u) .lt. Lt_roof &
                                    .and. abs(y_u) .lt. Wt_roof &
                                    .and. (abs(x_u) .gt. recirc_frac*(bld%Lt(ibuild)-court_frac) &
                                    .or. abs(y_u) .gt. recirc_frac*(bld%Wt(ibuild)-court_frac)) &
                                    .and. bld%icellflag(i,j,k) .ne. 0 .and. bld%icellflag(i-1,j,k) .ne. 0 &
                                    .and. xu_rel .lt. 0)then
                                 quwinds%uo(i,j,k)=indoor_vel_frac*uo_h*(4*abs(zb-bld%zfo_roof)/bld%zfo_roof)
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in partial roof courtyard',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              if(abs(x_v) .lt. Lt_roof &
                                    .and. abs(y_v) .lt. Wt_roof &
                                    .and. (abs(x_v) .gt. recirc_frac*(bld%Lt(ibuild)-court_frac) &
                                    .or. abs(y_v) .gt. recirc_frac*(bld%Wt(ibuild)-court_frac)) &
                                    .and. bld%icellflag(i,j,k) .ne. 0  &
                                    .and. xv_rel .lt. 0)then !.and. bld%icellflag(i,j-1,k) .ne. 0
                                 quwinds%vo(i,j,k)=indoor_vel_frac*vo_h*(4*abs(zb-bld%zfo_roof)/bld%zfo_roof)
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in partial roof courtyard',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               endif
            case(5) ! Elliptical stadium
               tol=0.01*min(qugrid%dx,qugrid%dy)
               if(abs(bld%aa(ibuild)-bld%bb(ibuild)) .lt. tol)then
                  circle_flag=1
               else
                  circle_flag=0
               endif
               if(bld%flag%wake .eq. 2)then
                  cav_fac=1.1
                  wake_fac=0.1
               else
                  cav_fac=1.
                  wake_fac=0.
               endif
               if(circle_flag .eq. 1)then
                  thetamin=-0.5*pi
                  thetamax=0.5*pi
                  ymax=bld%Lt(ibuild)
                  ymin=-bld%Lt(ibuild)
               else
                  y_u=0.
                  y_v=0.
                  do i=1,180
                     thetai=real(180-i)*pi/180.
                     y_u=radius(bld%aa(ibuild),bld%bb(ibuild),thetai,&
                         bld%gamma(ibuild)-upwind_dir)*sin(thetai)
                     if(y_u .lt. y_v)then
                        exit
                     endif
                     y_v=y_u
                  enddo
                  thetamax=thetai+pi/180.
                  thetamin=thetamax-pi
                  ymax=y_v
                  ymin=-y_v
               endif
               ynorm=max(ymax,-ymin)
               courtyard_Lr=7.6*bld%Ht(ibuild)
               if(bld%roof(ibuild) .eq. 0)then
                  eff_height=0.9*bld%Ht(ibuild)
                  do k=bld%kstart(ibuild),qugrid%nz-1
                     k_top=k
                     if(eff_height .lt. zm(k+1))exit
                  enddo
                  do k=bld%kstart(ibuild),k_top
                     zb=zm(k)-bld%zfo(ibuild)
                     do j=bld%jstart(ibuild),bld%jend(ibuild)
                        do i=bld%istart(ibuild),bld%iend(ibuild)
                           x_u=((real(i)-1)*qugrid%dx-xco)*cos(upwind_dir)+		&
                                        ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                           y_u=-((real(i)-1)*qugrid%dx-xco)*sin(upwind_dir)+		&
                                        ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)
                           x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+	&
                                        ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)
                           y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&
                                        ((real(j)-1)*qugrid%dy-yco)*cos(upwind_dir)
! u values           
                           if(y_u .gt. ymin .and. y_u .lt. ymax &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i-1,j,k) .ne. 0
                              if(circle_flag .eq. 1)then 
                                 xmax=sqrt((bld%Lt(ibuild)**2.)-(y_u**2.))
                                 xnorm=-xmax
                              else
                                 xnorm=-xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                                       bld%gamma(ibuild)-upwind_dir,-y_u,thetamin,thetamax,qugrid%dxy)
                                 xmax=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                                       bld%gamma(ibuild)-upwind_dir,y_u,thetamin,thetamax,qugrid%dxy)
                              endif
                              if(abs(ynorm) > epsilon .and. eff_height > epsilon)then
                                 dNu=sqrt(exp(-10.*(y_u/ynorm)**2.)*   &
                                          (1.-((zb)/eff_height)**2.)*courtyard_Lr**2)
                              else
                              	dNu=0.
                              endif
! far wake                    
                              if(x_u .gt. xnorm+dNu .and. x_u-xnorm .le. &
                                    farwake_factor*dNu .and. x_u .lt. xmax-qugrid%dxy)then
                                 farwake_velocity=quwinds%uo(i,j,k)*&
                                    (1.-(dNu/(x_u-xnorm+wake_fac*dNu))**(farwake_exponent))
                                 if(abs(farwake_velocity).lt.abs(quwinds%uo(i,j,k)))then
                                    quwinds%uo(i,j,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       if(bld%icellflag(i,j,k).ne.0)then
                                          bld%icellflag(i,j,k)=5
                                       endif
                                    endif
                                 endif
                              endif
! cavity                      
                              if(x_u .gt. xnorm +qugrid%dxy .and. x_u .le. dNu+xnorm .and. x_u .lt. xmax-qugrid%dxy)then
                                 quwinds%uo(i,j,k)=-uo_h*min(sqrt(1.-abs(y_u/ynorm)),1.)*(1.-(x_u-xnorm)/(cav_fac*dNu))**2.
                                 quwinds%wo(i,j,k)=0.
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in courtyard',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           endif
! v values           
                           if(y_v .gt. ymin .and. y_v .lt. ymax &
                                 .and. bld%icellflag(i,j,k) .ne. 0 )then !.and. bld%icellflag(i,j-1,k) .ne. 0
                              if(circle_flag .eq. 1)then
                                 xmax=sqrt((bld%Lt(ibuild)**2.)-(y_v**2.))
                                 xnorm=-xmax
                              else
                                 xnorm=-xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                                       bld%gamma(ibuild)-upwind_dir,-y_v,thetamin,thetamax,qugrid%dxy)
                                 xmax=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                                       bld%gamma(ibuild)-upwind_dir,y_v,thetamin,thetamax,qugrid%dxy)
                              endif
                              if(abs(ynorm) > epsilon .and. eff_height > epsilon)then
                                 dNv=sqrt(exp(-10.*(y_v/ynorm)**2.)*   &
                                          (1.-((zb)/eff_height)**2.)*courtyard_Lr**2)
                              else
                              	dNv=0.
                              endif
! far wake                    
                              if(x_v .gt. xnorm+dNv .and. x_v-xnorm .le. &
                                    farwake_factor*dNv .and. x_v .lt. xmax-qugrid%dxy)then
                                 farwake_velocity=quwinds%vo(i,j,k)*&
                                    (1.-(dNv/(x_v-xnorm+wake_fac*dNv))**(farwake_exponent))
                                 if(abs(farwake_velocity).lt.abs(quwinds%vo(i,j,k)))then
                                    quwinds%vo(i,j,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                    if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                       if(bld%icellflag(i,j,k).ne.0)then
                                          bld%icellflag(i,j,k)=5
                                       endif
                                    endif
                                 endif
                              endif
! cavity                      
                              if(x_v .gt. xnorm +qugrid%dxy .and. x_v .le. dNv+xnorm .and. x_v .lt. xmax-qugrid%dxy)then
                                 quwinds%vo(i,j,k)=-vo_h*min(sqrt(1.-abs(y_v/ynorm)),1.)*(1.-(x_v-xnorm)/(cav_fac*dNv))**2.
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in courtyard',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 quwinds%wo(i,j,k)=0.
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           endif
                        enddo        
                     enddo         
                  enddo
               else
                  roof_ratio=0.8
                  bld%zfo_roof=(bld%Ht(ibuild)-bld%zfo_actual(ibuild))*roof_ratio
                  court_frac=bld%wall(ibuild)
                  do j=bld%jstart(ibuild),bld%jend(ibuild)
                     do i=bld%istart(ibuild),bld%iend(ibuild)
                        x_u=((real(i)-1)*qugrid%dx-xco)*cos(upwind_dir)+		&
                                     ((real(j)-0.5)*qugrid%dy-yco)*sin(upwind_dir)
                        y_u=-((real(i)-1)*qugrid%dx-xco)*sin(upwind_dir)+		&
                                     ((real(j)-0.5)*qugrid%dy-yco)*cos(upwind_dir)
                        x_v=((real(i)-0.5)*qugrid%dx-xco)*cos(upwind_dir)+	&
                                     ((real(j)-1)*qugrid%dy-yco)*sin(upwind_dir)
                        y_v=-((real(i)-0.5)*qugrid%dx-xco)*sin(upwind_dir)+	&
                                     ((real(j)-1)*qugrid%dy-yco)*cos(upwind_dir)
                        r_u=sqrt((x_u**2.)+(y_u**2.))
                        r_v=sqrt((x_v**2.)+(y_v**2.))
                        theta_u=atan2(y_u,x_u)
                        theta_v=atan2(y_v,x_v)
                        ro_u=radius(bld%aa(ibuild),bld%bb(ibuild),theta_u,bld%gamma(ibuild))
                        ro_v=radius(bld%aa(ibuild),bld%bb(ibuild),theta_v,bld%gamma(ibuild))
                        ri_u=radius(bld%aa(ibuild)-court_frac,bld%bb(ibuild)-court_frac,theta_u,bld%gamma(ibuild))
                        ri_v=radius(bld%aa(ibuild)-court_frac,bld%bb(ibuild)-court_frac,theta_v,bld%gamma(ibuild))
                        do k=bld%kstart(ibuild),bld%kend(ibuild)
                           zb=zm(k)-bld%zfo_actual(ibuild)
                           if(zb .le. bld%zfo_roof)then
                              if(r_u .lt. ro_u-0.75*qugrid%dxy &
                                 .and. bld%icellflag(i,j,k) .ne. 0  &
                                 .and. x_u .lt. 0. .and. r_u .gt. recirc_frac*ri_u)then !.and. bld%icellflag(i-1,j,k) .ne. 0 
                                 quwinds%uo(i,j,k)=-indoor_vel_frac*uo_h*(1-zb/bld%zfo_roof)
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in courtyard',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              if(r_v .lt. ro_v-0.75*qugrid%dxy &
                                 .and. bld%icellflag(i,j,k) .ne. 0  &
                                 .and. x_v .lt. 0. .and. r_v .gt. recirc_frac*ri_v)then ! .and. bld%icellflag(i,j-1,k) .ne. 0
                                 quwinds%vo(i,j,k)=-indoor_vel_frac*vo_h*(1-zb/bld%zfo_roof)
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in courtyard',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           else
                              rr_u=ro_u-(ro_u-ri_u)*((zb-bld%zfo_roof)/(0.25*bld%zfo_roof))**2.
                              rr_v=ro_v-(ro_v-ri_v)*((zb-bld%zfo_roof)/(0.25*bld%zfo_roof))**2.
                              if(r_u .lt. rr_u-0.75*qugrid%dxy &
                                 .and. bld%icellflag(i,j,k) .ne. 0 &
                                 .and. x_u .lt. 0..and. r_u .gt. recirc_frac*ri_u)then ! .and. bld%icellflag(i-1,j,k) .ne. 0 
                                 quwinds%uo(i,j,k)=indoor_vel_frac*uo_h*(4*abs(zb-bld%zfo_roof)/bld%zfo_roof)
                                 if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized U exceeds max in courtyard',&
                                       quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                              if(r_v .lt. rr_v-0.75*qugrid%dxy &
                                 .and. bld%icellflag(i,j,k) .ne. 0  &
                                 .and. x_v .lt. 0. .and. r_v .gt. recirc_frac*ri_v)then ! .and. bld%icellflag(i,j-1,k) .ne. 0
                                 quwinds%vo(i,j,k)=indoor_vel_frac*vo_h*(4*abs(zb-bld%zfo_roof)/bld%zfo_roof)
                                 if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized V exceeds max in courtyard',&
                                       quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                 endif
                                 if(i .lt. qugrid%nx .and. j .lt. qugrid%ny .and. k .lt. qugrid%nz)then
                                    if(bld%icellflag(i,j,k) .ne. 0)then
                                       bld%icellflag(i,j,k)=4
                                    endif
                                 endif
                              endif
                           endif
                        enddo
                     enddo
                  enddo
               endif
            case(6)
               farwake_exponent=1.5
               farwake_factor=3
               if(bld%flag%wake .eq. 2)then
                  cav_fac=1.1
                  wake_fac=0.1
               else
                  cav_fac=1.
                  wake_fac=0.
               endif
               ynormp=0.
               ynormm=0.
               tol=0.01*pi/180.
               do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
                  y1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
                  if(y1 .lt. ynormm)ynormm=y1
                  if(y1 .gt. ynormp)ynormp=y1
                  if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                        .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))exit
               enddo
               do k=2,bld%kstart(ibuild)
                  kbottom=k
                  if(bld%zfo_actual(ibuild) .le. zm(k))exit
               enddo
               do k=bld%kstart(ibuild),qugrid%nz-1
                  ktop=k
                  if(bld%Ht(ibuild) .lt. zm(k+1))exit
               enddo
               eff_height=bld%Ht(ibuild)-bld%zfo_actual(ibuild)
               ivert=bld%startidx(ibuild)
               ipoly=1
               do while(ivert .lt. bld%stopidx(ibuild))
                  if(ipoly .gt. 1)then
                     xw1=(bld%x(ivert)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*sin(upwind_dir)
                     yw1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
                     xw3=(bld%x(ivert+1)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*sin(upwind_dir)
                     yw3=-(bld%x(ivert+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*cos(upwind_dir)
                     upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
                     if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
                     if(abs(upwind_rel) .lt. 0.5*pi)then
                        if(abs(upwind_rel) .lt. tol)then
                           perpendicular_flag=1
                        else
                           perpendicular_flag=0
                        endif
                        yw2=min(yw1,yw3)
                        do k=ktop,kbottom,-1
                           zb=zm(k)-bld%zfo(ibuild)
                           do y_idx=0,2*int(abs(yw1-yw3)/qugrid%dxy)+1
                              yc=0.5*real(y_idx)*qugrid%dxy+yw2
                              if(perpendicular_flag .gt. 0)then
                                 xwall=xw1
                              else
                                 xwall=((xw3-xw1)/(yw3-yw1))*(yc-yw1)+xw1
                              endif
                              if(yc .ge. 0.)then
                                 ynorm=ynormp
                              else
                                 ynorm=ynormm
                              endif
                              x_idx_min=-1
                              do x_idx=1,2*int(farwake_factor*bld%Lr(ibuild)/qugrid%dxy)+1
                                 uwakeflag=1
                                 vwakeflag=1
                                 xc=0.5*real(x_idx)*qugrid%dxy !
                                 i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                 j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                                 if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                                    exit
                                 endif
                                 if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                                    x_idx_min=x_idx
                                 endif
                                 if(bld%icellflag(i,j,k) .eq. 0)then
                                    exit
                                 endif
! u values
! Far wake
                                 if(bld%icellflag(i,j,k) .ne. 0 )then
                                    iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                    ju=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                                    xp=real(iu-1)*qugrid%dx-xco
                                    yp=(real(ju)-0.5)*qugrid%dy-yco
                                    xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                    yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                    if(perpendicular_flag .gt. 0)then
                                       xwallu=xw1
                                    else
                                       xwallu=((xw3-xw1)/(yw3-yw1))*(yu-yw1)+xw1
                                    endif
                                    xu=xu-xwallu
                                    if(ynorm > epsilon .and. eff_height > epsilon &
                                          .and. abs(yu) < abs(ynorm) .and. zb < eff_height)then
                                       dNu=sqrt((1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(bld%Lr(ibuild))**2)
                                    else
                                       dNu = 0.0
                                    endif
                                    if(xu .gt. farwake_factor*dNu)uwakeflag=0
                                    if(dNu > 0.0 .and. uwakeflag .eq. 1)then
                                       if(xu .gt. dNu)then
                                          farwake_velocity=quwinds%uo(iu,ju,k)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                                          quwinds%uo(iu,ju,k)=farwake_velocity
                                          quwinds%wo(i,j,k)=0.
                                          if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=5
! Cavity                   
                                       else
                                          quwinds%uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                                          if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in courtyard',&
                                                quwinds%uo(iu,ju,k),quwinds%max_velmag,iu,ju,k
                                          endif
                                          quwinds%wo(i,j,k)=0.
                                          if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=4
                                       endif
                                    endif
! v values
! Far wake
                                    iv=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                    jv=nint(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                                    xp=(real(iv)-0.5)*qugrid%dx-xco
                                    yp=real(jv-1)*qugrid%dy-yco
                                    xv=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                    yv=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                    if(perpendicular_flag .gt. 0)then
                                       xwallv=xw1
                                    else
                                       xwallv=((xw3-xw1)/(yw3-yw1))*(yv-yw1)+xw1
                                    endif
                                    xv=xv-xwallv
                                    if(ynorm > epsilon .and. eff_height > epsilon &
                                          .and. abs(yv) < abs(ynorm) .and. zb < eff_height)then
                                       dNv=sqrt((1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(bld%Lr(ibuild))**2)
                                    else
                                       dNv = 0.0
                                    endif
                                    if(xv .gt. farwake_factor*dNv)vwakeflag=0
                                    if(dNv > 0.0 .and. vwakeflag .eq. 1)then
                                       if(xv .gt. dNv)then
                                          farwake_velocity=quwinds%vo(iv,jv,k)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                                          quwinds%vo(iv,jv,k)=farwake_velocity
                                          quwinds%wo(i,j,k)=0.
                                          if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=5
! Cavity                   
                                       else
                                          quwinds%vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                                          if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in courtyard',&
                                                quwinds%vo(iv,jv,k),quwinds%max_velmag,iv,jv,k
                                          endif
                                          quwinds%wo(iv,jv,k)=0.
                                          if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=4
                                       endif  
                                    endif
                                 endif
                                 if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0)exit
                              enddo     
                           enddo  
                        enddo
                     endif
                  endif
                  ivert=ivert+1
                  if(bld%x(ivert) .eq. bld%x(bld%startidx(ibuild)) &
                        .and. bld%y(ivert) .eq. bld%y(bld%startidx(ibuild)))then
                     ivert=ivert+1
                     ipoly=ipoly+1
                  endif
               enddo
         endselect
         return
      end