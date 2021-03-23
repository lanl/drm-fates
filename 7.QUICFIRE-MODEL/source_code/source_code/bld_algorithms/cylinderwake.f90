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
      subroutine cylinderwake(ibuild)
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
         integer circle_flag,uwakeflag,vwakeflag,wwakeflag
         real uo_h,vo_h,upwind_dir,xco,yco,tol,zb
         real farwake_exponent,farwake_factor,farwake_velocity
         real thetamax,thetamin,thetai,LoverH,WoverH
         real ynorm,radius,xnorm_bisect,cav_fac,wake_fac,eff_height
         real canyon_factor,xc,yc,yw1,yw3,y1,y2,dNu,dNv,xwall
         real xu,yu,xv,yv,xp,yp,xwallu,xwallv
         real xw,yw,dNw,xwallw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,iw,jw
         integer ktop,kbottom,i,j,k
         real epsilon,bridge_thickness,cavity_velocity,wakesign
         
         epsilon = 10e-10
         
         if(bld%geometry(ibuild) .eq. 5 .and. bld%roof(ibuild) .gt. 0)then
            eff_height=0.8*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild)
         else
            eff_height=bld%Ht(ibuild)
         endif
         farwake_exponent=1.5
         farwake_factor=3
         xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
         yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1) ! +1
         vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1) ! +1
         do k=1,nz
            bld%uProfile(k)=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
            bld%vProfile(k)=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
            bld%speedProfile(k)=sqrt(bld%uProfile(k)*bld%uProfile(k)+bld%vProfile(k)*bld%vProfile(k))
         enddo
! find upwind direction and determine the type of flow regime
         upwind_dir=atan2(vo_h,uo_h)
         tol=0.01*qugrid%dxy
         if(abs(bld%aa(ibuild)-bld%bb(ibuild)) .lt. tol)then
            circle_flag=1
         else
            circle_flag=0
         endif
         if(circle_flag .eq. 1)then
            thetamin=-0.5*pi
            thetamax=0.5*pi
            yw1=bld%Lt(ibuild)
            yw3=-bld%Lt(ibuild)
            bld%Weff(ibuild)=bld%Lti(ibuild)
            bld%Leff(ibuild)=bld%Lti(ibuild)
         else
            y1=0.
            y2=0.
            do i=1,180
               thetai=real(180-i)*pi/180.
               y1=radius(bld%aa(ibuild),bld%bb(ibuild),thetai,&
                   bld%gamma(ibuild)-upwind_dir)*sin(thetai)
               if(y1 .lt. y2)then
                  exit
               endif
               y2=y1
            enddo
            thetamax=thetai+pi/180.
            thetamin=thetamax-pi
            yw1=y2
            yw3=-y2
            bld%Weff(ibuild)=2*yw1
            bld%Leff(ibuild)=2*radius(bld%aa(ibuild),bld%bb(ibuild),upwind_dir,bld%gamma(ibuild))
         endif
         if(bld%flag%wake .eq. 2)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         
         do k=2,bld%kstart(ibuild)
            kbottom=k
            if(bld%zfo(ibuild) .le. zm(k))exit
         enddo
         do k=bld%kstart(ibuild),nz-1
            ktop=k
            if(eff_height .lt. zm(k+1))exit
         enddo
         if(bld%btype(ibuild) .eq. 4)then
            bridge_thickness=0.5*(bld%Ht(ibuild)-bld%zfo(ibuild))
            eff_height=bridge_thickness   
         else
            eff_height=eff_height-bld%zfo(ibuild)
         endif
         LoverH=bld%Leff(ibuild)/eff_height
         WoverH=bld%Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         bld%Lr(ibuild)=0.9*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         ynorm=yw1
         do k=bld%kstart(ibuild),bld%kend(ibuild)
            kk=k
            if(0.75*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild) .le. zm(k))exit
         enddo
lp003:   do k=ktop,kbottom,-1
            if(bld%btype(ibuild) .eq. 4)then
               zb=zm(k)-(bridge_thickness+bld%zfo(ibuild))
            else
               zb=zm(k)-bld%zfo(ibuild)
            endif
            bld%ufarwake(:,:)=quwinds%uo(:,:,k)
            bld%vfarwake(:,:)=quwinds%vo(:,:,k)
            bld%icellwake(:,:)=bld%icellflag(:,:,k)
            !$omp parallel do private(yc,xwall,canyon_factor,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,xwallv,dNv,iw,jw,xw,yw,xwallw,dNw,wwakeflag, &
            !$omp cavity_velocity,wakesign)
lp002:      do y_idx=1,ceiling(2.*abs(yw1-yw3)/qugrid%dxy)
               yc=0.5*real(y_idx)*qugrid%dxy+yw3
               if(yc .gt. yw1) cycle !yc=yw1
               if(circle_flag .eq. 1)then
                  if(abs(yc) .gt. bld%Lt(ibuild))then
                     cycle
                  else
                     xwall=sqrt((bld%Lt(ibuild)**2.)-(yc**2.))
                  endif
               else
                  xwall=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                        bld%gamma(ibuild)-upwind_dir,yc,thetamin,thetamax,qugrid%dxy)
               endif
               ! check for building that will disrupt the wake
               canyon_factor=1.
               x_idx_min=-1
               do x_idx=1,ceiling(bld%Lr(ibuild)/qugrid%dxy)
                  xc=real(x_idx)*qugrid%dxy
                  i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                  j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                  if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(bld%icellflag(i,j,kk) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(bld%icellflag(i,j,kk) .eq. 0 .and. ibldflag(i,j,kk) .ne. bld%num(ibuild) .and. x_idx_min .gt. 0)then
                     ! canyon_factor=xc/bld%Lr(ibuild) ! turn off shortening of the wake and cavity in streetcanyons
                     exit
                  endif
               enddo
               x_idx_min=-1
lp001:         do x_idx=0,2*ceiling(farwake_factor*bld%Lr(ibuild)/qugrid%dxy)
                  uwakeflag=1
                  vwakeflag=1
                  wwakeflag=1
                  xc=0.5*real(x_idx)*qugrid%dxy
                  i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                  j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                  if(i .ge. nx-1 .or. i .le. 1 .or. j .ge. ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(bld%icellflag(i,j,k) .eq. 0)then
                     if(x_idx_min .ge. 0)then
                        if(ibldflag(i,j,k) .eq. bld%num(ibuild))then
                           x_idx_min=-1
                        elseif(canyon_factor .lt. 1.)then
                           exit
                        endif
                     endif
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
                     if(circle_flag .eq. 1)then
                        if(abs(yu) .gt. bld%Lt(ibuild))then
                           cycle
                        else
                           xwallu=sqrt((bld%Lt(ibuild)**2.)-(yu**2.))
                        endif
                     else
                        xwallu=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                              bld%gamma(ibuild)-upwind_dir,yu,thetamin,thetamax,qugrid%dxy)
                     endif
                     xu=xu-xwallu

                     
                     if(ynorm > epsilon .and. abs(yu) < ynorm &
                           .and. eff_height > epsilon  .and. zb < eff_height)then
                        dNu = (1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*bld%Lr(ibuild))**2
                        dNu=sqrt(dNu)
                     else
                        dNu = 0.0
                     endif
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. bld%icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then                           
                           farwake_velocity=bld%uProfile(k)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then ! .and. bld%icellflag(iu,ju,k) .ne. 4
                              select case(bld%flag%blending)
                                 case(0)
                                    quwinds%uo(iu,ju,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                 case default
                                    if(bld%icellwake(iu,ju) .eq. 1)then
                                       quwinds%uo(iu,ju,k)=farwake_velocity
                                       quwinds%wo(i,j,k)=0.
                                    elseif(bld%icellwake(iu,ju) .eq. 31)then
                                       if(farwake_velocity .ge. 0.)then
                                          wakesign=1.0
                                       else
                                          wakesign=-1.0
                                       endif
                                       if(wakesign*farwake_velocity .lt. wakesign*quwinds%uo(iu,ju,k))then
                                          quwinds%uo(iu,ju,k)=farwake_velocity
                                          quwinds%wo(i,j,k)=0.
                                       endif
                                    endif
                              end select
                              ! if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           cavity_velocity=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           select case(bld%flag%blending)
                              case(0)
                                 quwinds%uo(iu,ju,k)=cavity_velocity
                              case(1)
                                 if(bld%icellwake(iu,ju) .ne. 30)then
                                    quwinds%uo(iu,ju,k)=cavity_velocity
                                 endif
                              case(2)
                                 if(bld%icellwake(iu,ju) .ne. 30)then
                                    quwinds%uo(iu,ju,k)=cavity_velocity
                                 else
                                    quwinds%uo(iu,ju,k)=0.5*(quwinds%uo(iu,ju,k)+cavity_velocity)
                                 endif
                              case(3)
                                 if(cavity_velocity .ge. 0.)then
                                    wakesign=1.0
                                 else
                                    wakesign=-1.0
                                 endif
                                 if(bld%icellwake(iu,ju) .ne. 30)then
                                    quwinds%uo(iu,ju,k)=cavity_velocity
                                 elseif(wakesign*cavity_velocity .gt. wakesign*quwinds%uo(iu,ju,k))then
                                    quwinds%uo(iu,ju,k)=cavity_velocity
                                 endif
                           end select
                           if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Parameterized U exceeds max in cylinder wake',&
                                 quwinds%uo(iu,ju,k),quwinds%max_velmag,iu,ju,k
                           endif
                           quwinds%wo(i,j,k)=0.
                           ! if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=4
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
                     if(circle_flag .eq. 1)then
                        if(abs(yv) .gt. bld%Lt(ibuild))then
                           cycle
                        else
                           xwallv=sqrt((bld%Lt(ibuild)**2.)-(yv**2.))
                        endif
                     else
                        xwallv=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                              bld%gamma(ibuild)-upwind_dir,yv,thetamin,thetamax,qugrid%dxy)
                     endif
                     xv=xv-xwallv
                     
                     if(ynorm > epsilon .and. abs(yv) < ynorm &
                           .and. eff_height > epsilon .and. zb < eff_height)then
                        dNv = (1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*bld%Lr(ibuild))**2
                        dNv=sqrt(dNv)
                     else
                        dNv = 0.0
                     endif
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. bld%icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           farwake_velocity=bld%vProfile(k)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                           if(canyon_factor .eq. 1.)then !  .and. bld%icellflag(iv,jv,k) .ne. 4
                              select case(bld%flag%blending)
                                 case(0)
                                    quwinds%vo(iv,jv,k)=farwake_velocity
                                    quwinds%wo(i,j,k)=0.
                                 case default
                                    if(bld%icellwake(iv,jv) .eq. 1)then
                                       quwinds%vo(iv,jv,k)=farwake_velocity
                                       quwinds%wo(i,j,k)=0.
                                    elseif(bld%icellwake(iv,jv) .eq. 31)then
                                       if(farwake_velocity .ge. 0.)then
                                          wakesign=1.0
                                       else
                                          wakesign=-1.0
                                       endif
                                       if(wakesign*farwake_velocity .lt. wakesign*quwinds%vo(iv,jv,k))then
                                          quwinds%vo(iv,jv,k)=farwake_velocity
                                          quwinds%wo(i,j,k)=0.
                                       endif
                                    endif
                              end select
                              ! if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           cavity_velocity=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           select case(bld%flag%blending)
                              case(0)
                                 quwinds%vo(iv,jv,k)=cavity_velocity
                              case(1)
                                 if(bld%icellwake(iv,jv) .ne. 30)then
                                    quwinds%vo(iv,jv,k)=cavity_velocity
                                 endif
                              case(2)
                                 if(bld%icellwake(iv,jv) .ne. 30)then
                                    quwinds%vo(iv,jv,k)=cavity_velocity
                                 else
                                    quwinds%vo(iv,jv,k)=0.5*(quwinds%vo(iv,jv,k)+cavity_velocity)
                                 endif
                              case(3)
                                 if(cavity_velocity .ge. 0.)then
                                    wakesign=1.0
                                 else
                                    wakesign=-1.0
                                 endif
                                 if(bld%icellwake(iv,jv) .ne. 30)then
                                    quwinds%vo(iv,jv,k)=cavity_velocity
                                 elseif(wakesign*cavity_velocity .gt. wakesign*quwinds%vo(iv,jv,k))then
                                    quwinds%vo(iv,jv,k)=cavity_velocity
                                 endif
                           end select
                           if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Parameterized V exceeds max in cylinder wake',&
                                 quwinds%vo(iv,jv,k),quwinds%max_velmag,iv,jv,k
                           endif
                           quwinds%wo(iv,jv,k)=0.
                           ! if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=4
                        endif  
                     endif
! celltype values
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                     xp=(real(iw)-0.5)*qugrid%dx-xco
                     yp=(real(jw)-0.5)*qugrid%dy-yco
                     xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                     yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                     if(circle_flag .eq. 1)then
                        if(abs(yw) .gt. bld%Lt(ibuild))then
                           cycle
                        else
                           xwallw=sqrt((bld%Lt(ibuild)**2.)-(yw**2.))
                        endif
                     else
                        xwallw=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                              bld%gamma(ibuild)-upwind_dir,yw,thetamin,thetamax,qugrid%dxy)
                     endif
                     xw=xw-xwallw
                     if(ynorm > epsilon .and. abs(yw) < ynorm &
                           .and. eff_height > epsilon .and. zb < eff_height)then
                        dNw = (1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*bld%Lr(ibuild))**2
                        dNw=sqrt(dNw)
                     else
                        dNw = 0.0
                     endif
                     if(xw .gt. farwake_factor*dNw)wwakeflag=0
                     if(dNw .gt. 0. .and. wwakeflag .eq. 1 .and. bld%icellflag(iw,jw,k) .ne. 0)then
                        if(xw .gt. dNw)then
                           if(canyon_factor .eq. 1.)then
                              ! if(bld%icellflag(iw,jw,k) .eq. 4)then
                              !    bld%icellflag(iw,jw,k)=12
                              ! else
                              !    bld%icellflag(iw,jw,k)=5
                              ! endif
                              if(bld%icellflag(iw,jw,k) .eq. 1)then
                                 bld%icellflag(iw,jw,k)=31 ! 5 is the old value
                              endif
                           endif
! Cavity                   
                        else
                           bld%icellflag(iw,jw,k)=30   ! 4 is the old value
                        endif  
                     endif
                  elseif(x_idx_min .ge. 0 .and. Ht(bld%num(ibldflag(i,j,k))) .ge. 0.9*eff_height)then
                     exit
                  endif
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit
               enddo   lp001      
            enddo   lp002
            !$omp end parallel do
         enddo   lp003
         return
      end
      
      
      real function radius(a,b,bld%gamma,theta)
         implicit none
         real a,b,bld%gamma,theta
         radius=a*b/sqrt( (a*sin(theta-bld%gamma))**2. + (b*cos(theta-bld%gamma))**2. )
         return
      end
      
      
      real function xnorm_bisect(a,b,bld%gamma,y,thetamin,thetamax,qugrid%dxy)
         implicit none
         integer i
         real, intent(in) :: a,b,bld%gamma,y,qugrid%dxy,thetamin,thetamax
         real yguess,yguess_low,eps,prod,theta,rad,radius,thetalow,thetahigh
         i=0
         eps=a+b
         thetalow=thetamin
         thetahigh=thetamax
         do while (i .lt. 100 .and. eps .gt. 0.1*qugrid%dxy)
            i=i+1
            theta=0.5*(thetalow+thetahigh)
            rad=radius(a,b,bld%gamma,theta)
            yguess=rad*sin(theta)-y
            yguess_low=radius(a,b,bld%gamma,thetalow)*sin(thetalow)-y
            eps=abs(yguess)
            prod=yguess*yguess_low
            if(prod .lt. 0)then
               thetahigh=theta
            elseif(prod .gt. 0)then
               thetalow=theta
            else
               eps=0.
            endif
         enddo
         xnorm_bisect=rad*cos(theta)
         return
      end