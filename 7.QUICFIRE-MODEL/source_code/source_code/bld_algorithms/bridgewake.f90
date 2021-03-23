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
      subroutine bridgewake
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
         use constants
         use bld_module
         use flags_module
         use grid_module
         use winds_module

         implicit none
         
         integer perpendicular_flag,uwakeflag,vwakeflag,wwakeflag
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real xw1,yw1,xw2,yw2,xw3,yw3,xf2,yf2,tol,zb,ynorm
         real farwake_exponent,farwake_factor,farwake_velocity
         real cav_fac,wake_fac,beta,LoverH,upwind_rel_norm
         real bridge_thickness,xc,yc,dNu,dNv,dNw,xwall
         real xu,yu,xv,yv,xp,yp,xwallu,xwallv,xwallw,xw,yw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,iw,jw
         integer ktop,kbottom,i,j,k,ibuild
         
         xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
         yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         ! find upwind direction and determine the type of flow regime
         uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         upwind_rel=upwind_dir-bld%gamma(ibuild)
         if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
         if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
         upwind_rel_norm=upwind_rel+0.5*pi
         if(upwind_rel_norm .gt. pi)upwind_rel_norm=upwind_rel_norm-2*pi

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
            xw2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yf2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            perpendicular_flag=0
         elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
            xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
            xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yf2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            perpendicular_flag=0
         elseif(abs(upwind_rel) .le. tol)then
            xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
            xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xf2=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yf2=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            perpendicular_flag=0
         elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
            xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            perpendicular_flag=1
         elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
            xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yf2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            perpendicular_flag=0
         else
            xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
         endif
         if(bld%flag%wake .gt. 1)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         select case(bld%flag%wake)
            case(1)
               bld%Weff(ibuild)=abs(yw3-yw1)
               if(perpendicular_flag .eq. 1)then
                  bld%Leff(ibuild)=abs(xf2-xw1)
               else
                  bld%Leff(ibuild)=abs(xf2-xw2)
               endif
            case(2)
               beta=abs(atan2(bld%Lti(ibuild),bld%Wti(ibuild)))
               if(abs(upwind_rel) .gt. 0.5*pi-beta .and. &
                     abs(upwind_rel) .lt. 0.5*pi+beta)then
                  bld%Leff(ibuild)=abs(bld%Wti(ibuild)/sin(upwind_rel))
               else
                  bld%Leff(ibuild)=abs(bld%Lti(ibuild)/cos(upwind_rel))
               endif
               if(abs(upwind_rel_norm) .gt. 0.5*pi-beta .and. &
                     abs(upwind_rel_norm) .lt. 0.5*pi+beta)then
                  bld%Weff(ibuild)=abs(bld%Wti(ibuild)/sin(upwind_rel_norm))
               else
                  bld%Weff(ibuild)=abs(bld%Lti(ibuild)/cos(upwind_rel_norm))
               endif
            case(3)
               bld%Leff(ibuild)=bld%Wti(ibuild)*bld%Lti(ibuild)/abs(yw3-yw1)
               if(perpendicular_flag .eq. 1)then
                  bld%Weff(ibuild)=bld%Wti(ibuild)*bld%Lti(ibuild)/abs(xf2-xw1)
               else
                  bld%Weff(ibuild)=bld%Wti(ibuild)*bld%Lti(ibuild)/abs(xf2-xw2)
               endif
            case(4)
               bld%Weff(ibuild)=abs(yw3-yw1)
               bld%Leff(ibuild)=bld%Wti(ibuild)*bld%Lti(ibuild)/bld%Weff(ibuild) 
         endselect
         bridge_thickness=0.5*(bld%Ht(ibuild)-bld%zfo(ibuild))
         LoverH=bld%Leff(ibuild)/bridge_thickness
         if(LoverH.gt.3.)LoverH=3.
         if(LoverH.lt.0.3)LoverH=0.3
         bld%Lr(ibuild)=0.9*bld%Weff(ibuild)/((LoverH**(0.3))*   &
                        (1+0.24*bld%Weff(ibuild)/bridge_thickness))
         do k=2,bld%kstart(ibuild)
            kbottom=k
            if(bld%zfo(ibuild) .le. zm(k))exit
         enddo
         do k=bld%kstart(ibuild),nz-1
            ktop=k
            if(bld%Ht(ibuild) .lt. zm(k+1))exit
         enddo
lp003:   do k=ktop,kbottom,-1
            zb=zm(k)-(bridge_thickness+bld%zfo(ibuild))
            bld%ufarwake(:,:)=0
            bld%vfarwake(:,:)=0
            !$omp parallel do private(yc,xwall,ynorm,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,xwallv,dNv,iw,jw,xw,yw,xwallw,dNw,wwakeflag)
lp002:      do y_idx=1,2*int((yw1-yw3)/qugrid%dxy)-1
               yc=0.5*real(y_idx)*qugrid%dxy+yw3
               if(perpendicular_flag .gt. 0)then
                  xwall=xw1
               elseif(yc.ge.yw2)then
                  xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
               else
                  xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
               endif
               if(yc .ge. 0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif    
               x_idx_min=-1
lp001:         do x_idx=0,2*ceiling(farwake_factor*bld%Lr(ibuild)/qugrid%dxy)
                  uwakeflag=1
                  vwakeflag=1
                  wwakeflag=1
                  xc=0.5*real(x_idx)*qugrid%dxy
                  i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                  j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                  if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(bld%icellflag(i,j,k) .eq. 0)then
                     if(x_idx_min .ge. 0)then
                        if(bld%flag%isbld(i,j,k) .eq. bld%num(ibuild))then
                           x_idx_min=-1
                        else
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
                     if(perpendicular_flag .gt. 0)then
                        xwallu=xw1
                     elseif(yu.ge.yw2)then
                        xwallu=((xw2-xw1)/(yw2-yw1))*(yu-yw1)+xw1
                     else
                        xwallu=((xw3-xw2)/(yw3-yw2))*(yu-yw2)+xw2
                     endif
                     xu=xu-xwallu
                     dNu=sqrt((1.-(yu/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(bld%Lr(ibuild))**2)
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .eq. dNu .and. uwakeflag .eq. 1 .and. bld%icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then
                           farwake_velocity=bld%ufarwake(iu,ju)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent))
                           ! if(bld%icellflag(iu,ju,k) .ne. 4)then
                              quwinds%uo(iu,ju,k)=farwake_velocity
                              quwinds%wo(i,j,k)=0.
                           ! endif
! Cavity                   
                        else
                           quwinds%uo(iu,ju,k)=-uo_h*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Parameterized U exceeds max in rectangle wake',&
                                 quwinds%uo(iu,ju,k),quwinds%max_velmag,iu,ju,k
                           endif
                           quwinds%wo(i,j,k)=0.
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
                     elseif(yv.ge.yw2)then
                        xwallv=((xw2-xw1)/(yw2-yw1))*(yv-yw1)+xw1
                     else
                        xwallv=((xw3-xw2)/(yw3-yw2))*(yv-yw2)+xw2
                     endif
                     xv=xv-xwallv
                     dNv=sqrt((1.-(yv/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(bld%Lr(ibuild))**2)
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .eq. dNv .and. vwakeflag .eq. 1 .and. bld%icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           farwake_velocity=bld%vfarwake(iv,jv)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent))
                           ! if(bld%icellflag(iv,jv,k) .ne. 4)then
                              quwinds%vo(iv,jv,k)=farwake_velocity
                              quwinds%wo(i,j,k)=0.
                           ! endif
! Cavity                   
                        else
                           quwinds%vo(iv,jv,k)=-vo_h*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Parameterized V exceeds max in rectangle wake',&
                                 quwinds%vo(iv,jv,k),quwinds%max_velmag,iv,jv,k
                           endif
                           quwinds%wo(iv,jv,k)=0.
                        endif  
                     endif
! check cell centers to mark cell flags
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                     xp=(real(iw)-0.5)*qugrid%dx-xco
                     yp=(real(jw)-0.5)*qugrid%dy-yco
                     xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                     yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                     if(perpendicular_flag .gt. 0)then
                        xwallw=xw1
                     elseif(yw.ge.yw2)then
                        xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                     else
                        xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                     endif
                     xw=xw-xwallw
                     dNw=sqrt((1.-(yw/ynorm)**2.)*(1.-((zb)/bridge_thickness)**2.)*(bld%Lr(ibuild))**2)
                     if(xw .gt. farwake_factor*dNw)wwakeflag=0
                     if(dNw .eq. dNw .and. wwakeflag .eq. 1 .and. bld%icellflag(iw,jw,k) .ne. 0)then
                        if(xw .gt. dNw)then
                           if(bld%icellflag(iw,jw,k) .eq. 4)then
                              bld%icellflag(iw,jw,k)=12
                           else  
                              bld%icellflag(iw,jw,k)=5
                           endif
! Cavity                   
                        else
                            bld%icellflag(iw,jw,k)=4
                        endif  
                     endif
                  endif
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit
               enddo   lp001      
            enddo   lp002 
            !$omp end parallel do
         enddo   lp003
         return
      end
