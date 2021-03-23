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
      subroutine rectanglewake(ibuild)
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
         integer perpendicular_flag,uwakeflag,vwakeflag,wwakeflag,i,j,k
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real xw1,yw1,xw2,yw2,xw3,yw3,xf2,yf2,tol,zb,ynorm
         real farwake_exponent,farwake_factor,farwake_velocity
         real cav_fac,wake_fac,beta,LoverH,WoverH,upwind_rel_norm,eff_height
         real canyon_factor,xc,yc,dNu,dNv,dNc,xwall,xu,yu,xv,yv,xp,yp,xwallu,xwallv,xwallw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,iw,jw
         real vd,hd,Bs,BL,shell_height,xw,yw,dNw
         integer roof_perpendicular_flag,ns_flag
         integer ktop,kbottom,nupwind
         real LrRect(3),LrLocal,LrLocalu,LrLocalv,LrLocalw
         real epsilon,bridge_thickness,cavity_velocity,wakesign
         real cavityDir,cavityDir1a,cavityDir1b,cavityDir2a,cavityDir2b,referenceSpeed
         real cavityDir1,canyonEffectiveHeight,nonCanyonEffectiveHeight
         real wakeDir1a,wakeDir2a,wakeDir1,wakeDir2,wakeDir,deltaDir1,deltaDir2,wakeAngleWeight
         real cavityCrossFactor,cavityAlongFactor,cosCavityDir,sinCavityDir,cosWindDir,sinWindDir
         real rotationFlag
         
         epsilon = 10e-10
         
         if(bld%geometry(ibuild) .eq. 4 .and. bld%roof(ibuild) .gt. 0)then
            eff_height=0.8*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild)
         else
            eff_height=bld%Ht(ibuild)
         endif
         
         xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
         yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         ! find upwind direction and determine the type of flow regime
         uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1) ! +1
         vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1) ! +1
         referenceSpeed=sqrt(uo_h*uo_h+vo_h*vo_h)
         do k=1,qugrid%nz
            bld%uProfile(k)=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
            bld%vProfile(k)=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),k)
            bld%speedProfile(k)=sqrt(bld%uProfile(k)*bld%uProfile(k)+bld%vProfile(k)*bld%vProfile(k))
         enddo
         upwind_dir=atan2(vo_h,uo_h)
         cosWindDir=cos(upwind_dir)
         sinWindDir=sin(upwind_dir)
         
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
            cavityDir1a=bld%gamma(ibuild)-0.5*pi
            cavityDir2a=bld%gamma(ibuild)
            wakeDir1a=bld%gamma(ibuild)+pi
            wakeDir2a=bld%gamma(ibuild)+0.5*pi
         elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
            xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
            cavityDir1=bld%gamma(ibuild)-0.5*pi
            wakeDir1=upwind_dir
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
            cavityDir1a=bld%gamma(ibuild)-pi
            cavityDir2a=bld%gamma(ibuild)-0.5*pi
            wakeDir1a=bld%gamma(ibuild)+0.5*pi
            wakeDir2a=bld%gamma(ibuild)
         elseif(abs(upwind_rel) .le. tol)then
            xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
            xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xf2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            perpendicular_flag=1
            cavityDir1=bld%gamma(ibuild)-pi
            wakeDir1=upwind_dir
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
            cavityDir1a=bld%gamma(ibuild)+0.5*pi
            cavityDir2a=bld%gamma(ibuild)-pi
            wakeDir1a=bld%gamma(ibuild)
            wakeDir2a=bld%gamma(ibuild)-0.5*pi
         elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
            xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
            xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xf2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
            perpendicular_flag=1
            cavityDir1=bld%gamma(ibuild)+0.5*pi
            wakeDir1=upwind_dir
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
            cavityDir1a=bld%gamma(ibuild)
            cavityDir2a=bld%gamma(ibuild)+0.5*pi
            wakeDir1a=bld%gamma(ibuild)-0.5*pi
            wakeDir2a=bld%gamma(ibuild)+pi
         else
            xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
            yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
            xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
            yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
            xf2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
            perpendicular_flag=1
            cavityDir1=bld%gamma(ibuild)
            wakeDir1=upwind_dir
         endif
         if(perpendicular_flag .eq. 0)then
            if(cavityDir1a .gt. pi)cavityDir1a=cavityDir1a-2.*pi
            if(cavityDir2a .gt. pi)cavityDir2a=cavityDir2a-2.*pi
            if(cavityDir1a .le. -pi)cavityDir1a=cavityDir1a+2.*pi
            if(cavityDir2a .le. -pi)cavityDir2a=cavityDir2a+2.*pi
            if(wakeDir1a .gt. pi)wakeDir1a=wakeDir1a-2.*pi
            if(wakeDir2a .gt. pi)wakeDir2a=wakeDir2a-2.*pi
            if(wakeDir1a .le. -pi)wakeDir1a=wakeDir1a+2.*pi
            if(wakeDir2a .le. -pi)wakeDir2a=wakeDir2a+2.*pi
            deltaDir1=abs(wakeDir1a-upwind_dir)
            deltaDir2=abs(wakeDir2a-upwind_dir)
            if(deltaDir1 .gt. 0.5*pi)then
               !print*,'wrap1',deltaDir1*180/pi,wakeDir1a*180/pi
               if(upwind_dir .ge. 0.)then
                  wakeDir1a=wakeDir1a+2.*pi
               else
                  wakeDir1a=wakeDir1a-2.*pi
               endif
               deltaDir1=abs(wakeDir1a-upwind_dir)
            endif
            if(deltaDir2 .gt. 0.5*pi)then
               !print*,'wrap2',deltaDir2*180/pi,wakeDir2a*180/pi
               if(upwind_dir .ge. 0.)then
                  wakeDir2a=wakeDir2a+2.*pi
               else
                  wakeDir2a=wakeDir2a-2.*pi
               endif
               deltaDir2=abs(wakeDir2a-upwind_dir)
            endif
            wakeAngleWeight=1-(4*abs(deltaDir1-0.25*pi)/pi)**2.
            wakeDir1=wakeDir1a*wakeAngleWeight+upwind_dir*(1-wakeAngleWeight)
            wakeAngleWeight=1-(4*abs(deltaDir2-0.25*pi)/pi)**2.
            wakeDir2=wakeDir2a*wakeAngleWeight+upwind_dir*(1-wakeAngleWeight)
            !wakeDir1=wakeDir1a
            !wakeDir2=wakeDir2a
            ! if(wakeDir1 .gt. pi)wakeDir1=wakeDir1-2.*pi
            ! if(wakeDir2 .gt. pi)wakeDir2=wakeDir2-2.*pi
            ! if(wakeDir1 .le. -pi)wakeDir1=wakeDir1+2.*pi
            ! if(wakeDir2 .le. -pi)wakeDir2=wakeDir2+2.*pi
            cavityDir1b=upwind_dir+0.5*pi
            cavityDir2b=upwind_dir-0.5*pi
            if(cavityDir1b .gt. pi)cavityDir1b=cavityDir1b-2.*pi
            if(cavityDir2b .gt. pi)cavityDir2b=cavityDir2b-2.*pi
            if(cavityDir1b .le. -pi)cavityDir1b=cavityDir1b+2.*pi
            if(cavityDir2b .le. -pi)cavityDir2b=cavityDir2b+2.*pi
            ! cavityDir1b=0.5*(cavityDir1b+cavityDir1a)
            ! cavityDir2b=0.5*(cavityDir2b+cavityDir2a)
            !print*,upwind_dir*180/pi,wakeDir1*180/pi,wakeDir2*180/pi
            !print*,wakeDir1a*180/pi,wakeDir2a*180/pi
            !print*,deltaDir1*180/pi,deltaDir2*180/pi       
         else
            if(cavityDir1 .gt. pi)cavityDir1=cavityDir1-2.*pi
            if(cavityDir1 .le. -pi)cavityDir1=cavityDir1+2.*pi
            if(wakeDir1 .gt. pi)wakeDir1=wakeDir1-2.*pi
            if(wakeDir1 .le. -pi)wakeDir1=wakeDir1+2.*pi
            !print*,upwind_dir*180/pi,wakeDir1*180/pi
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
         
         do k=2,bld%kstart(ibuild)
            kbottom=k
            if(bld%zfo(ibuild) .le. zm(k))exit
         enddo
         do k=bld%kstart(ibuild),qugrid%nz-1
            ktop=k
            if(eff_height .lt. zm(k+1))exit
         enddo
         if(bld%btype(ibuild) .eq. 4)then
            bridge_thickness=0.5*(bld%Ht(ibuild)-bld%zfo(ibuild))
            eff_height=bridge_thickness
         else
            eff_height=eff_height-bld%zfo(ibuild)
         endif
         canyonEffectiveHeight=eff_height
         LoverH=bld%Leff(ibuild)/eff_height
         WoverH=bld%Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         bld%Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         if((bld%btype(ibuild) .eq. 1 .or. bld%btype(ibuild) .eq. 10) .and. bld%flag%roof .eq. 2)then
            tol=30*pi/180.
            roof_perpendicular_flag=0
            ns_flag=-1
            if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
               roof_perpendicular_flag=0
            elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
               roof_perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
               roof_perpendicular_flag=0
            elseif(abs(upwind_rel) .le. tol)then
               roof_perpendicular_flag=1
               ns_flag=0
            elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
               roof_perpendicular_flag=0
            elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
               roof_perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
               roof_perpendicular_flag=0
            else
               roof_perpendicular_flag=1
               ns_flag=0
            endif
            bld%rooftop_flag(ibuild)=1
            k=bld%kend(ibuild)
            nupwind=0
            do y_idx=1,ceiling(2.*(yw1-yw3)/qugrid%dxy)
               yc=0.5*real(y_idx)*qugrid%dxy+yw3
               if(yc .gt. yw1)cycle !yc=yw1
               if(perpendicular_flag .gt. 0)then
                  xwall=xf2
               elseif(yc.ge.yf2)then
                  xwall=((xf2-xw1)/(yf2-yw1))*(yc-yw1)+xw1
               else
                  xwall=((xw3-xf2)/(yw3-yf2))*(yc-yf2)+xf2
               endif
               if(yc.ge.0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif
               do x_idx=int(bld%Lr(ibuild)/qugrid%dxy)+1,1,-1
                  xc=-real(x_idx)*qugrid%dxy
                  i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                  j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                  if(i .ge. 1 .and. i .le. qugrid%nx-1 .and. j .ge. 1 .and. j .le. qugrid%ny-1)then
                     ! print*,i,j,k,bld%icellflag(i,j,k)
                     if(bld%icellflag(i,j,k) .eq. 0)then
                        nupwind=nupwind+1
                        exit
                     endif
                  endif
               enddo
            enddo
            if(nupwind .ge. int((yw1-yw3)/qugrid%dxy))bld%rooftop_flag(ibuild)=0
            ! print*,ibuild,bld%xfo(ibuild),bld%rooftop_flag(ibuild),nupwind,int((yw1-yw3)/qugrid%dxy),bld%Lr(ibuild),xwall
            if(roof_perpendicular_flag .eq. 1 .and. bld%rooftop_flag(ibuild) .eq. 1)then
               if(ns_flag .eq. 1)then
                  hd=bld%Wti(ibuild)
               else
                  hd=bld%Lti(ibuild)
               endif
               Bs=min(bld%Weff(ibuild),bld%Ht(ibuild))
               BL=max(bld%Weff(ibuild),bld%Ht(ibuild))
               bld%Rscale(ibuild) = ((Bs**(2./3.))*(BL**(1./3.)))
               bld%Rcx(ibuild)=(0.9*bld%Rscale(ibuild))
               vd= 0.5*0.22*bld%Rscale(ibuild)
               if(hd .lt. bld%Rcx(ibuild))then
                  shell_height=vd*sqrt(1-((0.5*bld%Rcx(ibuild)-hd)/(0.5*bld%Rcx(ibuild)))**2.)
                  if(shell_height .gt. 0)then
                     eff_height=eff_height+0.5*shell_height
                  endif
               endif
               do k=bld%kstart(ibuild),qugrid%nz-1
                  ktop=k
                  if(eff_height+bld%zfo(ibuild) .lt. zm(k+1))exit
               enddo
            endif
         endif
         nonCanyonEffectiveHeight=eff_height
         LoverH=bld%Leff(ibuild)/eff_height
         WoverH=bld%Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         bld%Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         tol=0.01*pi/180.
         if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
            LrRect(1)=bld%Lr(ibuild)*abs(cos(upwind_rel))
            LrRect(3)=bld%Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
         elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
            LrRect(1)=bld%Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
            LrRect(3)=bld%Lr(ibuild)*abs(cos(upwind_rel))
         elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
            LrRect(1)=bld%Lr(ibuild)*abs(cos(upwind_rel))
            LrRect(3)=bld%Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
         elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
            LrRect(1)=bld%Lr(ibuild)*abs(cos(upwind_rel+0.5*pi))
            LrRect(3)=bld%Lr(ibuild)*abs(cos(upwind_rel))
         endif
         if(perpendicular_flag .eq. 0)then
            LrRect(2)=((yw1-yw2)*LrRect(1)+(yw2-yw3)*LrRect(3))/(yw1-yw3)
            bld%Lr(ibuild)=((yw1-yw2)*0.5*(LrRect(1)+LrRect(2))+(yw2-yw3)*0.5*(LrRect(2)+LrRect(3)))/(yw1-yw3)
         endif
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
            !$omp parallel do private(yc,LrLocal,xwall,ynorm,canyon_factor,x_idx_min,x_idx, &
            !$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,LrLocalu,xwallu,dNu,farwake_velocity, &
            !$omp iv,jv,xv,yv,LrLocalv,xwallv,dNv,iw,jw,xw,yw,LrLocalw,xwallw,dNw,wwakeflag, &
            !$omp cavity_velocity,wakesign,cavityDir,eff_height,dNc, &
            !$omp wakeDir,wakeAngleWeight,cavityCrossFactor,cavityAlongFactor, &
            !$omp cosCavityDir,sinCavityDir)
lp002:      do y_idx=1,ceiling(2.*(yw1-yw3)/qugrid%dxy)
               yc=0.5*real(y_idx)*qugrid%dxy+yw3
               if(yc .ge. yw1)cycle !yc=yw1
               eff_height=nonCanyonEffectiveHeight
               if(perpendicular_flag .gt. 0)then
                  xwall=xw1
                  LrLocal=bld%Lr(ibuild)
                  cavityDir=cavityDir1
                  ! uo_h=referenceSpeed*cos(cavityDir)
                  ! vo_h=referenceSpeed*sin(cavityDir)
                  wakeDir=wakeDir1
               elseif(yc .ge. yw2)then
                  xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
                  LrLocal=LrRect(1)+(yc-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                  cavityDir=cavityDir2a
                  !uo_h=referenceSpeed*cos(cavityDir2)
                  !vo_h=referenceSpeed*sin(cavityDir2)
                  wakeDir=wakeDir2
               else
                  xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
                  LrLocal=LrRect(2)+(yc-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                  cavityDir=cavityDir1a
                  !if(abs(cavityDir2a-cavityDir2b) .gt. pi)then
                  !   cavityDir1=cavityDir1a
                  !else
                  !   cavityDir1=((cavityDir1a-cavityDir1b)/(yw3-yw2))*(yc-yw2)+cavityDir1b
                  !endif
                  !uo_h=referenceSpeed*cos(cavityDir1)
                  !vo_h=referenceSpeed*sin(cavityDir1)
                  wakeDir=wakeDir1
               endif
               cosCavityDir=cos(cavityDir)
               sinCavityDir=sin(cavityDir)
               if(yc .ge. 0.)then
                  ynorm=yw1
               else
                  ynorm=yw3
               endif
               !check for building that will disrupt the wake
               canyon_factor=1.
               rotationFlag=1.
               x_idx_min=-1
               if(abs(yc) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                     .and. zb < eff_height .and. eff_height > epsilon) then 
               	dNc=sqrt((1.-(yc/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocal)**2)
               else
               	dNc = 0.
               endif
               do x_idx=1,int(max(bld%Lr(ibuild),farwake_factor*dNc)/qugrid%dxy)+1
                  xc=real(x_idx)*qugrid%dxy
                  i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                  j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                  if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                     exit
                  endif
                  if(bld%icellflag(i,j,kk) .ne. 0 .and. x_idx_min .lt. 0)then
                     x_idx_min=x_idx
                  endif
                  if(bld%icellflag(i,j,kk) .eq. 0 .and. bld%flag%isbld(i,j,kk) .ne. bld%num(ibuild) .and. x_idx_min .gt. 0)then
                     if(xc .le. bld%Lr(ibuild))then
                        eff_height=canyonEffectiveHeight
                        !canyon_factor=xc/bld%Lr(ibuild)  ! turn off shortening of the wake and cavity in streetcanyons
                     endif
                     rotationFlag=0.
                     exit
                  endif
               enddo
               x_idx_min=-1
lp001:         do x_idx=0,2*ceiling(farwake_factor*bld%Lr(ibuild)/qugrid%dxy)
                  uwakeflag=1
                  vwakeflag=1
                  wwakeflag=1
                  xc=0.5*real(x_idx)*qugrid%dxy !
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
                     if(perpendicular_flag .gt. 0)then
                        xwallu=xw1
                        LrLocalu=bld%Lr(ibuild)
                     elseif(yu.ge.yw2)then
                        xwallu=((xw2-xw1)/(yw2-yw1))*(yu-yw1)+xw1
                        LrLocalu=LrRect(1)+(yu-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallu=((xw3-xw2)/(yw3-yw2))*(yu-yw2)+xw2
                        LrLocalu=LrRect(2)+(yu-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xu=xu-xwallu
                     if(abs(yu) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNu=sqrt((1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalu)**2)
                     else
                     	dNu = 0.
                     endif
                     if(xu .gt. farwake_factor*dNu)uwakeflag=0
                     if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. bld%icellflag(iu,ju,k) .ne. 0)then
                        if(xu .gt. dNu)then
                           cavityCrossFactor=(1-xu/(farwake_factor*dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                           farwake_velocity=bld%uProfile(k)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent)) & 
                              - rotationFlag*referenceSpeed*cavityCrossFactor*sinWindDir*(sinCavityDir*cosWindDir-cosCavityDir*&
                              &sinWindDir)
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
                              ! bld%icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           !cavity_velocity=referenceSpeed*cosCavityDir*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           !cavityCrossFactor=min(sqrt(1.-abs(yu/ynorm)),1.)
                           if(rotationFlag .eq. 1.)then
                              cavityCrossFactor=(1-xu/(farwake_factor*dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                           else
                              cavityCrossFactor=(1-xu/(dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                           endif
                           cavityAlongFactor=min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                           cavity_velocity=referenceSpeed*(cavityAlongFactor*cosWindDir*(cosWindDir*cosCavityDir + sinCavityDir*&
                              &sinWindDir) - cavityCrossFactor*sinWindDir*(sinCavityDir*cosWindDir-cosCavityDir*sinWindDir))
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
                              write(325,*)'Parameterized U exceeds max in rectangle wake',&
							      	   quwinds%uo(iu,ju,k),quwinds%max_velmag,iu,ju,k
						         endif
						         quwinds%wo(i,j,k)=0.
						         ! bld%icellflag(i,j,k)=4
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
                        LrLocalv=bld%Lr(ibuild)
                     elseif(yv.ge.yw2)then
                        xwallv=((xw2-xw1)/(yw2-yw1))*(yv-yw1)+xw1
                        LrLocalv=LrRect(1)+(yv-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallv=((xw3-xw2)/(yw3-yw2))*(yv-yw2)+xw2
                        LrLocalv=LrRect(2)+(yv-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xv=xv-xwallv
                     if(abs(yv) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNv=sqrt((1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalv)**2)
                     else 
                     	dNv = 0.
                     endif
                     if(xv .gt. farwake_factor*dNv)vwakeflag=0
                     if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. bld%icellflag(iv,jv,k) .ne. 0)then
                        if(xv .gt. dNv)then
                           cavityCrossFactor=(1-xv/(farwake_factor*dNv))*min(sqrt(1.-abs(yv/ynorm)),1.)
                           farwake_velocity=bld%vProfile(k)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent)) &
                              + rotationFlag*referenceSpeed*cavityCrossFactor*cosWindDir*(sinCavityDir*cosWindDir-cosCavityDir*&
                               &sinWindDir)
                           if(canyon_factor .eq. 1.)then ! .and. bld%icellflag(iv,jv,k) .ne. 4
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
                              ! bld%icellflag(i,j,k)=5
                           endif
! Cavity                   
                        else
                           !cavity_velocity=referenceSpeed*sinCavityDir*min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.) ! -
                           !cavityCrossFactor=min(sqrt(1.-abs(yv/ynorm)),1.)
                           if(rotationFlag .eq. 1.)then
                              cavityCrossFactor=(1-xv/(farwake_factor*dNv))*min(sqrt(1.-abs(yv/ynorm)),1.)
                           else
                              cavityCrossFactor=(1-xv/(dNv))*min(sqrt(1.-abs(yv/ynorm)),1.)
                           endif
                           cavityAlongFactor=min((1.-xv/(cav_fac*dNv))**2.,1.)*min(sqrt(1.-abs(yv/ynorm)),1.)
                           cavity_velocity=referenceSpeed*(cavityAlongFactor*sinWindDir*(cosWindDir*cosCavityDir + sinCavityDir*&
                               sinWindDir) + cavityCrossFactor*cosWindDir*(sinCavityDir*cosWindDir-cosCavityDir*sinWindDir))
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
                              write(325,*)'Parameterized V exceeds max in rectangle wake',&
                                 quwinds%vo(iv,jv,k),quwinds%max_velmag,iv,jv,k
                           endif
                           quwinds%wo(iv,jv,k)=0.
                           ! bld%icellflag(i,j,k)=4
                        endif  
                     endif
! Far wake
                     iw=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                     jw=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                     xp=(real(iw)-0.5)*qugrid%dx-xco
                     yp=(real(jw)-0.5)*qugrid%dy-yco
                     xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                     yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                     if(perpendicular_flag .gt. 0)then
                        xwallw=xw1
                        LrLocalw=bld%Lr(ibuild)
                     elseif(yw .ge. yw2)then
                        xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                        LrLocalw=LrRect(1)+(yv-yw1)*(LrRect(2)-LrRect(1))/(yw2-yw1)
                     else
                        xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                        LrLocalw=LrRect(2)+(yv-yw2)*(LrRect(3)-LrRect(2))/(yw3-yw2)
                     endif
                     xw=xw-xwallw
                     if(abs(yw) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon) then 
                     	dNw=sqrt((1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalw)**2)
                     else 
                     	dNw = 0.
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
                        elseif(bld%icellflag(iw,jw,bld%kend(ibuild)) .ne. 0)then
                           bld%icellflag(iw,jw,k)=30 ! 4 is the old value
                        endif  
                     endif
                  elseif(x_idx_min .ge. 0 .and. Ht(bld%num(bld%flag%isbld(i,j,k))) .ge. 0.9*eff_height)then
                     exit
                  endif
                  if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)then
                     exit !  .and. wwakeflag .eq. 0
                  endif
               enddo   lp001      
            enddo   lp002
            !$omp end parallel do      
         enddo   lp003
         return
      end
