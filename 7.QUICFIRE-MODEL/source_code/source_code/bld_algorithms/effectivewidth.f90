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
      subroutine effectivewidth
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
			use constants
			use bld_module
			
         implicit none
			
         integer perpendicular_flag,circle_flag
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4,xp,yp
         real xw1,yw1,xw2,yw2,xw3,yw3,xf2,yf2,tol
         real beta,upwind_rel_norm,eff_height
         real thetamax,thetamin,thetai,radius
         real polyarea
			integer :: i,j,ibuild
         
         do ibuild=1,bld%number
            if(bld%damage(ibuild) .eq. 2)cycle
            if(bld%btype(ibuild) .eq. 2 .or. bld%btype(ibuild) .eq. 5)cycle
            if((bld%geometry(ibuild) .eq. 4 .or. bld%geometry(ibuild) .eq. 5) .and. bld%roof(ibuild) .gt. 0)then
               eff_height=0.8*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild)
            else
               eff_height=bld%Ht(ibuild)
            endif
            select case(bld%geometry(ibuild))
               case(1,4)
                  xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
                  yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
                  ! find upwind direction and determine the type of flow regime
                  uo_h=quwinds%uo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
                  vo_h=quwinds%vo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
                  upwind_dir=atan2(vo_h,uo_h)
                  upwind_rel=upwind_dir-bld%gamma(ibuild)
                  if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
                  if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
                  upwind_rel_norm=upwind_rel+0.5*pi
                  if(upwind_rel_norm .gt. pi)upwind_rel_norm=upwind_rel_norm-2*pi
                  tol=0.01*pi/180.
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
               case(2,5)
                  xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
                  yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
                  uo_h=quwinds%uo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
                  vo_h=quwinds%vo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
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
               case(6)
                  xco = bld%cx(ibuild)!CENTER of building in QUIC domain coordinates
                  yco = bld%cy(ibuild)
                  ! find upwind direction and determine the type of flow regime
                  uo_h=quwinds%uo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
                  vo_h=quwinds%vo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
                  upwind_dir=atan2(vo_h,uo_h)
                  tol=0.01*pi/180.
                  x1=0.
                  x2=0.
                  y1=0.
                  y2=0.
                  polyarea=0.
                  do j=bld%startidx(ibuild),bld%stopidx(ibuild)
                     polyarea=polyarea+bld%x(j)*bld%y(j+1)-bld%x(j+1)*bld%y(j)
                     xp=(bld%x(j)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(j)-bld%cy(ibuild))*sin(upwind_dir)
                     yp=-(bld%x(j)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(j)-bld%cy(ibuild))*cos(upwind_dir)
                     if(xp .lt. x1)x1=xp
                     if(xp .gt. x2)x2=xp
                     if(yp .lt. y1)y1=yp
                     if(yp .gt. y2)y2=yp
                     if(bld%x(j+1) .eq. bld%x(bld%startidx(ibuild)) .and. bld%y(j+1) .eq. bld%y(bld%startidx(ibuild)))exit
                  enddo
                  polyarea=0.5*abs(polyarea)
                  if(bld%flag%wake .eq. 4)then
                     bld%Weff(ibuild)=(y2-y1)
                     bld%Leff(ibuild)=polyarea/bld%Weff(ibuild)
                  else
                     bld%Weff(ibuild)=polyarea/(x2-x1)
                     bld%Leff(ibuild)=polyarea/(y2-y1)
                  endif
            endselect
            if(bld%zfo_actual(ibuild) .gt. 0)then
               call building_connect(ibuild)
            endif
         enddo
         return
      end