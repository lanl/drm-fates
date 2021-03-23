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
      subroutine defbuild
!************************************************************************
! defbuild- define building geometries    
!    - called by bcsetup.f       
!    - calls pentagon.f                   
! this is a new subroutine as of 7/26/03 that allows buildings that are
! none orthagonal to the main coordinate system to be used.
! bld%gamma - building angle, is valid +/- 45 degrees
!
!ERP 2003         
!************************************************************************
         
         use constants
         use bld_module
         use canopy_module
         use landuse_module
         use grid_module
         
         implicit none
         
         real x1,x2,x3,x4,y1,y2,y3,y4,xL1,xL2,yL3,yL4,x_c,y_c         
         real slope,xmin,ymin,xmax,ymax         
         real yco,xco
         real radius,thetacell
         real rayintersect
         integer ilevel,ivert,startpoly,numcrossing
         integer :: i,j,k,ibuild
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! this is just a building generation part
         if(landuse%canopy_flag .gt. 0 .or. inumcanopy.gt.0)then
            canopy%top = 0.
            canopy%atten  = 0.
         endif
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate buildings
! this now includes a building generation loop that allows
! for multiple buildings
! calculate building spacing s and set wake flags
         !$omp parallel do
         do ibuild=1,bld%number
            select case(bld%geometry(ibuild))
               case(1,4)
                  if(bld%gamma(ibuild) .eq. 0)then
                     bld%bld%jstart(ibuild)=int(bld%xfo(ibuild)/qugrid%dx)     ! front of the building
                     bld%iend(ibuild)=int((bld%xfo(ibuild)+bld%Lti(ibuild))/qugrid%dx)+1 ! back of the bld
                     bld%jend(ibuild)=int((bld%yfo(ibuild)+bld%Wt(ibuild))/qugrid%dy)+1  ! far side of bld
                     jstart(ibuild)=int((bld%yfo(ibuild)-bld%Wt(ibuild))/qugrid%dy)! close side of bld
                     if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                     if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                     if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                     if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
                  else
                     x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))
                     y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))
                     x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
                     y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
                     x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))
                     y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))
                     x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
                     y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
                     if(bld%gamma(ibuild).gt.0)then
                        xmin=x4
                        xmax=x2
                        ymin=y1
                        ymax=y3
                     endif
                     if(bld%gamma(ibuild).lt.0)then
                        xmin=x1
                        xmax=x3
                        ymin=y2
                        ymax=y4
                     endif
                     bld%bld%jstart(ibuild)=int(xmin/qugrid%dx)
                     bld%iend(ibuild)=int(xmax/qugrid%dx)+1
                     jstart(ibuild)=int(ymin/qugrid%dy)
                     bld%jend(ibuild)=int(ymax/qugrid%dy)+1
                     if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                     if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                     if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                     if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
                  endif
               case(2,5)
                  if(bld%gamma(ibuild) .eq. 0)then
                     bld%bld%jstart(ibuild)=int(bld%xfo(ibuild)/qugrid%dx)     ! front of the building
                     bld%iend(ibuild)=int((bld%xfo(ibuild)+bld%Lti(ibuild))/qugrid%dx)+1 ! back of the bld
                     bld%jend(ibuild)=int((bld%yfo(ibuild)+bld%Wt(ibuild))/qugrid%dy)+1  ! far side of bld
                     jstart(ibuild)=int((bld%yfo(ibuild)-bld%Wt(ibuild))/qugrid%dy)! close side of bld
                     if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                     if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                     if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                     if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
                  else
                     xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))
                     yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
                     bld%bld%jstart(ibuild)=int((xco-max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dx)
                     bld%iend(ibuild)=int((xco+max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dx)+1
                     jstart(ibuild)=int((yco-max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dy)
                     bld%jend(ibuild)=int((yco+max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dy)+1
                  endif     
               case(3)
                  bld%bld%jstart(ibuild)=int(bld%xfo(ibuild)/qugrid%dx)     ! front of the building
                  bld%iend(ibuild)=int((bld%xfo(ibuild)+bld%Lti(ibuild))/qugrid%dx)+1 ! back of the bld
                  bld%jend(ibuild)=int((bld%yfo(ibuild)+bld%Wt(ibuild))/qugrid%dy)+1  ! far side of bld
                  jstart(ibuild)=int((bld%yfo(ibuild)-bld%Wt(ibuild))/qugrid%dy)! close side of bld
                  if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                  if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                  if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                  if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
               case(6)
                  x1=bld%x(bld%startidx(ibuild))
                  x2=x1
                  y1=bld%y(bld%startidx(ibuild))
                  y2=y1
                  do ivert=bld%startidx(ibuild)+1,bld%stopidx(ibuild)
                     if(bld%x(ivert) .eq. bld%x(bld%startidx(ibuild)) .and. &
                           bld%y(ivert) .eq. bld%y(bld%startidx(ibuild)))exit
                     if(bld%x(ivert) .lt. x1)x1=bld%x(ivert)
                     if(bld%x(ivert) .gt. x2)x2=bld%x(ivert)
                     if(bld%y(ivert) .lt. y1)y1=bld%y(ivert)
                     if(bld%y(ivert) .gt. y2)y2=bld%y(ivert)
                  enddo
                  bld%bld%jstart(ibuild)=int(x1/qugrid%dx)      !front of the building  
                  bld%iend(ibuild)=int(x2/qugrid%dx)+1  !back of the bld  
                  bld%jend(ibuild)=int(y2/qugrid%dy)  !far side of bld  
                  jstart(ibuild)=int(y1/qugrid%dy)+1!close side of bld
                  if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                  if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                  if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                  if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
            endselect
         enddo
         !$omp end parallel do
         if(landuse%canopy_flag .gt. 0)then
            !$omp parallel do private(i,j,k)
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  canopy%top(i,j)=landuse%height(i,j)
                  if(canopy%top(i,j) .gt. 0.)then
                     do k=2,qugrid%nz-1
                        if(k .eq. 2 .and. canopy%top(i,j) .lt. zm(k))then
                           canopy%top(i,j)=0.
                           exit
                        endif
                        canopy%atten(i,j,k)=landuse%atten(i,j)
                        if(canopy%top(i,j) .le. zm(k))exit
                     enddo
                  endif
               enddo
            enddo
            !$omp end parallel do
         endif
!building Loop2Loop2Loop2!Loop2Loop2Loop2!Loop2Loop2Loop2!Loop2Loop2 begin
         
!ccccc need an even number of cells in building
         !!$omp parallel do private(ibuild,i,j,k,ilevel,x1,x2,x3,x4,y1,y2,y3,y4, &
         !!$omp xmin,xmax,ymin,ymax,x_c,y_c,slope,xL1,xL2,yL3,yL4,xco,yco, &
         !!$omp thetacell,roof_ratio,roof_bld%zfo,kroof,z_c,court_frac, &
         !!$omp bld%xfoin,bld%yfoin,x1in,x2in,x3in,x4in,y1in,y2in,y3in,y4in, &
         !!$omp xL1in,xL2in,yL3in,yL4in,z_c_x,z_c_y,wall_thickness, &
         !!$omp radius_out,radius_in,r_c,ivert,numcrossing,startpoly,rayintersect)
lp002:   do ibuild=1,bld%number         !begin building gen loop 2
! set non-fluid cell type flags
! icellflag = 0 is a solid cell
            if(bld%damage(ibuild) .eq. 2 .or. (bld%btype(ibuild) .eq. 1 .and. bld%geometry(ibuild) .ne. 3))then
               cycle
            endif
            ilevel=0
            if(landuse%canopy_flag .gt. 0 .and. bld%btype(ibuild) .eq. 2 &
                  .and. bld%kend(ibuild) .lt. 2)then
               bld%kend(ibuild)=2
            endif
! if the building is orthoganol to the main coordinate system
! just set the cellflags as follows
            !print*,bld%num(ibuild),bld%Ht(ibuild),atten(ibuild)
            select case(bld%geometry(ibuild))
               case(1) !Rectangular buildings
                  if(bld%gamma(ibuild) .eq. 0)then
!erp  do k=int(bld%zfo(ibuild)),bld%kend(ibuild)
! int changed to nint on next line 8-14-06   
                     do j=jstart(ibuild),bld%jend(ibuild)
                        do i=bld%bld%jstart(ibuild),bld%iend(ibuild)
                           select case(bld%btype(ibuild))
                              case(0)
                                 icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=1
                                 ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                              case(2)
                                 do k=bld%kstart(ibuild),bld%kend(ibuild)
                                    if(icellflag(i,j,k) .ne. 0)then
                                       if(landuse%canopy_flag .gt. 0)then
                                          if(canopy%top(i,j) .eq. landuse%height(i,j))then
                                             if(k .eq. bld%kstart(ibuild))canopy%atten(i,j,bld%kstart(ibuild):qugrid%nz-1)=0
                                             if(bld%Ht(ibuild) .lt. 0.5*qugrid%dz_array(1))then
                                                canopy%top(i,j)=0.
                                             else
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                          elseif(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                             canopy%top(i,j)=bld%Ht(ibuild)
                                          endif
                                          if(zm(k) .lt. canopy%top(i,j))canopy%atten(i,j,k)=atten(ibuild)
                                       else
                                          if(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                             canopy%top(i,j)=bld%Ht(ibuild)
                                          endif
                                          canopy%atten(i,j,k)=atten(ibuild)
                                       endif
                                    endif
                                 enddo
                              case(3)
                                 ilevel=0
                                 do k=bld%kstart(ibuild),bld%kend(ibuild)
                                    ilevel=ilevel+1
                                    if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                    icellflag(i,j,k)=0
                                    ! ibldflag(i,j,k)=ibuild
                                 enddo
                              case default
                                 ! icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=0
                                 ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                           endselect
                        enddo
                     enddo
! if the building is NON-orthoganol to the main coordinate system
! use the following algorithm
                  else
! calculate corner coordinates of the building
                     x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))
                     y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))
                     x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
                     y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
                     x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))
                     y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))
                     x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
                     y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
 
                     if(bld%gamma(ibuild).gt.0)then
                        xmin=x4
                        xmax=x2
                        ymin=y1
                        ymax=y3
                     endif
                     if(bld%gamma(ibuild).lt.0)then
                        xmin=x1
                        xmax=x3
                        ymin=y2
                        ymax=y4
                     endif
                     bld%bld%jstart(ibuild)=int(xmin/qugrid%dx)
                     bld%iend(ibuild)=int(xmax/qugrid%dx)+1
                     jstart(ibuild)=int(ymin/qugrid%dy)
                     bld%jend(ibuild)=int(ymax/qugrid%dy)+1
                     if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                     if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                     if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                     if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
!erp  do k=int(bld%zfo(ibuild)),bld%kend(ibuild)  
!erp        do j=int(ymin),int(ymax)
!erp     do i=int(xmin),int(xmax)
!erp     x_c=real(i) !x coordinate to be checked
!erp     y_c=real(j) !y coordinate to be checked
! changed int to nint in next three lines 8-14-06
                     do j=nint(ymin/qugrid%dy)+1,nint(ymax/qugrid%dy)+1   !convert back to real world unit, TZ 10/29/04
                        y_c=(real(j)-0.5)*qugrid%dy !y coordinate to be checked   !convert back to real world unit, TZ 10/29/04
                        do i=nint(xmin/qugrid%dx)+1,nint(xmax/qugrid%dx)+1   !convert back to real world unit, TZ 10/29/04
                           x_c=(real(i)-0.5)*qugrid%dx !x coordinate to be checked   !convert back to real world unit, TZ 10/29/04
!calculate the equations of the lines making up the 4 walls of the
!building
						   if( x4 .eq. x1)x4=x4+.0001
                           slope = (y4-y1)/(x4-x1) !slope of L1
                           xL1 = x4 + (y_c-y4)/slope
                           if( x3 .eq. x2)x3=x3+.0001
                           slope = (y3-y2)/(x3-x2) !slope of L2
                           xL2 = x3 + (y_c-y3)/slope
                           if( x2 .eq. x1)x2=x2+.0001
                           slope = (y2-y1)/(x2-x1) !slope of L3
                           yL3 = y1 + slope*(x_c-x1)
                           if( x3 .eq. x4)x3=x3+.0001
                           slope = (y3-y4)/(x3-x4) !slope of L4
                           yL4 = y4 + slope*(x_c-x4)
                           if(x_c.gt.xL1.and.x_c.lt.xL2.and.y_c.gt.yL3.and.y_c.lt.yL4)then
                              select case(bld%btype(ibuild))
                                 case(0)
                                    icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=1
                                    ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                                 case(2)
                                    do k=bld%kstart(ibuild),bld%kend(ibuild)
                                       if(icellflag(i,j,k) .ne. 0)then
                                          if(landuse%canopy_flag .gt. 0)then
                                             if(canopy%top(i,j) .eq. landuse%height(i,j))then
                                                if(k .eq. bld%kstart(ibuild))canopy%atten(i,j,bld%kstart(ibuild):qugrid%nz-1)=0
                                                if(bld%Ht(ibuild) .lt. 0.5*qugrid%dz_array(1))then
                                                   canopy%top(i,j)=0.
                                                else
                                                   canopy%top(i,j)=bld%Ht(ibuild)
                                                endif
                                             elseif(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                             if(zm(k) .lt. canopy%top(i,j))canopy%atten(i,j,k)=atten(ibuild)
                                          else
                                             if(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                             canopy%atten(i,j,k)=atten(ibuild)
                                          endif
                                       endif
                                    enddo
                                 case(3)
                                    ilevel=0
                                    do k=bld%kstart(ibuild),bld%kend(ibuild)
                                       ilevel=ilevel+1
                                       if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                       icellflag(i,j,k)=0
                                       ! ibldflag(i,j,k)=ibuild
                                    enddo
                                 case default
                                    ! icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=0
                                    ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                              endselect
                           endif
                        enddo
                     enddo
                  endif
! generate cylindrical buildings
! need to specify a and b as the major and minor axis of
! the ellipse
! xco and yco are the coordinates of the center of the ellipse
               case(2)
                  if(bld%aa(ibuild) .gt. 0. .and. bld%bb(ibuild) .gt. 0.)then
                     if(bld%gamma(ibuild) .ne. 0.)then
                        xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))
                        yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
                        bld%bld%jstart(ibuild)=int((xco-max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dx)
                        bld%iend(ibuild)=int((xco+max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dx)+1
                        jstart(ibuild)=int((yco-max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dy)
                        bld%jend(ibuild)=int((yco+max(bld%Lt(ibuild),bld%Wt(ibuild)))/qugrid%dy)+1
                     else
                        xco = bld%xfo(ibuild) + bld%Lt(ibuild)
                        yco = bld%yfo(ibuild)
                     endif
                     if(bld%bld%jstart(ibuild) .lt. 1)bld%bld%jstart(ibuild)=1
                     if(bld%iend(ibuild) .gt. qugrid%nx-1)bld%iend(ibuild)=qugrid%nx-1
                     if(jstart(ibuild) .lt. 1)jstart(ibuild)=1
                     if(bld%jend(ibuild) .gt. qugrid%ny-1)bld%jend(ibuild)=qugrid%ny-1
!erp 7/23/03 do k=1,bld%kend(ibuild)
!erp  do k=int(bld%zfo(ibuild)),bld%kend(ibuild)  !erp 7/23/03
! int changed to nint in next line 8-14-06
                     do j=jstart(ibuild),bld%jend(ibuild)
                        do i=bld%bld%jstart(ibuild),bld%iend(ibuild)
                           x_c=(real(i)-0.5)*qugrid%dx-xco
                           y_c=(real(j)-0.5)*qugrid%dy-yco
                           thetacell=atan2(y_c,x_c)
                           if(sqrt(x_c**2.+y_c**2.) .le. radius(bld%aa(ibuild),bld%bb(ibuild),&
                                 thetacell,bld%gamma(ibuild)))then
                              select case(bld%btype(ibuild))
                                 case(0)
                                    icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=1
                                    ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                                 case(2)
                                    do k=bld%kstart(ibuild),bld%kend(ibuild)
                                       if(icellflag(i,j,k) .ne. 0)then
                                          if(landuse%canopy_flag .gt. 0)then
                                             if(canopy%top(i,j) .eq. landuse%height(i,j))then
                                                if(k .eq. bld%kstart(ibuild))canopy%atten(i,j,bld%kstart(ibuild):qugrid%nz-1)=0
                                                if(bld%Ht(ibuild) .lt. 0.5*qugrid%dz_array(1))then
                                                   canopy%top(i,j)=0.
                                                else
                                                   canopy%top(i,j)=bld%Ht(ibuild)
                                                endif
                                             elseif(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                             if(zm(k) .lt. canopy%top(i,j))canopy%atten(i,j,k)=atten(ibuild)
                                          else
                                             if(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                             canopy%atten(i,j,k)=atten(ibuild)
                                          endif
                                       endif
                                    enddo
                                 case(3)
                                    ilevel=0
                                    do k=bld%kstart(ibuild),bld%kend(ibuild)
                                       ilevel=ilevel+1
                                       if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                       icellflag(i,j,k)=0
                                       ! ibldflag(i,j,k)=ibuild
                                    enddo
                                 case default
                                    ! icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=0
                                    ! ibldflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=ibuild
                              endselect
                           endif
                        enddo
                     enddo
                  endif
! build a Pentagon shaped building
               case(3)
               case(6)
                  do j=jstart(ibuild),bld%jend(ibuild)
                     y_c=(real(j)-0.5)*qugrid%dy
                     do i=bld%bld%jstart(ibuild),bld%iend(ibuild)
                        x_c=(real(i)-0.5)*qugrid%dx
                        ivert=bld%startidx(ibuild)
                        startpoly=ivert
                        numcrossing=0
                        !Based on Wm. Randolph Franklin, "PNPOLY - Point Inclusion in Polygon Test" Web Page (2000)
                        do while(ivert .lt. bld%stopidx(ibuild))
                           if(((bld%y(ivert) .le. y_c) .and. (bld%y(ivert+1) .gt. y_c)) &
                                 .or. ((bld%y(ivert) .gt. y_c) .and. (bld%y(ivert+1) .le. y_c)))then
                              rayintersect=(y_c-bld%y(ivert))/(bld%y(ivert+1)-bld%y(ivert))
                              if(x_c .lt. bld%x(ivert)+rayintersect*(bld%x(ivert+1)-bld%x(ivert)))then
                                 numcrossing=numcrossing+1
                              endif
                           endif
                           ivert=ivert+1
                           if(bld%x(ivert) .eq. bld%x(startpoly) .and. &
                                 bld%y(ivert) .eq. bld%y(startpoly))then
                              ivert=ivert+1
                              startpoly=ivert
                           endif
                        enddo
                        if(numcrossing/2 .ne. ceiling(0.5*real(numcrossing)))then
                        
                        
                           select case(bld%btype(ibuild))
                              case(0)
                                 icellflag(i,j,bld%kstart(ibuild):bld%kend(ibuild))=1                                 
                              case(2)
                                 do k=bld%kstart(ibuild),bld%kend(ibuild)
                                    if(icellflag(i,j,k) .ne. 0)then
                                       if(landuse%canopy_flag .gt. 0)then
                                          if(canopy%top(i,j) .eq. landuse%height(i,j))then
                                             if(k .eq. bld%kstart(ibuild))canopy%atten(i,j,bld%kstart(ibuild):qugrid%nz-1)=0
                                             if(bld%Ht(ibuild) .lt. 0.5*qugrid%dz_array(1))then
                                                canopy%top(i,j)=0.
                                             else
                                                canopy%top(i,j)=bld%Ht(ibuild)
                                             endif
                                          elseif(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                             canopy%top(i,j)=bld%Ht(ibuild)
                                          endif
                                          if(zm(k) .lt. canopy%top(i,j))canopy%atten(i,j,k)=atten(ibuild)
                                       else
                                          if(bld%Ht(ibuild) .gt. canopy%top(i,j))then
                                             canopy%top(i,j)=bld%Ht(ibuild)
                                          endif
                                          canopy%atten(i,j,k)=atten(ibuild)
                                       endif
                                    endif
                                 enddo
                              case(3)
                                 ilevel=0
                                 do k=bld%kstart(ibuild),bld%kend(ibuild)
                                    ilevel=ilevel+1
                                    if(ilevel/2 .ne. ceiling(0.5*real(ilevel)))cycle
                                    icellflag(i,j,k)=0                                    
                                 enddo
                              case default                                 
                           endselect
                        endif
                     enddo
                  enddo
            endselect
! erp 1/31/2003
         enddo   lp002
         !!$omp end parallel do
         do ibuild=1,bld%number
            if(bld%geometry(ibuild) .eq. 3)call pentagon
         enddo

! end building generation 2a
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         return
      end


!      real function det(Pxi,Pyi,Pxip1,Pyip1,Rx,Ry)
!         implicit none
!         real Pxi,Pxip1,Pyi,Pyip1,Rx,Ry
!         det=(Pxi-Rx)*(Pyip1-Ry)-(Pxip1-Rx)*(Pyi-Ry)
!         return
!      end
