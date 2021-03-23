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
      subroutine streetcanyon
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! STREETCANYON applies the canyon vortices for building separations
! that are less than the critical value for skimming flow.
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
			use bld_module
			use constants
         use flags_module
         use grid_module
         use winds_module
         
         implicit none
         
         integer perpendicular_flag,x_idx,y_idx,x_idx_max,canyon_flag
			integer :: i,j,k
         integer k_ref,top_flag,dbuild,ic,jc,circle_flag,iu,ju,iv,jv
         integer x_idx_min,reverse_flag
         real uo_h,vo_h,upwind_dir,upwind_rel,downwind_rel,xco,yco
         real beta,xd,yd,thetad,xcd,ycd,rd,along_dir,cross_dir
         real velmag,canyon_dir,along_mag,cross_mag,usign,xc,yc
         real x1,y1,x2,y2,x3,y3,x4,y4,xwall,xpos,S
         real xu,yu,xv,yv,xw,yw,xp,yp,xwallu,xwallv,xwallw
         real xw1,yw1,xw2,yw2,xw3,yw3,tol
         real thetamax,thetamin,thetai,xnorm_bisect,radius
         real upwind_norm,upwind_norm1,upwind_norm2,angle_tol
         real ucomponent,vcomponent,numu,numv
         integer ivert,jvert,numCanyonSlices,y_idx_min,y_idx_max
         real segmentLength,canyonWidth
         integer, allocatable :: canyonFlagArray(:,:)
         integer kref2,kcanyonref
         real fractionslope,localfraction,canyonheight
         real crossmagreflow,crossmagrefhigh,verticalreflow
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Begin Street Canyon vortex subsection for skimming 
! flow regime 
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
         angle_tol=3*pi/4
         do ibuild=1,bld%number
            ! print*,ibuild,bld%geometry(ibuild)
            if(bld%damage(ibuild) .eq. 2)cycle
            if(bld%btype(ibuild) .eq. 2 .or. bld%btype(ibuild) .eq. 5)cycle
            if(bld%geometry(ibuild) .eq. 6)then
               xco = bld%cx(ibuild)
               yco = bld%cy(ibuild)
            else
               xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))!CENTER of building in QUIC domain coordinates
               yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
            endif
            ! find upwind direction and determine the type of flow regime
            uo_h=quwinds%uo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
            vo_h=quwinds%vo(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1)
            ! if(bld%icellflag(nint(xco/qugrid%dx),nint(yco/qugrid%dy),bld%kend(ibuild)+1) .eq. 4)exit ! eliminate step down canyons
            upwind_norm=0.
            upwind_norm1=0.
            upwind_norm2=0.
            upwind_dir=atan2(vo_h,uo_h)
            upwind_rel=upwind_dir-bld%gamma(ibuild)
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
            select case(bld%flag%streetcanyon)
               case(5)
                  canyonWidth=bld%Ht(ibuild)/0.65
               CASE DEFAULT
                  canyonWidth=bld%Lr(ibuild)
            end select
            select case(bld%geometry(ibuild))
               case(1,4)
                  tol=0.01*pi/180.
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
                     xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     upwind_norm1=pi+bld%gamma(ibuild)
                     upwind_norm2=0.5*pi+bld%gamma(ibuild)
                     perpendicular_flag=0
                     usign=1
                  elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
                     xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     xw3=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw3=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     upwind_norm=0.5*pi+bld%gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
                     xw1=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw1=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     xw2=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw2=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     upwind_norm1=0.5*pi+bld%gamma(ibuild)
                     upwind_norm2=bld%gamma(ibuild)
                     perpendicular_flag=0
                     usign=-1
                  elseif(abs(upwind_rel) .le. tol)then
                     xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw3=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw3=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     upwind_norm=bld%gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
                     xw1=x3*cos(upwind_dir)+y3*sin(upwind_dir)
                     yw1=-x3*sin(upwind_dir)+y3*cos(upwind_dir)
                     xw2=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw2=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     upwind_norm1=bld%gamma(ibuild)
                     upwind_norm2=-0.5*pi+bld%gamma(ibuild)
                     perpendicular_flag=0
                     usign=1
                  elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
                     xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw3=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw3=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     upwind_norm=-0.5*pi+bld%gamma(ibuild)
                     perpendicular_flag=1
                  elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
                     xw1=x2*cos(upwind_dir)+y2*sin(upwind_dir)
                     yw1=-x2*sin(upwind_dir)+y2*cos(upwind_dir)
                     xw2=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw2=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     upwind_norm1=-0.5*pi+bld%gamma(ibuild)
                     upwind_norm2=-pi+bld%gamma(ibuild)
                     perpendicular_flag=0
                     usign=-1
                  else
                     xw1=x1*cos(upwind_dir)+y1*sin(upwind_dir)
                     yw1=-x1*sin(upwind_dir)+y1*cos(upwind_dir)
                     xw3=x4*cos(upwind_dir)+y4*sin(upwind_dir)
                     yw3=-x4*sin(upwind_dir)+y4*cos(upwind_dir)
                     upwind_norm=-pi+bld%gamma(ibuild)
                     perpendicular_flag=1
                  endif
                  if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                  if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                  if(upwind_norm1 .gt. pi)upwind_norm1=upwind_norm1-2*pi
                  if(upwind_norm1 .le. -pi)upwind_norm1=upwind_norm1+2*pi
                  if(upwind_norm2 .gt. pi)upwind_norm2=upwind_norm2-2*pi
                  if(upwind_norm2 .le. -pi)upwind_norm2=upwind_norm2+2*pi
                  numCanyonSlices=ceiling(2.*abs(yw1-yw3)/qugrid%dxy)
                  allocate(canyonFlagArray(numCanyonSlices,bld%kend(ibuild)))
                  canyonFlagArray(:,:)=0
                  ! print*,numCanyonSlices
                  do y_idx=1,numCanyonSlices
                     ! print*,y_idx
                     yc=0.5*real(y_idx)*qugrid%dxy+yw3
                     if(yc .gt. yw1)cycle !yc=yw1
                     top_flag=0
                     if(perpendicular_flag .gt. 0)then
                        xwall=xw1
                     elseif(yc.ge.yw2)then
                        xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
                        upwind_norm=upwind_norm1
                     else
                        xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
                        upwind_norm=upwind_norm2
                     endif
                     ! stepupcanyonflag=0 ! apply street canyons only for step up canyons 
                     do k=bld%kend(ibuild),bld%kstart(ibuild),-1
                     ! k=bld%kend(ibuild)+1
                     ! do while(k .ge. bld%kstart(ibuild))
                     !    k=k-1
                        !if(y_idx .eq. int(0.5*real(numCanyonSlices)))print*,1,k
                        canyon_flag=0
                        S=0.
                        x_idx_min=-1
                        reverse_flag=0
                        do x_idx=1,ceiling(2.*canyonWidth/qugrid%dxy)
                           xc=0.5*real(x_idx)*qugrid%dxy
                           i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                           j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                           if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                              exit
                           endif
                           if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                              x_idx_min=x_idx
                           endif
                           if(bld%icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                              canyon_flag=1
                              x_idx_max=x_idx !-1
                              S=0.5*real(x_idx_max-x_idx_min)*qugrid%dxy
                              if(top_flag .eq. 0 .and. S .gt. 0.)then
                                 k_ref=k+1
                                 kref2=k
                                 ic=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                 jc=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                 do while(bld%icellflag(ic,jc,k_ref) .eq. 30 .or. bld%icellflag(ic,jc,k_ref) .eq. 50) ! 4 is the old value
                                    k_ref=k_ref+1
                                 enddo
                                 if(zm(k_ref) .gt. 1.5*zm(bld%kend(ibuild)))k_ref=bld%kend(ibuild)+1
                                 if(bld%icellflag(ic,jc,k_ref) .ne. 0)then
                                    ! k=k_ref-1
                                    !if(y_idx .eq. int(0.5*real(numCanyonSlices)))print*,3,k
                                    numu=0.
                                    numv=0.
                                    ucomponent=0.
                                    vcomponent=0.
                                    if(bld%icellflag(ic-1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+quwinds%uo(ic,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic+1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+quwinds%uo(ic+1,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic,jc-1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+quwinds%vo(ic,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic,jc+1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+quwinds%vo(ic,jc+1,k_ref)
                                    endif
                                    if( ucomponent .ne. 0. .and. numu > 0.) then
                                       ucomponent=ucomponent/numu
                                    else
                                       ucomponent=0.
                                    endif
                                    if( vcomponent .ne. 0. .and. numv > 0.) then
                                       vcomponent=vcomponent/numv
                                    else
                                       vcomponent=0.
                                    endif
                                    if(numu .eq. 0. .and. numv .eq. 0.)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       ! k=kref2
                                       exit
                                    elseif(numu .gt. 0 .and. numv .gt. 0.)then
                                       velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                       canyon_dir=atan2(vcomponent,ucomponent)
                                    elseif(numu .gt. 0)then
                                       velmag=abs(ucomponent)
                                       if(ucomponent .gt. 0.)then
                                          canyon_dir=0.
                                       else
                                          canyon_dir=pi
                                       endif
                                    else
                                       velmag=abs(vcomponent)
                                       if(vcomponent .gt. 0.)then
                                          canyon_dir=0.5*pi
                                       else
                                          canyon_dir=-0.5*pi
                                       endif
                                    endif
                                    top_flag=1
                                    dbuild=bld%invnum(bld%flag%isbld(i,j,kref2))
                                    if(Ht(bld%invnum(bld%flag%isbld(i,j,kref2))) .lt. bld%Ht(ibuild) .and. z(kref2)/S .lt. 0.65)then
                                       ! k=kref2
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    endif
                                 else
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    top_flag=0
                                    S=0.
                                    exit
                                 endif
                                 if(velmag .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized velocity exceeds max in street canyon',&
                                       velmag,quwinds%max_velmag,i,j,k,yc
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                              else
                                 dbuild=bld%invnum(bld%flag%isbld(i,j,k))
                              endif
! Find the along canyon and cross canyon directions
                              i=ceiling(((xc-0.5*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                              j=ceiling(((xc-0.5*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                              if(bld%btype(dbuild) .eq. 4)then
                                 !if(top_flag .gt. 0)k=kref2
                                 canyon_flag=0
                                 S=0.
                                 top_flag=0
                                 exit
                              endif
                              select case(bld%geometry(dbuild))
                                 case(1,4)
                                    beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                    downwind_rel=canyon_dir-bld%gamma(dbuild)
                                    if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                    if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                    xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                    ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                    xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                    yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                       if(downwind_rel .le. 0.)then
                                          if(downwind_rel .le. -0.5*pi)then
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                       if(abs(downwind_rel) .ge. 0.5*pi)then
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                       if(downwind_rel .ge. 0.)then
                                          if(downwind_rel .le. 0.5*pi)then
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    else
                                       if(abs(downwind_rel) .lt. 0.5*pi)then
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    endif
                                 case(2,5)
                                    downwind_rel=canyon_dir-bld%gamma(dbuild)
                                    if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                    if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                    xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                    ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                    xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                    yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                       ((Wt(dbuild)*cos(thetad))**2.))
                                    along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                       (Lt(dbuild)**2.)*rd*sin(thetad))
                                    if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+bld%gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       else
                                          along_dir=along_dir-pi+bld%gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       endif
                                    else
                                       reverse_flag=1
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+bld%gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       else
                                          along_dir=along_dir-pi+bld%gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       endif
                                    endif
                                 case(6)
                                    do jvert=bld%startidx(dbuild),bld%stopidx(dbuild)
                                       cross_dir=atan2(bld%y(jvert+1)-bld%y(jvert),bld%x(jvert+1)-bld%x(jvert))+0.5*pi
                                       if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                       xcd=0.5*(bld%x(jvert+1)+bld%x(jvert))
                                       ycd=0.5*(bld%y(jvert+1)+bld%y(jvert))
                                       xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(cross_dir)+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*sin(cross_dir)
                                       yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(cross_dir)+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*cos(cross_dir)
                                       if(abs(xd) .lt. qugrid%dxy)then
                                          segmentLength=sqrt(((bld%y(jvert+1)-bld%y(jvert))**2.) &
                                             +((bld%x(jvert+1)-bld%x(jvert))**2.))
                                          if(abs(yd) .le. 0.5*segmentLength)then
                                             downwind_rel=canyon_dir-cross_dir
                                             if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                             if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir-0.5*pi
                                                else
                                                   along_dir=cross_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=0
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir+0.5*pi
                                                else
                                                   along_dir=cross_dir-0.5*pi
                                                endif
                                             endif
                                             if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                             if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                             exit
                                          endif
                                       endif
                                       if(bld%x(jvert+1) .eq. bld%x(bld%startidx(dbuild)) &
                                             .and. bld%y(jvert+1) .eq. bld%y(bld%startidx(dbuild)))exit
                                    enddo
                              end select
                              if(along_dir .gt. pi)along_dir=along_dir-2*pi
                              if(along_dir .le. -pi)along_dir=along_dir+2*pi
                              if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                              if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                              if(reverse_flag .eq. 1)then
                                 if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              else
                                 if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              endif
                              exit
                           endif
                        enddo
                        if(k .eq. bld%kend(ibuild))then
                           canyonFlagArray(y_idx,k)=canyon_flag
                        elseif(k .lt. bld%kend(ibuild) .and. canyonFlagArray(y_idx,bld%kend(ibuild)) .gt. 0)then
                           canyonFlagArray(y_idx,k)=canyon_flag
                        endif
                     enddo
                  enddo
                  do y_idx=1,numCanyonSlices
                     yc=0.5*real(y_idx)*qugrid%dxy+yw3
                     if(yc .gt. yw1)cycle !yc=yw1
                     top_flag=0
                     if(perpendicular_flag .gt. 0)then
                        xwall=xw1
                     elseif(yc.ge.yw2)then
                        xwall=((xw2-xw1)/(yw2-yw1))*(yc-yw1)+xw1
                        upwind_norm=upwind_norm1
                     else
                        xwall=((xw3-xw2)/(yw3-yw2))*(yc-yw2)+xw2
                        upwind_norm=upwind_norm2
                     endif
                     ! stepupcanyonflag=0 ! apply street canyons only for step up canyons
                     do k=bld%kend(ibuild),2,-1
                     ! k=bld%kend(ibuild)+1
                     ! do while(k .ge. 2) !bld%kstart(ibuild))
                     !    k=k-1
                        !if(y_idx .eq. int(0.5*real(numCanyonSlices)))print*,2,k
                        if(top_flag .gt. 0 .and. k .gt. bld%kend(ibuild))then
                           canyon_flag=1
                        else
                           canyon_flag=0
                           S=0.
                           x_idx_min=-1
                           reverse_flag=0
                           do x_idx=1,ceiling(2.*canyonWidth/qugrid%dxy)
                              xc=0.5*real(x_idx)*qugrid%dxy
                              i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                              j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                              if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                                 exit
                              endif
                              if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                                 x_idx_min=x_idx
                              endif
                              if(bld%icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                                 canyon_flag=1
                                 x_idx_max=x_idx !-1
                                 S=0.5*real(x_idx_max-x_idx_min)*qugrid%dxy
                                 if(top_flag .eq. 0 .and. S .gt. 0.)then
                                    k_ref=k+1
                                    kref2=k
                                    ic=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                    jc=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                    do while(bld%icellflag(ic,jc,k_ref) .eq. 30 .or. bld%icellflag(ic,jc,k_ref) .eq. 50 .and. k_ref .lt. qugrid%nz) ! 4 is the old value
                                       k_ref=k_ref+1
                                    enddo
                                    if(zm(k_ref) .gt. 1.5*zm(bld%kend(ibuild)))k_ref=bld%kend(ibuild)+1
                                    !print*,ibuild,y_idx,k_ref
                                    if(bld%icellflag(ic,jc,k_ref) .ne. 0)then
                                       ! k=k_ref-1
                                       !if(y_idx .eq. int(0.5*real(numCanyonSlices)))print*,4,k
                                       numu=0.
                                       numv=0.
                                       ucomponent=0.
                                       vcomponent=0.
                                       if(bld%icellflag(ic-1,jc,k_ref) .ne. 0)then
                                          numu=numu+1.
                                          ucomponent=ucomponent+quwinds%uo(ic,jc,k_ref)
                                       endif
                                       if(bld%icellflag(ic+1,jc,k_ref) .ne. 0)then
                                          numu=numu+1.
                                          ucomponent=ucomponent+quwinds%uo(ic+1,jc,k_ref)
                                       endif
                                       if(bld%icellflag(ic,jc-1,k_ref) .ne. 0)then
                                          numv=numv+1.
                                          vcomponent=vcomponent+quwinds%vo(ic,jc,k_ref)
                                       endif
                                       if(bld%icellflag(ic,jc+1,k_ref) .ne. 0)then
                                          numv=numv+1.
                                          vcomponent=vcomponent+quwinds%vo(ic,jc+1,k_ref)
                                       endif
                                       if( ucomponent .ne. 0. .and. numu > 0.) then
                                          ucomponent=ucomponent/numu
                                       else
                                          ucomponent=0.
                                       endif
                                       if( vcomponent .ne. 0. .and. numv > 0.) then
                                          vcomponent=vcomponent/numv
                                       else
                                          vcomponent=0.
                                       endif
                                       if(numu .eq. 0. .and. numv .eq. 0.)then
                                          canyon_flag=0
                                          top_flag=0
                                          S=0.
                                          ! k=kref2
                                          exit
                                       elseif(numu .gt. 0 .and. numv .gt. 0.)then
                                          velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                          canyon_dir=atan2(vcomponent,ucomponent)
                                       elseif(numu .gt. 0)then
                                          velmag=abs(ucomponent)
                                          if(ucomponent .gt. 0.)then
                                             canyon_dir=0.
                                          else
                                             canyon_dir=pi
                                          endif
                                       else
                                          velmag=abs(vcomponent)
                                          if(vcomponent .gt. 0.)then
                                             canyon_dir=0.5*pi
                                          else
                                             canyon_dir=-0.5*pi
                                          endif
                                       endif
                                       top_flag=1
                                       canyonheight=z(k)
                                       dbuild=bld%invnum(bld%flag%isbld(i,j,kref2))
                                       if(Ht(bld%invnum(bld%flag%isbld(i,j,kref2))) .lt. bld%Ht(ibuild) .and. z(kref2)/S .lt. 0.65)then
                                          ! if(top_flag .gt. 0)k=kref2
                                          canyon_flag=0
                                          top_flag=0
                                          S=0.
                                          exit
                                       endif
                                       canyonheight=z(k)
                                    else
                                       ! if(top_flag .gt. 0)k=kref2
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       exit
                                    endif
                                    if(velmag .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       write(325,*)'Parameterized velocity exceeds max in street canyon',&
                                          velmag,quwinds%max_velmag,i,j,k,yc
                                       ! if(top_flag .gt. 0)k=kref2
                                       canyon_flag=0
                                       S=0.
                                       top_flag=0
                                       exit
                                    endif
                                 else
                                    dbuild=bld%invnum(bld%flag%isbld(i,j,k))
                                 endif
! Find the along canyon    and cross canyon directions
                                 i=ceiling(((xc-0.5*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                 j=ceiling(((xc-0.5*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                 if(bld%btype(dbuild) .eq. 4)then
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                                 select case(bld%geometry(dbuild))
                                    case(1,4)
                                       beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                       downwind_rel=canyon_dir-bld%gamma(dbuild)
                                       if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                       if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                       xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                       ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                       xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                       yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                       thetad=atan2(yd,xd)
                                       if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                          if(downwind_rel .le. 0.)then
                                             if(downwind_rel .le. -0.5*pi)then
                                                along_dir=-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(downwind_rel .ge. 0.5*pi)then
                                                along_dir=-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          endif
                                       elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                          if(abs(downwind_rel) .ge. 0.5*pi)then
                                             if(downwind_rel .lt. 0)then
                                                along_dir=-0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(downwind_rel .lt. 0)then
                                                along_dir=-0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          endif
                                       elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                          if(downwind_rel .ge. 0.)then
                                             if(downwind_rel .le. 0.5*pi)then
                                                along_dir=bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(downwind_rel .ge. 0.5*pi)then
                                                along_dir=bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          endif
                                       else
                                          if(abs(downwind_rel) .lt. 0.5*pi)then
                                             if(downwind_rel .ge. 0.)then
                                                along_dir=0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=-0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(downwind_rel .ge. 0.)then
                                                along_dir=0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=-0.5*pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          endif
                                       endif
                                    case(2,5)
                                       downwind_rel=canyon_dir-bld%gamma(dbuild)
                                       if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                       if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                       xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                       ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                       xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                       yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                       thetad=atan2(yd,xd)
                                       rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                          ((Wt(dbuild)*cos(thetad))**2.))
                                       along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                          (Lt(dbuild)**2.)*rd*sin(thetad))
                                       if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                          if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                             along_dir=along_dir+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=along_dir-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                             along_dir=along_dir+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=along_dir-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    case(6)
                                       do jvert=bld%startidx(dbuild),bld%stopidx(dbuild)
                                          cross_dir=atan2(bld%y(jvert+1)-bld%y(jvert),bld%x(jvert+1)-bld%x(jvert))+0.5*pi
                                          if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                          xcd=0.5*(bld%x(jvert+1)+bld%x(jvert))
                                          ycd=0.5*(bld%y(jvert+1)+bld%y(jvert))
                                          xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(cross_dir)+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*sin(cross_dir)
                                          yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(cross_dir)+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*cos(cross_dir)
                                          if(abs(xd) .lt. qugrid%dxy)then
                                             segmentLength=sqrt(((bld%y(jvert+1)-bld%y(jvert))**2.) &
                                                +((bld%x(jvert+1)-bld%x(jvert))**2.))
                                             if(abs(yd) .le. 0.5*segmentLength)then
                                                downwind_rel=canyon_dir-cross_dir
                                                if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                                if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                                if(abs(downwind_rel) .lt. 0.5*pi)then
                                                   reverse_flag=1
                                                   if(downwind_rel .ge. 0.)then
                                                      along_dir=cross_dir-0.5*pi
                                                   else
                                                      along_dir=cross_dir+0.5*pi
                                                   endif
                                                else
                                                   reverse_flag=0
                                                   if(downwind_rel .ge. 0.)then
                                                      along_dir=cross_dir+0.5*pi
                                                   else
                                                      along_dir=cross_dir-0.5*pi
                                                   endif
                                                endif
                                                if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                                if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                                exit
                                             endif
                                          endif
                                          if(bld%x(jvert+1) .eq. bld%x(bld%startidx(dbuild)) &
                                                .and. bld%y(jvert+1) .eq. bld%y(bld%startidx(dbuild)))exit
                                       enddo
                                 end select
                                 if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                 if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                 if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                                 if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                                 if(reverse_flag .eq. 1)then
                                    if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                       ! if(top_flag .gt. 0)k=kref2
                                       canyon_flag=0
                                       S=0.
                                       top_flag=0
                                    endif
                                 else
                                    if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                       ! if(top_flag .gt. 0)k=kref2
                                       canyon_flag=0
                                       S=0.
                                       top_flag=0
                                    endif
                                 endif
                                 exit
                              endif
                           enddo  
                        endif
                        if(canyonFlagArray(y_idx,k) .gt. 0)then
                           y_idx_min=y_idx-2
                           if(y_idx_min .lt. 1)then
                              y_idx_min=1
                           endif
                           y_idx_max=y_idx+2
                           if(y_idx_max .gt. numCanyonSlices)then
                              y_idx_max=numCanyonSlices
                           endif
                           if(sum(canyonFlagArray(y_idx_min:y_idx_max,k)) .lt. 3)then
                              canyon_flag=0
                           endif
                        endif
                        
                        ! if(k .eq. bld%kend(ibuild))then ! apply street canyons only for step up canyons
                        !    stepupcanyonflag=canyon_flag
                        ! elseif(stepupcanyonflag .eq. 0)then
                        !    canyon_flag=0
                        !    S=0.
                        !    top_flag=0
                        ! endif
                        ! if(top_flag .eq. 1)print*,yc,along_dir,cross_dir
                        if(canyon_flag .eq. 1 .and. S .gt. 0.9*qugrid%dxy)then
                           
                           fractionslope=(bld%params%downwindfraction-bld%params%upwindfraction)/S
                           along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                              log(zm(k)/zo)/log(zm(k_ref)/zo)
                           cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                           if(abs(along_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Along canyon exceeds max in street canyon',&
                                 along_mag,quwinds%max_velmag,i,j,k
                           endif
                           if(abs(cross_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Cross canyon exceeds max in street canyon',&
                                 cross_mag,quwinds%max_velmag,i,j,k
                           endif
                           do x_idx=x_idx_min,x_idx_max
                              xc=0.5*real(x_idx)*qugrid%dxy
                              i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                              j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                              
                              if(bld%icellflag(i,j,k) .ne. 0)then ! bld%icellflag(i,j,k) .ne. 6 .and. 
! u component
                                 iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                 ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
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
                                 xpos=xu-xwallu
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then
                                          quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+bld%params%canyonwidthfactor* &
															cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                          ! if(k .eq. bld%kend(ibuild))print*,ju,quwinds%uo(iu,ju,k),cross_mag,cross_dir
                                          if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in street canyon',&
                                                quwinds%uo(iu,ju,k),velmag,quwinds%max_velmag,iu,ju,k
                                          endif
                                       else
                                          crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                          crossmagrefhigh=-bld%params%canyonwidthfactor*bld%params%canyonwidthfactor*cross_mag*cos(cross_dir)
                                          quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                             ((1.-localfraction)*canyonheight))*(zm(k)-localfraction*canyonheight)+crossmagreflow
                                       endif
                                    endif
                                 endif
! v component
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
                                 xpos=xv-xwallv
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then ! 4 is the old value
                                          quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*&
                                          (1.-xpos/S)*sin(cross_dir)
                                          if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in street canyon',&
                                                quwinds%vo(iv,jv,k),velmag,quwinds%max_velmag,iv,jv,k
                                          endif
                                       else
                                          crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                          crossmagrefhigh=-bld%params%canyonwidthfactor*bld%params%canyonwidthfactor*cross_mag*sin(cross_dir)
                                          quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                             ((1.-localfraction)*canyonheight))* (zm(k)-localfraction*canyonheight)+crossmagreflow
                                       endif
                                    endif
                                 endif
! w component
                                 xp=(real(i)-0.5)*qugrid%dx-xco
                                 yp=(real(j)-0.5)*qugrid%dy-yco
                                 xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                 yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                 if(perpendicular_flag .gt. 0)then
                                    xwallw=xw1
                                 elseif(yw.ge.yw2)then
                                    xwallw=((xw2-xw1)/(yw2-yw1))*(yw-yw1)+xw1
                                 else
                                    xwallw=((xw3-xw2)/(yw3-yw2))*(yw-yw2)+xw2
                                 endif
                                 xpos=xw-xwallw
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then
                                          if(bld%icellflag(i,j,k-1) .ne. 0)then
                                             if(reverse_flag .eq. 0)then
                                                quwinds%wo(i,j,k)=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             else
                                                quwinds%wo(i,j,k)=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             endif
                                             if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized W exceeds max in street canyon',&
                                                   quwinds%wo(i,j,k),velmag,quwinds%max_velmag,i,j,k
                                             endif
                                          endif
                                          bld%icellflag(i,j,k)=50 ! 6 is the old value
                                       else
                                          if(reverse_flag .eq. 0)then
                                                verticalreflow=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*&
                                                   (S-xpos)/S)
                                          else
                                                verticalreflow=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*&
                                                   (S-xpos)/S)
                                          endif
                                          quwinds%wo(i,j,k)=((-verticalreflow)/((1.-localfraction)*canyonheight))* &
                                             (z(k-1)-localfraction*canyonheight)+verticalreflow
                                          bld%icellflag(i,j,k)=51 ! 6 is the old value
                                       endif
                                    endif
                                 endif
                              endif
                           enddo
                        endif
                     enddo
                  enddo
                  deallocate(canyonFlagArray)
               case(2,5)
                  tol=0.01*min(qugrid%dx,qugrid%dy)
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
                  endif
                  do y_idx=1,ceiling(2.*abs(yw1-yw3)/qugrid%dxy)
                     yc=0.5*real(y_idx)*qugrid%dxy+yw3
                     if(yc .gt. yw1)cycle !yc=yw1
                     top_flag=0
                     if(circle_flag .eq. 1)then
                        xwall=sqrt((bld%Lt(ibuild)**2.)-(yc**2.))
                        upwind_norm=atan2(xwall*sin(upwind_dir)+yc*cos(upwind_dir),&
                              xwall*cos(upwind_dir)-yc*sin(upwind_dir))
                     else
                        xwall=xnorm_bisect(bld%aa(ibuild),bld%bb(ibuild),&
                              bld%gamma(ibuild)-upwind_dir,yc,thetamin,thetamax,qugrid%dxy)
                        xd=xwall*cos(bld%gamma(ibuild))+yc*sin(bld%gamma(ibuild))
                        yd=-xwall*sin(bld%gamma(ibuild))+yc*cos(bld%gamma(ibuild))
                        thetad=atan2(yd,xd)
                        rd=bld%Lt(ibuild)*bld%Wt(ibuild)/sqrt(((bld%Lt(ibuild)*sin(thetad))**2.)+&
                           ((bld%Wt(ibuild)*cos(thetad))**2.))
                        upwind_norm=atan2(-(bld%Wt(ibuild)**2.)*rd*cos(thetad),&
                           (bld%Lt(ibuild)**2.)*rd*sin(thetad))+bld%gamma(ibuild)+0.5*pi
                        if(cos(upwind_norm-upwind_dir) .le. 0.)upwind_norm=upwind_norm-pi
                        if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                        if(upwind_norm .le. -pi)upwind_norm=upwind_norm+2*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                     endif
                     ! stepupcanyonflag=0 ! apply street canyons only for step up canyons
                     do k=bld%kend(ibuild),2,-1
                     ! k=bld%kend(ibuild)+1
                     ! do while(k .ge. 2) !bld%kstart(ibuild))
                     !    k=k-1
                        canyon_flag=0
                        S=0.
                        x_idx_min=-1
                        reverse_flag=0
                        do x_idx=0,ceiling(2.*canyonWidth/qugrid%dxy)
                           xc=0.5*real(x_idx)*qugrid%dxy
                           i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                           j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                           if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                              exit
                           endif
                           if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                              x_idx_min=x_idx
                           endif
                           if(bld%icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                              canyon_flag=1
                              x_idx_max=x_idx !-1
                              S=0.5*real(x_idx_max-x_idx_min)*qugrid%dxy
                              if(top_flag .eq. 0)then
                                 k_ref=k+1
                                 ! kref2=k
                                 ic=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                 jc=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                 ! do while(bld%icellflag(ic,jc,k_ref) .eq. 30 .or. bld%icellflag(ic,jc,k_ref) .eq. 50) ! 4 is the old value
                                 !    k_ref=k_ref+1
                                 ! enddo
                                 ! if(zm(k_ref) .gt. 1.5*zm(bld%kend(ibuild)))k_ref=bld%kend(ibuild)+1
                                 ! k=k_ref-1
                                 if(bld%icellflag(ic,jc,k_ref) .ne. 0)then
                                    numu=0.
                                    numv=0.
                                    ucomponent=0.
                                    vcomponent=0.
                                    if(bld%icellflag(ic-1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+quwinds%uo(ic,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic+1,jc,k_ref) .ne. 0)then
                                       numu=numu+1.
                                       ucomponent=ucomponent+quwinds%uo(ic+1,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic,jc-1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+quwinds%vo(ic,jc,k_ref)
                                    endif
                                    if(bld%icellflag(ic,jc+1,k_ref) .ne. 0)then
                                       numv=numv+1.
                                       vcomponent=vcomponent+quwinds%vo(ic,jc+1,k_ref)
                                    endif
                                    if(numu > 0.)ucomponent=ucomponent/numu
                                    if(numv > 0.)vcomponent=vcomponent/numv
                                    if(numu == 0. .and. numv == 0.)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       ! k=kref2
                                       exit
                                    elseif(numu > 0. .and. numv > 0.)then
                                       velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                       canyon_dir=atan2(vcomponent,ucomponent)
                                    elseif(numu > 0.)then
                                       velmag=abs(ucomponent)
                                       if(ucomponent .gt. 0.)then
                                          canyon_dir=0.
                                       else
                                          canyon_dir=pi
                                       endif
                                    else
                                       velmag=abs(vcomponent)
                                       if(vcomponent .gt. 0.)then
                                          canyon_dir=0.5*pi
                                       else
                                          canyon_dir=-0.5*pi
                                       endif
                                    endif
                                    top_flag=1
                                    dbuild=bld%invnum(bld%flag%isbld(i,j,k)) ! k vs kref2
                                    if(S <= 0.)then
                                    	canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       ! k=kref2
                                       exit
                                    elseif(Ht(bld%invnum(bld%flag%isbld(i,j,k))) .lt. bld%Ht(ibuild) .and. z(k)/S .lt. 0.65)then
                                       canyon_flag=0
                                       top_flag=0
                                       S=0.
                                       ! k=kref2
                                       exit
                                    endif
                                    canyonheight=z(k)
                                 else
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                                 if(velmag .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Parameterized velocity exceeds max in street canyon',&
                                       velmag,quwinds%max_velmag,i,j,k
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                    exit
                                 endif
                              else
                                 dbuild=bld%invnum(bld%flag%isbld(i,j,k))
                              endif
! Find the along canyon and cross canyon directions
                              i=ceiling(((xc-0.5*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                              j=ceiling(((xc-0.5*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                              if(bld%btype(dbuild) .eq. 4)then
                                 ! if(top_flag .gt. 0)k=kref2
                                 canyon_flag=0
                                 S=0.
                                 top_flag=0
                                 exit
                              endif
                              select case(bld%geometry(dbuild))
                                 case(1,4)
                                    beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                    downwind_rel=canyon_dir-bld%gamma(dbuild)
                                    xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                    ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                    xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                    yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                       if(downwind_rel .le. 0.)then
                                          if(downwind_rel .le. -0.5*pi)then
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                       if(abs(downwind_rel) .ge. 0.5*pi)then
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .lt. 0)then
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                       if(downwind_rel .ge. 0.)then
                                          if(downwind_rel .le. 0.5*pi)then
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.5*pi)then
                                             along_dir=bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       endif
                                    else
                                       if(abs(downwind_rel) .lt. 0.5*pi)then
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          else
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          endif
                                       else
                                          reverse_flag=1
                                          if(downwind_rel .ge. 0.)then
                                             along_dir=0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir-0.5*pi
                                          else
                                             along_dir=-0.5*pi+bld%gamma(dbuild)
                                             cross_dir=along_dir+0.5*pi
                                          endif
                                       endif
                                    endif
                                 case(2,5)
                                    downwind_rel=canyon_dir-bld%gamma(dbuild)
                                    xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                    ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                    xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                    yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                       ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                    thetad=atan2(yd,xd)
                                    rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                       ((Wt(dbuild)*cos(thetad))**2.))
                                    along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                       (Lt(dbuild)**2.)*rd*sin(thetad))
                                    if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+bld%gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       else
                                          along_dir=along_dir-pi+bld%gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       endif
                                    else
                                       reverse_flag=1
                                       if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                          along_dir=along_dir+bld%gamma(dbuild)
                                          cross_dir=along_dir-0.5*pi
                                       else
                                          along_dir=along_dir-pi+bld%gamma(dbuild)
                                          cross_dir=along_dir+0.5*pi
                                       endif
                                    endif
                                 case(6)
                                    do jvert=bld%startidx(dbuild),bld%stopidx(dbuild)
                                       cross_dir=atan2(bld%y(jvert+1)-bld%y(jvert),bld%x(jvert+1)-bld%x(jvert))+0.5*pi
                                       if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                       xcd=0.5*(bld%x(jvert+1)+bld%x(jvert))
                                       ycd=0.5*(bld%y(jvert+1)+bld%y(jvert))
                                       xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(cross_dir)+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*sin(cross_dir)
                                       yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(cross_dir)+&
                                          ((real(j)-0.5)*qugrid%dy-ycd)*cos(cross_dir)
                                       if(abs(xd) .lt. qugrid%dxy)then
                                          segmentLength=sqrt(((bld%y(jvert+1)-bld%y(jvert))**2.) &
                                             +((bld%x(jvert+1)-bld%x(jvert))**2.))
                                          if(abs(yd) .le. 0.5*segmentLength)then
                                             downwind_rel=canyon_dir-cross_dir
                                             if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                             if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir-0.5*pi
                                                else
                                                   along_dir=cross_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=0
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=cross_dir+0.5*pi
                                                else
                                                   along_dir=cross_dir-0.5*pi
                                                endif
                                             endif
                                             if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                             if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                             exit
                                          endif
                                       endif
                                       if(bld%x(jvert+1) .eq. bld%x(bld%startidx(dbuild)) &
                                             .and. bld%y(jvert+1) .eq. bld%y(bld%startidx(dbuild)))exit
                                    enddo
                              end select
                              if(along_dir .gt. pi)along_dir=along_dir-2*pi
                              if(along_dir .le. -pi)along_dir=along_dir+2*pi
                              if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                              if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                              if(reverse_flag .eq. 1)then
                                 if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              else
                                 if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                    ! if(top_flag .gt. 0)k=kref2
                                    canyon_flag=0
                                    S=0.
                                    top_flag=0
                                 endif
                              endif
                              exit
                           endif
                        enddo
                        ! if(k .eq. bld%kend(ibuild))then ! apply street canyons only for step up canyons
                        !    stepupcanyonflag=canyon_flag
                        ! elseif(stepupcanyonflag .eq. 0)then
                        !    canyon_flag=0
                        ! endif
                        if(canyon_flag .eq. 1 .and. S .gt. 0.9*qugrid%dxy)then
                           fractionslope=(bld%params%downwindfraction-bld%params%upwindfraction)/S
                           along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                              log(zm(k)/zo)/log(zm(k_ref)/zo)
                           cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                           if(abs(along_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Along canyon exceeds max in street canyon',&
                                 along_mag,quwinds%max_velmag,i,j,k
                           endif
                           if(abs(cross_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                              write(325,*)'Cross canyon exceeds max in street canyon',&
                                 cross_mag,quwinds%max_velmag,i,j,k
                           endif
                           do x_idx=x_idx_min,x_idx_max
                              xc=0.5*real(x_idx)*qugrid%dxy
                              i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                              j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                              if(bld%icellflag(i,j,k) .ne. 0)then ! bld%icellflag(i,j,k) .ne. 6 .and. 
! u component
                                 iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                 ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
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
                                 xpos=xu-xwallu
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then 
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then
                                          quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*&
                                             (1.-xpos/S)*cos(cross_dir)
                                          ! if(k .eq. bld%kend(ibuild))print*,ju,quwinds%uo(iu,ju,k),cross_mag,cross_dir
                                          if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in street canyon',&
                                                quwinds%uo(iu,ju,k),velmag,quwinds%max_velmag,iu,ju,k
                                          endif
                                       else
                                          crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                          crossmagrefhigh=-bld%params%canyonwidthfactor*bld%params%canyonwidthfactor*cross_mag*cos(cross_dir)
                                          quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                             ((1.-localfraction)*canyonheight))*(zm(k)-localfraction*canyonheight)+crossmagreflow
                                       endif
                                    endif
                                 endif
                                 ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy .and. &
                                 !       (bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                 !    quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                 !    if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 !       write(325,*)'Parameterized U exceeds max in street canyon',&
                                 !          quwinds%uo(iu,ju,k),quwinds%max_velmag,iu,ju,k
                                 !    endif
                                 ! endif
! v component
                                 iv=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
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
                                 xpos=xv-xwallv
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then ! 4 is the old value
                                          quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*&
                                             (1.-xpos/S)*sin(cross_dir)
                                          if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in street canyon',&
                                                quwinds%vo(iv,jv,k),velmag,quwinds%max_velmag,iv,jv,k
                                          endif
                                       else
                                          crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                          crossmagrefhigh=-bld%params%canyonwidthfactor*bld%params%canyonwidthfactor*cross_mag*sin(cross_dir)
                                          quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                             ((1.-localfraction)*canyonheight))*(zm(k)-localfraction*canyonheight)+crossmagreflow
                                       endif
                                    endif
                                 endif
                                 ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy .and. &
                                 !       (bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                 !    quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                 !    if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 !       write(325,*)'Parameterized V exceeds max in street canyon',&
                                 !          quwinds%vo(iv,jv,k),quwinds%max_velmag,iv,jv,k
                                 !    endif
                                 ! endif
! w component
                                 xp=(real(i)-0.5)*qugrid%dx-xco
                                 yp=(real(j)-0.5)*qugrid%dy-yco
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
                                 xpos=xw-xwallw
                                 if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                    !if((bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                    localfraction=fractionslope*xpos+bld%params%upwindfraction
                                    do kcanyonref=k,1,-1
                                       if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                          bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                          exit
                                       endif
                                    enddo
                                    if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                       if(zm(k) .le. localfraction*canyonheight)then
                                          if(bld%icellflag(i,j,k-1) .ne. 0)then
                                             if(reverse_flag .eq. 0)then
                                                quwinds%wo(i,j,k)=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             else
                                                quwinds%wo(i,j,k)=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                             endif
                                             if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized W exceeds max in street canyon',&
                                                   quwinds%wo(i,j,k),velmag,quwinds%max_velmag,i,j,k
                                             endif
                                          endif
                                          bld%icellflag(i,j,k)=50 ! 6 is the old value
                                       else
                                          if(reverse_flag .eq. 0)then
                                                verticalreflow=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                   (1.-2.*(S-xpos)/S)
                                          else
                                                verticalreflow=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                (1.-2.*(S-xpos)/S)
                                          endif
                                          quwinds%wo(i,j,k)=((-verticalreflow)/((1.-localfraction)*canyonheight))* &
                                             (z(k-1)-localfraction*canyonheight)+verticalreflow
                                          bld%icellflag(i,j,k)=51 ! 6 is the old value
                                       endif
                                    endif
                                 endif
                                 ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                 !    if((bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                 !       if(bld%icellflag(i,j,k-1) .ne. 0)then
                                 !          if(reverse_flag .eq. 0)then
                                 !             quwinds%wo(i,j,k)=-abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                 !          else
                                 !             quwinds%wo(i,j,k)=abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                 !          endif
                                 !          if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                 !             write(325,*)'Parameterized W exceeds max in street canyon',&
                                 !                quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                 !          endif
                                 !       endif
                                 !       bld%icellflag(i,j,k)=50 ! 6 is the old value
                                 !    else
                                 !       bld%icellflag(i,j,k)=51 ! 6 is the old value
                                 !    endif
                                 ! endif
                              endif
                           enddo
                        endif
                     enddo
                  enddo
               case(6)
                  tol=0.01*pi/180.
                  do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
                     xw1=(bld%x(ivert)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*sin(upwind_dir)
                     yw1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
                     xw3=(bld%x(ivert+1)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*sin(upwind_dir)
                     yw3=-(bld%x(ivert+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*cos(upwind_dir)
                     upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
                     if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
                     if(abs(upwind_rel) .lt. 0.5*pi)then
                        if(abs(upwind_rel) .gt. pi-tol .or. abs(upwind_rel) .lt. tol)then
                           perpendicular_flag=1
                        else
                           perpendicular_flag=0
                        endif
                        yw2=min(yw1,yw3)
                        upwind_norm=atan2(bld%y(ivert+1)-bld%y(ivert),bld%x(ivert+1)-bld%x(ivert))+0.5*pi
                        if(upwind_norm .gt. pi)upwind_norm=upwind_norm-2*pi
                        do y_idx=1,ceiling(2.*abs(yw1-yw3)/qugrid%dxy)
                           yc=0.5*real(y_idx)*qugrid%dxy+yw2
                           if(yc .gt. yw1)cycle !yc=yw1
                           top_flag=0
                           if(perpendicular_flag .gt. 0)then
                              xwall=xw1
                           else
                              xwall=((xw3-xw1)/(yw3-yw1))*(yc-yw1)+xw1
                           endif
                           ! stepupcanyonflag=0 ! apply street canyons only for step up canyons
                           do k=bld%kend(ibuild),2,-1
                           ! k=bld%kend(ibuild)+1
                           ! do while(k .ge. 2) !bld%kstart(ibuild))
                           !    k=k-1
                           !    if(k .eq. qugrid%nz-1)exit
                              canyon_flag=0
                              S=0.
                              x_idx_min=-1
                              reverse_flag=0
                              do x_idx=1,ceiling(2.*canyonWidth/qugrid%dxy)
                                 xc=0.5*real(x_idx)*qugrid%dxy
                                 i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                 j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                 if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                                    exit
                                 endif
                                 if(bld%icellflag(i,j,k) .ne. 0 .and. x_idx_min .lt. 0)then
                                    x_idx_min=x_idx
                                 endif
                                 if(bld%icellflag(i,j,k) .eq. 0 .and. x_idx_min .ge. 0)then
                                    canyon_flag=1
                                    x_idx_max=x_idx !-1
                                    S=0.5*real(x_idx_max-x_idx_min)*qugrid%dxy
                                    if(top_flag .eq. 0)then
                                       k_ref=k+1
                                       ! kref2=k
                                       ic=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                       jc=ceiling(((0.25*real(x_idx_max)*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                       ! do while(bld%icellflag(ic,jc,k_ref) .eq. 30 .or. bld%icellflag(ic,jc,k_ref) .eq. 50) ! 4 is the old value
                                       !    k_ref=k_ref+1
                                       ! enddo
                                       ! if(zm(k_ref) .gt. 1.5*zm(bld%kend(ibuild)))k_ref=bld%kend(ibuild)+1
                                       ! print*,ibuild,k,k_ref,kref2
                                       ! k=k_ref-1
                                       if(bld%icellflag(ic,jc,k_ref) .ne. 0)then
                                          numu=0.
                                          numv=0.
                                          ucomponent=0.
                                          vcomponent=0.
                                          if(bld%icellflag(ic-1,jc,k_ref) .ne. 0)then
                                             numu=numu+1.
                                             ucomponent=ucomponent+quwinds%uo(ic,jc,k_ref)
                                          endif
                                          if(bld%icellflag(ic+1,jc,k_ref) .ne. 0)then
                                             numu=numu+1.
                                             ucomponent=ucomponent+quwinds%uo(ic+1,jc,k_ref)
                                          endif
                                          if(bld%icellflag(ic,jc-1,k_ref) .ne. 0)then
                                             numv=numv+1.
                                             vcomponent=vcomponent+quwinds%vo(ic,jc,k_ref)
                                          endif
                                          if(bld%icellflag(ic,jc+1,k_ref) .ne. 0)then
                                             numv=numv+1.
                                             vcomponent=vcomponent+quwinds%vo(ic,jc+1,k_ref)
                                          endif
                                          if( ucomponent .ne. 0. .and. numu > 0.) then
                                             ucomponent=ucomponent/numu
                                          else
                                             ucomponent=0.
                                          endif
                                          if( vcomponent .ne. 0. .and. numv > 0.) then
                                             vcomponent=vcomponent/numv
                                          else
                                             vcomponent=0.
                                          endif
                                          if(numu == 0 .and. numv == 0)then
                                             canyon_flag=0
                                             top_flag=0
                                             S=0.
                                             ! k=kref2
                                             exit
                                          elseif(numu > 0 .and. numv > 0)then
                                             velmag=sqrt((ucomponent*ucomponent)+(vcomponent*vcomponent))
                                             canyon_dir=atan2(vcomponent,ucomponent)
                                          elseif(numu > 0)then
                                             velmag=abs(ucomponent)
                                             if(ucomponent .gt. 0.)then
                                                canyon_dir=0.
                                             else
                                                canyon_dir=pi
                                             endif
                                          else
                                             velmag=abs(vcomponent)
                                             if(vcomponent .gt. 0.)then
                                                canyon_dir=0.5*pi
                                             else
                                                canyon_dir=-0.5*pi
                                             endif
                                          endif
                                          top_flag=1
                                          dbuild=bld%invnum(bld%flag%isbld(i,j,k))   ! k vs kref2
                                          if(abs(S) .gt. 0.) then
                                             if(Ht(bld%invnum(bld%flag%isbld(i,j,k))) .lt. bld%Ht(ibuild) .and. z(k)/S .lt. 0.65)then
                                                canyon_flag=0
                                                top_flag=0
                                                S=0.
                                                ! k=kref2
                                                exit
                                             endif
                                          endif
                                          canyonheight=z(k)
                                       else
                                          canyon_flag=0
                                          top_flag=0
                                          S=0.
                                          ! k=kref2
                                          exit
                                       endif
                                       if(velmag .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                          write(325,*)'Parameterized velocity exceeds max in street canyon',&
                                             velmag,quwinds%max_velmag,i,j,k
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                          ! k=kref2
                                          exit
                                       endif
                                    else
                                       dbuild=bld%invnum(bld%flag%isbld(i,j,k))
                                    endif
! Find the along canyon and cross canyon directions
                                    i=ceiling(((xc-0.5*qugrid%dxy+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                    j=ceiling(((xc-0.5*qugrid%dxy+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                    if(bld%btype(dbuild) .eq. 4)then
                                       ! if(top_flag .gt. 0)k=kref2
                                       canyon_flag=0
                                       S=0.
                                       top_flag=0
                                       exit
                                    endif
                                    select case(bld%geometry(dbuild))
                                       case(1,4)
                                          beta=abs(atan2(Lti(dbuild),Wti(dbuild)))
                                          downwind_rel=canyon_dir-bld%gamma(dbuild)
                                          if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                          if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                          xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                          ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                          xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                          yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                          thetad=atan2(yd,xd)
                                          if(thetad .le. 0.5*pi+beta .and. thetad .ge. 0.5*pi-beta)then
                                             if(downwind_rel .le. 0.)then
                                                if(downwind_rel .le. -0.5*pi)then
                                                   along_dir=-pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.5*pi)then
                                                   along_dir=-pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             endif
                                          elseif(thetad .lt. 0.5*pi-beta .and. thetad .gt. -0.5*pi+beta)then
                                             if(abs(downwind_rel) .ge. 0.5*pi)then
                                                if(downwind_rel .lt. 0)then
                                                   along_dir=-0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .lt. 0)then
                                                   along_dir=-0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             endif
                                          elseif(thetad .le. -0.5*pi+beta .and. thetad .ge. -0.5*pi-beta)then
                                             if(downwind_rel .ge. 0.)then
                                                if(downwind_rel .le. 0.5*pi)then
                                                   along_dir=bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=-pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.5*pi)then
                                                   along_dir=bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=-pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             endif
                                          else
                                             if(abs(downwind_rel) .lt. 0.5*pi)then
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                else
                                                   along_dir=-0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                endif
                                             else
                                                reverse_flag=1
                                                if(downwind_rel .ge. 0.)then
                                                   along_dir=0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir-0.5*pi
                                                else
                                                   along_dir=-0.5*pi+bld%gamma(dbuild)
                                                   cross_dir=along_dir+0.5*pi
                                                endif
                                             endif
                                          endif
                                       case(2,5)
                                          downwind_rel=canyon_dir-bld%gamma(dbuild)
                                          if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                          if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                          xcd=bld%xfo(dbuild)+Lt(dbuild)*cos(bld%gamma(dbuild))
                                          ycd=bld%yfo(dbuild)+Lt(dbuild)*sin(bld%gamma(dbuild))
                                          xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(bld%gamma(dbuild))+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*sin(bld%gamma(dbuild))
                                          yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(bld%gamma(dbuild))+&
                                             ((real(j)-0.5)*qugrid%dy-ycd)*cos(bld%gamma(dbuild))
                                          thetad=atan2(yd,xd)
                                          rd=Lt(dbuild)*Wt(dbuild)/sqrt(((Lt(dbuild)*sin(thetad))**2.)+&
                                             ((Wt(dbuild)*cos(thetad))**2.))
                                          along_dir=atan2(-(Wt(dbuild)**2.)*rd*cos(thetad),&
                                             (Lt(dbuild)**2.)*rd*sin(thetad))
                                          if(cos(upwind_dir-canyon_dir) .ge. 0.)then
                                             if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                                along_dir=along_dir+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             else
                                                along_dir=along_dir-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             endif
                                          else
                                             reverse_flag=1
                                             if(abs(downwind_rel-along_dir) .le. 0.5*pi)then
                                                along_dir=along_dir+bld%gamma(dbuild)
                                                cross_dir=along_dir-0.5*pi
                                             else
                                                along_dir=along_dir-pi+bld%gamma(dbuild)
                                                cross_dir=along_dir+0.5*pi
                                             endif
                                          endif
                                       case(6)
                                          do jvert=bld%startidx(dbuild),bld%stopidx(dbuild)
                                             cross_dir=atan2(bld%y(jvert+1)-bld%y(jvert),bld%x(jvert+1)-bld%x(jvert))+0.5*pi
                                             if(cross_dir .gt. pi)cross_dir=cross_dir-2.*pi
                                             xcd=0.5*(bld%x(jvert+1)+bld%x(jvert))
                                             ycd=0.5*(bld%y(jvert+1)+bld%y(jvert))
                                             xd=((real(i)-0.5)*qugrid%dx-xcd)*cos(cross_dir)+&
                                                ((real(j)-0.5)*qugrid%dy-ycd)*sin(cross_dir)
                                             yd=-((real(i)-0.5)*qugrid%dx-xcd)*sin(cross_dir)+&
                                                ((real(j)-0.5)*qugrid%dy-ycd)*cos(cross_dir)
                                             if(abs(xd) .lt. 0.75*qugrid%dxy)then
                                                segmentLength=sqrt(((bld%y(jvert+1)-bld%y(jvert))**2.) &
                                                   +((bld%x(jvert+1)-bld%x(jvert))**2.))
                                                if(abs(yd) .le. 0.5*segmentLength)then
                                                   downwind_rel=canyon_dir-cross_dir
                                                   if(downwind_rel .gt. pi)downwind_rel=downwind_rel-2*pi
                                                   if(downwind_rel .le. -pi)downwind_rel=downwind_rel+2*pi
                                                   if(abs(downwind_rel) .lt. 0.5*pi)then
                                                      reverse_flag=1
                                                      if(downwind_rel .ge. 0.)then
                                                         along_dir=cross_dir-0.5*pi
                                                      else
                                                         along_dir=cross_dir+0.5*pi
                                                      endif
                                                   else
                                                      reverse_flag=0
                                                      if(downwind_rel .ge. 0.)then
                                                         along_dir=cross_dir+0.5*pi
                                                      else
                                                         along_dir=cross_dir-0.5*pi
                                                      endif
                                                   endif
                                                   if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                                   if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                                   exit
                                                endif
                                             endif
                                             if(bld%x(jvert+1) .eq. bld%x(bld%startidx(dbuild)) &
                                                   .and. bld%y(jvert+1) .eq. bld%y(bld%startidx(dbuild)))exit
                                          enddo
                                    end select
                                    if(along_dir .gt. pi)along_dir=along_dir-2*pi
                                    if(along_dir .le. -pi)along_dir=along_dir+2*pi
                                    if(cross_dir .gt. pi)cross_dir=cross_dir-2*pi
                                    if(cross_dir .le. -pi)cross_dir=cross_dir+2*pi
                                    if(reverse_flag .eq. 1)then
                                       if(cos(cross_dir-upwind_norm) .lt. -cos(angle_tol))then
                                          ! if(top_flag .gt. 0)k=kref2
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                       endif
                                    else
                                       if(cos(cross_dir-upwind_norm) .gt. cos(angle_tol))then
                                          ! if(top_flag .gt. 0)k=kref2
                                          canyon_flag=0
                                          S=0.
                                          top_flag=0
                                       endif
                                    endif
                                    exit
                                 endif
                              enddo
                              ! if(k .eq. bld%kend(ibuild))then ! apply street canyons only for step up canyons
                              !    stepupcanyonflag=canyon_flag
                              ! elseif(stepupcanyonflag .eq. 0)then
                              !    canyon_flag=0
                              ! endif
                              if(canyon_flag .eq. 1 .and. S .gt. 0.9*qugrid%dxy)then
                                 fractionslope=(bld%params%downwindfraction-bld%params%upwindfraction)/S
                                 along_mag=abs(velmag*cos(canyon_dir-along_dir))*&
                                    log(zm(k)/zo)/log(zm(k_ref)/zo)
                                 cross_mag=abs(velmag*cos(canyon_dir-cross_dir))
                                 if(abs(along_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Along canyon exceeds max in street canyon',&
                                       along_mag,quwinds%max_velmag,i,j,k
                                 endif
                                 if(abs(cross_mag) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                    write(325,*)'Cross canyon exceeds max in street canyon',&
                                       cross_mag,quwinds%max_velmag,i,j,k
                                 endif
                                 do x_idx=x_idx_min,x_idx_max
                                    xc=0.5*real(x_idx)*qugrid%dxy
                                    i=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
                                    j=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                    
                                    if(bld%icellflag(i,j,k) .ne. 0)then ! bld%icellflag(i,j,k) .ne. 6 .and. 
! u component
                                       iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                                       ju=ceiling(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)
                                       xp=real(iu-1)*qugrid%dx-xco
                                       yp=(real(ju)-0.5)*qugrid%dy-yco
                                       xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                       yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                       if(perpendicular_flag .gt. 0)then
                                          xwallu=xw1
                                       else
                                          xwallu=((xw3-xw1)/(yw3-yw1))*(yu-yw1)+xw1
                                       endif
                                       xpos=xu-xwallu
                                       if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then 
                                          localfraction=fractionslope*xpos+bld%params%upwindfraction
                                          do kcanyonref=k,1,-1
                                             if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                                bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                                exit
                                             endif
                                          enddo
                                          if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                             if(zm(k) .le. localfraction*canyonheight)then
                                                quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*&
                                                   (1.-xpos/S)*cos(cross_dir)
                                                ! if(k .eq. bld%kend(ibuild))print*,ju,quwinds%uo(iu,ju,k),cross_mag,cross_dir
                                                if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                   write(325,*)'Parameterized U exceeds max in street canyon',&
                                                      quwinds%uo(iu,ju,k),velmag,quwinds%max_velmag,iu,ju,k
                                                endif
                                             else
                                                crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                                crossmagrefhigh=-bld%params%canyonwidthfactor*bld%params%canyonwidthfactor*cross_mag*cos(cross_dir)
                                                quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                                   ((1.-localfraction)*canyonheight))*(zm(k)-localfraction*canyonheight)+&
                                                   &crossmagreflow
                                             endif
                                          endif
                                       endif
                                       ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy .and. &
                                       !       (bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                       !    quwinds%uo(iu,ju,k)=along_mag*cos(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*cos(cross_dir)
                                       !    if(abs(quwinds%uo(iu,ju,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       !       write(325,*)'Parameterized U exceeds max in street canyon',&
                                       !          quwinds%uo(iu,ju,k),velmag,quwinds%max_velmag,iu,ju,k
                                       !    endif
                                       ! endif
! v component
                                       iv=ceiling(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)
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
                                       xpos=xv-xwallv
                                       if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                          localfraction=fractionslope*xpos+bld%params%upwindfraction
                                          do kcanyonref=k,1,-1
                                             if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                                bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                                exit
                                             endif
                                          enddo
                                          if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                             if(zm(k) .le. localfraction*canyonheight)then ! 4 is the old value
                                                quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*&
                                                   (1.-xpos/S)*sin(cross_dir)
                                                if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                   write(325,*)'Parameterized V exceeds max in street canyon',&
                                                      quwinds%vo(iv,jv,k),velmag,quwinds%max_velmag,iv,jv,k
                                                endif
                                             else
                                                crossmagreflow=bld%params%canyonwidthfactor*cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                                crossmagrefhigh=-bld%params%canyonwidthfactor*cross_mag*sin(cross_dir)
                                                quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+((crossmagrefhigh-crossmagreflow)/&
                                                   ((1.-localfraction)*canyonheight))*(zm(k)-localfraction*canyonheight)+&
                                                   &crossmagreflow
                                             endif
                                          endif
                                       endif
                                       ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy .and. &
                                       !       (bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                       !    quwinds%vo(iv,jv,k)=along_mag*sin(along_dir)+cross_mag*(2*xpos/S)*2.*(1.-xpos/S)*sin(cross_dir)
                                       !    if(abs(quwinds%vo(iv,jv,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       !       write(325,*)'Parameterized V exceeds max in street canyon',&
                                       !          quwinds%vo(iv,jv,k),velmag,quwinds%max_velmag,iv,jv,k
                                       !    endif
                                       ! endif
! w component
                                       xp=(real(i)-0.5)*qugrid%dx-xco
                                       yp=(real(j)-0.5)*qugrid%dy-yco
                                       xw=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                                       yw=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                                       if(perpendicular_flag .gt. 0)then
                                          xwallw=xw1
                                       else
                                          xwallw=((xw3-xw1)/(yw3-yw1))*(yw-yw1)+xw1
                                       endif
                                       xpos=xw-xwallw
                                       if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                          localfraction=fractionslope*xpos+bld%params%upwindfraction
                                          do kcanyonref=k,1,-1
                                             if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50 .or.&
                                                bld%icellflag(i,j,kcanyonref) .eq. 0)then
                                                exit
                                             endif
                                          enddo
                                          if(bld%icellflag(i,j,kcanyonref) .eq. 30 .or. bld%icellflag(i,j,kcanyonref) .eq. 50)then
                                             if(zm(k) .le. localfraction*canyonheight)then
                                                if(bld%icellflag(i,j,k-1) .ne. 0)then
                                                   if(reverse_flag .eq. 0)then
                                                      quwinds%wo(i,j,k)=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                          (1.-2.*(S-xpos)/S)
                                                   else
                                                      quwinds%wo(i,j,k)=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                          (1.-2.*(S-xpos)/S)
                                                   endif
                                                   if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                      write(325,*)'Parameterized W exceeds max in street canyon',&
                                                         quwinds%wo(i,j,k),velmag,quwinds%max_velmag,i,j,k
                                                   endif
                                                endif
                                                bld%icellflag(i,j,k)=50 ! 6 is the old value
                                             else
                                                if(reverse_flag .eq. 0)then
                                                      verticalreflow=-bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                         (1.-2.*(S-xpos)/S)
                                                else
                                                      verticalreflow=bld%params%canyonwidthfactor*abs(0.5*cross_mag*(1.-2.*xpos/S))*&
                                                          (1.-2.*(S-xpos)/S)
                                                endif
                                                quwinds%wo(i,j,k)=((-verticalreflow)/((1.-localfraction)*canyonheight))* &
                                                   (z(k-1)-localfraction*canyonheight)+verticalreflow
                                                bld%icellflag(i,j,k)=51 ! 6 is the old value
                                             endif
                                          endif
                                       endif
                                       ! if(xpos .le. S .and. xpos .gt. -0.5*qugrid%dxy)then
                                       !    if((bld%icellflag(i,j,k) .eq. 30 .or. bld%icellflag(i,j,k) .eq. 50))then ! 4 is the old value
                                       !       if(bld%icellflag(i,j,k-1) .ne. 0)then
                                       !          if(reverse_flag .eq. 0)then
                                       !             quwinds%wo(i,j,k)=-abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       !          else
                                       !             quwinds%wo(i,j,k)=abs(0.5*cross_mag*(1.-2.*xpos/S))*(1.-2.*(S-xpos)/S)
                                       !          endif
                                       !          if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                       !             write(325,*)'Parameterized W exceeds max in street canyon',&
                                       !                quwinds%wo(i,j,k),velmag,quwinds%max_velmag,i,j,k
                                       !          endif
                                       !       endif
                                       !       bld%icellflag(i,j,k)=50 ! 6 is the old value
                                       !    else
                                       !       bld%icellflag(i,j,k)=51 ! 6 is the old value
                                       !    endif
                                       ! endif
                                    endif
                                 enddo
                              endif
                           enddo
                        enddo
                     endif
                     if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                        .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))exit
                  enddo
            end select
         enddo
         return
      end
