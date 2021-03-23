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
      subroutine sidewall(ibuild)
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
         
         
         use constants
         use bld_module
         use grid_module
			
         implicit none
			
         integer,intent(IN) :: ibuild
         real xco,yco,uo_h,vo_h,upwind_dir,cosUpwindDir,sinUpwindDir,tol
         real upwind_rel,x1,y1,x2,y2,x3,y3,x4,y4
         real xStartRight,yStartRight,xStartLeft,yStartLeft
         real xEndRight,yEndRight,xEndLeft,yEndLeft
         real eff_height,Bs,Bl,RscaleSide,RcxSide,vd,ypref,invvd
         real wallLength,wallDir,sinWallDir,cosWallDir,invzo
         real x_u,y_u,x_v,y_v,shellWidth,internalBLWidth,shellWidthCalc
         real x,y,xp,yp,xp_u,yp_u,xp_v,yp_v,xp_c,yp_c,dxp
         real uo_right,vo_right,uo_left,vo_left
         integer perpendicular_flag,right_flag,left_flag,i,j,k
         integer iStartRight,jstartRight,iStartLeft,jStartLeft
         integer x_idx,y_idx,ic,jc,ivert,endFace,sidewallPolyFlag,upwindFace
         integer iu,ju,iv,jv
         if(bld%flag%sidewall .eq. 0)then
            return
         endif
         if(bld%geometry(ibuild) .eq. 6)then
            xco=bld%cx(ibuild)
            yco=bld%cy(ibuild)
         else
            xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))! CENTER of building in QUIC domain coordinates
            yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         endif
         ic=ceiling(xco/qugrid%dx)
         jc=ceiling(yco/qugrid%dy)
         uo_h=quwinds%uo(ic,jc,bld%kend(ibuild)+1)
         vo_h=quwinds%vo(ic,jc,bld%kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         cosUpwindDir=cos(upwind_dir)
         sinUpwindDir=sin(upwind_dir)
         tol=bld%params%sidewallAngleRange*pi/180.
         dxp=0.5*min(qugrid%dx,qugrid%dy)
         invzo=1./zo
         if(bld%geometry(ibuild) .eq. 6)then
            eff_height=bld%Ht(ibuild)-bld%zfo_actual(ibuild)
            sidewallPolyFlag=0
            do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
               x1=(bld%x(ivert)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*sin(upwind_dir)
               y1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
               x2=(bld%x(ivert+1)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*sin(upwind_dir)
               y2=-(bld%x(ivert+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*cos(upwind_dir)
               bld%FaceRelWindDir(ivert)=atan2(y2-y1,x2-x1)+0.5*pi
               if(bld%FaceRelWindDir(ivert) .gt. pi)bld%FaceRelWindDir(ivert)=bld%FaceRelWindDir(ivert)-2.*pi
               if(abs(bld%FaceRelWindDir(ivert)) .ge. 0.5*pi-tol .and. &
                     abs(bld%FaceRelWindDir(ivert)) .le. 0.5*pi+tol)then
                  sidewallPolyFlag=1 ! indicates that there are faces that are nominally parallel with the wind
               endif
               if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                     .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))then
                  endFace=ivert
                  exit
               endif
            enddo
            if(sidewallPolyFlag .eq. 1)then
               Bs=min(bld%Weff(ibuild),eff_height)
               BL=max(bld%Weff(ibuild),eff_height)
               RscaleSide = ((Bs**(2./3.))*(BL**(1./3.)))
               RcxSide=(0.9*RscaleSide)
               vd= 0.5*0.22*RscaleSide
               ypref=(vd/sqrt(0.5*RcxSide))
               invvd=1./vd
               do ivert=bld%startidx(ibuild),endFace
                  ! check if the face is perpendicular to the wind
                  if(abs(bld%FaceRelWindDir(ivert)) .ge. 0.5*pi-tol .and. &
                        abs(bld%FaceRelWindDir(ivert)) .le. 0.5*pi+tol)then
                     right_flag=0
                     left_flag=0
                     if(bld%FaceRelWindDir(ivert) .gt. 0.)then ! Left side
                        if(ivert .eq. bld%startidx(ibuild))then
                           upwindFace=endFace
                        else
                           upwindFace=ivert-1
                        endif
                        if(abs(bld%FaceRelWindDir(upwindFace)) .ge. pi-tol)then
                           left_flag=1
                           xStartLeft=bld%x(ivert)
                           yStartLeft=bld%y(ivert)
                           xEndLeft=bld%x(ivert+1)
                           yEndLeft=bld%y(ivert+1)
                           wallLength=sqrt(((xStartLeft-xEndLeft)**2.)+((yStartLeft-yEndLeft)**2.))
                           wallDir=atan2(yEndLeft-yStartLeft,xEndLeft-xStartLeft)
                           cosWallDir=cos(wallDir)
                           sinWallDir=sin(wallDir)
                        endif
                     else
                        if(ivert .eq. endFace)then
                           upwindFace=bld%startidx(ibuild)
                        else
                           upwindFace=ivert+1
                        endif
                        if(abs(bld%FaceRelWindDir(upwindFace)) .ge. pi-tol)then
                           right_flag=1
                           xStartRight=bld%x(ivert+1)
                           yStartRight=bld%y(ivert+1)
                           xEndRight=bld%x(ivert)
                           yEndRight=bld%y(ivert)
                           wallLength=sqrt(((xStartRight-xEndRight)**2.)+((yStartRight-yEndRight)**2.))
                           wallDir=atan2(yEndRight-yStartRight,xEndRight-xStartRight)
                           cosWallDir=cos(wallDir)
                           sinWallDir=sin(wallDir)
                        endif
                     endif
                     do k=bld%kend(ibuild),bld%kstart(ibuild),-1
                        if(right_flag .eq. 1)then
                           iStartRight=ceiling(xStartRight/qugrid%dx)
                           jStartRight=ceiling(yStartRight/qugrid%dy)
                           do j=max(1,jStartRight-1),min(qugrid%ny-1,jStartRight+1)
                              do i=max(1,iStartRight-1),min(qugrid%nx-1,iStartRight+1)
                                 select case(bld%icellflag(i,j,k))
                                    case(0) ! solid cell use to define undisturbed reference velocities
                                       uo_right=quwinds%uo(i,j,k)
                                       vo_right=quwinds%vo(i,j,k)
                                    case(1,8,20,21) ! undisturbed, vegetation, or either of the two upwind cell types
                                       ! do nothing as these are acceptable values
                                    case DEFAULT ! detection of any other flow algorithms cause the algorithm to stop searchin on this side
                                       right_flag=0
                                 end select
                              enddo
                           enddo
                           if(right_flag .eq. 1)then
                              do x_idx=1,ceiling(wallLength/dxp) !min(wallLength,RcxSide)
                                 xp=x_idx*dxp
                                 shellWidth=ypref*sqrt(xp)
                                 do y_idx=1,ceiling(shellWidth/dxp)+2
                                    yp=-dxp*y_idx
                                    x=xStartRight+cosWallDir*xp-sinWallDir*yp
                                    y=yStartRight+sinWallDir*xp+cosWallDir*yp
                                    i=ceiling(x/qugrid%dx)
                                    j=ceiling(y/qugrid%dy)
                                    if(bld%icellflag(i,j,k) .gt. 0)then
                                       ! U values
                                       x_u=real(i-1)*qugrid%dx
                                       y_u=(real(j)-0.5)*qugrid%dy
                                       iu=ceiling(x_u/qugrid%dx)
                                       ju=ceiling(y_u/qugrid%dy)
                                       xp_u=(x_u-xStartRight)*cosWallDir+(y_u-yStartRight)*sinWallDir
                                       yp_u=-(x_u-xStartRight)*sinWallDir+(y_u-yStartRight)*cosWallDir
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_u)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_u .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_u)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_u) .le. shellWidth)then
                                          quwinds%uo(iu,ju,k)=-uo_right*abs((shellWidth-abs(yp_u))*invvd)
                                       elseif(abs(yp_u) .le. internalBLWidth)then
                                          quwinds%uo(iu,ju,k)=uo_right*log((abs(yp_u)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                       endif
                                       ! V values
                                       x_v=(real(i)-0.5)*qugrid%dx
                                       y_v=real(j-1)*qugrid%dy
                                       iv=ceiling(x_v/qugrid%dx)
                                       jv=ceiling(y_v/qugrid%dy)
                                       xp_v=(x_v-xStartRight)*cosWallDir+(y_v-yStartRight)*sinWallDir
                                       yp_v=-(x_v-xStartRight)*sinWallDir+(y_v-yStartRight)*cosWallDir
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_v)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_v .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_v)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_v) .le. shellWidth)then
                                          quwinds%vo(iv,jv,k)=-vo_right*abs((shellWidth-abs(yp_v))*invvd)
                                       elseif(abs(yp_v) .le. internalBLWidth)then
                                          quwinds%vo(iv,jv,k)=vo_right*log((abs(yp_v)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                       endif
                                       ! check for change of celltype
                                       xp_c=(x_v-xStartRight)*cosWallDir+(y_u-yStartRight)*sinWallDir
                                       yp_c=-(x_v-xStartRight)*sinWallDir+(y_u-yStartRight)*cosWallDir
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_c)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_c .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_c)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_c) .le. shellWidth)then
                                          bld%icellflag(i,j,k)=60
                                       elseif(abs(yp_c) .le. internalBLWidth)then
                                          bld%icellflag(i,j,k)=61
                                       endif
                                       ! if(abs(yp_c) .le. shellWidth .or. abs(yp_c) .le. internalBLWidth)then
                                       !    bld%icellflag(i,j,k)=61 ! 11 is the old value
                                       ! endif
                                    endif
                                 enddo
                              enddo
                           endif
                        endif
                        if(left_flag .eq. 1)then
                           iStartLeft=ceiling(xStartLeft/qugrid%dx)
                           jStartLeft=ceiling(yStartLeft/qugrid%dy)
                           do j=max(1,jStartLeft-1),min(qugrid%ny-1,jStartLeft+1)
                              do i=max(1,iStartLeft-1),min(qugrid%nx-1,iStartLeft+1)
                                 select case(bld%icellflag(i,j,k))
                                    case(0) ! solid cell use to define undisturbed reference velocities
                                       uo_left=quwinds%uo(i,j,k)
                                       vo_left=quwinds%vo(i,j,k)
                                    case(1,8,20,21)  ! undisturbed, vegetation, or either of the two upwind cell types
                                       ! do nothing as these are acceptable values
                                    case DEFAULT ! detection of any other flow algorithms cause the algorithm to stop searchin on this side
                                       left_flag=0
                                 end select
                              enddo
                           enddo
                           if(left_flag .eq. 1)then
                              do x_idx=1,ceiling(wallLength/dxp) ! min(wallLength,RcxSide)
                                 xp=x_idx*dxp
                                 shellWidth=ypref*sqrt(xp)
                                 do y_idx=1,ceiling(shellWidth/dxp)+2
                                    yp=dxp*y_idx
                                    x=xStartLeft+cosWallDir*xp-sinWallDir*yp
                                    y=yStartLeft+sinWallDir*xp+cosWallDir*yp
                                    i=ceiling(x/qugrid%dx)
                                    j=ceiling(y/qugrid%dy)
                                    if(bld%icellflag(i,j,k) .gt. 0)then
                                       ! U values
                                       x_u=real(i-1)*qugrid%dx
                                       y_u=(real(j)-0.5)*qugrid%dy
                                       iu=ceiling(x_u/qugrid%dx)
                                       ju=ceiling(y_u/qugrid%dy)
                                       xp_u=(x_u-xStartLeft)*cosWallDir+(y_u-yStartLeft)*sinWallDir
                                       yp_u=-(x_u-xStartLeft)*sinWallDir+(y_u-yStartLeft)*cosWallDir
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_u)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_u .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_u)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_u) .le. shellWidth)then
                                          quwinds%uo(iu,ju,k)=-uo_left*abs((shellWidth-abs(yp_u))*invvd)
                                       elseif(abs(yp_u) .le. internalBLWidth)then
                                          quwinds%uo(iu,ju,k)=uo_left*log((abs(yp_u)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                       endif
                                       ! V values
                                       x_v=(real(i)-0.5)*qugrid%dx
                                       y_v=real(j-1)*qugrid%dy
                                       iv=ceiling(x_v/qugrid%dx)
                                       jv=ceiling(y_v/qugrid%dy)
                                       xp_v=(x_v-xStartLeft)*cosWallDir+(y_v-yStartLeft)*sinWallDir
                                       yp_v=-(x_v-xStartLeft)*sinWallDir+(y_v-yStartLeft)*cosWallDir
                                       ! if(k .eq. 10)print*,xp_u,yp_u,xp_v,yp_v
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_v)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_v .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_v)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_v) .le. shellWidth)then
                                          quwinds%vo(iv,jv,k)=-vo_left*abs((shellWidth-abs(yp_v))*invvd)
                                       elseif(abs(yp_v) .le. internalBLWidth)then
                                          quwinds%vo(iv,jv,k)=vo_left*log((abs(yp_v)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                       endif
                                       ! check for change of celltype
                                       xp_c=(x_v-xStartLeft)*cosWallDir+(y_u-yStartLeft)*sinWallDir
                                       yp_c=-(x_v-xStartLeft)*sinWallDir+(y_u-yStartLeft)*cosWallDir
                                       shellWidthCalc = 1-((0.5*RcxSide-xp_c)/(0.5*RcxSide))**2.
                                       if(shellWidthCalc .gt. 0.)then
                                          shellWidth=vd*sqrt(shellWidthCalc)
                                       else
                                          shellWidth=0.0
                                       endif
                                       if(xp_c .gt. 0.)then
                                          internalBLWidth=ypref*sqrt(xp_c)
                                       else
                                          internalBLWidth=0.
                                       endif
                                       if(abs(yp_c) .le. shellWidth)then
                                          bld%icellflag(i,j,k)=60
                                       elseif(abs(yp_c) .le. internalBLWidth)then
                                          bld%icellflag(i,j,k)=61
                                       endif
                                       ! if(abs(yp_c) .le. shellWidth .or. abs(yp_c) .le. internalBLWidth)then
                                       !    bld%icellflag(i,j,k)= 61 ! 11 is the old value
                                       ! endif
                                    endif
                                 enddo
                              enddo
                           endif
                        endif
                     enddo
                  endif
               enddo
            endif
         else
            if(bld%geometry(ibuild) .eq. 4)then
               eff_height=0.8*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))
            else
               eff_height=bld%Ht(ibuild)-bld%zfo_actual(ibuild)
            endif
            right_flag=1
            left_flag=1
            upwind_rel=upwind_dir-bld%gamma(ibuild)
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
            ! Location of corners relative to the center of the building
            x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))
            y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))
            x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))
            y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))
            x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            perpendicular_flag=0
            if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
               ! sidewalls are ommited
            elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
               perpendicular_flag=1
               xStartRight=x2
               yStartRight=y2
               xEndRight=x3
               yEndRight=y3
               xStartLeft=x1
               yStartLeft=y1
               xEndLeft=x4
               yEndLeft=y4
               wallLength=bld%Wti(ibuild)
            elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
               ! sidewalls are ommited
            elseif(abs(upwind_rel) .le. tol)then
               perpendicular_flag=1
               xStartRight=x1
               yStartRight=y1
               xEndRight=x2
               yEndRight=y2
               xStartLeft=x4
               yStartLeft=y4
               xEndLeft=x3
               yEndLeft=y3
               wallLength=bld%Lti(ibuild)
            elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
               ! sidewalls are ommited
            elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
               perpendicular_flag=1
               xStartRight=x4
               yStartRight=y4
               xEndRight=x1
               yEndRight=y1
               xStartLeft=x3
               yStartLeft=y3
               xEndLeft=x2
               yEndLeft=y2
               wallLength=bld%Wti(ibuild)
            elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
               ! sidewalls are ommited
            else
               perpendicular_flag=1
               xStartRight=x3
               yStartRight=y3
               xEndRight=x4
               yEndRight=y4
               xStartLeft=x2
               yStartLeft=y2
               xEndLeft=x1
               yEndLeft=y1
               wallLength=bld%Lti(ibuild)
            endif
            if(perpendicular_flag .eq. 1)then
               wallDir=atan2(yEndRight-yStartRight,xEndRight-xStartRight)
               cosWallDir=cos(wallDir)
               sinWallDir=sin(wallDir)
               Bs=min(bld%Weff(ibuild),eff_height)
               BL=max(bld%Weff(ibuild),eff_height)
               RscaleSide = ((Bs**(2./3.))*(BL**(1./3.)))
               RcxSide=(0.9*RscaleSide)
               vd= 0.5*0.22*RscaleSide
               ypref=(vd/sqrt(0.5*RcxSide))
               invvd=1./vd
               ! search from the top down on the building and stop searching if anything
               ! other than a front recirculation or undisturbed flow is found
               do k=bld%kend(ibuild),bld%kstart(ibuild),-1
                  if(right_flag .eq. 1)then
                     iStartRight=ceiling(xStartRight/qugrid%dx)
                     jStartRight=ceiling(yStartRight/qugrid%dy)
                     do j=max(1,jStartRight-1),min(qugrid%ny-1,jStartRight+1)
                        do i=max(1,iStartRight-1),min(qugrid%nx-1,iStartRight+1)
                           select case(bld%icellflag(i,j,k))
                              case(0) ! solid cell use to define undisturbed reference velocities
                                 uo_right=quwinds%uo(i,j,k)
                                 vo_right=quwinds%vo(i,j,k)
                              case(1,8,20,21)  ! undisturbed, vegetation, or either of the two upwind cell types
                                 ! do nothing as these are acceptable values
                              case DEFAULT ! detection of any other flow algorithms cause the algorithm to stop searchin on this side
                                 right_flag=0
                           end select
                        enddo
                     enddo
                     if(right_flag .eq. 1)then
                        do x_idx=1,ceiling(wallLength/dxp) ! min(wallLength,RcxSide)
                           xp=x_idx*dxp
                           shellWidth=ypref*sqrt(xp)
                           do y_idx=1,ceiling(shellWidth/dxp)+2
                              yp=-dxp*y_idx
                              x=xStartRight+cosWallDir*xp-sinWallDir*yp
                              y=yStartRight+sinWallDir*xp+cosWallDir*yp
                              i=ceiling(x/qugrid%dx)
                              j=ceiling(y/qugrid%dy)
                              if(bld%icellflag(i,j,k) .gt. 0)then
                                 ! U values
                                 x_u=real(i-1)*qugrid%dx
                                 y_u=(real(j)-0.5)*qugrid%dy
                                 iu=ceiling(x_u/qugrid%dx)
                                 ju=ceiling(y_u/qugrid%dy)
                                 xp_u=(x_u-xStartRight)*cosWallDir+(y_u-yStartRight)*sinWallDir
                                 yp_u=-(x_u-xStartRight)*sinWallDir+(y_u-yStartRight)*cosWallDir
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_u)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_u .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_u)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_u) .le. shellWidth)then
                                    quwinds%uo(iu,ju,k)=-uo_right*abs((shellWidth-abs(yp_u))*invvd)
                                 elseif(abs(yp_u) .le. internalBLWidth)then
                                    quwinds%uo(iu,ju,k)=uo_right*log((abs(yp_u)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                 endif
                                 ! V values
                                 x_v=(real(i)-0.5)*qugrid%dx
                                 y_v=real(j-1)*qugrid%dy
                                 iv=ceiling(x_v/qugrid%dx)
                                 jv=ceiling(y_v/qugrid%dy)
                                 xp_v=(x_v-xStartRight)*cosWallDir+(y_v-yStartRight)*sinWallDir
                                 yp_v=-(x_v-xStartRight)*sinWallDir+(y_v-yStartRight)*cosWallDir
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_v)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_v .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_v)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_v) .le. shellWidth)then
                                    quwinds%vo(iv,jv,k)=-vo_right*abs((shellWidth-abs(yp_v))*invvd)
                                 elseif(abs(yp_v) .le. internalBLWidth)then
                                    quwinds%vo(iv,jv,k)=vo_right*log((abs(yp_v)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                 endif
                                 ! check for change of celltype
                                 xp_c=(x_v-xStartRight)*cosWallDir+(y_u-yStartRight)*sinWallDir
                                 yp_c=-(x_v-xStartRight)*sinWallDir+(y_u-yStartRight)*cosWallDir
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_c)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_c .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_c)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_c) .le. shellWidth)then
                                    bld%icellflag(i,j,k)=60
                                 elseif(abs(yp_c) .le. internalBLWidth)then
                                    bld%icellflag(i,j,k)=61
                                 endif
                                 ! if(abs(yp_c) .le. shellWidth .or. abs(yp_c) .le. internalBLWidth)then
                                 !    bld%icellflag(i,j,k)=61 ! 11 is the old value
                                 ! endif
                              endif
                           enddo
                        enddo
                     endif
                  endif
                  if(left_flag .eq. 1)then
                     iStartLeft=ceiling(xStartLeft/qugrid%dx)
                     jStartLeft=ceiling(yStartLeft/qugrid%dy)
                     do j=max(1,jStartLeft-1),min(qugrid%ny-1,jStartLeft+1)
                        do i=max(1,iStartLeft-1),min(qugrid%nx-1,iStartLeft+1)
                           select case(bld%icellflag(i,j,k))
                              case(0) ! solid cell use to define undisturbed reference velocities
                                 uo_left=quwinds%uo(i,j,k)
                                 vo_left=quwinds%vo(i,j,k)
                              case(1,8,20,21)  ! undisturbed, vegetation, or either of the two upwind cell types
                                 ! do nothing as these are acceptable values
                              case DEFAULT ! detection of any other flow algorithms cause the algorithm to stop searchin on this side
                                 left_flag=0
                           end select
                        enddo
                     enddo
                     if(left_flag .eq. 1)then
                        do x_idx=1,ceiling(wallLength/dxp) ! min(wallLength,RcxSide)
                           xp=x_idx*dxp
                           shellWidth=ypref*sqrt(xp)
                           do y_idx=1,ceiling(shellWidth/dxp)+2
                              yp=dxp*y_idx
                              x=xStartLeft+cosWallDir*xp-sinWallDir*yp
                              y=yStartLeft+sinWallDir*xp+cosWallDir*yp
                              i=ceiling(x/qugrid%dx)
                              j=ceiling(y/qugrid%dy)
                              if(bld%icellflag(i,j,k) .gt. 0)then
                                 ! U values
                                 x_u=real(i-1)*qugrid%dx
                                 y_u=(real(j)-0.5)*qugrid%dy
                                 iu=ceiling(x_u/qugrid%dx)
                                 ju=ceiling(y_u/qugrid%dy)
                                 xp_u=(x_u-xStartLeft)*cosWallDir+(y_u-yStartLeft)*sinWallDir
                                 yp_u=-(x_u-xStartLeft)*sinWallDir+(y_u-yStartLeft)*cosWallDir
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_u)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_u .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_u)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_u) .le. shellWidth)then
                                    quwinds%uo(iu,ju,k)=-uo_left*abs((shellWidth-abs(yp_u))*invvd)
                                 elseif(abs(yp_u) .le. internalBLWidth)then
                                    quwinds%uo(iu,ju,k)=uo_left*log((abs(yp_u)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                 endif
                                 ! V values
                                 x_v=(real(i)-0.5)*qugrid%dx
                                 y_v=real(j-1)*qugrid%dy
                                 iv=ceiling(x_v/qugrid%dx)
                                 jv=ceiling(y_v/qugrid%dy)
                                 xp_v=(x_v-xStartLeft)*cosWallDir+(y_v-yStartLeft)*sinWallDir
                                 yp_v=-(x_v-xStartLeft)*sinWallDir+(y_v-yStartLeft)*cosWallDir
                                 ! if(k .eq. 10)print*,xp_u,yp_u,xp_v,yp_v
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_v)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_v .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_v)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_v) .le. shellWidth)then
                                    quwinds%vo(iv,jv,k)=-vo_left*abs((shellWidth-abs(yp_v))*invvd)
                                 elseif(abs(yp_v) .le. internalBLWidth)then
                                    quwinds%vo(iv,jv,k)=vo_left*log((abs(yp_v)+zo)*invzo)/log((internalBLWidth+zo)*invzo)
                                 endif
                                 ! check for change of celltype
                                 xp_c=(x_v-xStartLeft)*cosWallDir+(y_u-yStartLeft)*sinWallDir
                                 yp_c=-(x_v-xStartLeft)*sinWallDir+(y_u-yStartLeft)*cosWallDir
                                 shellWidthCalc = 1-((0.5*RcxSide-xp_c)/(0.5*RcxSide))**2.
                                 if(shellWidthCalc .gt. 0.)then
                                    shellWidth=vd*sqrt(shellWidthCalc)
                                 else
                                    shellWidth=0.0
                                 endif
                                 if(xp_c .gt. 0.)then
                                    internalBLWidth=ypref*sqrt(xp_c)
                                 else
                                    internalBLWidth=0.
                                 endif
                                 if(abs(yp_c) .le. shellWidth)then
                                    bld%icellflag(i,j,k)=60
                                 elseif(abs(yp_c) .le. internalBLWidth)then
                                    bld%icellflag(i,j,k)=61
                                 endif
                                 ! if(abs(yp_c) .le. shellWidth .or. abs(yp_c) .le. internalBLWidth)then
                                 !    bld%icellflag(i,j,k)=61 ! 11 is the old value
                                 ! endif
                              endif
                           enddo
                        enddo
                     endif
                  endif
               enddo
            endif
         endif
         return
      end