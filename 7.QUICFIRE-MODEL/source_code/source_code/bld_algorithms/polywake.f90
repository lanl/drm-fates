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
      subroutine polywake(ibuild)
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
         real x1,x2,polyarea
         real y1,y2,y3,y4,xw1,yw1,xw3,yw3,tol,zb,ynorm,ynormp,ynormm
         real farwake_exponent,farwake_factor,farwake_velocity,LrAve,seglength,totseglength
         real cav_fac,wake_fac,LoverH,WoverH,eff_height,LrLocal,LrLocalu,LrLocalv
         real canyon_factor,xc,yc,dNu,dNv,dNc,xwall,xu,yu,xv,yv,xp,yp,xwallu,xwallv
         real xw,yw,dNw,xwallw,LrLocalw
         integer x_idx,y_idx,x_idx_min,iu,ju,iv,jv,kk,stop_idx,iw,jw
         integer ktop,kbottom,ivert,endFace,lastFace,nextFace, test1, test2
         real epsilon,bridge_thickness,cavity_velocity,wakesign
         real referenceSpeed,wakeDirA,wakeDir,deltaDir,wakeAngleWeight
         real cavityCrossFactor,cavityAlongFactor,cosCavityDir,sinCavityDir,cosWindDir,sinWindDir
         real rotationFlag
         
         epsilon = 10e-10
         
         
         eff_height=bld%Ht(ibuild)
         xco = bld%cx(ibuild)!CENTER of building in QUIC domain coordinates
         yco = bld%cy(ibuild)
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
         tol=0.01*pi/180.
         farwake_exponent=1.5
         farwake_factor=3
         x1=0.
         x2=0.
         y1=0.
         y2=0.
         polyarea=0.
         stop_idx = 0
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
         if(bld%flag%wake .gt. 1)then
            cav_fac=1.1
            wake_fac=0.1
         else
            cav_fac=1.
            wake_fac=0.
         endif
         bld%LrNode(:)=0.
         bld%LrFace(:)=-1.
         ynormp=0.
         ynormm=0.
         do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
            y1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
            if(y1 .lt. ynormm)ynormm=y1
            if(y1 .gt. ynormp)ynormp=y1
            if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                  .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))exit
         enddo
         ! if(bld%zfo_actual(ibuild) .gt. 0)then
         !    call building_connect
         ! endif
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
         LoverH=bld%Leff(ibuild)/eff_height
         WoverH=bld%Weff(ibuild)/eff_height
         if(LoverH .gt. 3.)LoverH=3.
         if(LoverH .lt. 0.3)LoverH=0.3
         if(WoverH .gt. 10.)WoverH=10.
         bld%Lr(ibuild)=1.8*eff_height*WoverH/((LoverH**(0.3))*(1+0.24*WoverH))
         do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
            xw1=(bld%x(ivert)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*sin(upwind_dir)
            yw1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
            xw3=(bld%x(ivert+1)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*sin(upwind_dir)
            yw3=-(bld%x(ivert+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*cos(upwind_dir)
            upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi  
            bld%FaceWakeDir(ivert)=atan2(bld%y(ivert+1)-bld%y(ivert),bld%x(ivert+1)-bld%x(ivert))-0.5*pi
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(bld%FaceWakeDir(ivert) .le. -pi)bld%FaceWakeDir(ivert)=bld%FaceWakeDir(ivert)+2*pi
            if(abs(upwind_rel) .lt. 0.5*pi)then
               bld%LrFace(ivert)=bld%Lr(ibuild)*cos(upwind_rel)
            endif
            if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                  .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))then
               endFace=ivert
               exit
            endif
         enddo
         LrAve=0.
         totseglength=0.
         do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
            if(bld%LrFace(ivert) .gt. 0.)then
               if(ivert .eq. 1)then
                  lastFace=endFace
                  nextFace=ivert+1
               elseif(ivert .eq. endFace)then
                  lastFace=ivert-1
                  nextFace=bld%startidx(ibuild)
               else
                  lastFace=ivert-1
                  nextFace=ivert+1
               endif
               y1=-(bld%x(lastFace)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(lastFace)-bld%cy(ibuild))*cos(upwind_dir)
               y2=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
               y3=-(bld%x(nextFace)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(nextFace)-bld%cy(ibuild))*cos(upwind_dir)
               y4=-(bld%x(nextFace+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(nextFace+1)-bld%cy(ibuild))*cos(upwind_dir)
               if(bld%LrFace(lastFace) .lt. 0. .and. bld%LrFace(nextFace) .lt. 0.)then
                  bld%LrNode(ivert)=bld%LrFace(ivert)
                  bld%LrNode(ivert+1)=bld%LrFace(ivert)
               elseif(bld%LrFace(lastFace) .lt. 0.)then
                  bld%LrNode(ivert)=bld%LrFace(ivert)
                  bld%LrNode(ivert+1)=((y3-y4)*bld%LrFace(nextFace)+(y2-y3)*bld%LrFace(ivert))/(y2-y4)
               elseif(bld%LrFace(nextFace) .lt. 0.)then
                  bld%LrNode(ivert)=((y2-y3)*bld%LrFace(ivert)+(y1-y2)*bld%LrFace(lastFace))/(y1-y3)
                  bld%LrNode(ivert+1)=bld%LrFace(ivert)
               else
                  bld%LrNode(ivert)=((y2-y3)*bld%LrFace(ivert)+(y1-y2)*bld%LrFace(lastFace))/(y1-y3)
                  bld%LrNode(ivert+1)=((y3-y4)*bld%LrFace(nextFace)+(y2-y3)*bld%LrFace(ivert))/(y2-y4)
               endif
               seglength=y2-y3
               LrAve=LrAve+bld%LrFace(ivert)*seglength
               totseglength=totseglength+seglength
            endif
            
            if(ivert .eq. bld%stopidx(ibuild) .or. (ivert .gt. bld%startidx(ibuild) .and. &
                  bld%x(ivert+1) .gt. bld%x(bld%startidx(ibuild)) - .01*qugrid%dx .and. &
                  bld%x(ivert+1) .lt. bld%x(bld%startidx(ibuild)) + .01*qugrid%dx .and. &
                  bld%y(ivert+1) .gt. bld%y(bld%startidx(ibuild)) - .01*qugrid%dy .and. &
                  bld%y(ivert+1) .lt. bld%y(bld%startidx(ibuild)) + .01*qugrid%dy))then
               stop_idx=ivert
               exit
            endif
         enddo
         bld%Lr(ibuild)=LrAve/totseglength
         do k=bld%kstart(ibuild),bld%kend(ibuild)
            kk=k
            if(0.75*(bld%Ht(ibuild)-bld%zfo_actual(ibuild))+bld%zfo_actual(ibuild) .le. zm(k))exit
         enddo
         xco=bld%cx(ibuild)
         yco=bld%cy(ibuild)
lp004:   do k=ktop,kbottom,-1
            if(bld%btype(ibuild) .eq. 4)then
               zb=zm(k)-(bridge_thickness+bld%zfo(ibuild))
            else
               zb=zm(k)-bld%zfo(ibuild)
            endif
            bld%ufarwake=quwinds%uo(:,:,k)
            bld%vfarwake=quwinds%vo(:,:,k)
            bld%icellwake=bld%icellflag(:,:,k)
            test1 = ibuild
            test2 = bld%startidx(ibuild)
            !!$omp parallel do private(xw1,yw1,xw3,yw3,upwind_rel,perpendicular_flag, &
            !!$omp y_idx,yc,LrLocal,xwall,ynorm,canyon_factor,x_idx_min,x_idx, &
            !!$omp xc,i,j,uwakeflag,vwakeflag,iu,ju,xp,yp,xu,yu,LrLocalu,xwallu,dNu,farwake_velocity, &
            !!$omp iv,jv,xv,yv,LrLocalv,xwallv,dNv,cavity_velocity,wakesign) &
            !!$omp firstprivate(ibuild)
            
            
            
lp003:      do ivert=bld%startidx(ibuild),stop_idx
               xw1=(bld%x(ivert)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*sin(upwind_dir)
               yw1=-(bld%x(ivert)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert)-bld%cy(ibuild))*cos(upwind_dir)
               xw3=(bld%x(ivert+1)-bld%cx(ibuild))*cos(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*sin(upwind_dir)
               yw3=-(bld%x(ivert+1)-bld%cx(ibuild))*sin(upwind_dir)+(bld%y(ivert+1)-bld%cy(ibuild))*cos(upwind_dir)
               upwind_rel=atan2(yw3-yw1,xw3-xw1)+0.5*pi
               if(upwind_rel .gt. pi)upwind_rel=upwind_rel-2*pi
               if(abs(upwind_rel) .lt. 0.5*pi)then
                  cosCavityDir=cos(bld%FaceWakeDir(ivert))
                  sinCavityDir=sin(bld%FaceWakeDir(ivert))
                  !uo_h=referenceSpeed*cos(bld%FaceWakeDir(ivert))
                  !vo_h=referenceSpeed*sin(bld%FaceWakeDir(ivert))
                  if(abs(upwind_rel) .lt. tol)then
                     perpendicular_flag=1
                     wakeDir=upwind_dir
                  else
                     perpendicular_flag=0
                     if(upwind_rel .gt. 0.)then
                        wakeDirA=bld%FaceWakeDir(ivert)+0.5*pi
                     else
                        wakeDirA=bld%FaceWakeDir(ivert)-0.5*pi
                     endif
                     if(wakeDirA .gt. pi)wakeDirA=wakeDirA-2.*pi
                     if(wakeDirA .gt. pi)wakeDirA=wakeDirA-2.*pi
                     deltaDir=abs(wakeDirA-upwind_dir)
                     if(deltaDir .gt. 0.5*pi)then
                        if(upwind_dir .ge. 0.)then
                           wakeDirA=wakeDirA+2.*pi
                        else
                           wakeDirA=wakeDirA-2.*pi
                        endif
                        deltaDir=abs(wakeDirA-upwind_dir)
                     endif
                     wakeAngleWeight=1-(4*abs(deltaDir-0.25*pi)/pi)**2.
                     wakeDir=wakeDirA*wakeAngleWeight+upwind_dir*(1-wakeAngleWeight)
                  endif
                  
lp002:            do y_idx=0,ceiling(2.*abs(yw1-yw3)/qugrid%dxy)+1
                     yc=yw1-0.5*real(y_idx)*qugrid%dxy
                     if(yc .lt. yw3)cycle !yc=yw3
                     LrLocal=bld%LrNode(ivert)+(yc-yw1)*(bld%LrNode(ivert+1)-bld%LrNode(ivert))/(yw3-yw1)
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
                     canyon_factor=1.
                     rotationFlag=1.
                     x_idx_min=-1
                     if(abs(yc) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                           .and. zb < eff_height .and. eff_height > epsilon)then
                        dNc = (1.-(yc/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocal)**2
                        dNc=sqrt(dNc)
                     else
                     	dNc=0.
                     endif
                     do x_idx=1,ceiling(max(bld%Lr(ibuild),farwake_factor*dNc)/qugrid%dxy)
                        xc=real(x_idx)*qugrid%dxy
                        i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                        j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                        if(i .ge. qugrid%nx-1 .or. i .le. 1 .or. j .ge. qugrid%ny-1 .or. j .le. 1)then
                           exit
                        endif
                        if(bld%icellflag(i,j,kk) .ne. 0 .and. x_idx_min .lt. 0)then
                           x_idx_min=x_idx
                        endif
                        if(bld%icellflag(i,j,kk) .eq. 0  .and. bld%flag%isbld(i,j,kk) .ne. &
                           bld%num(ibuild) .and. x_idx_min .gt. 0)then
                           !if(xc .le. bld%Lr(ibuild))then
                              ! canyon_factor=xc/bld%Lr(ibuild)  ! turn off shortening of the wake and cavity in streetcanyons
                           !endif
                           rotationFlag=0.
                           exit
                        endif
                     enddo
                     x_idx_min=-1
lp001:               do x_idx=1,2*ceiling(farwake_factor*LrLocal/qugrid%dxy)
                        uwakeflag=1
                        vwakeflag=1
                        wwakeflag=1
                        xc=0.5*real(x_idx)*qugrid%dxy
                        i=int(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                        j=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
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
                              elseif(bld%icellflag(i,j,kk) .eq. 0)then
                                 exit
                              endif
                           endif
                        endif
! u values
! Far wake
                        if(bld%icellflag(i,j,k) .ne. 0)then
                           iu=nint(((xc+xwall)*cos(upwind_dir)-yc*sin(upwind_dir)+xco)/qugrid%dx)+1
                           ju=int(((xc+xwall)*sin(upwind_dir)+yc*cos(upwind_dir)+yco)/qugrid%dy)+1
                           xp=real(iu-1)*qugrid%dx-xco
                           yp=(real(ju)-0.5)*qugrid%dy-yco
                           xu=xp*cos(upwind_dir)+yp*sin(upwind_dir)
                           yu=-xp*sin(upwind_dir)+yp*cos(upwind_dir)
                           LrLocalu=bld%LrNode(ivert)+(yu-yw1)*(bld%LrNode(ivert+1)-bld%LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallu=xw1
                           else
                              xwallu=((xw3-xw1)/(yw3-yw1))*(yu-yw1)+xw1
                           endif
                           xu=xu-xwallu
                           
                           if(abs(yu) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                          	   dNu = (1.-(yu/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalu)**2
                          	   dNu=sqrt(dNu)
                           else
                           	 dNu=0.
                           endif
                           if(xu .gt. farwake_factor*dNu)uwakeflag=0
                           if(dNu .gt. 0. .and. uwakeflag .eq. 1 .and. yu .le. yw1 .and. yu .ge. yw3 .and. &
                                 bld%icellflag(iu,ju,k) .ne. 0)then
                              if(xu .gt. dNu)then
                                 cavityCrossFactor=(1-xu/(farwake_factor*dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 farwake_velocity=uProfile(k)*(1.-(dNu/(xu+wake_fac*dNu))**(farwake_exponent)) & 
                                    - rotationFlag*referenceSpeed*cavityCrossFactor*sinWindDir*(sinCavityDir*cosWindDir&
                                    &-cosCavityDir*sinWindDir)
                                 if(canyon_factor .eq. 1.)then ! bld%icellflag(iu,ju,k) .ne. 4
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
                                 !cavity_velocity=referenceSpeed*cosCavityDir*min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 !cavityCrossFactor=min(sqrt(1.-abs(yu/ynorm)),1.)
                                 if(rotationFlag .eq. 1.)then
                                    cavityCrossFactor=(1-xu/(farwake_factor*dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 else
                                    cavityCrossFactor=(1-xu/(dNu))*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 endif
                                 cavityAlongFactor=min((1.-xu/(cav_fac*dNu))**2.,1.)*min(sqrt(1.-abs(yu/ynorm)),1.)
                                 cavity_velocity=referenceSpeed*(cavityAlongFactor*cosWindDir*(cosWindDir*cosCavityDir +&
                                       sinCavityDir*sinWindDir) - cavityCrossFactor*sinWindDir*(sinCavityDir*cosWindDir-&
                                       &cosCavityDir*sinWindDir))
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
                                    write(325,*)'Parameterized U exceeds max in polygon wake',&
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
                           LrLocalv=bld%LrNode(ivert)+(yv-yw1)*(bld%LrNode(ivert+1)-bld%LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallv=xw1
                           else
                              xwallv=((xw3-xw1)/(yw3-yw1))*(yv-yw1)+xw1
                           endif
                           xv=xv-xwallv
                           
                           if(abs(yv) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                           	dNv = (1.-(yv/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalv)**2
                           	dNv=sqrt(dNv)
                           else
                           	dNv=0.
                           endif
                           if(xv .gt. farwake_factor*dNv)vwakeflag=0
                           if(dNv .gt. 0. .and. vwakeflag .eq. 1 .and. yv .le. yw1 .and. yv .ge. yw3 .and. &
                                 bld%icellflag(iv,jv,k) .ne. 0)then
                              if(xv .gt. dNv)then
                                 cavityCrossFactor=(1-xv/(farwake_factor*dNv))*min(sqrt(1.-abs(yv/ynorm)),1.)
                                 farwake_velocity=bld%vProfile(k)*(1.-(dNv/(xv+wake_fac*dNv))**(farwake_exponent)) &
                                    + rotationFlag*referenceSpeed*cavityCrossFactor*cosWindDir*(sinCavityDir*cosWindDir-&
                                     &cosCavityDir*sinWindDir)
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
                                    !if(bld%icellflag(i,j,k) .ne. 0)bld%icellflag(i,j,k)=5
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
                                 cavity_velocity=referenceSpeed*(cavityAlongFactor*sinWindDir*(cosWindDir*cosCavityDir + &
                                       sinCavityDir*sinWindDir) + cavityCrossFactor*cosWindDir*(sinCavityDir*cosWindDir-&
                                       &cosCavityDir*sinWindDir))
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
                                    write(325,*)'Parameterized V exceeds max in polygon wake',&
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
                           LrLocalw=bld%LrNode(ivert)+(yw-yw1)*(bld%LrNode(ivert+1)-bld%LrNode(ivert))/(yw3-yw1)
                           if(perpendicular_flag .gt. 0)then
                              xwallw=xw1
                           else
                              xwallw=((xw3-xw1)/(yw3-yw1))*(yw-yw1)+xw1
                           endif
                           xw=xw-xwallw
                           
                           if(abs(yw) < abs(ynorm)  .and. abs(ynorm) > epsilon & 
                                 .and. zb < eff_height .and. eff_height > epsilon)then
                           	dNw = sqrt((1.-(yw/ynorm)**2.)*(1.-((zb)/eff_height)**2.)*(canyon_factor*LrLocalw)**2)
                           else
                           	dNw=0.
                           endif
                           if(xw .gt. farwake_factor*dNw)wwakeflag=0
                           if(dNw .gt. 0. .and. wwakeflag .eq. 1 .and. yw .le. yw1 .and. yw .ge. yw3 .and. &
                                 bld%icellflag(iw,jw,k) .ne. 0)then
                              if(xw .gt. dNw)then
                                 if(canyon_factor .eq. 1.)then
                                    ! if(bld%icellflag(iv,jv,k) .eq. 4)then
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
                                 bld%icellflag(iw,jw,k)=30 ! 4 is the old value
                              endif  
                           endif
                        elseif(x_idx_min .ge. 0  .and. Ht(bld%num(bld%flag%isbld(i,j,k))) .ge. 0.9*eff_height)then
                           exit
                        endif
                        if(uwakeflag .eq. 0 .and. vwakeflag .eq. 0 .and. wwakeflag .eq. 0)exit
                     enddo   lp001      
                  enddo   lp002
               endif
            enddo   lp003
            !!$omp end parallel do
         enddo   lp004
         return
      end
