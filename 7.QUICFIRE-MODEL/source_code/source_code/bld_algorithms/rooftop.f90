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
      subroutine rooftop(ibuild)
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
         use constants
         use bld_module
         use flags_module
         use winds_module
         
         implicit none
			
			integer,intent(IN) :: ibuild
         real perpendicular_flag,ns_flag
         integer roofflag_temp,uflag,vflag,wflag,k_ref,kendv,i,j,k
         integer ivert,k_shellu,k_shellv,k_shellw
         real uo_h,vo_h,upwind_dir,upwind_rel,xco,yco
         real x1,y1,x2,y2,x3,y3,x4,y4
         real tol,xfront,yfront
         real zr,x_u,y_u,x_v,y_v,x_w,y_w
         real vd,Bs,BL,roofangle,hx,hy,hdu,hdv,hdwx,hdwy
         real xnorm,ynorm,vel_mag,zref2,zref2u,zref2v
         real shell_heightu,shell_heightv,rotationAngle
         real sinRotationAngle,cosRotationAngle,vel_magu,vel_magv
         real cosUpwindDir,sinUpwindDir,denomu,denomv,invzo,invvd
         real hxu,hyu,hxv,hyv,hxw,hyw,tanRoofangle
         real cosXnorm,sinXnorm,cosXnormpPi,sinXnormpPi
         real cosYnorm,sinYnorm,cosYnormpPi,sinYnormpPi
         real shell_heightu_part, shell_heightv_part,vortexCenterFactor
         
         if(roofflag .eq. 0)then
            bld%rooftop_flag(ibuild)=0
            return
         endif
         ! vortexCenterFactor=0.5 ! place the center of the vortex at the half height
         vortexCenterFactor=1.0 ! old version
         if(bld%geometry(ibuild) .eq. 6)then
            xco=bld%cx(ibuild)
            yco=bld%cy(ibuild)
         else
            xco = bld%xfo(ibuild) + bld%Lt(ibuild)*cos(bld%gamma(ibuild))! CENTER of building in QUIC domain coordinates
            yco = bld%yfo(ibuild) + bld%Lt(ibuild)*sin(bld%gamma(ibuild))
         endif
         ! find upwind direction and determine the type of flow regime
         uo_h=quwinds%uo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
         vo_h=quwinds%vo(nint(xco/dx),nint(yco/dy),bld%kend(ibuild)+1)
         upwind_dir=atan2(vo_h,uo_h)
         cosUpwindDir=cos(upwind_dir)
         sinUpwindDir=sin(upwind_dir)
         tol=30*pi/180.
         invzo=1./zo
         ns_flag=0
         if(bld%geometry(ibuild) .eq. 6)then
            perpendicular_flag=0
            xfront=0.
            yfront=0.
            do ivert=bld%startidx(ibuild),bld%stopidx(ibuild)
               x1=(bld%x(ivert)-bld%cx(ibuild))*cosUpwindDir+(bld%y(ivert)-bld%cy(ibuild))*sinUpwindDir
               y1=-(bld%x(ivert)-bld%cx(ibuild))*sinUpwindDir+(bld%y(ivert)-bld%cy(ibuild))*cosUpwindDir
               if(x1 .lt. xfront)then
                  xfront=x1
                  yfront=y1
               endif
               if(bld%x(ivert+1) .eq. bld%x(bld%startidx(ibuild)) &
                     .and. bld%y(ivert+1) .eq. bld%y(bld%startidx(ibuild)))exit
            enddo
            rotationAngle=upwind_dir
         else
            rotationAngle=bld%gamma(ibuild)
            upwind_rel=upwind_dir-bld%gamma(ibuild)
            if(upwind_rel.gt.pi)upwind_rel=upwind_rel-2*pi
            if(upwind_rel.le.-pi)upwind_rel=upwind_rel+2*pi
            ! Location of corners relative to the center of the building
            x1=bld%xfo(ibuild)+bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
            y1=bld%yfo(ibuild)-bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
            x2=x1+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y2=y1+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            x4=bld%xfo(ibuild)-bld%Wt(ibuild)*sin(bld%gamma(ibuild))-xco
            y4=bld%yfo(ibuild)+bld%Wt(ibuild)*cos(bld%gamma(ibuild))-yco
            x3=x4+bld%Lti(ibuild)*cos(bld%gamma(ibuild))
            y3=y4+bld%Lti(ibuild)*sin(bld%gamma(ibuild))
            perpendicular_flag=0
            if(upwind_rel .gt. 0.5*pi+tol .and. upwind_rel .lt. pi-tol)then
               xfront=bld%Lt(ibuild)
               yfront=-bld%Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(0.5*pi-upwind_rel)-2*abs(0.75*pi-upwind_rel)))
               xnorm=bld%gamma(ibuild)!+roofangle
               ynorm=bld%gamma(ibuild)-0.5*pi!-roofangle
            elseif(upwind_rel .ge. 0.5*pi-tol .and. upwind_rel .le. 0.5*pi+tol)then
               xfront=bld%Lt(ibuild)
               yfront=-bld%Wt(ibuild)
               perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .gt. tol .and. upwind_rel .lt. 0.5*pi-tol)then
               xfront=-bld%Lt(ibuild)
               yfront=-bld%Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(upwind_rel)-2*abs(0.25*pi-upwind_rel)))
               xnorm=bld%gamma(ibuild)+pi!-roofangle
               ynorm=bld%gamma(ibuild)-0.5*pi!+roofangle
            elseif(abs(upwind_rel) .le. tol)then
               xfront=-bld%Lt(ibuild)
               yfront=-bld%Wt(ibuild)
               perpendicular_flag=1
               ns_flag=0
            elseif(upwind_rel .lt. -tol .and. upwind_rel .gt. -0.5*pi+tol)then
               xfront=-bld%Lt(ibuild)
               yfront=bld%Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(upwind_rel)-2.0*abs(-0.25*pi-upwind_rel)))
               xnorm=bld%gamma(ibuild)-pi!+roofangle
               ynorm=bld%gamma(ibuild)+0.5*pi!-roofangle
            elseif(upwind_rel .lt. -0.5*pi+tol .and. upwind_rel .gt. -0.5*pi-tol)then
               xfront=-bld%Lt(ibuild)
               yfront=bld%Wt(ibuild)
               perpendicular_flag=1
               ns_flag=1
            elseif(upwind_rel .lt. -0.5*pi-tol .and. upwind_rel .gt. -pi+tol)then
               xfront=bld%Lt(ibuild)
               yfront=bld%Wt(ibuild)
               perpendicular_flag=0
               roofangle=0.0513*exp(1.7017*(abs(-0.5*pi-upwind_rel)-2.0*abs(-0.75*pi-upwind_rel)))
               xnorm=bld%gamma(ibuild)!-roofangle
               ynorm=bld%gamma(ibuild)+0.5*pi!+roofangle
            else
               xfront=bld%Lt(ibuild)
               yfront=bld%Wt(ibuild)
               perpendicular_flag=1
               ns_flag=0
            endif
            if(perpendicular_flag .lt. 1)then
               tanRoofangle=tan(roofangle)
               cosXnorm=cos(xnorm)
               sinXnorm=sin(xnorm)
               cosXnormpPi=cos(xnorm+pi)
               sinXnormpPi=sin(xnorm+pi)
               cosYnorm=cos(ynorm)
               sinYnorm=sin(ynorm)
               cosYnormpPi=cos(ynorm+pi)
               sinYnormpPi=sin(ynorm+pi)
            endif
         endif
         sinRotationAngle=sin(rotationAngle)
         cosRotationAngle=cos(rotationAngle)
         ! MAN 07/25/2008 stretched vertical grid
         do k=bld%kend(ibuild)+1,nz-1
            k_ref=k
            if(bld%params%vortexHeightFactor*bld%Ht(ibuild) .lt. z(k))exit
         enddo
         if(k_ref .lt. nz)then
            Bs=min(bld%Weff(ibuild),bld%Ht(ibuild))
            BL=max(bld%Weff(ibuild),bld%Ht(ibuild))
            bld%Rscale(ibuild) = ((Bs**(2./3.))*(BL**(1./3.)))
            bld%Rcx(ibuild)=(0.9*bld%Rscale(ibuild))
            vd= 0.5*0.22*bld%Rscale(ibuild)
            zref2=(vd/sqrt(0.5*bld%Rcx(ibuild)))
            invvd=1./(vd*vortexCenterFactor)
            ! MAN 07/25/2008 stretched vertical grid
            do k=bld%kend(ibuild)+1,k_ref
               kendv=k
               if(bld%Ht(ibuild)+bld%params%vortexHeightFactor*vd .lt. zm(k))exit
            enddo
            if(bld%flag%roof .eq. 2 .and. (bld%geometry(ibuild) .eq. 1) &
                  .and. bld%rooftop_flag(ibuild) .eq. 1)then
               roofflag_temp=2
            else
               roofflag_temp=1
               bld%rooftop_flag(ibuild)=0
            endif
            select case(roofflag_temp)
               case(1)
                  !$omp parallel do private(i,j,k,uflag,vflag,wflag,x_u,y_u,x_v,y_v,hx,hy,hdu,hdv, &
                  !$omp zref2u,zref2v,k_shellu,k_shellv,denomu,denomv,vel_magu,vel_magv,zr)
lp002:            do j=bld%jstart(ibuild),bld%jend(ibuild)
lp001:               do i=bld%istart(ibuild),bld%iend(ibuild)
                        uflag=0
                        vflag=0
                        wflag=0
                        !check to see if velocity vector is above the building or in a street canyon cell
                        if(bld%icellflag(i,j,bld%kend(ibuild)) .eq. 0)then
                           uflag=1
                           vflag=1
                           wflag=1
                        else
                           if(bld%icellflag(i-1,j,bld%kend(ibuild)) .eq. 0)uflag=1
                           if(bld%icellflag(i,j-1,bld%kend(ibuild)) .eq. 0)vflag=1
                        endif
                        if(bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 50 .or. bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 30)then ! 6 and 4 are the old values for Streetcanyon and Cavity regions, respectively
                           uflag=0
                           vflag=0
                           wflag=0
                        endif
                        if(uflag+vflag+wflag .gt. 0 .and. bld%icellflag(i,j,bld%kend(ibuild)+1) .gt. 0)then
                           x_u=((real(i)-1)*dx-xco)*cosRotationAngle+ &
                                 ((real(j)-0.5)*dy-yco)*sinRotationAngle
                           y_u=-((real(i)-1)*dx-xco)*sinRotationAngle+ &
                                 ((real(j)-0.5)*dy-yco)*cosRotationAngle
                           x_v=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                 ((real(j)-1)*dy-yco)*sinRotationAngle
                           y_v=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                 ((real(j)-1)*dy-yco)*cosRotationAngle
                           hx=abs(x_u-xfront)
                           hy=abs(y_u-yfront)
                           hdu=(1-perpendicular_flag)*min(hx,hy)+perpendicular_flag*(ns_flag*hy+(1-ns_flag)*hx)
                           hx=abs(x_v-xfront)
                           hy=abs(y_v-yfront)
                           hdv=(1.-perpendicular_flag)*min(hx,hy)+perpendicular_flag*(ns_flag*hy+(1-ns_flag)*hx)
                           zref2u=zref2*sqrt(hdu)
                           zref2v=zref2*sqrt(hdv)
                           k_shellu=1
                           k_shellv=1
                           do k=bld%kend(ibuild),nz-1
                              if(zref2u+bld%Ht(ibuild) .lt. zm(k+1) .and. k_shellu .lt. 1)then
                                 k_shellu=k
                              endif
                              if(zref2v+bld%Ht(ibuild) .lt. zm(k+1) .and. k_shellv .lt. 1)then
                                 k_shellv=k
                              endif
                              if(k_shellu .gt. 0 .and. k_shellv .gt. 0)exit
                           enddo
                           if(k_shellu .le. bld%kend(ibuild))then
                              uflag=0
                              denomu=1.
                           else
                              denomu=1./log((zm(k_shellu)-bld%Ht(ibuild))*invzo)
                           endif
                           if(k_shellv .le. bld%kend(ibuild))then
                              vflag=0
                              denomv=1.
                           else
                              denomv=1./log((zm(k_shellv)-bld%Ht(ibuild))*invzo)
                           endif
                           vel_magu=quwinds%uo_roof(i,j,k_shellu)
                           vel_magv=quwinds%vo_roof(i,j,k_shellv)
                           do k=bld%kend(ibuild)+1,k_ref
                              if(bld%icellflag(i,j,k) .lt. 1)then
                                 exit
                              else
                                 k_shellu=0
                                 k_shellv=0
                                 zr=zm(k)-bld%Ht(ibuild)
                                 if(uflag .gt. 0)then
                                    if(zr .le. zref2u)then
                                       quwinds%uo(i,j,k)=vel_magu*log(zr*invzo)*denomu
                                       if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                          write(325,*)'Parameterized U exceeds max in rooftop',&
                                             quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                       endif
                                       k_shellu=1
                                       if(wflag .eq. 1)then
                                          bld%icellflag(i,j,k)=41 ! 3 is the old value
                                       endif
                                    endif
                                 endif
                                 if(vflag .gt. 0)then
                                    if(zr .le. zref2v)then
                                       quwinds%vo(i,j,k)=vel_magv*log(zr*invzo)*denomv
                                       if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                          write(325,*)'Parameterized V exceeds max in rooftop',&
                                             quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                       endif
                                       k_shellv=1
                                       if(wflag .eq. 1)then
                                          bld%icellflag(i,j,k)=41 ! 3 is the old value
                                       endif
                                    endif
                                 endif
                              endif
                              if(k_shellu+k_shellv .lt. 1)exit
                           enddo
                        endif
                     enddo lp001
                  enddo lp002
                  !$omp end parallel do
               case(2)
                  if(perpendicular_flag .gt. 0)then
                     !$omp parallel do private(i,j,k,uflag,vflag,wflag,x_u,y_u,x_v,y_v,hx,hy,hdu,hdv, &
                     !$omp zref2u,zref2v,k_shellu,k_shellv,shell_heightu,shell_heightv,denomu,denomv, &
                     !$omp vel_magu,vel_magv,zr, shell_heightu_part, shell_heightv_part)
lp005:               do j=bld%jstart(ibuild),bld%jend(ibuild)
lp004:                  do i=bld%istart(ibuild),bld%iend(ibuild)
                           uflag=0
                           vflag=0
                           wflag=0
                           !check to see if velocity vector is above the building or in a street canyon cell
                           if(bld%icellflag(i,j,bld%kend(ibuild)) .eq. 0)then
                              uflag=1
                              vflag=1
                              wflag=1
                           else
                              if(bld%icellflag(i-1,j,bld%kend(ibuild)) .eq. 0)uflag=1
                              if(bld%icellflag(i,j-1,bld%kend(ibuild)) .eq. 0)vflag=1
                           endif
                           if(bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 50 .or. bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 30)then ! 6 and 4 are the old values for Streetcanyon and Cavity regions, respectively
                              uflag=0
                              vflag=0
                              wflag=0
                           endif
                           if(uflag+vflag+wflag .gt. 0 .and. bld%icellflag(i,j,bld%kend(ibuild)+1) .gt. 0)then
                              x_u=((real(i)-1)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_u=-((real(i)-1)*dx-xco)*sinRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              x_v=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-1)*dy-yco)*sinRotationAngle
                              y_v=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-1)*dy-yco)*cosRotationAngle
                              hx=abs(x_u-xfront)
                              hy=abs(y_u-yfront)
                              hdu=(ns_flag*hy+(1-ns_flag)*hx)
                              hx=abs(x_v-xfront)
                              hy=abs(y_v-yfront)
                              hdv=(ns_flag*hy+(1-ns_flag)*hx)
                              zref2u=zref2*sqrt(hdu)
                              zref2v=zref2*sqrt(hdv)
                              k_shellu=0
                              k_shellv=0
                              do k=bld%kend(ibuild),nz-1
                                 if(bld%params%vortexHeightFactor*zref2u+bld%Ht(ibuild) .lt. zm(k+1) .and. k_shellu .lt. 1)then
                                    k_shellu=k
                                 endif
                                 if(bld%params%vortexHeightFactor*zref2v+bld%Ht(ibuild) .lt. zm(k+1) .and. k_shellv .lt. 1)then
                                    k_shellv=k
                                 endif
                                 if(k_shellu .gt. 0 .and. k_shellv .gt. 0)exit
                              enddo
                              shell_heightu_part = 1-((0.5*bld%Rcx(ibuild)-hdu)/(0.5*bld%Rcx(ibuild)))**2.
                              shell_heightv_part = 1-((0.5*bld%Rcx(ibuild)-hdv)/(0.5*bld%Rcx(ibuild)))**2.
                              if(shell_heightu_part .gt. 0)then
                              	shell_heightu=vd*sqrt(1-((0.5*bld%Rcx(ibuild)-hdu)/(0.5*bld%Rcx(ibuild)))**2.)
                              else
                              	shell_heightu=0.
                              endif
                              if(shell_heightv_part .gt. 0.) then
                              	shell_heightv=vd*sqrt(1-((0.5*bld%Rcx(ibuild)-hdv)/(0.5*bld%Rcx(ibuild)))**2.)
                              else
                              	shell_heightv=0.
                              endif
                              
                              if(k_shellu .le. bld%kend(ibuild))then
                                 uflag=0
                                 denomu=1.
                              else
                                 denomu=1./log((zm(k_shellu)-bld%Ht(ibuild))*invzo)
                              endif
                              if(k_shellv .le. bld%kend(ibuild))then
                                 vflag=0
                                 denomv=1.
                              else
                                 denomv=1./log((zm(k_shellv)-bld%Ht(ibuild))*invzo)
                              endif
                              vel_magu=quwinds%uo_roof(i,j,k_shellu)
                              vel_magv=quwinds%vo_roof(i,j,k_shellv)
                              ! if(j .eq. 50)print*,hdu,shell_heightu_part,shell_heightu,k_shellu
                              kendv=max(k_shellu,kendv)
                              kendv=max(kendv,k_shellv)
lp003:                        do k=bld%kend(ibuild)+1,kendv
                                 k_shellu=0
                                 k_shellv=0
                                 zr=zm(k)-bld%Ht(ibuild)
                                 if(bld%icellflag(i,j,k) .lt. 1)then
                                    exit
                                 else
                                    if(uflag .eq. 1)then
                                       if(zr .le. bld%params%vortexHeightFactor*zref2u)then
                                          quwinds%uo(i,j,k)=vel_magu*log(zr*invzo)*denomu
                                          if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in rooftop',&
                                                quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellu=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=41 ! 3 is the old value
                                          endif
                                       endif
                                       if(hdu .lt. 0.5*bld%Rcx(ibuild) .and. zr .le. bld%params%vortexHeightFactor*shell_heightu)then
                                          quwinds%uo(i,j,k)=-quwinds%uo_roof(i,j,k)*(vortexCenterFactor*vd-zr)*invvd
                                          if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in rooftop',&
                                                quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellu=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=40 ! 3 is the old value
                                          endif
                                       elseif(hdu .lt. bld%Rcx(ibuild) .and. zr .le. bld%params%vortexHeightFactor*shell_heightu)then
                                          quwinds%uo(i,j,k)=-quwinds%uo_roof(i,j,k)*(vortexCenterFactor*shell_heightu-zr)*invvd
                                          ! quwinds%uo(i,j,k)=-quwinds%uo_roof(i,j,k)*(vortexCenterFactor*vd-zr)*invvd
                                          if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized U exceeds max in rooftop',&
                                                quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellu=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=40 ! 3 is the old value
                                          endif
                                       endif
                                    endif
                                    if(vflag .eq. 1)then
                                       if(zr .le. bld%params%vortexHeightFactor*zref2v)then
                                          quwinds%vo(i,j,k)=vel_magv*log(zr*invzo)*denomv
                                          if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in rooftop',&
                                                quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellv=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=41 ! 3 is the old value
                                          endif
                                       endif
                                       if(hdv .lt. 0.5*bld%Rcx(ibuild) .and. zr .le. bld%params%vortexHeightFactor*shell_heightv)then
                                          quwinds%vo(i,j,k)=-quwinds%vo_roof(i,j,k)*(vortexCenterFactor*vd-zr)*invvd
                                          if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in rooftop',&
                                                quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellv=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=40 ! 3 is the old value
                                          endif
                                       elseif(hdv .lt. bld%Rcx(ibuild) .and. zr .le. bld%params%vortexHeightFactor*shell_heightv)then
                                          quwinds%vo(i,j,k)=-quwinds%vo_roof(i,j,k)*(vortexCenterFactor*shell_heightv-zr)*invvd
                                          !quwinds%vo(i,j,k)=-quwinds%vo_roof(i,j,k)*(vortexCenterFactor*vd-zr)*invvd
                                          if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized V exceeds max in rooftop',&
                                                quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellv=1
                                          if(wflag .eq. 1)then
                                             bld%icellflag(i,j,k)=40 ! 3 is the old value
                                          endif
                                       endif
                                    endif
                                 endif
                                 if(k_shellu+k_shellv .lt. 1)exit
                              enddo  lp003
                           endif
                        enddo   lp004      
                     enddo   lp005
                     !$omp end parallel do
                  else !delta wing vortex
                     !$omp parallel do private(i,j,k,uflag,vflag,wflag,x_u,y_u,x_v,y_v,x_w,y_w, &
                     !$omp hxu,hyu,hdu,hxv,hyv,hdv,hxw,hyw,hdwx,hdwy,k_shellu,k_shellv,k_shellw, &
                     !$omp vel_mag,zr)
lp008:               do j=bld%jstart(ibuild),bld%jend(ibuild)
lp007:                  do i=bld%istart(ibuild),bld%iend(ibuild)
                           uflag=0
                           vflag=0
                           wflag=0
                           !check to see if velocity vector is above the building or in a street canyon cell
                           if(bld%icellflag(i,j,bld%kend(ibuild)) .eq. 0)then
                              uflag=1
                              vflag=1
                              wflag=1
                           else
                              if(bld%icellflag(i-1,j,bld%kend(ibuild)) .eq. 0)uflag=1
                              if(bld%icellflag(i,j-1,bld%kend(ibuild)) .eq. 0)vflag=1
                           endif
                           if(bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 50 .or. bld%icellflag(i,j,bld%kend(ibuild)+1) .eq. 30)then ! 6 and 4 are the old values for Streetcanyon and Cavity regions, respectively
                              uflag=0
                              vflag=0
                              wflag=0
                           endif
                           if(uflag+vflag+wflag .gt. 0 .and. bld%icellflag(i,j,bld%kend(ibuild)+1) .gt. 0)then
                              x_u=((real(i)-1)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_u=-((real(i)-1)*dx-xco)*sinRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              x_v=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-1)*dy-yco)*sinRotationAngle
                              y_v=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-1)*dy-yco)*cosRotationAngle
                              x_w=((real(i)-0.5)*dx-xco)*cosRotationAngle+ &
                                    ((real(j)-0.5)*dy-yco)*sinRotationAngle
                              y_w=-((real(i)-0.5)*dx-xco)*sinRotationAngle+	&
                                    ((real(j)-0.5)*dy-yco)*cosRotationAngle
                              hxu=abs(x_u-xfront)
                              hyu=abs(y_u-yfront)
                              hdu=min(hxu,hyu)
                              hxv=abs(x_v-xfront)
                              hyv=abs(y_v-yfront)
                              hdv=min(hxv,hyv)
                              hxw=abs(x_w-xfront)
                              hyw=abs(y_w-yfront)
                              hdwx=hyw*tanRoofangle
                              hdwy=hxw*tanRoofangle
lp006:                        do k=bld%kend(ibuild)+1,k_ref
                                 k_shellu=0
                                 k_shellv=0
                                 k_shellw=0
                                 if(bld%icellflag(i,j,k) .lt. 1)then
                                    exit
                                 else
                                    zr=zm(k)-bld%Ht(ibuild)
                                    vel_mag=sqrt((quwinds%uo_roof(nint(xco/dx),nint(yco/dy),k)**2.)+&
                                       (quwinds%vo_roof(nint(xco/dx),nint(yco/dy),k)**2.))
                                    if(uflag .eq. 1)then
                                       if(hxu .le. min(bld%Rcx(ibuild),2*hyu*tanRoofangle))then
                                          if(zr .le. min(bld%Rcx(ibuild),hyu*tanRoofangle))then
                                             quwinds%uo(i,j,k)=vel_mag*cosXnorm
                                             if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized U exceeds max in rooftop',&
                                                   quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                          if(zr .le. min(bld%Rcx(ibuild),2*hyu*tanRoofangle))then
                                             quwinds%uo(i,j,k)=vel_mag*cosXnormpPi
                                             if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized U exceeds max in rooftop',&
                                                   quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                       endif
                                       if(hyu .le. min(bld%Rcx(ibuild),2*hxu*tanRoofangle))then
                                          if(zr .le. min(bld%Rcx(ibuild),hxu*tanRoofangle))then
                                             quwinds%uo(i,j,k)=vel_mag*cosYnorm
                                             if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized U exceeds max in rooftop',&
                                                   quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                          if(zr .le. min(bld%Rcx(ibuild),2*hxu*tanRoofangle))then
                                             quwinds%uo(i,j,k)=vel_mag*cosYnormpPi
                                             if(abs(quwinds%uo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized U exceeds max in rooftop',&
                                                   quwinds%uo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellu=1
                                          endif
                                       endif
                                    endif
                                    if(vflag .eq. 1)then
                                       if(hxv .le. min(bld%Rcx(ibuild),2*hyv*tanRoofangle))then
                                          if(zr .le. min(bld%Rcx(ibuild),hyv*tanRoofangle))then
                                             quwinds%vo(i,j,k)=vel_mag*sinXnorm
                                             if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized V exceeds max in rooftop',&
                                                   quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                          if(zr .le. min(bld%Rcx(ibuild),2*hyv*tanRoofangle))then
                                             quwinds%vo(i,j,k)=vel_mag*sinXnormpPi
                                             if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized V exceeds max in rooftop',&
                                                   quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                       endif
                                       if(hyv .le. min(bld%Rcx(ibuild),2*hxv*tanRoofangle))then
                                          if(zr .le. min(bld%Rcx(ibuild),hxv*tanRoofangle))then
                                             quwinds%vo(i,j,k)=vel_mag*sinYnorm
                                             if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized V exceeds max in rooftop',&
                                                   quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                          if(zr .le. min(bld%Rcx(ibuild),2*hxv*tanRoofangle))then
                                             quwinds%vo(i,j,k)=vel_mag*sinYnormpPi
                                             if(abs(quwinds%vo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                                write(325,*)'Parameterized V exceeds max in rooftop',&
                                                   quwinds%vo(i,j,k),quwinds%max_velmag,i,j,k
                                             endif
                                             k_shellv=1
                                          endif
                                       endif
                                    endif
                                    if(wflag .eq. 1)then
                                       if(hxw .le. min(bld%Rcx(ibuild),2*hdwx) .and. zr .le. min(bld%Rcx(ibuild),2*hdwx))then
                                          quwinds%wo(i,j,k)=0.1*vel_mag*((hdwx-hxw)/hdwx)*(1-abs((zr-hdwx)/hdwx))
                                          if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized W exceeds max in rooftop',&
                                                quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellw=1
                                          bld%icellflag(i,j,k)=40 ! 3 is the old value
                                       endif
                                       if(hyw .le. min(bld%Rcx(ibuild),2*hdwy) .and. zr .le. min(bld%Rcx(ibuild),2*hdwy))then
                                          quwinds%wo(i,j,k)=0.1*vel_mag*((hdwy-hyw)/hdwy)*(1-abs((zr-hdwy)/hdwy))
                                          if(abs(quwinds%wo(i,j,k)) .gt. quwinds%max_velmag .and. flag%errorWrite .gt. 0)then
                                             write(325,*)'Parameterized W exceeds max in rooftop',&
                                                quwinds%wo(i,j,k),quwinds%max_velmag,i,j,k
                                          endif
                                          k_shellw=1
                                          bld%icellflag(i,j,k)=40 ! 3 is the old value
                                       endif
                                    endif
                                 endif
                                 if(k_shellu+k_shellv+k_shellw .lt. 1)exit
                              enddo  lp006
                           endif
                        enddo   lp007      
                     enddo   lp008
                  endif
            endselect
         endif
         return
      end
