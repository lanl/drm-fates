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
       subroutine diffusion
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subroutine to calculate diffusive flux and update the velocity field with this flux  
! subroutine called by main.f90
! subroutine calls outfile.f90 to print out results   
! subroutine calls turbulence_model.f90 to calculate the turbulent viscosity        
! Akshay A. Gowardhan                  
! August 2006                    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
          use diffusion_module
          use grid_module
          use winds_module
          
          implicit none
          
			 integer ip,im,jp,jm,kp,km,i,j,k        
          real Tuuip, Tuuim, Tuvjp,Tuvjm, Tuwkp,Tuwkm
          real Tvuip,Tvuim, Tvvjp, Tvvjm, Tvwkp, Tvwkm
          real Twuip, Twuim, Twvjp, Twvjm, Twwkp, Twwkm, mzp, mz, mzm
			 real dt
                  
          quwinds%uo = quwinds%u
          quwinds%vo = quwinds%v
          quwinds%wo = quwinds%w
          
          call turbulence_model
          do k=2,qugrid%nz-1
               ! MAN 07/25/2008 stretched vertical grid
               mzp= qugrid%dz_array(k+1)/2.
		         mz = qugrid%dz_array(k)/2.
		         mzm= qugrid%dz_array(k-1)/2.
             do j=2,qugrid%ny-1
                do i=2,qugrid%nx-1
                   ip = i + 1
                   jp = j + 1
                   kp = k + 1
                   im = i - 1
                   jm = j - 1
                   km = k - 1
! X momentum      
                   Tuuip = diff%visc(i,j,k)*(2.*( qugrid%dxi * (quwinds%uo(ip,j,k)-quwinds%uo(i ,j,k) )))
                   Tuuim = diff%visc(i,j,k)*(2.*( qugrid%dxi * (quwinds%uo(i ,j,k)-quwinds%uo(im,j,k) )))
                   Tuvjp = diff%visc(i,j,k)*(qugrid%dyi*(quwinds%uo(i,jp,k)-quwinds%uo(i, j,k))+ &
                     qugrid%dxi*(quwinds%vo(ip,j, k)-quwinds%vo(i ,j ,k)))
                   Tuvjm = diff%visc(i,j,k)*(qugrid%dyi*(quwinds%uo(i,j ,k)-quwinds%uo(i,jm,k))+ &
                     qugrid%dxi*(quwinds%vo(ip,jm,k)-quwinds%vo(i ,jm,k)))
                   Tuwkp = diff%visc(i,j,k)*((1./(mzp+mz))*(quwinds%uo(i,j,kp)-quwinds%uo(i ,j,k))+ &
                     qugrid%dxi*(quwinds%wo(ip,j, k)-quwinds%wo(i ,j ,k)))
                   Tuwkm = diff%visc(i,j,k)*((1./(mzm+mz))*(quwinds%uo(i,j ,k)-quwinds%uo(i,j,km))+ &
                     qugrid%dxi*(quwinds%wo(ip,j,km)-quwinds%wo(i ,j,km)))
                   diff%Fxd(i,j,k) = qugrid%dxi*(Tuuip  - Tuuim )+ qugrid%dyi*( Tuvjp - Tuvjm )+ (0.5/mz)*( Tuwkp - Tuwkm )
! Y momentum
                   Tvuip = diff%visc(i,j,k)*(qugrid%dxi*(quwinds%vo(ip,j,k)-quwinds%vo(i ,j,k))+ &
                     qugrid%dyi*(quwinds%uo(i,jp ,k)-quwinds%uo(i ,j,k)))
                   Tvuim = diff%visc(i,j,k)*(qugrid%dxi*(quwinds%vo(i ,j,k)-quwinds%vo(im,j,k))+ &
                     qugrid%dyi*(quwinds%uo(im,jp,k)-quwinds%uo(im,j,k)))
                   Tvvjp = diff%visc(i,j,k)*(2.*qugrid%dyi*(quwinds%vo(i,jp,k)-quwinds%vo(i,j ,k)))
                   Tvvjm = diff%visc(i,j,k)*(2.*qugrid%dyi*(quwinds%vo(i,j ,k)-quwinds%vo(i,jm,k)))
                   Tvwkp = diff%visc(i,j,k)*((1./(mzp+mz))*(quwinds%vo(i,j,kp)-quwinds%vo(i, j,k))+ &
                     qugrid%dyi*(quwinds%wo(i,jp,k )-quwinds%wo(i ,j,k)))
                   Tvwkm = diff%visc(i,j,k)*((1./(mzm+mz))*(quwinds%vo(i,j,k )-quwinds%vo(i,j,km))+ &
                     qugrid%dyi*(quwinds%wo(i,jp,km)-quwinds%wo(i,j,km)))
                   diff%Fyd(i,j,k)= qugrid%dxi*( Tvuip - Tvuim)+  qugrid%dyi*(  Tvvjp - Tvvjm)+ (0.5/mz)*( Tvwkp - Tvwkm)
! Z momentum
                   Twuip = diff%visc(i,j,k)*(qugrid%dxi*(quwinds%wo(ip,j,k)-quwinds%wo(i ,j,k))+ &
                     (1./(mzp+mz))*(quwinds%uo(i ,j,kp)-quwinds%uo(i ,j,k)))
                   Twuim = diff%visc(i,j,k)*(qugrid%dxi*(quwinds%wo(i ,j,k)-quwinds%wo(im,j,k))+ &
                     (1./(mzp+mz))*(quwinds%uo(im,j,kp)-quwinds%uo(im,j,k)))
                   Twvjp = diff%visc(i,j,k)*(qugrid%dyi*(quwinds%wo(i,jp,k)-quwinds%wo(i, j,k))+ &
                     (1./(mzp+mz))*(quwinds%vo(i ,j,kp)-quwinds%vo(i ,j,k)))
                   Twvjm = diff%visc(i,j,k)*(qugrid%dyi*(quwinds%wo(i,j ,k)-quwinds%wo(i,jm,k))+ &
                   (1./(mzp+mz))*(quwinds%vo(i,jm,kp)-quwinds%vo(i,jm,k)))
                   Twwkp = diff%visc(i,j,k)*(2*(0.5/mzp)*(quwinds%wo(i,j,kp)-quwinds%wo(i,j ,k)))
                   Twwkm = diff%visc(i,j,k)*(2*(0.5/mzm)*(quwinds%wo(i,j,k )-quwinds%wo(i,j,km)))
                   diff%Fzd(i,j,k) = qugrid%dxi*( Twuip  - Twuim)+ qugrid%dyi*(Twvjp  - Twvjm)+ &
                     (1./(mzp+mz))*( Twwkp - Twwkm)
                   enddo
                enddo
             enddo
             dt = 0.25*(min(qugrid%dx,qugrid%dy,qugrid%dz))**2/maxval(diff%visc)
!Update velocity with diffusive fluxes  
             do k=2,qugrid%nz-1
                do j=2,qugrid%ny-1
                   do i=2,qugrid%nx-1
                      quwinds%uo(i,j,k)= quwinds%uo(i,j,k)+ (dt*(diff%Fxd(i,j,k)))
                      quwinds%vo(i,j,k)= quwinds%vo(i,j,k)+ (dt*(diff%Fyd(i,j,k)))
                      quwinds%wo(i,j,k)= quwinds%wo(i,j,k)+ (dt*(diff%Fzd(i,j,k)))
                   enddo
                enddo
             enddo
             return
          end
