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
       subroutine turbulence_model
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subroutine to turbulent viscosity (use smogorinsky model) 
! subroutine called by diffusion.f90      
! Akshay A. Gowardhan                  
! August 2006                    
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
          
			 use diffusion_module
          use bld_module
          use grid_module
          use winds_module
          
          implicit none
          
			 real  dzi, shear, cs_les, delta, mol_visc
			 integer :: i,j,k
			 
          cs_les = 0.2
          mol_visc = 1.8e-6
          delta = (qugrid%dx*qugrid%dy*qugrid%dz)**(1./3.)         
          diff%visc = 0.
!  Turbulence model (smagorinsky)
          shear=0.0
          !$omp parallel do private(i,j,k,shear,dzi)
          do k=2,qugrid%nz-1
             ! MAN 07/25/2008 stretched vertical grid
             dzi=1./(qugrid%dz_array(k))
             do j=2,qugrid%ny-1
                do i=2,qugrid%nx-1
                   shear= 2.0*( &
                         (qugrid%dxi*(quwinds%u(i,j,k)-quwinds%u(i-1,j,k)))**2+ &
                         (qugrid%dyi*(quwinds%v(i,j,k)-quwinds%v(i,j-1,k)))**2+ &
                         (dzi*(quwinds%w(i,j,k)-quwinds%w(i,j,k-1)))**2)
                   shear= shear +0.25*( &
                         ((qugrid%dyi*(quwinds%u(i,j+1,k)-quwinds%u(i,j,k)))+          &
                         (qugrid%dxi*(quwinds%v(i+1,j,k)-quwinds%v(i,j,k))))**2+			&
                         ((qugrid%dyi*(quwinds%u(i,j,k)-quwinds%u(i,j-1,k)))+          &
                         (qugrid%dxi*(quwinds%v(i+1,j-1,k)-quwinds%v(i,j-1,k))))**2+	&
                         ((qugrid%dyi*(quwinds%u(i-1,j+1,k)-quwinds%u(i-1,j,k)))+      &
                         (qugrid%dxi*(quwinds%v(i,j,k)-quwinds%v(i-1,j,k))))**2+	      &
                         ((qugrid%dyi*(quwinds%u(i-1,j,k)-quwinds%u(i-1,j-1,k)))+      &
                         (qugrid%dxi*(quwinds%v(i,j-1,k)-quwinds%v(i-1,j-1,k))))**2)
                   shear= shear +0.25*( &
                         ((dzi*(quwinds%u(i,j,k+1)-quwinds%u(i,j,k)))+                 &
                         (qugrid%dxi*(quwinds%w(i+1,j,k)-quwinds%w(i,j,k))))**2+		   &
                         ((dzi*(quwinds%u(i,j,k)-quwinds%u(i,j,k-1)))+                 &
                         (qugrid%dxi*(quwinds%w(i+1,j,k-1)-quwinds%w(i,j,k-1))))**2+   &
                         ((dzi*(quwinds%u(i-1,j,k+1)-quwinds%u(i-1,j,k)))+             &
                         (qugrid%dxi*(quwinds%w(i,j,k)-quwinds%w(i-1,j,k))))**2+	      &
                         ((dzi*(quwinds%u(i-1,j,k)-quwinds%u(i-1,j,k-1)))+             &
                         (qugrid%dxi*(quwinds%w(i,j,k-1)-quwinds%w(i-1,j,k-1))))**2)
                   shear= shear +0.25*( &
                         ((dzi*(quwinds%v(i,j,k+1)-quwinds%v(i,j,k)))+                 &
                         (qugrid%dyi*(quwinds%w(i,j+1,k)-quwinds%w(i,j,k))))**2+		   &
                         ((dzi*(quwinds%v(i,j,k)-quwinds%v(i,j,k-1)))+                 &
                         (qugrid%dyi*(quwinds%w(i,j+1,k-1)-quwinds%w(i,j,k-1))))**2+   &
                         ((dzi*(quwinds%v(i,j-1,k+1)-quwinds%v(i,j-1,k)))+             &
                         (qugrid%dyi*(quwinds%w(i,j,k)-quwinds%w(i,j-1,k))))**2+	      &
                         ((dzi*(quwinds%v(i,j-1,k)-quwinds%v(i,j-1,k-1)))+             &
                         (qugrid%dyi*(quwinds%w(i,j,k-1)-quwinds%w(i,j-1,k-1))))**2)
                   diff%visc(i,j,k)= (cs_les*delta)**2* sqrt(abs(shear)) + mol_visc
                enddo
             enddo
          enddo
          !$omp end parallel do
          diff%visc(1,:,:)=diff%visc(qugrid%nx-1,:,:)
          diff%visc(qugrid%nx,:,:)=diff%visc(2,:,:)
          diff%visc(:,:,1)=0.0
          diff%visc(:,:,qugrid%nz)=diff%visc(:,:,2)
          diff%visc(:,1,:)=diff%visc(:,qugrid%ny-1,:)
          diff%visc(:,qugrid%ny,:)=diff%visc(:,2,:)
          !$omp parallel do private(i,j)
          do k=2,qugrid%nz-2
             do j=2,qugrid%ny-2
                do i=2,qugrid%nx-2
                   if(bld%icellflag(i+1,j,k)==0 .or. bld%icellflag(i-1,j,k)==0 .or. &
                     bld%icellflag(i,j+1,k)==0 .or. bld%icellflag(i,j-1,k)==0       &  
                         .or. bld%icellflag(i,j,k+1)==0 .or. bld%icellflag(i,j,k-1)==0 )diff%visc(i,j,k)=0.0
                enddo
             enddo
          enddo
          !$omp end parallel do
          return
       end