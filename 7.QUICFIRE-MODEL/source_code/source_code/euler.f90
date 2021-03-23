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
   subroutine euler
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! subroutine to calculate mass consistent velocity field
! subroutine called by main.f90
! subroutine calls outfile.f90 to print out results
! alpha1 = horizontal gaussian precision moduli
! alpha2 = vertical gaussian precision moduli
! Eric Pardyjak
! December 2000
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      
		use sor_module
		use bld_module
		use grid_module
		use flags_module
		use winds_module
		
      implicit none
		
      real ovalph2, ovdz
		integer :: is,ie,js,je,ks,ke,i,j,k
			
		if (flag%isfire == 0) then			
! make sure outside cells still have uo value
			
			
			quwinds%u = quwinds%uo
			quwinds%v = quwinds%vo
			quwinds%w = quwinds%wo		
			
! using staggard grid arrangement calculate new mass consistent
! velocities. Note that we go from 2,nx-1 because we are using
! the p2's to the left and right. This is kind of a bruut force
! method.
			!$omp parallel do private(i,j,ovdz)
			do k=2,qugrid%nz-1
				! MAN 07/25/2008 stretched vertical grid				
				do j=2,qugrid%ny-1
					do i=2,qugrid%nx-1
						quwinds%u(i,j,k)=quwinds%uo(i,j,k)+sor%ovalph1*qugrid%dxi*(sor%p2(i,j,k)-sor%p2(i-1,j,k))
						quwinds%v(i,j,k)=quwinds%vo(i,j,k)+sor%ovalph1*qugrid%dyi*(sor%p2(i,j,k)-sor%p2(i,j-1,k))
						quwinds%w(i,j,k)=quwinds%wo(i,j,k)+sor%ovalph2*qugrid%dzmi(k)*(sor%p2(i,j,k)-sor%p2(i,j,k-1))
					enddo
				enddo
			enddo
			!$omp end parallel do
! all cells that are buildings have a zero velocity within them
			!$omp parallel do private(i,j,k)
			do k=1,qugrid%nz-1
				do j=1,qugrid%ny-1
					do i=1,qugrid%nx-1
						if(bld%icellflag(i,j,k).eq.0.)then ! MAN 7/8/2005 celltype definition change
							quwinds%u(i,j,k)  =0.
							quwinds%u(i+1,j,k)=0.
							quwinds%v(i,j,k)  =0.
							quwinds%v(i,j+1,k)=0.
							quwinds%w(i,j,k)  =0.
							quwinds%w(i,j,k+1)=0.
						endif
					enddo
				enddo
			enddo
			!$omp end parallel do
				 
		else
			
			is = max(sor%mc_is,1)
			ie = min(sor%mc_ie,qugrid%nx)
			js = max(sor%mc_js,1)
			je = min(sor%mc_je,qugrid%ny)
			ks = max(sor%mc_ks,1)
			ke = min(sor%mc_ke,qugrid%nz)
			
			! make sure outside cells still have uo value
			!$omp parallel do private(i,j,k)
			do k = ks,ke
				do j = js,je
					do i = is,ie
						quwinds%u(i,j,k)=quwinds%uo(i,j,k)
						quwinds%v(i,j,k)=quwinds%vo(i,j,k)
						quwinds%w(i,j,k)=quwinds%wo(i,j,k)
					enddo
				enddo
			enddo
			!$omp end parallel do
! using staggard grid arrangement calculate new mass consistent
! velocities. Note that we go from 2,nx-1 because we are using
! the p2's to the left and right. This is kind of a bruut force
! method.
			is = max(sor%mc_is,2)
			ie = min(sor%mc_ie,qugrid%nx-1)
			js = max(sor%mc_js,2)
			je = min(sor%mc_je,qugrid%ny-1)
			ks = max(sor%mc_ks,2)
			ke = min(sor%mc_ke,qugrid%nz-1)
			
			!$omp parallel do private(i,j,k,ovdz,ovalph2)
			do k = ks,ke
				! MAN 07/25/2008 stretched vertical grid
				ovdz=1./(0.5*(qugrid%dz_array(k)+qugrid%dz_array(k-1)))
				do j = js,je
					do i = is,ie
						ovalph2 = 0.5 / sor%alpha2_fire(i,j,k)**2
						quwinds%u(i,j,k)=quwinds%uo(i,j,k)+sor%ovalph1*qugrid%dxi*(sor%p2(i,j,k)-sor%p2(i-1,j,k))
						quwinds%v(i,j,k)=quwinds%vo(i,j,k)+sor%ovalph1*qugrid%dyi*(sor%p2(i,j,k)-sor%p2(i,j-1,k))
						quwinds%w(i,j,k)=quwinds%wo(i,j,k)+ovalph2*ovdz*(sor%p2(i,j,k)-sor%p2(i,j,k-1))
					enddo
				enddo
			enddo
			!$omp end parallel do
! all cells that are buildings have a zero velocity within them
			
			
		endif
		
      return
   end
