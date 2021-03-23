	! Notice
	! This program was prepared by the University of California (University)
	! under Contract W-7405-ENG-36 with the U.S. Department of Energy (DOE).
	! All rights in the program are reserved by DOE on behalf of the Government
	! and the University pursuant to the contract. You are authorized to use
	! this program for Government purposes but it is not to be released or
	! distributed to the public.
	! NEITHER THE UNITED STATES NOR THE UNITED STATES DEPARTMENT OF ENERGY,
	! NOR THE UNIVERSITY OF CALIFORNIA, NOR ANY OF THEIR EMPLOYEES,
	! MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY
	! OR RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS, OF
	! ANY INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
	! THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS.
	!
	subroutine divergence
	!cccccccccccccccccccccccccccccccccccccccccccccccc
	! This subroutine calculates the divergence of
	! the initial velocity field
	! Eric Pardyjak
	! December 2000
	! r(i,j,k)= (-2*alpha1**2)*divergence
	! where r is the right hand side of the functional being
	! minimized.
	!cccccccccccccccccccccccccccccccccccccccccccccccc
	
	use sor_module
	use time_module
	use grid_module
	use flags_module
	use winds_module

	implicit none
	
	real ovdz
	integer is,ie,js,je,ks,ke,i,j,k
		
	if (flag%isfire == 0) then
		!$omp parallel do private(i,j,ovdz)		
		do k=1,qugrid%nz-1
			! MAN 07/25/2008 stretched vertical grid
			ovdz=1./qugrid%dz_array(k)
			do j=1,qugrid%ny-1
				do i=1,qugrid%nx-1
					sor%r(i,j,k)=(-2.*sor%alpha1sq)*(qugrid%dxi*(quwinds%uo(i+1,j,k)-quwinds%uo(i,j,k)) + &
						qugrid%dyi*(quwinds%vo(i,j+1,k)-quwinds%vo(i,j,k)) + &
						ovdz*(quwinds%wo(i,j,k+1)-quwinds%wo(i,j,k)))					
				enddo
			enddo
		enddo
		!$omp end parallel do
		
	else		
		is = max(sor%mc_is,1)
		ie = min(sor%mc_ie,qugrid%nx-1)
		js = max(sor%mc_js,1)
		je = min(sor%mc_je,qugrid%ny-1)
		ks = max(sor%mc_ks,1)
		ke = min(sor%mc_ke,qugrid%nz-1)
		
		!$omp parallel do private(i,j,k,ovdz)
		do k = ks,ke
		! MAN 07/25/2008 stretched vertical grid
			ovdz=1./qugrid%dz_array(k)
			do j = js,je
				do i = is,ie
					sor%r(i,j,k)=(-2.*sor%alpha1sq)*(qugrid%dxi*(quwinds%uo(i+1,j,k)-quwinds%uo(i,j,k)) + &
						qugrid%dyi*(quwinds%vo(i,j+1,k)-quwinds%vo(i,j,k)) + &
						ovdz*(quwinds%wo(i,j,k+1)-quwinds%wo(i,j,k)))
				enddo
			enddo
		enddo
		!$omp end parallel do
	
	endif
		
	return
	
	end

