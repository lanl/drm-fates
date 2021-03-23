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
	subroutine sor3d
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! SOR3D is a 3 dimensional successive overrelaxation method solver.
	! err is an absolute error measurement
	! iter is the number of iterations spent in the solver
	! p1 & p2 are Lagrange multipliers
	! omegarelax is overelaxation coefficient specified in main.f
	! ERP
	! inlcudes PKKs f90 modifications 10/01/2002
	! includes twh's deallocation procedure 1/8 and 1/9 2003
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	
	use time_module
	use sor_module
	use grid_module
	use flags_module
	
	implicit none
	
	integer :: iter,is,ie,js,je,ks,ke,i,j,k
	real :: abse, eps
				
	iter = 0
	abse = 0
	
	call InitP1()

	if(flag%isfire == 0) then
		
		do iter = 1,sor%itermax
			abse=0.0
		
			sor%p2 = sor%p1		
		
			!$omp parallel do private(i,j)
			do k=2,qugrid%nz-2
				do j=2,qugrid%ny-2
					do i=2,qugrid%nx-2
						sor%p1(i,j,k)=sor%denom(i,j,k)*									&
							((sor%bc%e(i,j,k)*sor%p1(i+1,j,k)+							&
							sor%bc%f(i,j,k)*sor%p1(i-1,j,k))+							&
							(sor%bc%g(i,j,k)*sor%p1(i,j+1,k)+							&
							sor%bc%h(i,j,k)*sor%p1(i,j-1,k))+							&
							(sor%bc%m(i,j,k)*sor%p1(i,j,k+1)+							&
							sor%bc%n(i,j,k)*sor%p1(i,j,k-1)) &
							-sor%r(i,j,k)) + sor%one_minus_omegarelax*sor%p1(i,j,k)
					enddo
				enddo
			enddo
			!$OMP end parallel do
								
			! implementing the 'k' boundary condition, AAG & IS 07/03/06
			! dlambda / dz	= 0 at the interface
			sor%p1(:,:,1) = sor%p1(:,:,2)
			! calculating residual, AAG 07/03/06
			abse=sum(abs(sor%p1-sor%p2))*qugrid%one_over_nxnynz

			! MAN 09/26/2006 added residual reduction check instead of max error
			if(iter .eq. 1)then
				eps = abse*sor%res_red_factor
			endif

			! checking convergence, AAG 07/03/06
			if (abse < eps .or. abse .le. 1.e-9 .or. iter .eq. sor%itermax) exit
		
		enddo
		! MAN 06/04/2007 keeps the correct number of iterations for cases where iter equals itermax
		print*, '# Iterations = ',iter, ';   Average Residual = ',abse
			
	else
	
		is = max(sor%mc_is,2)
		ie = min(sor%mc_ie,qugrid%nx-2)
		js = max(sor%mc_js,2)
		je = min(sor%mc_je,qugrid%ny-2)
		ks = max(sor%mc_ks,2)
		ke = min(sor%mc_ke,qugrid%nz-2)			
		do iter = 1,sor%itermax
			sor%p2 = sor%p1
						
			!$omp parallel do private(i,j,k)
			do k = ks,ke
				do j = js,je
					do i = is,ie
						sor%p1(i,j,k)=sor%denom(i,j,k) * (		&
							sor%bc%e(i,j,k)*sor%p1(i+1,j,k)+		&
							sor%bc%f(i,j,k)*sor%p1(i-1,j,k)+		&
							sor%bc%g(i,j,k)*sor%p1(i,j+1,k)+		&
							sor%bc%h(i,j,k)*sor%p1(i,j-1,k)+		&
							(sor%bc%m(i,j,k)*sor%p1(i,j,k+1)+	&
							sor%bc%n(i,j,k)*sor%p1(i,j,k-1)) * sor%alpha2sq / sor%alpha2_fire(i,j,k)**2 &
							-sor%r(i,j,k)) + sor%one_minus_omegarelax*sor%p1(i,j,k)						
					enddo
				enddo
			enddo
			!$OMP end parallel do
				
			! implementing the 'k' boundary condition, AAG & IS 07/03/06
			sor%p1(:,:,1) = sor%p1(:,:,2)
			! calculating residual, AAG 07/03/06
			abse=sum(abs(sor%p1-sor%p2))*qugrid%one_over_nxnynz
					
			! MAN 09/26/2006 added residual reduction check instead of max error
			if(iter .eq. 1)then
				eps = abse*sor%res_red_factor
			endif

			! checking convergence, AAG 07/03/06
			if (abse < eps .or. abse .le. 1.e-9 .or. iter .eq. sor%itermax) exit
		
		enddo
	endif
	sor%p2 = sor%p1
	
	return
	
	end
!===========================================================================================	
!===========================================================================================
	SUBROUTINE InitP1()
	
	use sor_module
	use grid_module
	use flags_module

	implicit none
	
	if(sor%option == 0)then
		sor%p1 = 0.	
	endif
	
	END
!===========================================================================================	
!===========================================================================================
