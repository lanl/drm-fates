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
      subroutine wallbc

         
         use sor_module
         use bld_module
         use grid_module

         implicit none
			
			integer :: i,j,k
			
         ! non boundary cells
!erp 6/12/2006 the following includes the redefinition of b.c. based on 
!C. Bathke work
         sor%bc%e(:,:,:)=1.
         sor%bc%e(qugrid%nx-1,:,:)=0.
         sor%bc%f(:,:,:)=1.
         sor%bc%f(1,:,:)=0.
         sor%bc%g(:,:,:)=1.
         sor%bc%g(:,qugrid%ny-1,:)=0.
         sor%bc%h(:,:,:)=1.
         sor%bc%h(:,1,:)=0.
         sor%bc%m(:,:,:)=1.
         sor%bc%m(:,:,qugrid%nz-1)=0.
         sor%bc%n(:,:,:)=1.
         sor%bc%n(:,:,1)=0.
			
         !$omp parallel do private(i,j)
         do k=2,qugrid%nz-2
            do j=2,qugrid%ny-2
               do i=2,qugrid%nx-2
                  if(bld%icellflag(i,j,k) .ne. 0)then
                     if(bld%icellflag(i,j,k-1) +bld%icellflag(i,j,k+1)+bld%icellflag(i-1,j,k) &
                           +bld%icellflag(i+1,j,k)+bld%icellflag(i,j-1,k)+bld%icellflag(i,j+1,k) .lt. 1)then
                        bld%icellflag(i,j,k)=0
                     else
                        if(bld%icellflag(i,j,k-1) .eq. 0)sor%bc%n(i,j,k)=0.
                        if(bld%icellflag(i,j,k+1) .eq. 0)sor%bc%m(i,j,k)=0.
                        if(bld%icellflag(i-1,j,k) .eq. 0)sor%bc%f(i,j,k)=0.
                        if(bld%icellflag(i+1,j,k) .eq. 0)sor%bc%e(i,j,k)=0.
                        if(bld%icellflag(i,j-1,k) .eq. 0)sor%bc%h(i,j,k)=0.
                        if(bld%icellflag(i,j+1,k) .eq. 0)sor%bc%g(i,j,k)=0.
                     endif
                  endif
               enddo     
            enddo      
			enddo
         !$omp end parallel do
			
			sor%bc%e = sor%bc%e / qugrid%dx**2
         sor%bc%f = sor%bc%f / qugrid%dx**2
         sor%bc%g = sor%bc%g / qugrid%dy**2
         sor%bc%h = sor%bc%h / qugrid%dy**2
			
         !$omp parallel do
         do k=2,qugrid%nz-1
             sor%bc%m(:,:,k) = sor%bc%m(:,:,k)/ &
               (qugrid%dz_array(k)*0.5*(qugrid%dz_array(k)+qugrid%dz_array(k+1))) * sor%eta
             sor%bc%n(:,:,k) = sor%bc%n(:,:,k)/ &
               (qugrid%dz_array(k)*0.5*(qugrid%dz_array(k)+qugrid%dz_array(k-1))) * sor%eta
			enddo
			!$omp end parallel do
			sor%denom = sor%omegarelax/      &
				(sor%bc%e + sor%bc%f + sor%bc%g + &
				sor%bc%h + sor%bc%m + sor%bc%n)
			sor%denom0 = sor%denom
			
         return
      end
