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
      subroutine parking_garage(ibuild)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! This subroutine set the velocities within the parking garage
! ERP/AAG 8/2007
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
         use bld_module
         use winds_module
         
         implicit none
         
			integer, intent(IN) :: ibuild
			real K_garage
			integer :: i,j,k
			
         K_garage = 0.5	!multiplictive factor to reduce velocity in garage
         do k=bld%kstart(ibuild),bld%kend(ibuild)
            do j=bld%jstart(ibuild),bld%jend(ibuild)
               do i=bld%istart(ibuild),bld%iend(ibuild)
                  if(bld%icellflag(i,j,k) .ne. 0 .and. bld%icellflag(i,j,k+1) .eq. 0)then
                     quwinds%uo(i,j,k)  =quwinds%uo(i,j,k)*K_garage
                     if(bld%icellflag(i+1,j,k+1) .ne. 0)quwinds%uo(i+1,j,k)=quwinds%uo(i+1,j,k)*K_garage
                     quwinds%vo(i,j,k)  =quwinds%vo(i,j,k)*K_garage
                     if(bld%icellflag(i,j+1,k+1) .ne. 0)quwinds%vo(i,j+1,k)=quwinds%vo(i,j+1,k)*K_garage
                     quwinds%wo(i,j,k)  =quwinds%wo(i,j,k)*K_garage
                     quwinds%wo(i,j,k+1)=quwinds%wo(i,j,k+1)*K_garage
                     bld%icellflag(i,j,k) = 10
                  endif
               enddo
            enddo
         enddo
         return
      end
