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
      subroutine regress
!subroutine regress to do a linear regression to determine
!ustar and zocan
!ANU 8/30/2005
	 
         
			use canopy_module
         use constants
         use grid_module
         use winds_module
			
         implicit none

			integer :: i,j,k
         real sumx
         real sumy
         real sumxy
         real sumxsq
         integer no  !, index,c
         integer k_top
         real xm,ym,y
         real ao
         real local_mag
			
			do j=1,qugrid%ny-1
            do i=1,qugrid%nx-1
               if(canopy%top(i,j) .gt. 0.)then
                  do k=2,qugrid%nz-1
                     canopy%ktop(i,j)=k
                     if(canopy%top(i,j) .lt. qugrid%zm(k+1))exit
                  enddo
                  do k=canopy%ktop(i,j),qugrid%nz-1
                     k_top=k
                     if(2.*canopy%top(i,j) .lt. qugrid%zm(k+1))exit
                  enddo
                  if(k_top .eq. canopy%ktop(i,j))k_top=canopy%ktop(i,j)+1
                  if(k_top .gt. qugrid%nz-1)k_top=qugrid%nz-1
                  sumx = 0
                  sumy = 0
                  sumxy = 0
                  sumxsq = 0
                  no = 0
                  do k = canopy%ktop(i,j),k_top
                     no = no + 1
                     local_mag=sqrt((quwinds%uo(i,j,k)**2.)+(quwinds%vo(i,j,k)**2.))
                     ! MAN 07/25/2008 stretched vertical grid
                     y = log(qugrid%zm(k))
                     sumx = sumx + local_mag
                     sumy = sumy + y
                     sumxy = sumxy + local_mag * y
                     sumxsq = sumxsq + local_mag**2
                  enddo
                  canopy%ustar(i,j) = vk * (((no* sumxsq) - (sumx**2))/((no * sumxy) - (sumx * sumy)))
                  xm = sumx / no
                  ym = sumy / no
                  ao = ym - ((vk / canopy%ustar(i,j)) * xm)
                  canopy%zo(i,j) = exp(ao)
               endif
            enddo
         enddo
      end