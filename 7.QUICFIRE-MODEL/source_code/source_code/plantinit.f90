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
      subroutine plantinit

         !this function initializes the velocity profile in the specified
         !vegetative area area using the MacDonald (2000) approach.
         !ca - attenuation coefficient
         !vk - von Karmen constant

         !TMB 7/18/03 internal plant canopy variables
         !ANU 6/2005 implemented

         
			use canopy_module
			use constants
			use damage_module
         use bld_module
         use grid_module 
         use winds_module
         
         implicit none
   
         real bisect  !bisection function
         real vegvelfrac,avg_atten,num_atten, uh			
			integer :: i,j,k
			
         canopy%ktop(:,:)=0
         canopy%ustar(:,:)=0.
         canopy%zo(:,:)=0.
         canopy%d(:,:)=0.
         
!erp  add subroutine to calculate ustar and zo  using a least squares
! regression of the current initialization
			
			call regress
			
         do j=1,qugrid%ny-1
            do i=1,qugrid%nx-1					
               if(canopy%top(i,j) .gt. 0.)then
                  if(damage%flag .eq. 1)then
                     damage%Rbuild=sqrt((((real(i)-0.5)-damage%explosionx)**2.)+(((real(j)-0.5)-damage%explosiony)**2.))
                     if(damage%Rbuild .le. damage%Rdestroyed)then
                        canopy%top(i,j)=0.
                        do k=2,qugrid%nz
                           canopy%atten(i,j,k)=0.
                        enddo
                        cycle
                     endif
						endif						
                  canopy%d(i,j) = bisect(canopy%ustar(i,j),canopy%zo(i,j), &
                        canopy%top(i,j),canopy%atten(i,j,canopy%ktop(i,j)),0.)
                  if(canopy%d(i,j) .gt. 0.99*canopy%top(i,j))then
                     canopy%d(i,j)=0.7*canopy%top(i,j)
                     canopy%zo(i,j)=min(0.1*canopy%top(i,j),0.5*qugrid%dz_array(1))
						endif
                  uH = (canopy%ustar(i,j)/vk)*log((canopy%top(i,j)-canopy%d(i,j))/canopy%zo(i,j))
                  do k=2,qugrid%nz						
                     if(qugrid%zm(k) .le. canopy%top(i,j))then
                        if(canopy%atten(i,j,k) .gt. 0.)then
                           avg_atten = canopy%atten(i,j,k)
                           if(canopy%atten(i,j,k+1) .ne. canopy%atten(i,j,k) &
                                 .or. canopy%atten(i,j,k-1) .ne. canopy%atten(i,j,k))then
                              num_atten=1.
                              if(canopy%atten(i,j,k+1) .gt. 0.)then
                                 avg_atten = avg_atten + canopy%atten(i,j,k+1)
                                 num_atten=num_atten+1.
                              endif
                              if(canopy%atten(i,j,k-1) .gt. 0.)then
                                 avg_atten = avg_atten + canopy%atten(i,j,k-1)
                                 num_atten=num_atten+1.
                              endif
                              avg_atten=avg_atten/num_atten
                           endif
                           vegvelfrac=log((canopy%top(i,j)-canopy%d(i,j))/canopy%zo(i,j))*&
                                 exp(avg_atten*((qugrid%zm(k)/canopy%top(i,j))-1.))/&
                                 log(qugrid%zm(k)/canopy%zo(i,j))
                           ! if(abs(vegvelfrac) .gt. 1.)then
                           !    print*,i,j,k
                           !    print*,'vegvelfrac',vegvelfrac
                           !    print*,'atten',avg_atten
                           !    print*,'top',canopy%top(i,j)
                           !    print*,'d',canopy%d(i,j)
                           !    print*,'zo',canopy%zo(i,j)
                           !    print*,'z',qugrid%zm(k)
                           ! endif
                           if(vegvelfrac .gt. 1. .or. vegvelfrac .lt. 0.)vegvelfrac=1.
                           quwinds%uo(i,j,k)=quwinds%uo(i,j,k)*vegvelfrac
                           quwinds%vo(i,j,k)=quwinds%vo(i,j,k)*vegvelfrac
                           if(j .lt. qugrid%ny-1)then
                              if(canopy%atten(i,j+1,k) .eq. 0.)then
                                 quwinds%vo(i,j+1,k)=quwinds%vo(i,j+1,k)*vegvelfrac
                              endif
                           endif
                           if(i .lt. qugrid%nx-1)then
                              if(canopy%atten(i+1,j,k) .eq. 0.)then
                                 quwinds%uo(i+1,j,k)=quwinds%uo(i+1,j,k)*vegvelfrac
                              endif
                           endif
                           if(bld%icellflag(i,j,k) .gt. 0)then
                              bld%icellflag(i,j,k)=8
                           endif
                        endif
                     else
                        vegvelfrac=log((qugrid%zm(k)-canopy%d(i,j))/canopy%zo(i,j))/&
                                 log(qugrid%zm(k)/canopy%zo(i,j))
                        if(vegvelfrac .gt. 1. .or. vegvelfrac .lt. 0.)vegvelfrac=1
                        quwinds%uo(i,j,k)=quwinds%uo(i,j,k)*vegvelfrac
                        quwinds%vo(i,j,k)=quwinds%vo(i,j,k)*vegvelfrac
                        if(j .lt. qugrid%ny-1)then
                           if(canopy%atten(i,j+1,canopy%ktop(i,j)) .eq. 0.)then
                              quwinds%vo(i,j+1,k)=quwinds%vo(i,j+1,k)*vegvelfrac
                           endif
                        endif
                        if(i .lt. qugrid%nx-1)then
                           if(canopy%atten(i+1,j,canopy%ktop(i,j)) .eq. 0.)then
                              quwinds%uo(i+1,j,k)=quwinds%uo(i+1,j,k)*vegvelfrac
                           endif
                        endif
                     endif
                  enddo
               endif
            enddo
         enddo
         return
      end
