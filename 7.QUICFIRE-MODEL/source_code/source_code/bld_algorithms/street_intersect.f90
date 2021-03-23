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
!
!
! Subroutine to determine the location of a street intersection
! erp 1/2006
! NOTE: Need to make modification to handle multiple runs 1/24/2006

      subroutine street_intersect
         
         use bld_module
         use grid_module
         use winds_module
         
         implicit none

			integer :: i,j,k
         integer changeflag,intersect_flag,istart_intflag,jstart_intflag,iBuffer,iBufferMax !,EW_flag
         integer, allocatable :: intersect(:,:,:),intersect_1(:,:,:),intersect_2(:,:,:),intersect_1opp(:,:,:),intersect_2opp(:,:,:) 
         integer, allocatable :: E_W_flag(:,:,:),W_E_flag(:,:,:),N_S_flag(:,:,:),S_N_flag(:,:,:)  !SUP
         
         allocate(intersect(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1), &
            intersect_1(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),intersect_2(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
         allocate(intersect_1opp(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),intersect_2opp(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1))
         allocate(E_W_flag(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),W_E_flag(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),N_S_flag(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1),S_N_flag(qugrid%nx-1,qugrid%ny-1,qugrid%nz-1)) !SUP
         
         intersect_flag=0
         iBufferMax=4
! make sure that all cells that are buildings have a zero velocity within them
         do k=1,qugrid%nz-1
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  if(bld%icellflag(i,j,k) .eq. 0)then ! MAN 7/8/2005 celltype definition change
                     quwinds%uo(i,j,k)  =0.
                     quwinds%uo(i+1,j,k)=0.
                     quwinds%vo(i,j,k)  =0.
                     quwinds%vo(i,j+1,k)=0.
                     quwinds%wo(i,j,k)  =0.
                     quwinds%wo(i,j,k+1)=0.
                  endif
               enddo
            enddo
         enddo
         intersect(1:qugrid%nx-1,1:qugrid%ny-1,1:qugrid%nz-1) = 0
         intersect_1(1:qugrid%nx-1,1:qugrid%ny-1,1:qugrid%nz-1) = 0 !x-direction flag
         intersect_2(1:qugrid%nx-1,1:qugrid%ny-1,1:qugrid%nz-1) = 0 !y-direction flag
         intersect_1opp(1:qugrid%nx-1,1:qugrid%ny-1,1:qugrid%nz-1)=0
         intersect_2opp(1:qugrid%nx-1,1:qugrid%ny-1,1:qugrid%nz-1)=0
         E_W_flag=0
         W_E_flag=0
         N_S_flag=0
         S_N_flag=0
         changeflag = 0
! sweep through (x) to find intersections 
         do k=2,qugrid%nz-1
            do j=1,qugrid%ny-1
!SUP sweep through +x
               do i=2,qugrid%nx-1
!determine where the street interesection begins
                  if((bld%icellflag(i-1,j,k) .eq. 50 .or. bld%icellflag(i-1,j,k) .eq. 51) .and. &
                        bld%icellflag(i,j,k) .ne. 50 .and. bld%icellflag(i,j,k) .ne. 51 .and. bld%icellflag(i,j,k) .ne. 0)then !6 is the old value
                     changeflag=1
                     istart_intflag = i
                     iBuffer=0
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1)then
                     select case(bld%icellflag(i,j,k))
                        case(0) !run into another building 
                           changeflag=0
                        case(50) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(51) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(1) !run into free atm.
                           if(iBuffer .gt. iBufferMax)then
                              changeflag=0
                              intersect_1(istart_intflag:i,j,k) = 0
                           else
                              iBuffer=iBuffer+1
                           endif
                     end select
                     if(i .eq. qugrid%nx-1)intersect_1(istart_intflag:qugrid%nx-1,j,k) = 0
                  endif
                  intersect_1(i,j,k) = changeflag
               enddo
               changeflag = 0 !reset flag
!SUP sweep through -x
               do i=qugrid%nx-2,1,-1
!determine where the street interesection begins
                  if((bld%icellflag(i+1,j,k) .eq. 50 .or. bld%icellflag(i+1,j,k) .eq. 51) .and. &
                        bld%icellflag(i,j,k) .ne. 50 .and. bld%icellflag(i,j,k) .ne. 51 .and. bld%icellflag(i,j,k) .ne. 0)then ! 6 is the old value
                     changeflag=1
                     istart_intflag = i
                     iBuffer=0
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1)then
                     select case(bld%icellflag(i,j,k))
                        case(0) !run into another building 
                           changeflag=0
                        case(50) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(51) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(1) !run into free atm.
                           if(iBuffer .gt. iBufferMax)then
                              changeflag=0
                              intersect_1opp(i:istart_intflag,j,k) = 0
                           else
                              iBuffer=iBuffer+1
                           endif
                     end select
                     if(i .eq. 1)intersect_1opp(1:istart_intflag,j,k) = 0
                  endif
                  intersect_1opp(i,j,k) = changeflag
               enddo
               changeflag = 0 !reset flag
            enddo
         enddo
! now sweep in the j direction
         changeflag = 0
         do k=2,qugrid%nz-1
            do i=1,qugrid%nx-1
               do j=2,qugrid%ny-1
                  if((bld%icellflag(i,j-1,k) .eq. 50 .or. bld%icellflag(i,j-1,k) .eq. 51) .and. &
                        bld%icellflag(i,j,k) .ne. 50 .and. bld%icellflag(i,j,k) .ne. 51 .and. bld%icellflag(i,j,k) .ne. 0)then ! 6 is the old value
                     changeflag=1
                     jstart_intflag = j
                     iBuffer=0
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1)then
                     select case(bld%icellflag(i,j,k))
                        case(0) !run into another building 
                           changeflag=0
                        case(50) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(51) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(1) !run into free atm.
                           if(iBuffer .gt. iBufferMax)then
                              changeflag=0
                              intersect_2(i,jstart_intflag:j,k) = 0
                           else
                              iBuffer=iBuffer+1
                           endif
                     end select
                     if(j .eq. qugrid%ny-1)intersect_2(i,jstart_intflag:qugrid%ny-1,k) = 0
                  endif
                  intersect_2(i,j,k) = changeflag
               enddo
               changeflag = 0
!SUP sweep through -y
               do j=qugrid%ny-2,1,-1
                  if((bld%icellflag(i,j+1,k) .eq. 50 .or. bld%icellflag(i,j+1,k) .eq. 51) .and. & 
                        bld%icellflag(i,j,k) .ne. 50 .and. bld%icellflag(i,j,k) .ne. 51 .and. bld%icellflag(i,j,k) .ne. 0)then ! 6 is the old value
                     changeflag=1
                     jstart_intflag = j
                     iBuffer=0
                  endif
!determine where the street intersection ends
                  if(changeflag .eq. 1)then
                     select case(bld%icellflag(i,j,k))
                        case(0) !run into another building 
                           changeflag=0
                        case(50) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(51) !run into another street canyon ! 6 is the old value
                           changeflag=0
                        case(1) !run into free atm.
                           if(iBuffer .gt. iBufferMax)then
                              changeflag=0
                              intersect_2opp(i,j:jstart_intflag,k) = 0
                           else
                              iBuffer=iBuffer+1
                           endif
                     end select
                     if(j .eq. 1)intersect_2opp(i,1:jstart_intflag,k) = 0
                  endif
                  intersect_2opp(i,j,k) = changeflag
               enddo 
               changeflag = 0
            enddo
         enddo
         
         do k=2,qugrid%nz-1
            do j=2,qugrid%ny-2
               do i=2,qugrid%nx-2
                   if(intersect_1(i,j,k) + intersect_1opp(i,j,k) + intersect_2(i,j,k) + &
                         intersect_2opp(i,j,k) .gt. 1 .and. bld%icellflag(i,j,k) .gt. 0)intersect(i,j,k)=1
               enddo
            enddo
         enddo
!
         do k=2,qugrid%nz-1
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  if(intersect(i,j,k) .gt. 0)bld%icellflag(i,j,k)=52 ! 9 is the old value
               enddo
            enddo
         enddo
         deallocate(intersect_1,intersect_2,intersect)
         deallocate(E_W_flag,W_E_flag,N_S_flag,S_N_flag)
         return
      end
