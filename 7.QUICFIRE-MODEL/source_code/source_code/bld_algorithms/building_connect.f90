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
      subroutine building_connect(ibuild)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subroutine building_connect - this subroutine 
!	- called by bcsetup.f90
!	- calls none
! ERP 8/17/05
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
			use bld_module
			
         implicit none

			integer,intent(IN) :: ibuild
         integer jj,grpid_chk
         real Pc
         real debug1, debug2, debug3

!	Coeffeicent
         Pc = 1

         grpid_chk = bld%group_id(ibuild)
	
         do jj=1,bld%number
            if(jj.ne.ibuild)then
               if(bld%group_id(jj).eq.grpid_chk)then
                  if(Ht(jj).eq.bld%zfo_actual(ibuild)) then
!MAN 8/30/2005 stacked building fix
                     debug1 = bld%Weff(jj)
                     debug2 = bld%Weff(ibuild)
                     debug3 = Ht(jj)
                     if(bld%Weff(jj) > 0.0)then
                        bld%zfo(ibuild) = Ht(jj)*(1 - (bld%Weff(ibuild)/bld%Weff(jj))**Pc)
                     else
                        bld%zfo(ibuild) = 0.0
                     endif
                     bld%zfo(ibuild) = max(bld%zfo(ibuild),0.)
                  endif
               endif
            endif
         enddo

         return
      end
