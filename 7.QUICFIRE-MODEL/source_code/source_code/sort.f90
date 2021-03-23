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
      subroutine sort
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Subroutine sort - this subroutine sorts the buildings and reorders
! them in terms of height from smallest to tallest so that the 
! empirical algorithms in bcsetup.f90 apply them in that order
!	- called by main.f90
!	- calls none
! ERP 3/8/05
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

         
			use bld_module

         implicit none

			integer :: i,j
         integer, allocatable :: imax(:),bldnum_orig(:),bldtype_orig(:),group_id_orig(:)
         !MAN 08/10/2010 polygon buildings
         integer, allocatable :: bldgeometry_orig(:),bldstartidx_orig(:),bldstopidx_orig(:)
         integer, allocatable :: numpolygons_orig(:),bldroof_orig(:)
         real, allocatable :: bldwall_orig(:),bldcx_orig(:),bldcy_orig(:)
         !end MAN 08/10/2010
         real, allocatable :: Ht_orig(:),Wti_orig(:),Lti_orig(:),Htmax(:),aa_orig(:),bb_orig(:)
         real, allocatable :: xfo_orig(:),yfo_orig(:),zfo_orig(:),gamma_orig(:),atten_orig(:)
         allocate(Ht_orig(bld%number),Wti_orig(bld%number),Lti_orig(bld%number))
         allocate(bldnum_orig(bld%number),bldtype_orig(bld%number))
         allocate(bldgeometry_orig(bld%number),bldstartidx_orig(bld%number))
         allocate(bldstopidx_orig(bld%number),numpolygons_orig(bld%number))
         allocate(bldwall_orig(bld%number),bldroof_orig(bld%number))
         allocate(bldcx_orig(bld%number),bldcy_orig(bld%number))
         allocate(xfo_orig(bld%number),yfo_orig(bld%number),zfo_orig(bld%number))
         allocate(gamma_orig(bld%number),Htmax(bld%number),aa_orig(bld%number),bb_orig(bld%number))
         allocate(imax(bld%number),group_id_orig(bld%number),atten_orig(bld%number))

!define temporary arrays
         bldnum_orig=bld%num
         group_id_orig=bld%group_id ! Add group ID
         bldtype_orig = bld%btype
         bldgeometry_orig = bld%geometry
         bldstopidx_orig = bld%stopidx
         bldstartidx_orig = bld%startidx
         numpolygons_orig = bld%numpolygons
         bldwall_orig = bld%wall
         bldroof_orig = bld%roof
         bldcx_orig = bld%cx
         bldcy_orig = bld%cy
         Ht_orig=bld%Ht
         Wti_orig=bld%Wti
         Lti_orig=bld%Lti
         xfo_orig=bld%xfo
         yfo_orig=bld%yfo
         zfo_orig=bld%zfo
         gamma_orig=bld%gamma
         atten_orig=bld%atten
         aa_orig=bld%aa
         bb_orig=bld%bb

! initialize arrays
         do i=1,bld%number
            Htmax(i) = 0.
            imax(i) = i

!erp 1/6/2006 remove zero height buildings that were from vegetation
            if(bld%btype(i).eq.9)bld%Ht(i)=9999		!a very large number
            if(bld%btype(i).eq.0)bld%Ht(i)=9998
         enddo

         do j=1,bld%number
            do i=1,bld%number
               if(bld%Ht(i) .gt. Htmax(j))then
                  Htmax(j)=bld%Ht(i)
                  imax(j)=i
               endif
            enddo
            bld%Ht(imax(j))=-999
         enddo

         do i=1,bld%number
            !write(48,*)bldnum_orig(imax(bld%number+1-i))
            bld%invnum(imax(bld%number+1-i))=i
            bld%num(i)=bldnum_orig(imax(bld%number+1-i))
            bld%btype(i)=bldtype_orig(imax(bld%number+1-i))
            !MAN 08/10/2010 polygon buildings
            bld%geometry(i)=bldgeometry_orig(imax(bld%number+1-i))
            bld%startidx(i)=bldstartidx_orig(imax(bld%number+1-i))
            bld%stopidx(i)=bldstopidx_orig(imax(bld%number+1-i))
            bld%numpolygons(i)=numpolygons_orig(imax(bld%number+1-i))
            bld%wall(i)=bldwall_orig(imax(bld%number+1-i))
            bld%roof(i)=bldroof_orig(imax(bld%number+1-i))
            bld%cx(i)=bldcx_orig(imax(bld%number+1-i))
            bld%cy(i)=bldcy_orig(imax(bld%number+1-i))
            !end MAN 08/10/2010
            bld%group_id(i)=group_id_orig(imax(bld%number+1-i))
            bld%Ht(i)=Ht_orig(imax(bld%number+1-i))
            bld%Wti(i)=Wti_orig(imax(bld%number+1-i))
            bld%Lti(i)=Lti_orig(imax(bld%number+1-i))
            bld%xfo(i)=xfo_orig(imax(bld%number+1-i))
            bld%yfo(i)=yfo_orig(imax(bld%number+1-i))
            bld%zfo(i)=zfo_orig(imax(bld%number+1-i))
            bld%gamma(i)=gamma_orig(imax(bld%number+1-i))
            bld%atten(i)=atten_orig(imax(bld%number+1-i))
            bld%aa(i)=aa_orig(imax(bld%number+1-i))
            bld%bb(i)=bb_orig(imax(bld%number+1-i))
         enddo
         bld%zfo_actual=bld%zfo
         deallocate(Ht_orig,Wti_orig,Lti_orig,bldnum_orig,bldtype_orig,group_id_orig)
         deallocate(xfo_orig,yfo_orig,zfo_orig,gamma_orig,Htmax,aa_orig,bb_orig)
         deallocate(bldgeometry_orig,bldstartidx_orig,bldstopidx_orig,numpolygons_orig) !MAN 08/10/2010 polygon buildings
	      deallocate(bldwall_orig,bldroof_orig,bldcx_orig,bldcy_orig) !MAN 08/10/2010 polygon buildings
         return
      end


