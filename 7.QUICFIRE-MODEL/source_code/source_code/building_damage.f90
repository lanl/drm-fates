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
      subroutine building_damage
         
			use file_handling_module
			use bld_module
			use damage_module
			
         implicit none
         
         real xco,yco
         integer :: numberOfSources,j,ibuild
         
         
         open(unit=ID_FILE_QP_SOURCE,file=TRIM(workingDirectory)// &
				TRIM(fileSeparator)//'QP_source.inp',status='old')
         read(ID_FILE_QP_SOURCE,*) ! QUIC version
         read(ID_FILE_QP_SOURCE,*) numberOfSources ! number of sources
         if(numberOfSources .eq. 0)then
            damage%flag=0
            return
         endif
         read(ID_FILE_QP_SOURCE,*) ! number of source nodes
         read(ID_FILE_QP_SOURCE,*) ! start source comment line
         read(ID_FILE_QP_SOURCE,*) ! source name
         read(ID_FILE_QP_SOURCE,*) ! source strength unit flag
         read(ID_FILE_QP_SOURCE,*) ! source strength value
         read(ID_FILE_QP_SOURCE,*) ! source release type
         read(ID_FILE_QP_SOURCE,*) ! source start time
         read(ID_FILE_QP_SOURCE,*) ! source geometry
         read(ID_FILE_QP_SOURCE,*)damage%explosionx
         read(ID_FILE_QP_SOURCE,*)damage%explosiony
         read(ID_FILE_QP_SOURCE,*) ! explosion z value
         read(ID_FILE_QP_SOURCE,*)damage%hemass
         
         close(ID_FILE_QP_SOURCE)
         
         damage%Rdestroyed=255*((damage%hemass/907185)**(1./3.))
         damage%Rdamaged=1425*((damage%hemass/907185)**(1./3.))
         
         do ibuild=1,bld%number
            if(bld%btype(ibuild) .eq. 2)cycle
            if(bld%geometry(ibuild) .eq. 3)then
               xco = bld%xfo(ibuild)
               yco = bld%yfo(ibuild)
            else
               if(bld%gamma(ibuild) .ne. 0.)then
                  xco = bld%xfo(ibuild) + 0.5*bld%Lti(ibuild)*cos(bld%gamma(ibuild))
                  yco = bld%yfo(ibuild) + 0.5*bld%Lti(ibuild)*sin(bld%gamma(ibuild))
               else
                  xco = bld%xfo(ibuild) + 0.5*bld%Lti(ibuild)
                  yco = bld%yfo(ibuild)
               endif
            endif
            damage%Rbuild=sqrt(((xco-damage%explosionx)**2.)+((yco-damage%explosiony)**2.))
            if(damage%Rbuild .le. damage%Rdestroyed)then
               bld%damage(ibuild)=2
            elseif(damage%Rbuild .le. damage%Rdamaged)then
               bld%damage(ibuild)=1
            endif
         enddo
         do ibuild=1,bld%number
            if(bld%btype(ibuild) .eq. 2)cycle
            if(bld%damage(ibuild) .eq. 1 .and. bld%zfo_actual(ibuild) .gt. 0.)then
               do j=1,bld%number
                  if(bld%group_id(ibuild) .eq. bld%group_id(j) .and. ibuild .ne. j &
                        .and. bld%zfo_actual(j) .lt. bld%zfo_actual(ibuild) .and. bld%damage(j) .eq. 2 &
                        .and. bld%btype(j) .ne. 2)then
                     bld%damage(ibuild)=2
                  endif
               enddo
            endif
         enddo
         return
      end