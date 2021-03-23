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
      subroutine outfile
!********************************************************
!							
! this subroutine writes out the data files and closes all
! necessary files.
! ASCII files currently written:
!	1. uofield.dat - the initial non-mass consistent velocity field
!	2. uoutmat.dat - the final velocity field (Matlab Format)
!	
! ASCII files that can be uncommented to be written:
!	1. uoutfield.dat - the final velocity field (TECHPLOT format)
!	2. uoutmat.dat - final u-velocity field before cell cent. averaging
!	3. voutmat.dat - final u-velocity field before cell cent. averaging
!	4. woutmat.dat - final u-velocity field before cell cent. averaging
!
! 1/15/04 erp Unformatted binary writes have been added
! ERP
!********************************************************
         
			use file_handling_module
			use canopy_module
			use flags_module
         use bld_module
         use grid_module
         use flags_module
         use winds_module
         use time_module
         
         implicit none
         
			real x,y
			integer :: i,j,k,ibuild

 912     format(i9,  '       !Time Increment')
 913     format(i9,  '       !Number of time steps')
 914     format(i9,  '       !Unix Time')
 916     format(3i6,1x,3(f17.5,1x))

         if(flag%output_initialized .eq. 0)then
            if(flag%uofield .eq. 1)then !erp 3/2/05
               open(unit=ID_FILE_QU_UOFIELD,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
						"uofield.dat",status="unknown")
               write(ID_FILE_QU_UOFIELD,*)'% Inital velocity field i,j,k,uo,vo,wo'
               write(ID_FILE_QU_UOFIELD,913)qutime%nsteps
            endif
            if(flag%frm .eq. 1 .or. flag%frm .eq. 3)then  
               open(unit=ID_FILE_QU_VELOCITY,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
						"QU_velocity.dat",status="unknown")
               write(ID_FILE_QU_VELOCITY,*)'%matlab velocity output file'
               write(ID_FILE_QU_VELOCITY,913)qutime%nsteps
            endif
            if(flag%frm .eq. 1 .or. flag%frm .eq. 3)	& !erp 3/2/05
					open(unit=ID_FILE_QU_CELLTYPE_BIN,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_celltype.dat",status="unknown")
            if(flag%frm .eq. 2 .or. flag%frm .eq. 3) & !erp 3/2/05
					open(unit=ID_FILE_QU_VELOCITY_BIN,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_velocity.bin",form='unformatted',status="unknown")
            if(flag%frm .eq. 2 .or. flag%frm .eq. 3) & !erp 3/2/05   
               open(unit=ID_FILE_QU_CELLTYPE_BIN,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_celltype.bin",form='unformatted',status="unknown")
            open(unit=ID_FILE_QU_UNDIST_BIN,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_undisturbed.bin",form='unformatted',status="unknown")
            if(flag%staggered .eq. 1) &
               open(unit=ID_FILE_QU_STAGGERED,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_staggered_velocity.bin",form='unformatted',status="unknown")
            if(canopy%flag .gt. 0) &
               open(unit=ID_FILE_QU_VEG_BIN,file=TRIM(workingDirectory)//TRIM(fileSeparator)// &
					"QU_vegetation.bin",form='unformatted',status="unknown")
            call sufacecoords
            
            write(ID_FILE_QU_BUILDOUT,*)bld%number,' ! total number of buildings'
            write(ID_FILE_QU_BUILDOUT,*)canopy%number + bld%numberneg, ' ! total number of vegitative canopies'
            
            flag%output_initialized=1
         endif
         !MAN 07/28/2008 
         if(canopy%flag .gt. 0)then
            write(ID_FILE_QU_VEG_BIN)((canopy%ustar(i,j),i=1,qugrid%nx-1),j=1,qugrid%ny-1)
            write(ID_FILE_QU_VEG_BIN)(((canopy%atten(i,j,k),i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
         endif
         if(flag%uofield.eq.1) then !erp 3/02/05
            write(ID_FILE_QU_UOFIELD,*)'% Begin Output for new time step'
            write(ID_FILE_QU_UOFIELD,912)qutime%current_iter
            write(ID_FILE_QU_UOFIELD,914)qutime%unix(qutime%current_iter)
            do k=1,qugrid%nz
               do j=1,qugrid%ny
                  do i=1,qugrid%nx
                     write(ID_FILE_QU_UOFIELD,916)i,j,k,quwinds%uo(i,j,k),quwinds%vo(i,j,k),quwinds%wo(i,j,k)
                  enddo
               enddo
            enddo
         endif
	
         if(flag%frm.eq.1 .or. flag%frm.eq.3)  then
            write(ID_FILE_QU_VELOCITY,*)'% Begin Output for new time step'
            write(ID_FILE_QU_VELOCITY,912)qutime%current_iter
            write(ID_FILE_QU_VELOCITY,914)qutime%unix(qutime%current_iter)
         endif
         if(flag%staggered .eq. 1)then
            write(ID_FILE_QU_STAGGERED)                                              &
               (((quwinds%u(i,j,k),i=2,qugrid%nx),j=1,qugrid%ny-1),k=1,qugrid%nz-1), &
               (((quwinds%v(i,j,k),i=1,qugrid%nx-1),j=2,qugrid%ny),k=1,qugrid%nz-1), &
               (((quwinds%w(i,j,k),i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=2,qugrid%nz)
         endif
!erp modified to appropriately reflect physical locations
!erp 10/8/2003
!erp 1/18/2005 lines added to write out unformatted for binary read into Matlab
         if(flag%frm .eq. 1 .or. flag%frm .eq. 3)then
		      do k=1,qugrid%nz-1
               do j=1,qugrid%ny-1
                  do i=1,qugrid%nx-1
                     x=(0.5*(real(i+1)+real(i))-1.)*qugrid%dx
                     y=(0.5*(real(j+1)+real(j))-1.)*qugrid%dy
                     write(ID_FILE_QU_VELOCITY,101)x,y,qugrid%zm(k),	   &
								0.5*(quwinds%u(i,j,k)+quwinds%u(i+1,j,k)),		&
								0.5*(quwinds%v(i,j,k)+quwinds%v(i,j+1,k)),		&
								0.5*(quwinds%w(i,j,k)+quwinds%w(i,j,k+1))
                  enddo
               enddo
            enddo
         endif
			if(flag%frm .eq. 2 .or. flag%frm .eq. 3)then
            write(ID_FILE_QU_VELOCITY_BIN)										&
               (((0.5*(quwinds%u(i,j,k)+quwinds%u(i+1,j,k)),            &
                  i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1),	&
               (((0.5*(quwinds%v(i,j,k)+quwinds%v(i,j+1,k)),            &
                  i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1),	&
               (((0.5*(quwinds%w(i,j,k)+quwinds%w(i,j,k+1)),            &
                  i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
			endif
 101     format(6(f11.5,1x))
         
         if(flag%frm.eq.1 .or. flag%frm.eq.3)then
            do k=1,qugrid%nz-1
               do j=1,qugrid%ny-1
                  do i=1,qugrid%nx-1
                     write(ID_FILE_QU_CELLTYPE_BIN,71)		&
								(real(i)-.5)*qugrid%dx,					&
								(real(j)-.5)*qugrid%dy,					&
								qugrid%zm(k),bld%icellflag(i,j,k)
                  enddo
               enddo
            enddo
         endif
!erp 3/02/2005 lines added to write out unformatted for binary read into Matlab
         if(flag%frm .eq. 2 .or. flag%frm .eq. 3)then  !erp 3/2/05
            write(ID_FILE_QU_CELLTYPE_BIN)(((bld%icellflag(i,j,k),i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
         endif

 71      format(3(1x,f8.3),i5)
         
         
         write(ID_FILE_QU_UNDIST_BIN)(((quwinds%undisturbed(i,j,k),i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)

         do ibuild=1,bld%number
            if(bld%btype(ibuild) .ne. 2)then
! write to buildout.dat erp 1/15/2004
               write(ID_FILE_QU_BUILDOUT,*)bld%num(ibuild),' !number'
               write(ID_FILE_QU_BUILDOUT,*)bld%damage(ibuild),' !damage'
               if(bld%geometry(ibuild) .eq. 3)then
                  write(ID_FILE_QU_BUILDOUT,*)bld%Wti(ibuild),' !Weff'
                  write(ID_FILE_QU_BUILDOUT,*)bld%Lti(ibuild),' !Leff'
               else
                  write(ID_FILE_QU_BUILDOUT,*)bld%Weff(ibuild),' !Weff'
                  write(ID_FILE_QU_BUILDOUT,*)bld%Leff(ibuild),' !Leff'
               endif
               write(ID_FILE_QU_BUILDOUT,*)bld%Lf(ibuild),' !Lf'
               write(ID_FILE_QU_BUILDOUT,*)bld%Lr(ibuild),' !Lr'
            endif
         enddo
         return
      end
!===================================================================================
!===================================================================================
	subroutine output_vert_grid_qu()
	
	use grid_module
	use file_handling_module
	
	implicit none
	
	integer :: k
	
	open (unit=ID_FILE_Z_QU, file=TRIM(workingDirectory)//TRIM(fileSeparator)//'z_qu.bin',form='unformatted')
	write(ID_FILE_Z_QU)(qugrid%z(k),k=1,qugrid%nz)
	write(ID_FILE_Z_QU)(qugrid%zm(k),k=1,qugrid%nz)
	write(ID_FILE_Z_QU) firegrid%nz_en2atmos
	close(ID_FILE_Z_QU)
	
	end
!===================================================================================
!===================================================================================
	subroutine output_wind_qu_inst(simtime)

	use constants
	use file_handling_module
   use grid_module
   use winds_module

	implicit none
	
	integer,intent(IN) :: simtime	
	character(len=35) :: filename		
	integer i,j,k
	
	write(filename, "(A8,I0.5,A4)") "qu_windu", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_U, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_U) (((0.5*(quwinds%u(i,j,k)+quwinds%u(i+1,j,k)), &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
   close(ID_FILE_WIND_QU_U)

	write(filename, "(A8,I0.5,A4)") "qu_windv", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_V, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_V)(((0.5*(quwinds%v(i,j,k)+quwinds%v(i,j+1,k)),  &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
   close(ID_FILE_WIND_QU_V)
 
	write(filename, "(A8,I0.5,A4)") "qu_windw", int(simtime), ".bin"
   open (unit=ID_FILE_WIND_QU_W, file=TRIM(workingDirectory)//TRIM(fileSeparator)//filename, &
		form='unformatted',err=1101)
   write (ID_FILE_WIND_QU_W) (((0.5*(quwinds%w(i,j,k)+quwinds%w(i,j,k+1)),  &
      i=1,qugrid%nx-1),j=1,qugrid%ny-1),k=1,qugrid%nz-1)
	close(ID_FILE_WIND_QU_W)

	return

1101 write(msgoutfile,*)'Error while opening file '//TRIM(filename)//'.'
	 call TerminateProgram()

	end
!===================================================================================
!===================================================================================
