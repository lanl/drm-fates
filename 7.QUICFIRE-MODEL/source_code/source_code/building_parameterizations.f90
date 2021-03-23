	subroutine building_parameterizations_qu
!************************************************************************
! bcsetup - boundary condition setup program for qwic urb      
!    - called by main.f90                 
!    - calls defbuild.f90, upwind.f90, street_intersec.f90  
! The empirical parameterizations are applied in the following order:
!  1. uninterupted boundary layer
!  2. upwind vortex cavity (calls upwind.f90)
!  3. rooftop recirculation
!  4. wakes
! 
!                 
! THIS VERSION OF QWIC-URB CONTAINS MULTIPLE BUILDING CAPABILITY  
! all velocities are cell face quantities          
! bld%icellflag,celltype,Lagrange multipliers (p1,p2) and their    
! coeficients  (e,f,g,m,n,o,p,q,h,p) are a cell centered quantities. 
!                          
! Lfx is the length of the front eddy on a building in xdirection 
! Lfy is the length of the front eddy on a building in ydirection 
! Lr is the length of the rear vortex cavity behind  a building      
! theta is the angle of mean wind speed (meters/sec)        
! bld%number is the number of buildings in the city array    
! xfo is the x coord. left front center of bld
! yfo is the y coord. left front center of bld
! zfo is the z coord of building base 
! xlo is the x coord. lower front center of bld
! xuo is the x coord. upper front center of bld
!                          
! REFERENCE FRAME DEFINITIONS:                  
! below  is defined as  lesser  z value (ie z(k) is below z(k+1)  
! above  is defined as  greater z value (ie z(k+1) is above z(k)  
! right  is defined as  greater x value (ie x(i+1) is right of x(i)  
! left   is defined as  lesser  x value (ie x(i) is right of x(i+1)  
! front  is defined as  greater y value (ie y(j+1) is in front of y(j)  
! behind is defined as  lesser  y value (ie y(j) is behind y(j+1) 
! ERP dec/2000                      
! moving  over various wind angle routine from earlier version of code  
! ERP August/2002                   
! Most recent Modification:
! Changing of coordintate system to reflect true stagard grid   
! ERP Sept/Oct 2002
! erp 12/17 changes fixing north/south flows
! mdw 1/8/2003 changes include fix to jstartb
! erp 7/23/03 modifications to allow for building blocks to be stacked
!  this involves modifications to k loops,zfo and zbo
! erp 11/13/03 fixing Lfx bug
! erp 11/13/03 fixing rooftop vortex bug
! erp 1/29/04 added Petra Kastner-Klein's finalized street canyon technique
!  that she idependantly tested and verified
!  an option is included to use the CPB potential flow formulas for the velocity 
!  field initialization in the street canyon
!  the input file was modified, in line 9 the initialization method is chosen
!    streetcanyonflag=1 ==> original Roeckle approach
!    streetcanyonflag=2 ==> CPB approach
!  potential flow formulas are not applied near the lateral edges of the canyon
!  depth of lateral vortex zone equal to canyon heigth, u is set to zero, v=u0(z)
! erp 2/11/04 added Nilesh's upwind vortex 
! erp 03/09/04 added canyon flag check to remove wake for skimming flow?
! NLB 02/10/04 Added Upwind Vortex Parameterizations for Perpendicular and Varying Incident Wind Angles
! NLB 10/11/04 Added Rooftop Parameterizations for Perpendicular and Varying Incident Wind Angles 
! erp 03/02/05 This subroutine now can writeout the building grid locations in
!     both a binary (QU_celltype.bin) and ASCII format (celltype2.dat). If
!     format_flag = 1, ASCII only. If format_flag=2, binary only and if
!     format_flag = 3, both the ASCII and binaries are written out
!
! ERP 6/8/2006 this version includes rooftop fixes for both off angle and normal
!  angle calculations (Suhas Pols implementation of off angle fixes)
! Cellflag designations 
!  bld%icellflag = 0  building
!  bld%icellflag = 1  fluid cell BL parameterization
!  bld%icellflag = 2  upwind cavity parameterization
!  bld%icellflag = 3  rooftop parameterization
!  bld%icellflag = 4  near wake cavity parameterization
!  bld%icellflag = 5  far wake cavity parameterization
!  bld%icellflag = 6  street canyon parameterization
!  bld%icellflag = 8  vegetation parameterization
!  bld%icellflag = 9  street intersection parameterization
!  bld%icellflag = 10  parking garage parameterization
!
!
!************************************************************************
			
			use canopy_module
			use bld_module
         use landuse_module
         use grid_module
         use flags_module
         use winds_module         

         implicit none
			
         real vegvelfrac,avg_atten,num_atten
			integer :: i,j,k,ibuild
         
         if (flag%domain_initialized .eq. 0)then
! allocate arrays which are needed in bcsetup, PKK 10/01/02
! they are not passed to other subroutines and will be deallocated at the end of bcsetup 
            if(bld%number > 0) then
               allocate(bld%istart(bld%number),bld%iend(bld%number))
               allocate(bld%jstart(bld%number),bld%jend(bld%number))
               allocate(bld%kstart(bld%number),bld%kend(bld%number))
               allocate(bld%Lf(bld%number),bld%Lr(bld%number))
               allocate(bld%Weff(bld%number),bld%Leff(bld%number))
               allocate(bld%Wt(bld%number),bld%Lt(bld%number))
               allocate(bld%Rscale(bld%number),bld%Rcx(bld%number))  !NLB 10/11/04 For Rooftop
               allocate(quwinds%vo_roof(qugrid%nx,qugrid%ny,qugrid%nz),quwinds%uo_roof(qugrid%nx,qugrid%ny,qugrid%nz)) !NLB 10/10/05
            
   lp001:      do ibuild=1,bld%number
                  if(bld%damage(ibuild) .eq. 2)cycle
                  do k=2,qugrid%nz-1
                     bld%kstart(ibuild)=k
                     if(bld%zfo_actual(ibuild) .le. qugrid%zm(k))exit
                  enddo
                  do k=bld%kstart(ibuild),qugrid%nz-1
                     bld%kend(ibuild)=k
                     if(bld%Ht(ibuild) .lt. qugrid%zm(k+1))exit
                  enddo
               enddo   lp001
            endif
         endif    !end 1st time through if
! added AAG 09/20/06  for multiple time steps, quwinds%uo_roof was not atllocated  
        
!end veg 

lp002:   do i=1,bld%number
            bld%Lf(i)=-999.
            bld%Leff(i)=0.
            bld%Weff(i)=0.
            bld%Wt(i)=0.5*bld%Wti(i)
            bld%Lt(i)=0.5*bld%Lti(i)
         enddo    lp002

!erp 1/30/2003
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate buildings
! this now includes a building generation loop that allows
! for multiple buildings
! calculate building spacing s and set wake flags
! erp 3/9/05 Note that defbuild calls pentagon
         if(flag%domain_initialized .eq. 0)then
            print*,'Importing Building Data'
            call defbuild  !erp 1/05/05
            call wallbc
         else
            do k=2,qugrid%nz-1
               do j=1,qugrid%ny-1
                  do i=1,qugrid%nx-1
                     if(bld%icellflag(i,j,k) .gt. 1)bld%icellflag(i,j,k)=1
                  enddo
               enddo
            enddo
            do ibuild=1,bld%number
               if(bld%geometry(ibuild) .eq. 3)call pentagon(ibuild)
            enddo
         endif
         
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Generate Canopies         
         if(canopy%number .gt. 0. .or. (landuse%flag .eq. 1 .and. &
               (landuse%veg_flag .eq. 1 .or. landuse%urb_flag .eq. 1)))then
            print*,'Applying Vegetation Parameterizations'
            call plantinit
         endif
         
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc  
! Calculate effective width and adjust the effective base height
         if(bld%flag%wake .ne. 0 .or. bld%flag%upwind .gt. 0)then
            call effectivewidth
         endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Call upwind - routine to generate the upwind cavity on a building
         if(bld%flag%upwind .gt. 0)then
            print*,'Applying Upwind Rotor Parameterizations'
            do ibuild=1,bld%number
                if(bld%damage(ibuild) .eq. 2)cycle
                if(bld%btype(ibuild) .eq. 2 .or. bld%btype(ibuild) .eq. 5)cycle
                select case(bld%geometry(ibuild))
                  case(1,4,6)
                     if(bld%zfo(ibuild) .eq. 0.)call upwind(ibuild)                     
                end select
            enddo
         endif


!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! wake section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
         if(bld%flag%wake .ne. 0)then
            print*,'Applying Building Wake Parameterizations'
            do ibuild=1,bld%number               
               if(bld%damage(ibuild) .eq. 2)cycle
               if(bld%btype(ibuild) .eq. 2 .or. bld%btype(ibuild) .eq. 5)cycle               
               select case(bld%geometry(ibuild))
                  case(1,4)
                     call rectanglewake(ibuild)
                  case(2,5)
                     call cylinderwake(ibuild)
                  case(6)
                     call polywake(ibuild)
               endselect               
            enddo
         endif
         

         
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! street canyon section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
         if(bld%flag%streetcanyon .ne. 0)then
            print*,'Applying Street Canyon Parameterizations'
            call streetcanyon
         endif
         
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! sidewall section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(bld%flag%sidewall .ne. 0)then
            print*,'Applying Sidewall Rotor Parameterizations'
            do ibuild=1,bld%number
               if(bld%damage(ibuild) .eq. 2)cycle
               if(bld%btype(ibuild) .eq. 1 .or. bld%btype(ibuild) .eq. 3)then
                  select case(bld%geometry(ibuild))
                     case(1,4,6)
                        call sidewall(ibuild)
                  end select
               endif
            enddo               
         endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! rooftop or courtyard section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
         if(bld%flag%roof .gt. 0)then
            quwinds%uo_roof=quwinds%uo
            quwinds%vo_roof=quwinds%vo
            print*,'Applying Rooftop and/or Courtyard Parameterizations'
         endif
         do ibuild=1,bld%number
            if(bld%damage(ibuild) .eq. 2)cycle
            if(bld%btype(ibuild) .eq. 1 .or. bld%btype(ibuild) .eq. 3)then
               select case(bld%geometry(ibuild))
                  case(1,2)
                     call rooftop(ibuild)
                  case(4,5)
                     call courtyard(ibuild)
                  case(6)
                     call rooftop(ibuild)
                     if(bld%numpolygons(ibuild) .gt. 1)then
                        call courtyard(ibuild)
                     endif
               endselect
            endif
         enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! parking garage section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc 
         if(bld%number_garage .gt. 0)then
            print*,'Applying Garage Parameterizations' 
            do ibuild=1,bld%number
               if(bld%damage(ibuild) .eq. 2)cycle
               if(bld%btype(ibuild) .eq. 3)call parking_garage(ibuild)
            enddo
         endif
         
! end End Building parameterization section

!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
! street intersection section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc  
         if(bld%flag%streetcanyon .ne. 0 .and. bld%flag%intersection .eq. 1)then
            print*,'Applying Blended Region Parameterizations'
            if(bld%number.gt.0)then
               call street_intersect
               call poisson ! AG 04/06/2007 blends street intersection winds
            endif
         endif
! MAN 4/12/2007 Moved Poisson boundary conditions to new subroutine "wallbc"
!ANU 01/04/2006vegetation parameterization
         if(canopy%number.gt.0. .or. (landuse%flag .eq. 1 .and. &
               (landuse%veg_flag .eq. 1 .or. landuse%urb_flag .eq. 1)))then
            !$omp parallel do private(i,k,avg_atten,num_atten,vegvelfrac)
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  if(canopy%top(i,j) .gt. 0.)then
                     do k=2,canopy%ktop(i,j)
                        if(canopy%atten(i,j,k) .gt. 0. .and. bld%icellflag(i,j,k) .ne. 8)then
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
                           vegvelfrac=log(canopy%top(i,j)/canopy%zo(i,j))*&
                                 exp(avg_atten*((qugrid%zm(k)/canopy%top(i,j))-1.))/&
                                 log(qugrid%zm(k)/canopy%zo(i,j))
                           if(vegvelfrac .lt. 1. .and. vegvelfrac .ge. 0.)then
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
                           endif
                           if(bld%icellflag(i,j,k) .gt. 0)then
                              bld%icellflag(i,j,k)=8
                           endif
                        endif
                     enddo
                  endif
               enddo
            enddo
            !$omp end parallel do
         endif
            
            

!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!Celltype Coeficient Section
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
! begin defining celltypes using the cellflags based on boundary cells
         !$omp parallel do private(i,j)
         do k=1,qugrid%nz-1
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  if(bld%icellflag(i,j,k) .eq. 0)then ! MAN 7/8/2005 Celltype definition change
! all cells that are buildings have a zero velocity within them
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
         !$omp end parallel do
! MAN 7/8/2005 Celltype definition change



! new print out removed from upwind and now located here because of the
! removal of the upwind cavity under certain conidtions
! erp 12/06/04
         
         !$omp parallel do private(i,j)
         do k=2,qugrid%nz
            do j=1,qugrid%ny
               do i=1,qugrid%nx
                  if((quwinds%uo(i,j,k) .ne. quwinds%uo(i,j,k)) .or. &
                     (quwinds%vo(i,j,k) .ne. quwinds%vo(i,j,k)) .or. &
                     (quwinds%wo(i,j,k) .ne. quwinds%wo(i,j,k)))then
                     ! print*,'NaN found at ',i,j,k,' Cellflag = ',bld%icellflag(i,j,k)
                     if(quwinds%uo(i,j,k) .ne. quwinds%uo(i,j,k))then
                        quwinds%uo(i,j,k)=0.
                     endif
                     if(quwinds%vo(i,j,k) .ne. quwinds%vo(i,j,k))then
                        quwinds%vo(i,j,k)=0.
                     endif
                     if(quwinds%wo(i,j,k) .ne. quwinds%wo(i,j,k))then
                        quwinds%wo(i,j,k)=0.
                     endif
                  endif
               enddo
            enddo
         enddo
         !$omp end parallel do
!end change

         flag%domain_initialized=1
         
         return

                  
		end subroutine building_parameterizations_qu
!======================================================================================
!======================================================================================
		subroutine building_parameterizations_fire
		
      
      use grid_module
      use bld_module
      use flags_module
      use winds_module

      implicit none
         
		integer :: i,j,k
        
		if (flag%domain_initialized .eq. 0)then
			call wallbc
		else
			do k=2,qugrid%nz-1
            do j=1,qugrid%ny-1
               do i=1,qugrid%nx-1
                  if(bld%icellflag(i,j,k) .gt. 1)bld%icellflag(i,j,k)=1
               enddo
            enddo
			enddo
		endif
				
		!$omp parallel do private(i,j)
      do k=1,qugrid%nz-1
         do j=1,qugrid%ny-1
            do i=1,qugrid%nx-1
               if(bld%icellflag(i,j,k) .eq. 0)then ! MAN 7/8/2005 Celltype definition change
! all cells that are buildings have a zero velocity within them
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
      !$omp end parallel do
		flag%domain_initialized=1
		
		end subroutine building_parameterizations_fire
!======================================================================================
!======================================================================================

		subroutine defbuild()
		end subroutine defbuild
		
		subroutine effectivewidth
		end subroutine effectivewidth
		
		subroutine polywake(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine polywake
		
		subroutine cylinderwake(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine cylinderwake
		
		subroutine rectanglewake(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine rectanglewake
		
		subroutine streetcanyon
		end subroutine streetcanyon
		
		subroutine street_intersect
		end subroutine street_intersect
		
		subroutine parking_garage(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine parking_garage
		
		subroutine courtyard(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine courtyard
		
		subroutine rooftop(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine rooftop
		
		subroutine sidewall(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine sidewall
		
		subroutine upwind(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine upwind
		
		subroutine pentagon(ibuild)
		implicit none
		integer,intent(IN) :: ibuild
		end subroutine pentagon
		