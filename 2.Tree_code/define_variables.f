!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! define_variables definess all variables, both constant and 
! user-defined, used throughout the program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine define_constant_variables
      !-----------------------------------------------------------------
      ! Constant variables
      !-----------------------------------------------------------------
      use constant_variables
      implicit none

      PI = 3.14159

      end subroutine define_constant_variables

      subroutine define_grid_variables
      !-----------------------------------------------------------------
      ! Grid variables
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      
      implicit none

      integer i,j,k
      real,external:: zcart
           
      nfuel = ngrass+ntspecies*(ilitter+tfuelbins)
      allocate(rhof(nfuel,nx,ny,nz))
      allocate(sizescale(nfuel,nx,ny,nz))
      allocate(moist(nfuel,nx,ny,nz))
      allocate(fueldepth(nfuel,nx,ny,nz))
      
      !-----------------------------------------------------------------
      ! Create topo layer (Should be adjusted for non-flat topo)
      !-----------------------------------------------------------------
      
      allocate(zs(nx,ny))
      allocate(zheight(nx,ny,nz))
      
      if (topofile.eq.'flat'.or.topofile.eq.'') then
        zs(:,:)=0.0
        print *,'Not using topo'      
      else
        open (1,file=topofile,form='unformatted',status='old')
        read (1) zs
        close (1)
        print *,'Reading topo file = ',topofile
      endif
      
      do i=1,nx
        do j=1,ny
          do k=1,nz
            zheight(i,j,k) = zcart((k-1)*dz,zs(i,j))
            if(i.eq.1.and.j.eq.1) print*,'cell',k,'bottom height',zheight(i,j,k)
          enddo
        enddo
      enddo

      end subroutine define_grid_variables

      subroutine define_baseline_variables
      !-----------------------------------------------------------------
      ! Variables unique to the baseline
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables

      implicit none

      integer i,j,itree
      real,dimension(7+3*tfuelbins):: temp_array
      
      if (igrass.ne.0) then
        allocate(grhof(ngrass,nx,ny,nz))
        allocate(gsizescale(ngrass,nx,ny,nz))
        allocate(gmoist(ngrass,nx,ny,nz))
        allocate(gfueldepth(ngrass,nx,ny))
      endif
      if (itrees.ne.0) then
        allocate(trhof(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tsizescale(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tmoist(ntspecies*tfuelbins,nx,ny,nz))
        allocate(tfueldepth(ntspecies*tfuelbins,nx,ny))
      endif
      if (ilitter.ne.0) then
        allocate(lrhof(ntspecies,nx,ny,nz))
        allocate(lsizescale(ntspecies,nx,ny,nz))
        allocate(lmoist(ntspecies,nx,ny,nz))
        allocate(lfueldepth(ntspecies,nx,ny))
      endif

      !-----------------------------------------------------------------
      ! Groundfuel variables unique to the ground fuels baseline
      !-----------------------------------------------------------------
      if (igrass.eq.1) then  ! EJ modified it, before it was ne.0, which did not work for igrass=2
        allocate(grho(ngrass))
        allocate(gmoisture(ngrass))
        allocate(gss(ngrass))
        allocate(gdepth(ngrass))
        
        open (1,file=grassfile)
        read (1,*) grho       ! bulk density of grass [kg/m3]
        read (1,*) gmoisture  ! moisture content of grass
        read (1,*) gss        ! size scale of grass [m]
        read (1,*) gdepth     ! depth of grass [m]
        close (1)
      endif
      
      !-----------------------------------------------------------------
      ! Tree variables unique to the tree baseline
      !-----------------------------------------------------------------
      print*,"itrees:",itrees
      if (itrees.eq.1) then
        allocate(ntrees(ntspecies)) ! Number of trees for each species
        allocate(tcanopy(ntspecies)) ! Canopy closure [fraction]
        allocate(theight(2,ntspecies)) ! Tree heights [m]
        allocate(tcrownbotheight(2,ntspecies)) ! Height to live crown [m]
        allocate(tcrowndiameter(2,ntspecies)) ! Crown diameter [m]
        allocate(tcrownmaxheight(2,ntspecies)) ! Height to max crown diameter [m]
        allocate(tbulkdensity(tfuelbins,ntspecies)) ! Crown fuel bulk density [kg/m3]
        allocate(tmoisture(tfuelbins,ntspecies)) ! Crown fuel moisture content [fraction]
        allocate(tss(tfuelbins,ntspecies)) ! Crown fuel size scale [m]
        open (2,file=treefile)
        read (2,*) tcanopy
        read (2,*) theight(1,:)
        read (2,*) theight(2,:)
        read (2,*) tcrownbotheight(1,:)
        read (2,*) tcrownbotheight(2,:)
        read (2,*) tcrowndiameter(1,:)
        read (2,*) tcrowndiameter(2,:)
        read (2,*) tcrownmaxheight(1,:)
        read (2,*) tcrownmaxheight(2,:)
        do i=1,tfuelbins
          read(2,*) tbulkdensity(i,:)
          read(2,*) tmoisture(i,:)
          read(2,*) tss(i,:)
        enddo
        close (2)
      endif
      
      if (itrees.eq.2.or.itrees.eq.3) then
        itree = 0
        open (2,file=treefile)
        do
          read (2,*,end=10)
          itree = itree+1
        enddo
10      rewind(2)
        
        allocate(ntrees(1)) ! Total number of trees
        allocate(tspecies(itree)) ! Tree cartesian coordinates [m,m]
        allocate(tlocation(itree,2)) ! Tree cartesian coordinates [m,m]
        allocate(theight(itree,1))   ! Tree heights [m]
        allocate(tcrownbotheight(itree,1)) ! Height to live crown [m]
        allocate(tcrowndiameter(itree,1)) ! Crown diameter [m]
        allocate(tcrownmaxheight(itree,1)) ! Height to max crown diameter [m]
        allocate(tbulkdensity(tfuelbins,itree)) ! Crown fuel bulk density [kg/m3]
        allocate(tmoisture(tfuelbins,itree)) ! Crown fuel moisture content [fraction]
        allocate(tss(tfuelbins,itree)) ! Crown fuel size scale [m]
     
        ntrees(1) = itree 
        open (2,file=treefile)
        do i=1,itree
          read (2,*) temp_array(:)
          tspecies(i) = temp_array(1)
          tlocation(i,:) = temp_array(2:3)
          theight(i,1) = temp_array(4)
          tcrownbotheight(i,1) = temp_array(5)
          tcrowndiameter(i,1) = temp_array(6)
          tcrownmaxheight(i,1) = temp_array(7)
          do j=1,tfuelbins
            tbulkdensity(j,i) = temp_array(7+(j-1)*tfuelbins+1)
            tmoisture(j,i)    = temp_array(7+(j-1)*tfuelbins+2)
            tss(j,i)          = temp_array(7+(j-1)*tfuelbins+3)
          enddo
        enddo
        close (2)
      endif

      if (ilitter.ne.0) then ! Elchin changed it to not equal to escape segfault
        allocate(lrho(ntspecies))
        allocate(lmoisture(ntspecies))
        allocate(lss(ntspecies))
        allocate(ldepth(ntspecies))
      
        open (3,file=litterfile)
        read (3,*) lrho       ! bulk density of litter [kg/m3]
        read (3,*) lmoisture  ! moisture content of litter [fraction]
        read (3,*) lss        ! size scale of litter [m]
        read (3,*) ldepth     ! depth of litter [m]
        close (3)
      endif

      end subroutine define_baseline_variables
