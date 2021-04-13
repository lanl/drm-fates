!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! io contains the functions for reading input namelists and tree data 
! files and writing .dat files
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine namelist_input
      !-----------------------------------------------------------------
      ! namelist_input is a function reads in the user-defined variables 
      ! from the treelist then assigns those variables for use by the 
      ! trees softeware.
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      use treatment_variables
      implicit none

      namelist/fuellist/
     .   nx,ny,nz,dx,dy,dz,aa1,topofile,
     .   igrass,ngrass,grassconstant,grassfile,
     .   itrees,ntspecies,tfuelbins,tdnx,tdny,treefile,
     .   ilitter,litterconstant,litterfile,
     .   itreatment,sdnx,sdny,
     .   sdiameter,sheight,sprho,smoist,sdepth
      
      ! Area of interest arrays need to be allocated before calling namelist
      allocate(tdnx(2)) ! Array of the cell range (x)  where the trees are applied
      allocate(tdny(2)) ! Array of the cell range (x)  where the trees are applied
      allocate(sdnx(2)) ! Array of the cell range (x)  where the treatment is applied
      allocate(sdny(2)) ! Array of the cell range (x)  where the treatment is applied
       
      open(unit=15,file='Inputs/fuellist',form='formatted',status='old')
           read (15,nml=fuellist)
      close(15)

      ! Corrections for if variables not specifiedi on namelist
      if (tdnx(1).eq.0) then
        tdnx(1) = 1
        tdnx(2) = nx
        tdny(1) = 1
        tdny(2) = ny
      endif

      end subroutine namelist_input

      subroutine output
      !-----------------------------------------------------------------
      ! output is a function which writes the .dat files for use in 
      ! FIRETEC or QUIC-Fire
      !-----------------------------------------------------------------
      use grid_variables
      implicit none
      integer ift,i,j,k,lfuel
      real,dimension(nfuel):: nonzero
      
      nonzero(:) = 1
      lfuel = 1
      do ift=1,nfuel
        if (sum(rhof(ift,:,:,:)).le.0) then
          nonzero(ift) = 0
        else
          do k=1,nz
            if (sum(rhof(ift,:,:,k)).gt.0) lfuel = max(lfuel,k)
          enddo
        endif
      enddo

      print*,'Exporting data to .dat files'
      open (1,file='treesrhof.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) rhof(ift,:,:,1:lfuel)
      enddo
      close (1)

      open (1,file='treesfueldepth.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) fueldepth(ift,:,:,1)
      enddo
      close (1)

      open (1,file='treesss.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) sizescale(ift,:,:,1:lfuel)
      enddo
      close (1)

      open (1,file='treesmoist.dat',form='unformatted',status='unknown')
      do ift=1,nfuel
        if (nonzero(ift).ne.0)  write (1) moist(ift,:,:,1:lfuel)
      enddo
      close (1)

      print*,'Your nfuel is',int(sum(nonzero(:)))
      print*,'Your lfuel is',lfuel
      
      end subroutine output
