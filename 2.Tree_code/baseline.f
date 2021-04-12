!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! baseline contains the functions which construct the basic fuel map
! based off the forest and ground fuel dimensions defined in
! define_variables
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine baseline
      !-----------------------------------------------------------------
      ! baseline is a function which calls the grass and tree baselines
      ! and consolidates them to fill the rhof, sizescale, moisture, and
      ! fueldepth arrays.
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables

      implicit none

      integer i,j,k,ift
      real,external:: zcart

      call define_baseline_variables 
      
      ! Fill grass arrays
      if (igrass.ne.0) then     
        print*,'Filling Grass Baseline'
        call grass_baseline
        do ift=1,ngrass
          do i=1,nx
            do j=1,ny
              do k=1,nz
                rhof(ift,i,j,k)      = grhof(ift,i,j,k)
                sizescale(ift,i,j,k) = gsizescale(ift,i,j,k)
                moist(ift,i,j,k)     = gmoist(ift,i,j,k)
              enddo
              fueldepth(ift,i,j,1) = gfueldepth(ift,i,j)
            enddo
          enddo
        enddo 
      else
        grhof(:,:,:,:)      = 0
        gmoist(:,:,:,:)     = 0
        gsizescale(:,:,:,:) = 0
        gfueldepth(:,:,:)   = 0
      endif
      
      ! Fill tree arrays
      if (itrees.ne.0) then
        print*,'Filling Trees Baseline'
        call tree_baseline
        do ift=1,ntspecies*tfuelbins
          do i=1,nx
            do j=1,ny
              do k=1,nz
                rhof(ift+ngrass,i,j,k)      = trhof(ift,i,j,k)
                sizescale(ift+ngrass,i,j,k) = tsizescale(ift,i,j,k)
                moist(ift+ngrass,i,j,k)     = tmoist(ift,i,j,k)
              enddo
              fueldepth(ift+ngrass,i,j,1) = tfueldepth(ift,i,j)
            enddo
          enddo
        enddo
      else
        trhof(:,:,:,:)      = 0
        tmoist(:,:,:,:)     = 0
        tsizescale(:,:,:,:) = 0
        tfueldepth(:,:,:)   = 0
      endif

      ! Fill litter arrays
      if (ilitter.ne.0) then
        print*,'Filling Litter Baseline'
        call litter_baseline
        do ift=1,ntspecies
          if (sum(trhof(ift,:,:,1)).lt.sum(trhof(ift,:,:,:))*0.01
     +      .and.ldepth(ift).lt.zcart(1*dz,0)) then
            print*,'Little to no fuel from tree type',ift,'in first cell combining with litter'
            do i=1,nx
              do j=1,ny
                rhof(ngrass+(ift-1)*tfuelbins+1,i,j,1)      = lrhof(ift,i,j,1)
                sizescale(ngrass+(ift-1)*tfuelbins+1,i,j,1) = lsizescale(ift,i,j,1)
                moist(ngrass+(ift-1)*tfuelbins+1,i,j,1)     = lmoist(ift,i,j,1)
                fueldepth(ngrass+(ift-1)*tfuelbins+1,i,j,1) = lfueldepth(ift,i,j)
              enddo
            enddo
          else
            do i=1,nx
              do j=1,ny
                do k=1,nz
                  rhof(ift+ngrass+ntspecies*tfuelbins,i,j,k)      = lrhof(ift,i,j,k)
                  sizescale(ift+ngrass+ntspecies*tfuelbins,i,j,k) = lsizescale(ift,i,j,k)
                  moist(ift+ngrass+ntspecies*tfuelbins,i,j,k)     = lmoist(ift,i,j,k)
                enddo
                fueldepth(ift+ngrass+ntspecies*tfuelbins,i,j,1) = lfueldepth(ift,i,j)
              enddo
            enddo
          endif
        enddo
      else
        lrhof(:,:,:,:)      = 0
        lmoist(:,:,:,:)     = 0
        lsizescale(:,:,:,:) = 0
        lfueldepth(:,:,:)   = 0
      endif
      
      end subroutine baseline
      
      subroutine grass_baseline
      !-----------------------------------------------------------------
      ! grass_baseline is a function which computes the characteristics
      ! of a base grass field and fills rhof, sizescale, moisture, and
      ! fueldepth arrays.
      !-----------------------------------------------------------------
      use grid_variables
      use baseline_variables
      
      implicit none

      integer i,j,k,ift
      real target_mass,actual_mass
      real,allocatable:: rhofxy(:,:)

      if (igrass.eq.2) then
         print*,'Reading LLM WGlitter'
         allocate(rhofxy(nx,ny))
         open(1, file="LLM_litter_WG.txt")
           read(1,*) rhofxy
         close(1)
      endif

      do ift=1,ngrass
        do i=1,nx
          do j=1,ny
            ! updates grho values with the values from LLM WG litter file
            if(igrass.eq.2)grho(ift)=rhofxy(i,j)
            gfueldepth(ift,i,j) = gdepth(ift)
            do k=1,nz-1
              gmoist(ift,i,j,k) = gmoisture(ift)
              gsizescale(ift,i,j,k) = gss(ift)
              if (zheight(i,j,k+1).lt.gdepth(ift)) then
                grhof(ift,i,j,k) = grho(ift)
              else
                grhof(ift,i,j,k) = grho(ift)*(gdepth(ift)-zheight(i,j,k))/(zheight(i,j,k+1)-zheight(i,j,k))
                if (k.gt.zmax) zmax=k
                exit
              endif
            enddo
          enddo
        enddo
      enddo

      target_mass = 0
      actual_mass = 0
      do ift=1,ngrass
        target_mass = target_mass+gdepth(ift)*nx*dx*ny*dy*grho(ift)
        do i=1,nx
          do j=1,ny
            do k=1,zmax
              actual_mass = actual_mass+grhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Grass target fuel mass:',target_mass
      print*,'Grass actual fuel mass:',actual_mass
      print*,'Grass error:',actual_mass/target_mass*100,'%'
                         
      end subroutine grass_baseline

      subroutine tree_baseline 
      !-----------------------------------------------------------------
      ! tree_baseline is a function which computes the characteristics  
      ! of a forest from the variables designated in 
      ! define_variables.f. It then fills trhof, tsizescale, tmoisture, 
      ! and tfueldepth arrays.
      !-----------------------------------------------------------------
      use constant_variables
      use grid_variables
      use baseline_variables
      
      implicit none
      
      integer ift,i,j,k,ii,jj,kk,iii,jjj,kkk
      integer ii_real,jj_real
      integer ift_index
      real totarea
      real xtest,ytest
      integer xtop,xbot,ybot,ytop,zbot,ztop
      real canopytop,canopybot,canopydiameter,canopymaxh
      real atop,abot,test_height,bot_height,top_height
      integer cellcount
      real rhoftemp
      real target_mass,actual_mass
      real,external:: paraboloid,normal 

      !-----Determine the number of trees for each species
      if(itrees.eq.1) then
        totarea   = nx*dx*ny*dy
        do i=1,ntspecies
          ntrees(i) = ceiling(totarea*tcanopy(i)/(PI*(tcrowndiameter(1,i)/2.)**2.))
        enddo
        print*,'Number of trees of each species:',ntrees
      endif

      !----- Begin loop which fills arrays with information for each tree
      do i=1,ntspecies
        if (itrees.eq.1) print*,'Species',i,'with',ntrees(i),'trees'
        do j=1,ntrees(i)
          if (MOD(j,1000).eq.0) print*,'Placing tree',j,'of',ntrees(i)
          
          !----- Place tree location
          if (itrees.eq.1.or.itrees.eq.3) then
            ! Randomly place a tree
            call random_number(xtest)
            xtest = xtest*nx*dx
            call random_number(ytest)
            ytest = tdny(1)+ytest*(tdny(2)-tdny(1))*dy
          else
            ! Specific tree placement
            xtest = tlocation(j,1)
            ytest = tlocation(j,2)
          endif

          !----- Determine tree shape characteristics
          if (itrees.eq.1) then
            ! Sample shape from distributions
            canopytop = min(nz*dz,normal(theight(1,i),theight(2,i)))
            canopybot = max(0.,min(canopytop-0.02,normal(tcrownbotheight(1,i),tcrownbotheight(2,i))))
            canopydiameter = max(0.,normal(tcrowndiameter(1,i),tcrowndiameter(2,i)))
            canopymaxh = min(canopytop-0.01,max(canopybot+0.01,normal(tcrownmaxheight(1,i),tcrownmaxheight(2,i))))
          else
            ! Shape from tree file
            canopytop = theight(j,1)
            canopybot = tcrownbotheight(j,1)
            canopydiameter = tcrowndiameter(j,1)
            canopymaxh= tcrownmaxheight(j,1)
          endif

          !----- Translate tree shape to grid
          xbot = floor((xtest-canopydiameter/2.)/dx+1)
          xtop = floor((xtest+canopydiameter/2.)/dx+1)
          ybot = floor((ytest-canopydiameter/2.)/dy+1)
          ytop = floor((ytest+canopydiameter/2.)/dy+1)
          zbot = 0.
          ztop = 0.
          do k=1,nz-1
            if (canopybot.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),k+1)) then
              zbot = k
              exit
            endif
          enddo
          do kk=k,nz-1
            if (canopytop.le.zheight(nint(xtest/dx+1),nint(ytest/dy+1),kk+1)) then
              ztop = kk
              if (kk.gt.zmax) zmax=kk
              exit
            endif
          enddo
          
          ! Ellipitical paraboloid used to determine tree contours
          atop = 0.25*(canopydiameter)**2./(canopymaxh-canopytop)
          abot = 0.25*(canopydiameter)**2./(canopymaxh-canopybot)
          do ii=xbot,xtop
            if (ii.gt.nx) then
              ii_real = ii-nx 
            else if (ii.lt.1) then 
              ii_real = nx+ii
            else
              ii_real = ii
            endif
            do jj=ybot,ytop
              if (jj.gt.ny) then
                jj_real = jj-ny
              else if (jj.lt.1) then
                jj_real = ny+jj
              else
                jj_real = jj
              endif
              do kk=zbot,ztop
                ! Determine how many of subcells of a cell are within the paraboloid, the fraction of the subcells is equal to the fraction of the cell within the paraboloid
                cellcount = 0
                do iii=1,10
                  do jjj=1,10
                    do kkk=1,10
                      test_height= zheight(ii_real,jj_real,kk)+(2.*kkk-1.)/20.
     +                  *(zheight(ii_real,jj_real,kk+1)-zheight(ii_real,jj_real,kk))
                      bot_height = paraboloid(abot,((ii-1)+(2.*iii-1.)/20.)*dx,xtest,((jj-1)+(2.*jjj-1.)/20.)*dy,ytest,canopybot)
                      top_height = paraboloid(atop,((ii-1)+(2.*iii-1.)/20.)*dx,xtest,((jj-1)+(2.*jjj-1.)/20.)*dy,ytest,canopytop)
                      if (test_height.ge.bot_height.and.test_height.le.top_height) cellcount=cellcount+1
                    enddo
                  enddo
                enddo

                !----- Fill in the 3D arrays for a tree
                if (cellcount.gt.0) then
                  do ift=1,tfuelbins
                    if (itrees.eq.1) then
                      ift_index = (i-1)*tfuelbins+ift
                    else
                      ift_index = (tspecies(j)-1)*tfuelbins+ift
                    endif
                    rhoftemp = tbulkdensity(ift,i)*cellcount/1000.
                    tsizescale(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)*
     +                tsizescale(ift_index,ii_real,jj_real,kk)+rhoftemp*tss(ift,i))/(trhof(ift_index,ii_real,jj_real,kk)+rhoftemp)
                    tmoist(ift_index,ii_real,jj_real,kk) = (trhof(ift_index,ii_real,jj_real,kk)*
     +                tmoist(ift_index,ii_real,jj_real,kk)+rhoftemp*tmoisture(ift,i))/(trhof(ift_index,ii_real,jj_real,kk)+rhoftemp)
                    trhof(ift_index,ii_real,jj_real,kk) = trhof(ift_index,ii_real,jj_real,kk)+rhoftemp
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
        if (itrees.ne.1) exit
      enddo

      !----- Tree fuel depth is equal to height of the first cell
      do i=1,ntspecies
        do j=1,tfuelbins
          do ii=tdnx(1),tdnx(2)
            do jj=tdny(1),tdny(2)
              tfueldepth((i-1)*tfuelbins+j,ii,jj) = zheight(ii,jj,2)
            enddo
          enddo
        enddo
      enddo
      
      ! Print out the target and actual fuel masses for comparisons sake
      target_mass = 0
      if (itrees.eq.1) then
        do i=1,ntspecies
          do j=1,tfuelbins
            target_mass = target_mass + ntrees(i)*tbulkdensity(j,i)*PI*tcrowndiameter(1,i)**2./
     +        8.*(theight(1,i)-tcrownbotheight(1,i))
          enddo
        enddo
      endif
      if (itrees.ne.1) then
        do i=1,ntrees(1)
          do j=1,tfuelbins
            target_mass = target_mass + tbulkdensity(j,i)*PI*tcrowndiameter(i,1)**2./
     +        8.*(theight(i,1)-tcrownbotheight(i,1))
          enddo
        enddo
      endif
      print*,'Trees target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,ntspecies
        do i=tdnx(1),tdnx(2)
          do j=tdny(1),tdny(2)
            do k=1,zmax
              actual_mass = actual_mass+trhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      print*,'Trees actual fuel mass:',actual_mass
      print*,'Trees error:',actual_mass/target_mass*100,'%'
      
      end subroutine tree_baseline 

      subroutine litter_baseline
      !-----------------------------------------------------------------
      ! litter_baseline is a function which computes the characteristics
      ! of a base litter from trees and fills rhof, sizescale, moisture,
      ! and fueldepth arrays for litter and subracts mass of the grass 
      ! arrays for tree shading.
      !-----------------------------------------------------------------
      use constant_variables
      use grid_variables
      use baseline_variables
      
      implicit none

      integer i,j,k,ift,ift_grass
      real target_mass,actual_mass
      real rhoftemp,rhocolumn,coverfactor,shadefactor
      real,allocatable:: rhofxy(:,:)

      if (ilitter.eq.2) then
          print*,'Reading LLM tree litter'
          allocate(rhofxy(nx,ny))
          open(1, file="LLM_litter_trees.txt")
            read(1,*) rhofxy
          close(1)
      endif
      
      !----- Place litter on ground and remove grass to account for shading
      print*,'Placing litter and removing grass to account for shading'
      do ift = 1,ntspecies
        if (itrees.ne.1) theight(ift,1) = sum(theight(:,1),MASK=tspecies.eq.ift)/count(tspecies.eq.ift)
        do i=1,nx
          do j=1,ny
            if(ilitter.eq.2)lrho(ift)=rhofxy(i,j)
            ! Determine factors for placing litter and removing grass
            rhocolumn = 0
            do k=1,zmax
              rhocolumn = rhocolumn+sum(trhof((ift-1)*tfuelbins+1:ift*tfuelbins,i,j,k))*
     +          (zheight(i,j,k+1)-zheight(i,j,k))/theight(ift,1)
            enddo
            if (rhocolumn.gt.0) then
              shadefactor = exp(-grassconstant*rhocolumn/0.6)
              coverfactor = 1.-exp(-litterconstant*rhocolumn/0.6)

              ! Remove grass due to shadefactor
              if (igrass.eq.1) then
                do ift_grass=1,ngrass
                  do k=1,nz
                    if (zheight(i,j,k).gt.gdepth(ift_grass)) then
                      exit
                    else
                      rhof(ift_grass,i,j,k) = rhof(ift_grass,i,j,k)*shadefactor
                    endif
                  enddo
                enddo
              endif

              ! Add litter with dependence to coverfactor
              lfueldepth(ift,i,j) = coverfactor*ldepth(ift)
              do k=1,nz
                if (zheight(i,j,k).gt.lfueldepth(ift,i,j)) exit
                if (zheight(i,j,k+1).lt.lfueldepth(ift,i,j)) then
                  lrhof(ift,i,j,k) = lrho(ift)
                else 
                  lrhof(ift,i,j,k) = lrho(ift)*(lfueldepth(ift,i,j)-zheight(i,j,k))/(zheight(i,j,k+1)-zheight(i,j,k))
                endif
                lsizescale(ift,i,j,k) = lss(ift)
                lmoist(ift,i,j,k) = lmoisture(ift)
              enddo
            endif
          enddo
        enddo
        print*,'Finished litter for species',ift
      enddo
      
      ! Print out the target and actual fuel masses for comparisons sake
      target_mass = 0
      do i=1,ntspecies
        target_mass = target_mass + ntrees(i)*ldepth(i)*PI*tcrowndiameter(1,i)**2.*lrho(i)/4.
      enddo
      print*,'Litter target fuel mass:',target_mass
      actual_mass = 0
      do ift=1,ntspecies
        do i=1,nx
          do j=1,ny
            do k=1,zmax
              actual_mass = actual_mass+lrhof(ift,i,j,k)*dx*dy*(zheight(i,j,k+1)-zheight(i,j,k))
            enddo
          enddo
        enddo
      enddo
      !if (ilitter.eq.2)deallocate(rhofxy)
      print*,'Litter actual fuel mass:',actual_mass
      print*,'Litter error:',actual_mass/target_mass*100,'%'
      
      end subroutine litter_baseline 
