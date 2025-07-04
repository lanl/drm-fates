!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! variables declares all the constant variables used throughout the 
! program
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      module constant_variables
      !-----------------------------------------------------------------
      ! Constant variables and arrays
      !-----------------------------------------------------------------
      implicit none
            
      real  PI
      
      end module constant_variables
      
      module grid_variables
      !-----------------------------------------------------------------
      ! Grid and topo variables and arrays
      !-----------------------------------------------------------------
      implicit none
            
      integer nx,ny,nz
      real    dx,dy,dz
      real    aa1
      integer nfuel,zmax
      real,allocatable:: rhof(:,:,:,:),sizescale(:,:,:,:),moist(:,:,:,:),fueldepth(:,:,:,:)
      real,allocatable:: satarray(:,:,:)
      real,allocatable:: zs(:,:),zheight(:,:,:)
      character:: topofile*50
      
      end module grid_variables

      module baseline_variables
      !-----------------------------------------------------------------
      ! Tree and groundfuel variables unique to baseline establishment
      !-----------------------------------------------------------------
      implicit none
    
      integer:: igrass,itrees,ilitter 
      integer:: ngrass,ntspecies,tfuelbins
      real:: grassconstant,litterconstant
      real,allocatable:: grhof(:,:,:,:),gsizescale(:,:,:,:),gmoist(:,:,:,:),gfueldepth(:,:,:)
      real,allocatable:: trhof(:,:,:,:),tsizescale(:,:,:,:),tmoist(:,:,:,:),tfueldepth(:,:,:)
      real,allocatable:: lrhof(:,:,:,:),lsizescale(:,:,:,:),lmoist(:,:,:,:),lfueldepth(:,:,:)
      character:: grassfile*50,treefile*50,litterfile*50

      real,allocatable:: tcanopy(:),tlocation(:,:)
      real,allocatable:: tmoisture(:,:),tss(:,:),tbulkdensity(:,:)
      real,allocatable:: theight(:,:),tcrownbotheight(:,:),tcrownmaxheight(:,:),tcrowndiameter(:,:)
      integer,allocatable:: ntrees(:),tspecies(:)
      real,allocatable:: gdepth(:),grho(:),gss(:),gmoisture(:)
      real,allocatable:: ldepth(:),lrho(:),lss(:),lmoisture(:)
      integer,allocatable:: tdnx(:),tdny(:)

      end module baseline_variables

      module treatment_variables
      !-----------------------------------------------------------------
      ! Types of treatments which can occur
      !-----------------------------------------------------------------
      implicit none
      integer:: itreatment
      
      !-----------------------------------------------------------------
      ! Slash Treatment variables
      !-----------------------------------------------------------------
      
      integer,allocatable:: sdnx(:),sdny(:)
      real sdiameter,sheight
      real sdepth
      real sprho,smoist

      end module treatment_variables
