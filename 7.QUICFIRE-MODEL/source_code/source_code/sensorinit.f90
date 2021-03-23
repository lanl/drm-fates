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
	subroutine sensorinit
	! ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	! Subroutine to initialize the initial velocity field (uo,vo,wo) for each time step
	!
	! This program interpolates irregularly spaced data onto a
	! uniform grid using the Barnes Objective Map Analysis
	! scheme as implemented in Koch et al., 1983
	!
	! This program has been modified from Eric Pardyjak's original 2D code to work with
	! quicurb 3D
	!
	! this subroutine uses a quasi-3D barnes objective mapping technique
	! quasi-3D is just using sensor height (zref) and sensor wind speed (uref)
	! to extrapolate a vertical velocity profile at each sensor location
	! to get a velocity at every height at each location
	! from these hieght varying velocities, a regular 2D barnes objective map analysis
	! is done at each planar height throughout the whole domain.
	!
	! Called by met_init.f90
	!
	!
	!
	! Tom Booth 2/17/04
	!
	! Variable information:
	!  site_xcoord,site_ycoord,site_zcoord - the coordinates of each site (meters)
	!  t_site - is the time step for each site
	!  dir_site - is the direction of the wind for each site at each time step
	!  vel_site - is the magnitude of the wind for each site at each time step
	!  total_time_step - is the total number of time steps in a 24 hr period
	!  num_sites - is the total number of measurement sites
	!  sgamma - numerical convergence parameter
	!  lamda - weight parameter (ko)
	!  deln - average Radius (i.e. computed data spacing)
	! TMB/ERP 9/20/05
	!  Added a logarithmic interpolation below the the lowest input data point for
	!  input files that read in wind profile data

	!
	!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	use winds_module
	use grid_module
	
	implicit none

	quwinds%wo = 0.

	CALL INTERPOLATEWINDS()

	quwinds%max_velmag = 1.2 * maxval(sqrt(quwinds%uo(:,:,qugrid%nz)**2 + quwinds%vo(:,:,qugrid%nz)**2))
	
	return
	
	end
