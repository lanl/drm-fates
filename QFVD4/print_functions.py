# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:14:09 2022

@author: zcope
"""
import os
import numpy as np

def print_gridlist(ri):
    nx, ny, nz, dx, dy = ri['nx'], ri['ny'], ri['nz'], ri['dx'], ri['dy']
    fname = os.path.join(ri['RUN_PATH'], 'gridlist')
    with open(fname, 'w') as input_file:
        input_file.write("&compresslist\n")
        input_file.write("       irst=0 nt=60000 nts=10 ntp=10\n")
        input_file.write("       nprocx=10 nprocy=20\n")
        input_file.write("       n={} m={} l={} aa1=1.\n".format(nx, ny, nz))
        input_file.write("       nv=8\n")
        input_file.write("       dx={} dy={} dz=1.\n".format(dx, dy))
        input_file.write("       dts=0.001  ! small time step\n")
        input_file.write("       ih=3   ! ghost cells\n")
        input_file.write("       iord=2 idiv=1 ! advection order, divergent flows, old non oscilator is defautl\n")
        input_file.write("\n")
        input_file.write("       ! bc and dampers\n")
        input_file.write("       ibcx=0 ibcy=0 irlx=1 irly=1         !irlx>=1 and irly>=1 diff=erent relaxation options\n")
        input_file.write("       ibclatopen=0 ibctopopen=0          !openoutlet, freeslip at top (maybe not fully operational)\n")
        input_file.write("       nr=5                                ! number of cell for lateral relaxation (irlx,irly)\n")
        input_file.write("       iab=1 tow=10.                 ! different relaxation option at top and magnitude parameter\n")
        input_file.write("       zab=450. zabt=450.0               ! dampers bottom (m), zabt is zab for theta for iab=3\n")
        input_file.write("\n")
        input_file.write("       !topo params\n")
        input_file.write("       topofile=''\n")
        input_file.write("       ipotflow=0               ! potflow can be used for correcting xe in case of topography\n")
        input_file.write("       slopeangle=0.0                   ! slope angle in degree for rotated gravity\n")
        input_file.write("       slopeazimuth=0.0                 ! slope azimuth for rotated gravity: 0 slope is along x axis, 90 is along y\n")
        input_file.write("\n")
        input_file.write("       ! drag and turbulence parameters\n")
        input_file.write("       irod=1 iturb=2                      ! iturb=1 ka only, iturb=2 ka and kb\n")
        input_file.write("       idrag=2                             ! idrag=1 rrl drag, idrag=2 inra drag\n")
        input_file.write("       isa=2                               ! isa=1 sa fixed, isa=2, sa is dxdydz**0.33, isa=3 : saxy different from saz\n")
        input_file.write("       rturbprandtl = 2.0  !(2.0 in inra's version)    ! inverse of turb prandtl number for scalar diffusion (theta, O2, etc.)\n")
        input_file.write("       ilapdo=0                            ! oxygen diffusion\n")
        input_file.write("\n")
        input_file.write("       !radiation parameters\n")
        input_file.write("       irad=2 icallrad=10 irandseed=0 iseed=0     ! radiation mode, frequency, and seed (for montecarlo irad=2)\n")
        input_file.write("       isootmodel=0 ! 0 is rrl's soot model, 1 is inra's soot model\n")
        input_file.write("       crad=50.0 ! for isootmodel=0 (inra uses 20, lanl uses 50)\n")
        input_file.write("       iradeastflux=0 ! specific radiative flux ouput through a vertical hard coded target in the east (x>0) direction\n")
        input_file.write("\n")
        input_file.write("       ! ambient conditions (temp, pres, wind) and env forces (corio, lspgf)\n")
        input_file.write("       u0=5.36 uramp=0.0 uramptime=0.0 uswitch=2 zu=24.0\n")
        input_file.write("       v0=0.0 vramp=0.0 vramptime=0.0 vswitch=0\n")
        input_file.write("       ius=0 iue=0 jus=0 jue=0  !indices of reference zone for lai computation (when 0 whole domain is used) for uswitch=2\n")
        input_file.write("       tambient=300.      ! Ambient temperature (K)\n")
        input_file.write("       iperturb=0   ! perturbation for cyclic runs: 1 theta, 2 and 3 : pinwheel\n")
        input_file.write("       itheta = 0  ! stable layer at the domain top\n")
        input_file.write("       pressground=1.0e5  ! Pa\n")
        input_file.write("       zgroundref=0.   ! reference elevation for tambient and pressground\n")
        input_file.write("       icorio=0            ! coriolis 1  (xe have already coriolis effect), 2 (geostrophic wind)\n")
        input_file.write("       ilspgf=0 frqlspgf=1000  ! large scale pgf mode 1 to 3\n")
        input_file.write("       izlspgf = 1  ! reference for max flux is whole domaine (0) or a slice at height zu (1)\n")
        input_file.write("\n")
        input_file.write("       ! fuel and fire\n")
        input_file.write("       ifuel=1         ! flag for the pdf choice\n")
        input_file.write("       ivegread=1       ! flag for reading fuel files\n")
        input_file.write("       rhomicrovalue=500.0  ! (kg/m3)  inra uses 700 (except for grass, 500); lanl uses 500\n")
        input_file.write("       cpwood=2500.0        ! (J/kg/K) inra uses 1800 (albini and stocks); lanl uses 2500\n")
        input_file.write("       ifuelinra=0      ! inra default fuels for inra's test suite\n")
        input_file.write("       fuelinranumber=0\n")
        input_file.write("       ffparam = 1.0       ! reaction rate factor\n")
        input_file.write("       ignfile='ignite_atv.5line.dat'\n")
        input_file.write("\n")
        input_file.write("       ! io\n")
        input_file.write("       frqoutput=100\n")
        input_file.write("       outname='output/comp.out'\n")
        input_file.write("       restartfile='output/comp.out.40000'\n")
        input_file.write("       ioextra=0                           ! this entails to have extra io : rnetsol, convht, tempg..\n")
        input_file.write("       iwallclock=0                ! flag for walltimers\n")
        input_file.write("\n")
        input_file.write("       ! windfield in or out parameters\n")
        input_file.write("       icfmeflag=0\n")
        input_file.write("       windfieldstartfile='windfieldstart'\n")
        input_file.write("       xvbdataname='../../wind17/mid.dry.250x600.12mph/wf/xvbdata'\n")
        input_file.write("       iwindfieldout=0             ! Windfield output switch: 1->interpolate, 2-> one file\n")
        input_file.write("       iwindfieldin=1              ! Windfield input switch: 1->interpolate, 2-> one file\n")
        input_file.write("       windspeedupfactor=1         ! factor of speed up between wind run and fire run (should be an integer)\n")
        input_file.write("       itwindfield=60000              ! Timestep to start saving data\n")
        input_file.write("       itinterp=20                 ! # of interpolation timesteps\n")
        input_file.write("       ibcells=5                   ! # of x cells saved on boundary\n")
        input_file.write("       jbcells=5                   ! # of y cells saved on boundary\n")
        input_file.write("       is=0 ie=0 js=0 je=0 ! indices for production of xvdata file on a subdomain (when 0 the whole domain is used) for iwindfieldout=1\n")
        input_file.write("\n")
        input_file.write("       ! personal flags\n")
        input_file.write("       ifp=0                 ! francois 's personal flag\n")
        input_file.write("\n")
        input_file.write("       ! parameters below are not supposed to change (other choice in option removed)\n")
        input_file.write("       isoturb=1            ! flag for isotrope turbulence; unsotrope has not been used for a long time...\n")
        input_file.write("       ibctopbot=1          ! bc at the bottom and top of mesh (edges instead of cell center)\n")
        input_file.write("       islip=0\n")
        input_file.write("       iunstable=0\n")
        input_file.write("       idirt=0          ! flag for rhodirt\n")
        input_file.write("       isplit=0\n")
        input_file.write("       irhovapor=0\n")
        input_file.write("       isor=0\n")
        input_file.write("       nonos=1 nfct=1 nonosold=1\n")
        input_file.write("       inonlocal=0 isubgridgas=0\n")
        input_file.write("       !st=1.0e-05\n")
        input_file.write("/\n")


def print_QFire_Advanced_User_Inputs_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QFire_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("0.05			! N/A, fraction of the cells on fire that will launch a firebrand\n")
        input_file.write("40.			! N/A, scaling factor of the radius represented by the firebrands launched\n")
        input_file.write("1				! s, time step for the firebrands trajectory calculation\n")
        input_file.write("10				! s, how often to launch firebrands\n")
        input_file.write("500			! N/A, number of firebrands distributed over the landing area\n")
        input_file.write("20.			! N/A, FB_FRACTION_LAUNCHED_to_RT_ratio\n")
        input_file.write("50.			! N/A, min_b_value_coef\n")
        input_file.write("0.75			! N/A, fb_frac_of_max_size\n")
        input_file.write("180				! s, germination_delay\n")
        input_file.write("10.				! N/A, fraction of the cell on fire (to scale w)\n")
        input_file.write("50				! N/A, minimum number of ignitions via firebrands at a point\n")
        input_file.write("100			! N/A, maximum number of ignitions via firebrands at a point\n")
        input_file.write("0.523598		! rad, min_theta_value (pi/6)\n")


def print_QFire_Bldg_Advanced_User_Inputs_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QFire_Bldg_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("1				! N/A, flag to convert QUIC-URB buildings to fuel (0 = no, 1 = yes)\n")
        input_file.write("0.5			! kg/m3, thin fuel density within buildings (if no fuel is specified)\n")
        input_file.write("2.			! N/A, attenuation coefficient within buildings\n")
        input_file.write("0.1	  ! m, surface roughness within buildings\n")
        input_file.write("1				! N/A, flag to convert fuel to canopy for winds (0 = no, 1 = yes)\n")
        input_file.write("1				! N/A, update canopy winds when fuel is consumed\n")
        input_file.write("1.			! N/A, attenuation coefficient within fuel\n")
        input_file.write("0.1	  ! m, surface roughness within fuel\n")
        input_file.write("\n")


def print_QFire_Plume_Advanced_User_Inputs_inp(ri):
    ws = ri['windspeed']
    fname = os.path.join(ri['RUN_PATH'], 'QFire_Plume_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("150000			! N/A, max number of plume at each time step\n")
        input_file.write("0.1			! m/s, minimum vertical velocity of a plume. If wc is below minimum, the plume is eliminated\n")
        input_file.write("10			! m/s, maximum vertical velocity of a plume\n")
        #Per Rod convo on 10/20/2021 for how to calc MSR:
        if 0.1 * ws > 0.5:
            msr = round(0.5/ws, 3)
        else: msr = 0.1
        #Ig min speed ratio * wind speed <= 0.5
        input_file.write("{}			! N/A, minimum speed ratio (plume vertical velocity/wind speed). If below, the plume is eliminated\n".format(msr))
        input_file.write("0			! 1/s^2, brunt vaisala frequency squared\n")
        input_file.write("1			! N/A creeping flag: 0 = off, 1 = on\n")
        input_file.write("0			! N/A, plume time step flag (0 = fixed, 1 = adaptable)\n")
        input_file.write("1			! s, time step plume\n")
        input_file.write("1			! s, sor option\n")
        input_file.write("2			! N/A, alpha 2 for fire cells\n")
        input_file.write("0			! N/A, fraction of ignition energy that goes into en 2 atmos\n")
        input_file.write("30			! deg, max angle for merging two plumes\n")
        input_file.write("0.7			! N/A, fraction of a plume length that a point on a second plume can be for merging\n")
        input_file.write("0.0      !Plume cutoff [m] (<0 no cutoff)\n")
        input_file.write("\n")
        input_file.write("\n")


def print_QP_buildout_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QP_buildout.inp')
    with open(fname, 'w') as input_file:
        input_file.write("           0  ! total number of buildings\n")
        input_file.write("           0  ! total number of vegitative canopies\n")


def print_QUIC_fire_inp(ri):
    sim_time, print_times, nz = ri['SimTime'], ri['print_times'], ri['nz']
    custom_ig, xmin, ymin, xlen, ylen = ri['custom_ig'], ri['ig_xmin'], ri['ig_ymin'], ri['ig_xlen'], ri['ig_ylen']
    fname = os.path.join(ri['RUN_PATH'], 'QUIC_fire.inp')
    with open(fname, 'w') as input_file:
        input_file.write("1					! Fire flag: 1 = for fire; 0 = no fire\n")
        input_file.write("-1				! Random number generator: -1: use time and date, any other integer > 0 is used as the seed\n")
        input_file.write("! FIRE TIMES\n")
        input_file.write("1488794400		! When the fire is ignited in Unix Epoch time (integer seconds since 1970/1/1 00:00:00)\n")
        input_file.write("{}				! Total simulation time for the fire [s]\n".format(sim_time))
        input_file.write("1		   		! time step for the fire simulation [s]\n")
        input_file.write("1					! Number of fire time steps done before updating the quic wind field (integer, >= 1)\n")
        input_file.write("{}				! After how many fire time steps to print out fire-related files (excluding emissions and radiation)\n".format(print_times))
        input_file.write("{}					! After how many quic updates to print out wind-related files\n".format(print_times))
        input_file.write("{}					! After how many fire time steps to average emissions and radiation\n".format(print_times))
        input_file.write("{}					! After how many quic updates to print out averaged wind-related files\n".format(print_times))
        input_file.write("! FIRE GRID\n")
        input_file.write("{}					! Number of vertical layers of fire grid cells (integer)\n".format(nz))
        input_file.write("1					! x - fire grid ratio = (QUIC-URB cell size)/(fire cell size), integer, can only be >= 1\n")
        input_file.write("1					! y - fire grid ratio = (QUIC-URB cell size)/(fire cell size), integer, can only be >= 1\n")
        input_file.write("0					! Vertical stretching flag: 0 = uniform dz, 1 = custom\n")
        input_file.write("1					! height of of uniform cells\n")
        input_file.write("! FOLDER OF TREES AND IGNITION FILES (full path, empty line if none) -- USE FILE SEPARATOR AT THE END\n")
        input_file.write("\"\"\n")
        input_file.write("1			! 1 = all fuels in one file, 2 = separate files\n")
        input_file.write("1			! 1 = stream, 2 = with headers\n")
        input_file.write("! FUEL\n")
        input_file.write("5					! fuel density flag: 1 = uniform; 2 = provided thru QF_FuelDensity.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        input_file.write("5					! fuel moisture flag: 1 = uniform; 2 = provided thru QF_FuelMoisture.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        input_file.write("! IGNITION LOCATIONS\n")
        if custom_ig:
            input_file.write("7					! 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp), 5 = time-dependent ignitions (QF_IgnitionPattern.inp), 6 = ignite.dat (firetech)\n")
        else:
            input_file.write("1					! 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp), 5 = time-dependent ignitions (QF_IgnitionPattern.inp), 6 = ignite.dat (firetech)\n")
            input_file.write("{}\n".format(xmin))
            input_file.write("{}\n".format(ymin))
            input_file.write("{}\n".format(xlen))
            input_file.write("{}\n".format(ylen))
        input_file.write("100\n")
        input_file.write("! FIREBRANDS\n")
        input_file.write("0				! 0 = off, 1 = on\n")
        input_file.write("! OUTPUT FILES (formats depend on the grid type flag)\n")
        input_file.write("0					! Output gridded energy-to-atmosphere (fire grid)\n")
        input_file.write("0					! Output compressed array reaction rate (fire grid)\n")
        input_file.write("1					! Output compressed array fuel density (fire grid)\n")
        input_file.write("0					! Output gridded wind (u,v,w,sigma) (fire grid)\n")
        input_file.write("0					! Output gridded QU winds with fire effects, instantaneous (QUIC-URB grid)\n")
        input_file.write("0					! Output gridded QU winds with fire effects, averaged (QUIC-URB grid)\n")
        input_file.write("0					! Output plume trajectories\n")
        input_file.write("0					! Output compressed array fuel moisture (fire grid)\n")
        input_file.write("0					! Output vertically-integrated % mass burnt (fire grid)\n")
        input_file.write("0					! Output gridded file with plumes locations (QUIC-URB grid)\n")
        input_file.write("0					! Output compressed array emissions (fire grid)\n")
        input_file.write("0					! Output gridded thermal radiation (fire grid)\n")

def print_QU_buildings_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QU_buildings.inp')
    with open(fname, 'w') as input_file:
        input_file.write("!QUIC 6.26\n")
        input_file.write("0.1			!Wall roughness length (m)\n")
        input_file.write("0			!Number of Buildings\n")
        input_file.write("0			!Number of Polygon Building Nodes\n")


def print_QU_fileoptions_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QU_fileoptions.inp')
    with open(fname, 'w') as input_file:
        input_file.write("!QUIC 6.26\n")
        input_file.write("4   !output data file format flag (1=ascii, 2=binary, 3=both, 4=none)\n")
        input_file.write("0   !flag to write out non-mass conserved initial field (uofield.dat) (1=write,0=no write)\n")
        input_file.write("0   !flag to write out the file uosensorfield.dat, the initial sensor velocity field (1=write,0=no write)\n")
        input_file.write("0   !flag to write out the file QU_staggered_velocity.bin used by QUIC-Pressure(1=write,0=no write)\n")
        input_file.write("1   !Output fire energy per timestep\n")
        input_file.write("0   !flag for automatically killing simulation once fire behavior has quit (1=on,0=off)\n")
        input_file.write("0   !flag to output startup wind files for topo-influenced wind fields\n")


def print_QU_metparams_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QU_metparams.inp')
    with open(fname, 'w') as input_file:
        input_file.write("!QUIC 6.26\n")
        input_file.write("0 !Met input flag (0=QUIC,1=WRF,2=ITT MM5,3=HOTMAC)\n")
        input_file.write("1 !Number of measuring sites\n")
        input_file.write("1 !Maximum size of data points profiles\n")
        input_file.write("sensor1 !Site Name\n")
        input_file.write("!File name\n")
        input_file.write("sensor1.inp\n")


def print_QU_movingcoords_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QU_movingcoords.inp')
    with open(fname, 'w') as input_file:
        input_file.write("!QUIC 6.3\n")
        input_file.write("0   !Moving coordinates flag (0=no, 1=yes)\n")
        input_file.write("0   !Reference bearing of the ship relative to the non-rotated domain (degrees)\n")
        input_file.write("!Time (Unix Time), Ship Speed (m/s), Ship Bearing (deg), Ocean Current Speed (m/s), Ocean Current Direction (deg)\n")
        input_file.write("1488794400	0	0	0	0\n")
        input_file.write("1488794450	0	0	0	0\n")
        input_file.write("1488794500	0	0	0	0\n")
        input_file.write("1488794550	0	0	0	0\n")
        input_file.write("1488794600	0	0	0	0\n")
        input_file.write("1488794650	0	0	0	0\n")
        input_file.write("1488794700	0	0	0	0\n")
        input_file.write("1488794750	0	0	0	0\n")
        input_file.write("1488794800	0	0	0	0\n")
        input_file.write("1488794850	0	0	0	0\n")
        input_file.write("1488794900	0	0	0	0\n")
        input_file.write("1488794950	0	0	0	0\n")
        input_file.write("1488795000	0	0	0	0\n")


def print_QU_simparams_inp(ri):
    nx, ny, dx, dy = ri['nx'], ri['ny'], ri['dx'], ri['dy']
    max_topo = ri['max_topo']
    z_arr = calculate_QF_vertical_wind_grid(max_topo)
    nz = len(z_arr)
    fname = os.path.join(ri['RUN_PATH'], 'QU_simparams.inp')
    with open(fname, 'w') as input_file:
        input_file.write("!QUIC 6.26\n")
        input_file.write("{} !nx - Domain Length(X) Grid Cells\n".format(nx))
        input_file.write("{} !ny - Domain Width(Y) Grid Cells\n".format(ny))
        input_file.write("{} !nz - Domain Height(Z) Grid Cells\n".format(nz))
        input_file.write("{} !dx (meters)\n".format(dx))
        input_file.write("{} !dy (meters)\n".format(dy))
        input_file.write("3 !Vertical stretching flag(0=uniform,1=custom,2=parabolic Z,3=parabolic DZ,4=exponential)\n")
        input_file.write("1.500000 !Surface dz (meters)\n")
        input_file.write("4 !Number of uniform surface cells\n")
        input_file.write("!dz array (meters)\n")
        for i in z_arr:
            input_file.write("{}\n".format(i))
        input_file.write("1 !total time increments\n")
        input_file.write("0 !UTC conversion\n")
        input_file.write("!Begining of time step in Unix Epoch time (integer seconds since 1970/1/1 00:00:00)\n")
        input_file.write("1488794400\n")
        input_file.write("2 !rooftop flag (0-none, 1-log profile, 2-vortex)\n")
        input_file.write("3 !upwind cavity flag (0-none, 1-Rockle, 2-MVP, 3-HMVP)\n")
        input_file.write("4 !street canyon flag (0-none, 1-Roeckle, 2-CPB, 3-exp. param. PKK, 4-Roeckle w/ Fackrel)\n")
        input_file.write("1 !street intersection flag (0-off, 1-on)\n")
        input_file.write("3 !wake flag (0-none, 1-Rockle, 2-Modified Rockle, 3-Area Scaled)\n")
        input_file.write("1 !sidewall flag (0-off, 1-on)\n")
        input_file.write("2 !Canopy flag (1-Cionco w/o wakes, 2-Cionco w/ wakes)\n")
        input_file.write("1 !Season flag (1-Summer, 2-Winter, 3-Transition)\n")
        input_file.write("10 !Maximum number of iterations\n")
        # input_file.write("1.1 !omegarelax\n") #this line may be required with some versions of QFVD4
        input_file.write("3 !Residual Reduction (Orders of Magnitude)\n")
        input_file.write("0 !Use Diffusion Algorithm (1 = on)\n")
        input_file.write("20 !Number of Diffusion iterations\n")
        input_file.write("0 !Domain rotation relative to true north (cw = +)\n")
        input_file.write("0.0  !UTMX of domain origin (m)\n")
        input_file.write("0.   !UTMY of domain origin (m)\n")
        input_file.write("1 !UTM zone\n")
        input_file.write("17 !UTM zone leter (1=A,2=B,etc.)\n")
        input_file.write("0 !QUIC-CFD Flag\n")
        input_file.write("0 !Explosive building damage flag (1 = on)\n")
        input_file.write("0 !Building Array Flag (1 = on)\n")


def print_rasterorigin_txt(ri):
    fname = os.path.join(ri['RUN_PATH'], 'rasterorigin.txt')
    with open(fname, 'w') as input_file:
        input_file.write("0.\n")
        input_file.write("0.\n")
        input_file.write("752265.868913356\n")
        input_file.write("3752846.04249607\n")
        input_file.write("742265.868913356\n")
        input_file.write("3742846.04249607\n")
        input_file.write("10000\n")


def print_Runtime_Advanced_User_Inputs_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'Runtime_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("8\n")


def print_sensor1_inp(ri):
    ws, wd = ri['windspeed'], ri['winddir']
    fname = os.path.join(ri['RUN_PATH'], 'sensor1.inp')
    with open(fname, 'w') as input_file:
        input_file.write("sensor1 !Site Name\n")
        input_file.write("0 !Upper level flag\n")
        input_file.write("50 !Upper level height\n")
        input_file.write("1 !Site Coordinate Flag\n")
        input_file.write("1 !X coordinate\n")
        input_file.write("1 !Y coordinate\n")
        input_file.write("1488794400 !Time in unix Epoch Time\n")
        input_file.write("1 !site boundary layer\n")
        input_file.write("0.01 !site zo\n")
        input_file.write("0.0 !reciprocal Monin-Obukhov Length(1/m)\n")
        input_file.write("!height, speed, direction\n")
        input_file.write("6.1 {} {}\n".format(ws,wd))


def print_topo_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'topo.inp')
    if ri['topo_custom']:
        topo_flag = 5
        topo_fname = "topo.dat"
    else:
        topo_flag = 0
        topo_fname = '""'
    with open(fname, 'w') as input_file:
        input_file.write("!Relative filepath to topo .dat file (ex:: \"../path/to/topo.dat\")\n")
        input_file.write("{}\n".format(topo_fname))
        input_file.write("{}              !Topo flag 0:Flat 1:Gaussian Hill 3:Constant slope with flat section 5:Custom .dat\n".format(topo_flag))
        input_file.write("0              !Smoothing Flag\n")
        input_file.write("0              !Total startup iterations\n")
        input_file.write("0            !Startup Iterations Restart Period (set to zero for no restart)\n")
        input_file.write("0              !Preconditioning\n")

def calculate_QF_vertical_wind_grid(max_topo): 

    c = 1.41/2.0                # 1/sqrt(2)
    e = 2.71828/2.0             # Eulers Number / 2
    zo = e**c                   # Initial height for log profile      
    dist = zo                   # Keeps track of the max height
    zt = zo                     # Keeps track of previous cell's height
    zs = []                     # List storing the cell heights
    zs.append(zo)
    max_topo_adjusted = np.max([max_topo*3.,400.]) # Covering edge scenario for flat topo
    while(dist < max_topo_adjusted):
        zo = 1.125*zt*e**c
        dist+=zo
        zt = zo
        zs.append(zo)
        
    print (zs,len(zs))
    
    return zs