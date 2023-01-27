# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 10:08:25 2022

@author: zcope
"""
import os
import numpy as np

def print_gridlist(ri):
    nx, ny, nz, dx, dy = ri['nx'], ri['ny'], ri['nz'], ri['dx'], ri['dy']
    fname = os.path.join(ri['RUN_PATH'], 'gridlist')
    with open(fname, 'w') as input_file:
        input_file.write("       n={} m={} l={} aa1=1.\n".format(nx, ny, nz))
        input_file.write("       dx={} dy={} dz=1.\n".format(dx, dy))



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
        input_file.write("180			! s, germination_delay\n")
        input_file.write("5.				! N/A, fraction of the cell on fire (to scale w)\n")
        input_file.write("50				! N/A, minimum number of ignitions via firebrands at a point\n")
        input_file.write("100			! N/A, maximum number of ignitions via firebrands at a point\n")
        input_file.write("0.523598		! rad, min_theta_value (pi/6)\n")
        input_file.write("0.03        ! m, fb thickness\n")


def print_QFire_Bldg_Advanced_User_Inputs_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QFire_Bldg_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("1			! N/A, flag to convert QUIC-URB buildings to fuel (0 = no, 1 = yes)\n")
        input_file.write("0.5		! kg/m3, thin fuel density within buildings (if no fuel is specified)\n")
        input_file.write("2.			! N/A, attenuation coefficient within buildings\n")
        input_file.write("0.1	   ! m, surface roughness within buildings\n")
        input_file.write("1			! N/A, flag to convert fuel to canopy for winds (0 = no, 1 = yes)\n")
        input_file.write("1			! N/A, update canopy winds when fuel is consumed (0 = no, 1 = yes)\n")
        input_file.write("1			! N/A, attenuation coefficient within fuel for the Cionco profile (default = 1)\n")
        input_file.write("0.1	   ! m, surface roughness within fuel (default = 0.1 m)\n")
        input_file.write("\n")


def print_QFire_Plume_Advanced_User_Inputs_inp(ri):
    ws = ri['windspeed']
    fname = os.path.join(ri['RUN_PATH'], 'QFire_Plume_Advanced_User_Inputs.inp')
    with open(fname, 'w') as input_file:
        input_file.write("150000			! N/A, max number of plume at each time step\n")
        input_file.write("0.1			! m/s, minimum vertical velocity of a plume. If wc is below minimum, the plume is eliminated\n")
        input_file.write("100				! m/s, maximum vertical velocity of a plume (default 100 m/s)\n")
        #Per Rod convo on 10/20/2021 for how to calc MSR:
        if 0.1 * ws > 0.5:
            msr = round(0.5/ws, 3)
        else: msr = 0.1
        #Ig min speed ratio * wind speed <= 0.5
        input_file.write("{}			! N/A, minimum speed ratio (plume vertical velocity/wind speed). If below, the plume is eliminated\n".format(msr))
        input_file.write("0.					! 1/s2, brunt vaisala frequency squared\n")
        input_file.write("1					! N/A creeping flag: 0 = off, 1 = on\n")
        input_file.write("0					! Flag time step (0 = constant, 1 = adaptive)\n")
        input_file.write("1              ! s, flag = 0: plumes time step plume, flag = 1: minimum plumes time step\n")
        input_file.write("1              ! sor option (0 = reset lambda at each call, 1 = keep)\n")
        input_file.write("10             ! alpha 2 - plume centerline\n")
        input_file.write("2              ! alpha 2 - plume edges\n")
        input_file.write("30.            ! deg, max angle between plumes to merging\n")
        input_file.write("0.7            ! N/A, fraction of overlap of plumes for merging\n")
        input_file.write("1              ! N/A, which wplume-to-grid scheme: 0 = new, 1 = old\n")
        input_file.write("10             ! N/A, number of points per edge for new scheme (only new scheme)\n")
        input_file.write("1              ! N/A, w computation: 1 = max(w), 0 = sum(w^3)\n")


def print_QP_buildout_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QP_buildout.inp')
    with open(fname, 'w') as input_file:
        input_file.write("           0  ! total number of buildings\n")
        input_file.write("           0  ! total number of vegitative canopies\n")


def print_QUIC_fire_inp(ri):
    sim_time, print_times, nz, ig_method = ri['SimTime'], ri['print_times'], ri['nz'], ri['ig_method']
    custom_ig, xmin, ymin, xlen, ylen = ri['custom_ig'], ri['ig_xmin'], ri['ig_ymin'], ri['ig_xlen'], ri['ig_ylen']
    fname = os.path.join(ri['RUN_PATH'], 'QUIC_fire.inp')
    if ig_method == "aerial":
        energy_packets = 5
    elif ig_method == "drip":
        energy_packets = 1
    elif ig_method == "total":
        energy_packets = 100
    # ADD IG METHODS
    with open(fname, 'w') as input_file:
        input_file.write("1					! Fire flag: 1 = run fire; 0 = no fire\n")
        input_file.write("-1				! Random number generator: -1: use time and date, any other integer > 0 is used as the seed\n")
        input_file.write("! FIRE TIMES\n")
        input_file.write("0		! When the fire is ignited in Unix Epoch time (integer seconds since 1970/1/1 00:00:00). Must be greater or equal to the time of the first wind\n")
        input_file.write("{}				! Total simulation time for the fire [s]\n".format(sim_time))
        input_file.write("1					! time step for the fire simulation [s]\n")
        input_file.write("1					! Number of fire time steps done before updating the quic wind field (integer, >= 1)\n")
        input_file.write("{}					! After how many fire time steps to print out fire-related files (excluding emissions and radiation)\n".format(print_times))
        input_file.write("{}					! After how many quic updates to print out wind-related files\n".format(print_times))
        input_file.write("{}					! After how many fire time steps to average emissions and radiation\n".format(print_times))
        input_file.write("{}					! After how many quic updates to print out averaged wind-related files\n".format(print_times))
        input_file.write("! FIRE GRID\n")
        input_file.write("{}					! Number of vertical layers of fire grid cells (integer)\n".format(nz))
        input_file.write("0					! Vertical stretching flag: 0 = uniform dz, 1 = custom\n")
        input_file.write("1.             ! m, dz\n")
        input_file.write("! FILE PATH\n")
        input_file.write("\"\"\n")
        input_file.write("1              ! Fuel types are in separate files\n")
        input_file.write("1              ! File is stream (1) or headers (2)\n")
        input_file.write("! FUEL\n")
        input_file.write("4					! fuel density flag: 1 = uniform; 2 = provided thru QF_FuelDensity.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        input_file.write("4					! fuel moisture flag: 1 = uniform; 2 = provided thru QF_FuelMoisture.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        #Might want to change back
        # input_file.write("5					! fuel density flag: 1 = uniform; 2 = provided thru QF_FuelDensity.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        # input_file.write("5					! fuel moisture flag: 1 = uniform; 2 = provided thru QF_FuelMoisture.inp, 3 = Firetech files for quic grid, 4 = Firetech files for different grid (need interpolation)\n")
        input_file.write("! IGNITION LOCATIONS\n")
        if custom_ig:
            input_file.write("7					! 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp), 5 = time-dependent ignitions (QF_IgnitionPattern.inp), 6 = ignite.dat (firetech)\n")
        else:
            input_file.write("1					! 1 = rectangle, 2 = square ring, 3 = circular ring, 4 = file (QF_Ignitions.inp), 5 = time-dependent ignitions (QF_IgnitionPattern.inp), 6 = ignite.dat (firetech)\n")
            input_file.write("{}\n".format(xmin))
            input_file.write("{}\n".format(ymin))
            input_file.write("{}\n".format(xlen))
            input_file.write("{}\n".format(ylen))
        input_file.write("{}\n".format(energy_packets))
        input_file.write("! FIREBRANDS\n")
        input_file.write("0				! 0 = off, 1 = on\n")
        input_file.write("! OUTPUT FILES (formats depend on the grid type flag)\n")
        input_file.write("1					! Output gridded energy-to-atmosphere (3D fire grid + extra layers)\n")
        input_file.write("0					! Output compressed array reaction rate (fire grid)\n")
        input_file.write("1					! Output compressed array fuel density (fire grid)\n")
        input_file.write("0					! Output gridded wind (u,v,w,sigma) (3D fire grid)\n")
        input_file.write("1					! Output gridded QU winds with fire effects, instantaneous (QUIC-URB grid)\n")
        input_file.write("0					! Output gridded QU winds with fire effects, averaged (QUIC-URB grid)\n")
        input_file.write("0					! Output plume trajectories (ONLY FOR DEBUG)\n")
        input_file.write("1					! Output compressed array fuel moisture (fire grid)\n")
        input_file.write("0					! Output vertically-integrated % mass burnt (fire grid)\n")
        input_file.write("0					! Output trajectories firebrands\n")
        input_file.write("0					! Output compressed array emissions (fire grid)\n")
        input_file.write("0					! Output gridded thermal radiation (fire grid)\n")
        input_file.write("1              ! Output surface fire intensity at every fire time step\n")
        input_file.write("! AUTOKILL\n")
        input_file.write("0              ! Kill if the fire is out and there are no more ignitions or firebrands (0 = no, 1 = yes)\n")


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
        input_file.write("\n")
        input_file.write("2   !output data file format flag (1=ascii, 2=binary, 3=both)\n")
        input_file.write("0   !flag to write out non-mass conserved initial field (uofield.dat) (1=write,0=no write)\n")
        input_file.write("0   !flag to write out the file uosensorfield.dat, the initial sensor velocity field (1=write,0=no write)\n")
        input_file.write("0   !flag to write out the file QU_staggered_velocity.bin used by QUIC-Pressure(1=write,0=no write)\n")
        input_file.write("0   !flag to generate startup files\n")


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
        input_file.write("{} ! NZ - domain height(Z) grid cells\n".format(nz))
        input_file.write("{} !dx (meters)\n".format(dx))
        input_file.write("{} !dy (meters)\n".format(dy))
        input_file.write("3 !Vertical stretching flag (0 = uniform, 1 = custom, 2 = parabolic Z, 3 = parabolic DZ, 4 = exponential)\n")
        input_file.write("1.000000 ! Surface DZ [m]\n")
        input_file.write("5 ! Number of uniform surface cells\n")
        input_file.write("! DZ array [m]\n")
        for i in z_arr:
            input_file.write("{}\n".format(i))
        input_file.write("1 ! Number of time increments\n")
        input_file.write("0 ! UTC conversion [hours]\n")
        input_file.write("! Begining of time step in Unix epoch time [integer seconds since 1970/1/1 00:00:00]\n")
        input_file.write("0\n")
        input_file.write("2 ! Rooftop flag (0 = none, 1 = log profile, 2 = vortex)\n")
        input_file.write("3 ! Upwind cavity flag (0 = none, 1 = Rockle, 2 = MVP, 3 = HMVP)\n")
        input_file.write("4 ! Street canyon flag (0 = none, 1 = Roeckle, 2 = CPB, 3 = exp. param. PKK, 4 = Roeckle w/ Fackrel)\n")
        input_file.write("1 ! Street intersection flag (0 = off, 1 = on)\n")
        input_file.write("4 ! Wake flag (0 = none, 1 = Rockle, 2 = modified Rockle, 3 = area scaled)\n")
        input_file.write("1 ! Sidewall flag (0 = off, 1 = on)\n")
        input_file.write("1 ! Individual tree wake flag\n")
        input_file.write("2 ! Canopy flag (0 = off, 1 = Cionco w/o wakes, 2 = Cionco w/ wakes)\n")
        input_file.write("1 ! Season flag (1 = Summer, 2 = Winter, 3 = Transition)\n")
        input_file.write("10 ! Maximum number of iterations\n")
        input_file.write("3 ! Residual reduction (Orders of Magnitude)\n")
        input_file.write("0 ! Use diffusion algorithm (0 = off, 1 = on)\n")
        input_file.write("20 ! Number of diffusion iterations\n")
        input_file.write("0 ! Domain rotation relative to true north [degrees] (cw = +)\n")
        input_file.write("0.0 ! UTMX of domain origin [m]\n")
        input_file.write("0.0 ! UTMY of domain origin [m]\n")
        input_file.write("1 ! UTM zone\n")
        input_file.write("17 ! UTM zone leter (1 = A, 2 = B, etc.)\n")
        input_file.write("0 ! QUIC-CFD Flag\n")
        input_file.write("0 ! Explosive building damage flag (1 = on)\n")
        input_file.write("0 ! Building Array Flag (1 = on)\n")

def print_QU_TopoInputs_inp(ri):
    fname = os.path.join(ri['RUN_PATH'], 'QU_TopoInputs.inp')
    if ri['topo_custom']:
        topo_flag = 11
        topo_fname = "topo.dat"
    else:
        topo_flag = 0
        topo_fname = '""'
    with open(fname, 'w') as input_file:
        input_file.write("Specify file name for custom topo (full path)\n")
        input_file.write("{}\n".format(topo_fname))
        input_file.write("{}     ! N/A, topo flag: 0 = flat, 1 = Gaussian hill, 2 = hill pass, 3 = slope mesa, 4 = canyon, 5 = custom, 6 = half circle, 7 = sinusoid, 8 = cos hill, 9 = QP_elevation.inp, 10 = terrainOutput.txt (ARA), 11 = terrain.dat (firetec)\n".format(topo_flag))
        input_file.write("0     ! N/A, smoothing method (0 = none (default for idealized terrain), 1 = blur, 2 = David's)\n")
        input_file.write("0     ! N/A, number of smoothing passes\n")
        input_file.write("100   ! N/A, number of initial SOR iterations (only if fire is run)\n")
        input_file.write("4     ! N/A, sor cycles\n")
        input_file.write("1.3   ! N/A, sor relaxation parameter (default for flat is 1.78, ignored if there is no terrain)\n")
        input_file.write("0     ! N/A, add slope flow correction (0 = no, 1 = yes)\n")
        input_file.write("Command to run slope flow code (system specific)\n")
        input_file.write("\"\"\n")

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
        input_file.write("8  ! max number of cpu\n")


def print_sensor1_inp(ri):
    ws, wd = ri['windspeed'], ri['winddir']
    fname = os.path.join(ri['RUN_PATH'], 'sensor1.inp')
    with open(fname, 'w') as input_file:
        input_file.write("sensor1 !Site Name\n")
        input_file.write("0 !Upper level flag (1 = use this profile for upper level winds)\n")
        input_file.write("50 !Upper level height (meters)\n")
        input_file.write("1 !Site Coordinate Flag (1=QUIC, 2=UTM, 3=Lat/Lon)\n")
        input_file.write("1 !X coordinate (meters)\n")
        input_file.write("1 !Y coordinate (meters)\n")
        input_file.write("1653321600 !Begining of time step in Unix Epoch time (integer seconds since 1970/1/1 00:00:00)\n")
        input_file.write("1 !site boundary layer flag (1 = log, 2 = exp, 3 = urban canopy, 4 = discrete data points)\n")
        input_file.write("0.1 !site zo\n")
        input_file.write("0. ! 1/L\n")
        input_file.write("!Height (m),Speed	(m/s), Direction (deg relative to true N)\n")
        input_file.write("6.1 {} {}\n".format(ws,wd))
        input_file.write("\n")

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