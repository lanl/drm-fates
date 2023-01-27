# -*- coding: utf-8 -*-
"""
Created on Fri Jan 27 08:38:07 2023

@author: Niko Tutland
"""

import os
qf_options = {}
OG_PATH = os.getcwd()
PROJ_PATH = os.path.join(OG_PATH,"7.QUICFIRE-MODEL","projects")

###Settings to change###

# Simulation parameters
qf_options['QFVD'] = 4 #version 4 or 5
PROJ_FOLDER = "LandisTester" #folder containing QF inputs
qf_options['SimTime'] = 60
qf_options['print_times'] = 20

# Domain settings
qf_options['nx'] = 200 #starting values, may be modified in code
qf_options['ny'] = 200 #"
qf_options['nz'] = 22  #"

# Topo settings
qf_options['topo_custom'] = False
qf_options['max_topo'] = 0 #set for flat right now, may use custom topo in the future

# Wind settings
qf_options['windspeed'] = 6.0 #m/s
qf_options['winddir'] = 270

# Ignitions settings
qf_options['custom_ig'] = False
qf_options['ig_xmin'] = 395.0 #m (Adam's original pattern here. Modified in code for LANDIS runs)
qf_options['ig_ymin'] = 2.0
qf_options['ig_xlen'] = 5.0
qf_options['ig_ylen'] = 397.0
#below are options for potential future use of driptorch
qf_options['ig_method'] = "total" #options: "aerial", "drip", "total". Determines number of energy packets at ignitions
qf_options['ig_pattern'] = "strip" #options: "strip", "flank", "head", "back"
qf_options['ig_spacing'] = 20 #m

###DO NOT CHANGE###
qf_options['RUN_PATH'] = os.path.join(PROJ_PATH,PROJ_FOLDER)
qf_options['dx'] = qf_options['dy'] = 2 #m
qf_options['dz'] = 1