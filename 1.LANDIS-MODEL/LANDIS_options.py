# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 08:04:05 2022

@author: Niko Tutland
"""

opdict = {}

# Run info
opdict['states'] = ["CA","OR","WA","ID","NV"] #states containing and bordering study area
opdict['fia_spec'] = ["PSME","PILA","QUCH2","PIPO","ARME","CADE27","PICO",
                      "QUKE","LIDE3","ACMA3","ABGR","PIMO3","TABR2","ALRU2"] #fia species symbols for tree species in study area
opdict['landis_spec'] = ["PseuMenz","PinuLamb","QuerChry","PinuPond","ArbuMenz","CaloDecu","PinuCont", #landis designation for tree species in study area
                         "QuerKell","LithDens","AcerMacr","AbieGran","PinuMont","TaxuBrev","AlnuRubr"] #(input in same order as fia_spec)
opdict['region_flag'] = 1 # where is the AOI? 1 = California, 2 = other western state, 3 = midwest, eastern, or southern state
opdict['age_bin'] = 10
opdict['aoi_elev'] = 4000 #elevation of the study area in feet [should we calculate this in the code?]

# Fuels
opdict['bulk_density'] = 0.7 #constant canopy bulk density (kg/m^3)
opdict['cl_factor'] = 0.8 #percent of crown above point of maximum crown diameter (m)
opdict['moisture'] = 1 #constant canopy moisture content
opdict['sizescale'] = 0.0005 #canopy fuel radius (nominal fuel size - solids) (m)

# Spatial info
opdict['L2_res'] = 150 #landis resolution in m
opdict['QF_res'] = 2 # quic-fire resolution in m

## Make sure the community output biomass extension prints outputs at an interval that will capture both the spinup end year and the cycles' end year