# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 12:12:26 2022

@author: Niko Tutland
"""

import os
import sys
import subprocess
# import numpy as np
# import pandas as pd
import glob
import re
from time import sleep
import rasterio as rio
from LANDIS_options import opdict

def Landis(lp):
    
    batch_cmd = os.path.join(lp.landis_path,lp.batch_file)
    
    if lp.spinup == False:
        replace_IC(lp)
    
    replace_duration(lp)
    
    try:
        with subprocess.Popen(
            [batch_cmd], stdout=subprocess.PIPE
        ) as process:

            def poll_and_read():
                print(f"{process.stdout.read1().decode('utf-8')}")
            
            while process.poll() != 0:
                poll_and_read()
                sleep(1)
    except FileNotFoundError as exc:
        print(f"Batchfile not found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print(f"LANDIS run failed with return code {exc.returncode}\n{exc}")
        
    ### CROP DOMAIN HERE ###

class LandisParams:
    """
    Class containing parameters for controlling the LANDIS runs/loops
    """
    def __init__(self, OG_PATH, nyears, ncycyear, ncycle, cycle, spinup):  
        self.OG_PATH = OG_PATH               #path to the drm-fates directory
        
        self.nyears = int(nyears)            #number of years for spinup and transient runs
        self.ncycyear = int(ncycyear)        #number of cyclical year run
        self.ncycle = int(ncycle)            #number of loops
        self.cycle = int(cycle)              #current iteration
        self.spinup = spinup                 #is this the initial run?
        if spinup == True:
            year = int(nyears)
        else:
            year = int(ncycyear)
        self.year = int(year)
        self.states = opdict['states']
        self.fia_spec = opdict['fia_spec']
        self.landis_spec = opdict['landis_spec']
        self.region_flag = opdict['region_flag']
        self.age_bin = opdict['age_bin']
        self.aoi_elev = opdict['aoi_elev']
        self.bulk_density = opdict['bulk_density']
        self.cl_factor = opdict['cl_factor']
        self.moisture = opdict['moisture']
        self.sizescale = opdict['sizescale']
        self.QF_res = 2
        
        batch_file, scenario_file, necn_file, species_file, IC_file, IC_map = get_filenames(OG_PATH)
        landis_path = os.path.join(OG_PATH, "1.LANDIS-MODEL","LANDIS_run")
        
        self.landis_path = landis_path
        self.CIF_file = "community-input-file-"+str(year)+".csv" #name of community biomass output csv
        self.OC_file = "output-community-"+str(year)+".img" #name of community biomass output raster
        self.scenario_file = str(scenario_file)   #name of LANDIS scenario input file
        self.necn_file = str(necn_file)           #name of NECN input file
        self.batch_file = str(batch_file)     #name of LANDIS batchfile
        self.species_file = str(species_file) #name of species input file
        self.IC_file = str(IC_file)           #name of initial communities file
        self.IC_map = str(IC_map)             #name of initial communities raster
        
        with rio.open(os.path.join(landis_path, IC_map), 'r+') as IC_map :
            L2_res = IC_map.transform[0]
        self.L2_res = L2_res
        
        CIF_cropped = "community-input-file-"+str(year)+"_cropped.csv"
        IC_cropped = "IC_cycle"+str(cycle)+"_cropped.tif"
        OC_cropped = "IC_cycle"+str(cycle)+".tif"
        self.CIF_cropped = CIF_cropped
        self.IC_cropped = IC_cropped
        self.OC_cropped = OC_cropped

def get_filenames(path):
    os.chdir(os.path.join(path, "1.LANDIS-MODEL","LANDIS_run"))
    # Find batchfile and make sure it's the only one
    batchfile_list = glob.glob("*.bat") #case insensitive for windows, which should be what we want
    batch_file = batchfile_list[0]
    # Identify scenario file name from batchfile
    scenario_file = read_filename(batch_file,"call",1)
    # Identify NECN file name from scenario file
    necn_file = read_filename(scenario_file, '"NECN Succession"',1)
    # Identify Species file from scenario file
    species_file = read_filename(scenario_file, "Species",1)
    # Identify Initial Communities file and raster from NECN file
    IC_file = read_filename(necn_file, "InitialCommunities",1)
    IC_map = read_filename(necn_file, "InitialCommunities",2)
    return batch_file, scenario_file, necn_file, species_file, IC_file, IC_map

def read_filename(file,line_id,which_one):
    with open(file, 'r', encoding='utf-8') as file:
        filelist = file.readlines()
    lines = [match for match in filelist if line_id in match]
    line = lines[which_one-1]
    file_name = re.split(" |\t", line)[-1].strip()
    return file_name

def replace_duration(lp):
    if lp.spinup == True:
        runlen = lp.nyears
    else:
        runlen = lp.ncycyear
    
    with open(os.path.join(lp.landis_path,lp.scenario_file), 'r', encoding='utf-8') as file:
        filelist = file.readlines()

    matches = [match for match in filelist if "Duration" in match]
    s = matches[0]
    matched_indexes = []
    i = 0
    length = len(filelist)
    while i < length:
        if s == filelist[i]:
            matched_indexes.append(i)
        i += 1
    durationline = matched_indexes[0]
    
    filelist[durationline] = "Duration {}\n".format(runlen)
  
    with open(os.path.join(lp.landis_path,lp.scenario_file), 'w', encoding='utf-8') as file:
        file.writelines(filelist)

def replace_IC(lp):
    
    with open(os.path.join(lp.landis_path,lp.necn), 'r', encoding='utf-8') as file:
        filelist = file.readlines()
    
    matches = [match for match in filelist if "InitialCommunities" in match]
    matched_indexes = []
    print(matches)
    for j in range(0,len(matches)):
        i = 0
        length = len(filelist)
        while i < length:
            if matches[j] == filelist[i]:
                matched_indexes.append(i)
            i += 1
    print(matched_indexes)
    IC_txt = matched_indexes[0]
    IC_map = matched_indexes[1]
    
    filelist[IC_txt] = "InitialCommunities\t postfireIC-{}.txt\n".format(str(lp.year))
    filelist[IC_map] = "InitialCommunities\t output-community-{}.img\n".format(str(lp.year))
    
    with open(os.path.join(lp.landis_path,lp.necn), 'w', encoding='utf-8') as file:
        file.writelines(filelist)
 
# if __name__=="__main__":
#     main(sys.argv[1])









