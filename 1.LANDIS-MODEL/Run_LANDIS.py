# -*- coding: utf-8 -*-
"""
Created on Wed Nov  9 12:12:26 2022

@author: Niko Tutland
"""

"""
This script should control some inputs to the LANDIS runs, then do the initial model "spinup". 
After that, the LANDIS run should be cropped to the extent of the fire modeling domain,
(which for now is 400x400m) and buffered to help the accuracy of the fire behavior. The user
will have to somehow specify where in the LANDIS domain they want the fire domain to be. 
This script may also be where I put the functions to "continue" the LANDIS simulation after each
fire simulation (for looping forest change and fire behavior).

The run parameters that are controlled by the DRM_framework_coupling script are:
    nyears=10      # number of years for spinup and transient runs
    ncycyear=5    # number of cyclical year run
    ncycle=20      # number of loops

nyears will alter the Duration of the LANDIS scenario input file (line 3 for the scenario files I have)
for the FIRST run. This what we are calling the "spinup". I believe the minimum input for this would be
1, since we need at lease one year of simulation to produce a community-input-file used to generate a
treelist.

ncycyear will alter the same Duration input in the LANDIS scenario file for all subsequent runs.
I will have to make sure this is how the Community Biomass Output extension is supposed to work,
or if there is already a functionality within LANDIS for using community-input-files to continue
a run from where it left off. That would be useful since we wouldn't actually need the model to
do another climate spinup. If that exists, then ncycyear would change some other input file.

ncycle would not control LANDIS, it would simply tell this script when to stop looping LANDIS runs.

**NOTE: Since LANDIS runs, especially with a large domain, will take a lot of time to run, so if
possible we should try to run it in parallel

"""

import os
import sys
import subprocess
import glob
import re
import numpy as np
import pandas as pd
import Treelist_to_LANDIS as Treelist
import LANDIS_to_Treelist as Landis

def main():
    """
    1. Build LandisParams from DRM inputs and LANDIS files
    2. Run LANDIS "spinup"
    3. Landis_to_Treelist
    4. Fire modeling steps
    5. Treelist_to_Landis
    6. Alter LANDIS files
    7. Loop steps 3-6
    """
    # OG_PATH = os.getcwd()
    OG_PATH = os.path.abspath("C://Users/FireScience/Documents/2022_Projects/drm-fates-1")
    landis_path = os.path.join(OG_PATH,"1.LANDIS-MODEL")
    run_folder = "Klamath_BAU_Clipped"
    scenario_file = "Scenario1.txt"
    necn_file = "NECN_Succession.txt"
    batch_file = "RunIt.BAT"
    
    nyears=10      # number of years for spinup and transient runs
    ncycyear=5     # number of cyclical year run
    ncycle=20      # number of loops
    cycle = 1      # current iteration (will be looped through range(0,ncycle))
    
    spinup = True
    L2_params = LandisParams(nyears, ncycyear, ncycle, scenario_file, necn_file, batch_file)
    
    replace_duration(spinup,landis_path,L2_params)
    
    batch_cmd = os.path.join(landis_path,run_folder,batch_file)
    try:
        subprocess.run([batch_cmd])
    except FileNotFoundError as exc:
        print(f"Batchfile not found.\n{exc}")
    except subprocess.CalledProcessError as exc:
        print(f"LANDIS run failed with return code {exc.returncode}\n{exc}")
        
    
    Landis.toTreelist() # runs script to create a treelist from a landis run
    Treelist.toLandis() # runs script to create a landis input file from a treelist
    
    replace_IC(landis_path,L2_params)

class LandisParams:
    """
    Class containing parameters for controlling the LANDIS runs/loops
    """
    def __init__(self, nyears, ncycyear, ncycle, cycle, scenario_file, necn_file, batch_file):    
        self.nyears = int(nyears)            #number of years for spinup and transient runs
        self.ncycyear = int(ncycyear)        #number of cyclical year run
        self.ncycle = int(ncycle)            #number of loops
        self.cycle = int(cycle)              #current iteration
        self.scenario = str(scenario_file)   #name of LANDIS scenario input file
        self.necn = str(necn_file)           #name of NECN input file
        self.batchfile = str(batch_file)     #name of LANDIS batchfile 

def get_filenames():
    os.chdir("LANDIS_run")
    # Find batchfile and make sure it's the only one
    batchfile_list = glob.glob("*.bat") #case insensitive for windows, which should be what we want
    if len(batchfile_list) == 0:
        print("No batch file found.")
        return
    if len(batchfile_list) > 1:
        print("Multiple batch files in LANDIS directory. Please make sure there is only one.")
        return
    batch_file = batchfile_list[0]
    # Identify scenario file name from batchfile
    scenario_file = read_filename(batch_file,"call")
    # Identify NECN file name from scenario file
    necn_file = read_filename(scenario_file, '"NECN Succession"')
    # Identify Species file from scenario file
    species_file = read_filename(scenario_file, "Species")

def read_filename(file,line_id):
    with open(file, 'r', encoding='utf-8') as file:
        filelist = file.readlines()
    lines = [match for match in filelist if line_id in match]
    line = lines[0]
    file_name = re.split(" |\t", line)[-1].strip()
    return file_name
    
def replace_duration(spinup,path,lp):
    if spinup == True:
        runlen = lp.nyears
    else:
        runlen = lp.ncycyear
    
    with open(os.path.join(path,lp.scenario), 'r', encoding='utf-8') as file:
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
  
    with open(os.path.join(path,lp.scenario), 'w', encoding='utf-8') as file:
        file.writelines(filelist)

def replace_IC(path,lp):
    year = lp.nyears + lp.ncycyear*(lp.cycle-1)
    
    with open(os.path.join(path,lp.necn), 'r', encoding='utf-8') as file:
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
    
    filelist[IC_txt] = "InitialCommunities\t postfireIC-{}.txt\n".format(str(year))
    filelist[IC_map] = "InitialCommunities\t output-community-{}.img\n".format(str(year))
    
    with open(os.path.join(path,lp.necn), 'w', encoding='utf-8') as file:
        file.writelines(filelist)

# if __name__=="__main__":
#     main()

os.chdir("1.LANDIS-MODEL")








