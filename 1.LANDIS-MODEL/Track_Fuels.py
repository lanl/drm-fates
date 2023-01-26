# -*- coding: utf-8 -*-
"""
Created on Thu Jan 12 09:26:18 2023

@author: Niko Tutland
"""

"""
Writing a version of the "TreeOfLife" function from LLM_FT_utils.py that will work with LANDIS
Hopefully this will be more universal and can replace TreeOfLife
"""

import os
import numpy as np
import pandas as pd

class FuelsTracker:
    def __init__(self, VDM_folder, lp):
        os.chdir("..")
        OG_PATH = os.getcwd()
        os.chdir(VDM_folder)
        self.scorch_path = os.path.join(OG_PATH,"8.CROWN-SCORCH")
        
        self.tree_tracker = "TreeTracker.txt"
        self.percent_fuel = "PercentFuelChange.txt"
        self.treelist = "treelist_VDM.dat"
        self.grass_fuel = "VDM_litter_WG.dat"
        self.litter_fuel = "VDM_litter_trees.dat"
        
        self.VDM_folder = VDM_folder
        self.out_path = os.path.join(OG_PATH, VDM_folder, "FM2VDM")
        os.makedirs(self.out_path, exist_ok=True)
        
        self.trees_out = "AfterFireTrees.txt"
        self.grass_out = "AfterFireWG.txt"
        self.litter_out = "AfterFireLitter.txt"
        
        file_list = [
            self.tree_tracker,
            self.percent_fuel,
            self.treelist,
            self.grass_fuel,
            self.litter_fuel
            ]
        
        self.filelist = []
        for i in file_list:
            if os.path.exists(os.path.join(self.scorch_path,i)):
                self.filelist.append(i)
        
        # Decide threshold for tree death by species
        self.threshold_dict = get_thresholds(lp)
                

def postfire_fuels(ft,tp,lp):
    ## ft = FuelsTracker class
    ## tp = Treelist_params class
    ## lp = Landis_params class
    
    ## Set up empty arrays
        # tree and surface, 1-D and n-D
    # tree_arr = np.zeros((tp.nx,tp.ny,tp.nz))
    # tree_arr_1d = tree_arr.flatten()
    # surf_arr = np.zeros((tp.nx,tp.ny))
    # surf_arr_1d = surf_arr.flatten()
    
    # DEFINE DOMAIN
    Nx = tp.nx
    Ny = tp.ny
    Nz = tp.nz
    s = (Nx,Ny)
    array1d = Nx*Ny*Nz
    percentFuelChang1d = np.zeros(array1d)
    postfire_fuel = np.zeros(s)
    
    ## Import percent change in fuels for each cell
    with open(os.path.join(ft.scorch_path,ft.percent_fuel), 'r') as pfc:
        count = 0
        for line in pfc.readlines():
            percentFuelChang1d[count] = float(line)
            count = count + 1
            
    # ## Import percent change in fuels for each cell
    # pfc_2d = np.loadtxt(os.path.join(ft.scorch_path, ft.percent_fuel))
    # pfc_3d = pfc_2d.reshape(
    #     pfc_2d.shape[0],
    #     pfc_2d.shape[1] // tp.nz,
    #     tp.nz)
    # # pfc_3d.shape is in the format (y,x,z), so surface fuels are pfc_3d[:,:,0]
    # surface_change = pfc_3d[:,:0]
    # canopy_change = pfc_3d.copy()

    
    ## Calculate change in surface fuels
    for i in [ft.grass_fuel,ft.litter_fuel]:
        if i in ft.filelist:
            with open(os.path.join(ft.scorch_path,i), 'r') as pff:
                prefire_fuel = pff.read().split()
            for locy in range(Ny):
                for locx in range(Nx):
                    gloc = locx + (locy * Nx)
                    planarloc = (Nx*(Ny-1)-(locy*Nx)) + locx
                    postfire_fuel[locx,locy] = float(prefire_fuel[gloc]) * percentFuelChang1d[planarloc]
            if i==ft.grass_fuel:
                np.savetxt(os.path.join(ft.out_path,ft.grass_out), postfire_fuel, fmt='%10.2f')
            elif i==ft.litter_fuel:
                np.savetxt(os.path.join(ft.out_path,ft.litter_out), postfire_fuel, fmt='%10.2f')
    
    ## Which trees die?
    with open(os.path.join(ft.out_path,ft.trees_out), 'w') as aft:
        with open(os.path.join(ft.scorch_path,ft.tree_tracker), 'r') as tt:
            with open(os.path.join(ft.scorch_path,ft.treelist), 'r') as tl:
                treezip = zip(tt,tl)
                livetrees = 0
                deadtrees = 0
                for tt, tl in treezip:
                    line_tt = tt.split()
                    line_tl = tl.split()
                    cellnum = int(line_tt[1]) #number of cells with fuel from that tree
                    totfuel = 0.0
                    sppflag = int(line_tl[0]) #tree species identifier
                    threshold = ft.threshold_dict[sppflag]
                    for cell in range(cellnum):
                        cell_index = int(line_tt[2+cell]) #first two items in list are tree id and cellnum, so cound over starting on third item
                        fuel_conc = float(line_tt[2+cell+cellnum]) #next items correspond to inital fuel density in each cell
                        newfuel = percentFuelChang1d[cell_index] * fuel_conc
                        totfuel = totfuel + newfuel # sum the total amount of remaininf fuel for that tree
                    if totfuel > threshold:
                        aft.write(tl)
                        livetrees += 1
                    else:
                        deadtrees += 1
    LiveDeadList=[]
    LiveDeadList.append(livetrees)
    LiveDeadList.append(deadtrees)
    
    return LiveDeadList
                  

def get_thresholds(lp):
    # mortality threshold dict from LANDIS_options
    threshold_dict = dict(zip(lp.fia_spec,lp.mortality_thresholds))
    # Species IDs to species names from treelist
    treelist = pd.read_csv(os.path.join(lp.landis_path,"Treelist_alldata_cycle"+str(lp.cycle)+".csv"))
    spids = list(range(1,len(treelist["SPID"].unique())+1,1))
    sps = list(treelist["SPECIES_SYMBOL"].unique())
    spid_dict = dict(zip(sps,spids))
    # Match species ID to threshold by species name
    dicts = [threshold_dict,spid_dict]
    merge = merge_dictionary_list(dicts)
    mortality_dict = {}
    for i in merge.values():
        mortality_dict[i[0]] = i[1]

    return mortality_dict


def merge_dictionary_list(dict_list):
  return {
    k: [d.get(k) for d in dict_list if k in d] # explanation A
    for k in set().union(*dict_list) # explanation B
  }


