# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:14:55 2022

@author: Niko Tutland
"""
import pandas as pd
import os

def toLandis():
    ## User Inputs ##############
    # Run info
    aoi = "Klamath" #study area name
    year = 30 #year of landis run
    trt = "BAUnoharvest" #name of treatment. If no treatment, use None
    L2_res = 150 # LANDIS resolution in meters (100 = 100mx100m = 1ha)
    # File paths
    out_path = os.path.abspath("C://Users/FireScience/Documents/2022_Projects/landis_quicfire/Landis_to_Treelist") # where is the run path located? This should match the out_path specified in LANDIS_to_Treelist
    landis_path = os.path.abspath("C://Users/FireScience/Documents/2022_Projects/landis_quicfire/Klamath_BAU_Clipped") #where do you want to save the new LANDIS input file?
    ## End User Inputs ##########
    
    ## Name the run and run path
    if trt is None:
        run = "-".join([aoi,str(year)])
    else:
        run = "-".join([aoi,trt,str(year)])
    run_path = os.path.join(out_path,run)
    
    ## Read in post-fire treelist and assign necessary data to the remaining trees
    postfire = postfire_treelist(run_path,run)
    
    ## Group back into cohorts and sum the biomass
    community_input_file = treelist_to_cohorts(postfire,L2_res)
    
    ## Write LANDIS community input file
    community_input_file.to_csv(os.path.join(landis_path,"community-input-file"+str(year)+".csv"), index = False)
    
    ## End main function ##

def postfire_treelist(path,run):
    ## Read in treelist that has been altered by the simulated fire
    newtreelist = pd.read_csv(os.path.join(path,"Treelist_"+run+".csv")) # the filename will change depending on what we get from Adam's code
    newtreelist = newtreelist.sample(frac=0.75) # let's pretend some trees burned (temporary)
    ## Read in the pre-QF treelist file that contains addtional information about the trees, crucially biomass, species, and mapcode
    treelist_alldata = pd.read_csv(os.path.join(path,"Treelist_alldata_"+run+".csv"))
    newtreelist.columns = ["SPID","X","Y","HT_m","HTLC_m","CD_m","HTMCD_m","CBD","MOIST","SS"] # assign column names
    newtreelist = newtreelist[["SPID","X","Y"]] # use only columns needed to match to old treelist
    ## Assign attributes from pre-QF treelist to the remaining post-QF trees
    newtreelist_alldata = newtreelist.merge(treelist_alldata, how = "left", on = ["SPID","X","Y"])
    return newtreelist_alldata

def treelist_to_cohorts(x,L2_res):
    community_input_file = x.groupby(["MapCode","SPECIES_SYMBOL","AGE"]).sum("AGB_g")
    community_input_file = community_input_file.assign(CohortBiomass = lambda x: x.AGB_g/L2_res**2)
    community_input_file = community_input_file[["MapCode","SPECIES_SYMBOL","AGE","CohortBiomass"]]
    community_input_file = community_input_file.rename({"SPECIES_SYMBOL":"SpeciesName","AGE":"CohortAge"})
    return community_input_file

if __name__=="__main__":
    toLandis()
