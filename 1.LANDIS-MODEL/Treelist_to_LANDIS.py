# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:14:55 2022

@author: Niko Tutland
"""
import pandas as pd
import numpy as np
import rasterio as rio
import os
import re
import math

def toLandis(lp):
    print("Building LANDIS inputs from postfire treelist...")
    
    # File paths
    FM2VDM_path = os.path.join(lp.OG_PATH, "1.LANDIS-MODEL/FM2VDM")
    
    ## Read in post-fire treelist and assign necessary data to the remaining trees
    postfire = postfire_treelist(FM2VDM_path,lp.landis_path,lp.cycle)
    
    ## Group back into cohorts and sum the biomass
    spec_rename = dict(zip(lp.fia_spec,lp.landis_spec))
    postfire_cohorts = treelist_to_cohorts(postfire,lp.L2_res,spec_rename)
    
    ## Merge with uncropped treelist
    community_input_file = merge_cohorts(postfire_cohorts, lp.landis_path, 
                                         "community-input-file-"+str(lp.year_prev)+".csv", lp.cycle,
                                         "output-community-"+str(lp.year_prev)+".img")
    
    ## Replace fuels
    if lp.cycle == 1:
        IC_map = lp.IC_map
    else:
        IC_map = "output-community-"+str(lp.year_prev)+".tif"
    replace_fuels(lp.OG_PATH, lp.landis_path, lp.cycle, lp.coarseroots_map, IC_map, lp.L2_res, lp.year_prev, lp.crop_domain)
    
    ## Write new LANDIS community input file CSV
    # community_input_file.to_csv(os.path.join(lp.landis_path,"community-input-file-"+str(lp.year_prev)+".csv"), index = False)
    
    ## Create new Initial Communities file (not necessary I think?)
    write_IC(community_input_file,lp.landis_path,lp.cycle)
    
    print("Success!")
    ## End main function ##

def postfire_treelist(postfire_path,treelist_path,cycle):
    ## Read in treelist that has been altered by the simulated fire
    newtreelist = pd.read_csv(os.path.join(postfire_path,"AfterFireTrees."+str(cycle)+".txt"), sep=" ") 
    # newtreelist = newtreelist.sample(frac=0.75) # let's pretend some trees burned (temporary)
    ## Read in the pre-QF treelist file that contains addtional information about the trees, crucially biomass, species, and mapcode
    treelist_alldata = pd.read_csv(os.path.join(treelist_path,"Treelist_alldata_cycle"+str(cycle-1)+".csv"))
    newtreelist.columns = ["SPID","X","Y","HT_m","HTLC_m","CD_m","HTMCD_m","CBD","MOIST","SS"] # assign column names
    newtreelist = newtreelist[["SPID","X","Y"]] # use only columns needed to match to old treelist
    ## Assign attributes from pre-QF treelist to the remaining post-QF trees
    newtreelist_alldata = newtreelist.merge(treelist_alldata, how = "left", on = ["SPID","X","Y"])
    return newtreelist_alldata

def treelist_to_cohorts(x,L2_res,spec_rename):
    community_input_file = x.groupby(["MapCode","SPECIES_SYMBOL","AGE"], as_index=False).sum("AGB_g")
    community_input_file = community_input_file.assign(CohortBiomass = lambda x: x.AGB_g/L2_res**2)
    community_input_file = community_input_file[["MapCode","SPECIES_SYMBOL","AGE","CohortBiomass"]]
    community_input_file = community_input_file.rename({"SPECIES_SYMBOL":"SpeciesName","AGE":"CohortAge"}, axis = "columns")
    community_input_file = community_input_file.replace({"SpeciesName" : spec_rename})
    return community_input_file

def merge_cohorts(postfire,path,CIF_in,cycle,OC_in):
    prefire = pd.read_csv(os.path.join(path,"Treelist_alldata_cycle"+str(cycle-1)+".csv"))
    prefire_mc = prefire["MapCode"].unique()
    # print("prefire_mc:", prefire_mc)
    postfire_mc = postfire["MapCode"].unique()
    # print("postfire_mc:", postfire_mc)
    ## For any mapcodes with no fuels after fire, populate with zeros/None
    missing_mc = list(set(prefire_mc).difference(postfire_mc))
    if len(missing_mc) != 0:
        postfire_missing = pd.DataFrame({"MapCode":missing_mc,
                                         "SpeciesName":[None]*len(missing_mc),
                                         "CohortAge":[0]*len(missing_mc),
                                         "CohortBiomass":[0]*len(missing_mc)})
        postfire_all = pd.concat([postfire,postfire_missing])
    else:
        postfire_all = postfire
    ## Replace burn domain in landis run with updated fuels
    prefire_uncropped = pd.read_csv(os.path.join(path,CIF_in))
    burndomain_mc = postfire_all["MapCode"].unique()
    # print("burndomain_mc:", burndomain_mc)
    uncropped_mc = prefire_uncropped["MapCode"].unique()
    # print("prefire_uncropped:", uncropped_mc)
    if np.array_equal(burndomain_mc,uncropped_mc):
        postfire_landis = postfire_all
    else:
        outside_burndomain = prefire_uncropped[~prefire_uncropped["MapCode"].isin(burndomain_mc)]
        postfire_landis = pd.concat([outside_burndomain, postfire_all])
    postfire_landis = postfire_landis[["MapCode","SpeciesName","CohortAge","CohortBiomass"]]
    pf_landis_mc = postfire_landis["MapCode"].unique()
    with rio.open(os.path.join(path,OC_in)) as OC:
        OC_mc = OC.read(1).flatten()
    missing_mc = list(set(OC_mc).difference(pf_landis_mc))
    if len(missing_mc) != 0:
        postfire_missing = pd.DataFrame({"MapCode":missing_mc,
                                         "SpeciesName":[None]*len(missing_mc),
                                         "CohortAge":[0]*len(missing_mc),
                                         "CohortBiomass":[0]*len(missing_mc)})
        postfire_landis_all = pd.concat([postfire_landis,postfire_missing])
    else:
        postfire_landis_all = postfire_landis
    return postfire_landis_all

def write_IC(IC,path,cycle):
    with open(os.path.join(path,'postfireIC_cycle'+str(cycle)+".txt"), 'w') as file:
        file.write('LandisData "Initial Communities"\n')
        file.write("\n")
        for i in IC["MapCode"].unique():
            file.write("MapCode {}\n".format(i))
            IC_mc = IC[IC["MapCode"]==i]
            if len(IC_mc.index) == 0:
                file.write("\n")
            else:
                for j in IC_mc["SpeciesName"].unique():
                    if j==None:
                        file.write("\n")
                    else:
                        file.write("{} ".format(j))
                        IC_mc_sp = IC_mc[IC_mc["SpeciesName"]==j].reset_index()
                        for k in range(0,IC_mc_sp.shape[0]):
                            file.write("{} ({}) ".format(int(IC_mc_sp.loc[k,"CohortAge"]), int(math.ceil(IC_mc_sp.loc[k,"CohortBiomass"]))))
                        file.write("\n")
            file.write("\n")
        file.write("\n")

def replace_fuels(OG_PATH, landis_path, cycle, coarseroots_map, IC_map, L2_res, year_prev, crop_domain):
    litter_arr = np.loadtxt(os.path.join(OG_PATH,"1.LANDIS-MODEL","FM2VDM","AfterFireLitter."+str(cycle)+".txt"))
    litter_arr = litter_arr*1000
    ## Write georeferenced raster of postfire litter
    with rio.open(os.path.join(landis_path,"output-community-cycle"+str(cycle-1)+"_cropped.tif"), "r+") as IC:
        with rio.open(os.path.join(landis_path, "Postfire_litter_cycle"+str(cycle)+".tif"), 
                  mode="w",
                  height=IC.height,
                  width=IC.width,
                  count=1,
                  dtype=litter_arr.dtype,
                  crs="EPSG:5070",
                  transform=IC.transform) as pfl:
                pfl.write(litter_arr,1)
    ## Replace values of deadwood raster with the postfire litter values
    print("cycle =", cycle, "year_prev =", year_prev)
    with rio.open(os.path.join(landis_path,"NECN","SurfaceLitterBiomass-"+str(year_prev)+".img"), "r+") as slb:
        slb_arr = slb.read(1)
        with rio.open(os.path.join(landis_path,"NECN","ConiferNeedleBiomass-"+str(year_prev)+".img"), "r+") as cnb:
            cnb_arr = cnb.read(1)
            landis_arr = (slb_arr + cnb_arr)
            with rio.open(os.path.join(landis_path, "Postfire_litter_cycle"+str(cycle)+".tif"), "r+") as pfl:
                litter_arr = pfl.read(1)
                if crop_domain:
                    with rio.open(os.path.join(landis_path,IC_map), "r+") as IC:
                        print("IC_map:", IC_map) #this might be the problem
                        x_start = int((pfl.transform[2]-IC.transform[2])/L2_res)
                        y_start = int((IC.transform[5]-pfl.transform[5])/L2_res)
                        x_end = int(x_start + pfl.shape[1])
                        y_end = int(y_start + pfl.shape[0])
                        print("x_start:",x_start,"y_start:",y_start,"x_end:",x_end,"y_end:",y_end)
                        postfire_finefuel = landis_arr.copy()
                        print("litter_arr shape:", litter_arr.shape)
                        print("start/end shape:", postfire_finefuel[y_start:y_end,x_start:x_end].shape)
                        postfire_finefuel[y_start:y_end,x_start:x_end] = litter_arr
                        with rio.open(os.path.join(landis_path,"SurfaceFuels_cycle"+str(cycle)+".tif"),
                                      mode="w",
                                      height=IC.height,
                                      width=IC.width,
                                      count=1,
                                      dtype=postfire_finefuel.dtype,
                                      crs="EPSG:5070",
                                      transform=IC.transform) as out:
                            out.write(postfire_finefuel,1)
                else:
                    with rio.open(os.path.join(landis_path,"SurfaceFuels_cycle"+str(cycle)+".tif"),
                                  mode="w",
                                  height=pfl.height,
                                  width=pfl.width,
                                  count=1,
                                  dtype=litter_arr.dtype,
                                  crs="EPSG:5070",
                                  transform=pfl.transform) as out:
                        out.write(litter_arr,1)
    if cycle==1:
        ## Rename the original coarseroots raster so we can overwrite
        coarseroots_map_og = re.split("\.",coarseroots_map)[0] + "_original.tif"
        os.rename(os.path.join(landis_path,coarseroots_map), 
                  os.path.join(landis_path,coarseroots_map_og))
        ## Replace coarse roots values with zeros
        with rio.open(os.path.join(landis_path,coarseroots_map_og), mode="r") as crm:
            coarseroots_arr = crm.read(1)
            coarseroots_arr[:,:] = 0
            with rio.open(os.path.join(landis_path,coarseroots_map), 
                          mode="w",
                          height=crm.height,
                          width=crm.width,
                          count=1,
                          dtype=coarseroots_arr.dtype,
                          crs="EPSG:5070",
                          transform=crm.transform) as pfr:
                pfr.write(coarseroots_arr,1)
    

# if __name__=="__main__":
#     toLandis()

