# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 15:14:55 2022

@author: Niko Tutland
"""
import pandas as pd
import os

def toLandis(lp):
    ## User Inputs ##############
    # File paths
    landis_path = os.path.join(lp.OG_PATH, "1.LANDIS-MODEL/LANDIS_run")
    FM2VDM_path = os.path.join(lp.OG_PATH, "1.LANDIS-MODEL/FM2VDM")
    # VDM2FM_path = os.path.join(lp.OG_PATH, "1.LANDIS-MODEL/VDM2FM")
    ## End User Inputs ##########
    
    ## Read in post-fire treelist and assign necessary data to the remaining trees
    postfire = postfire_treelist(FM2VDM_path,landis_path,lp.year,lp.cycle)
    
    ## Group back into cohorts and sum the biomass
    postfire_cohorts = treelist_to_cohorts(postfire,lp.L2_res)
    
    ## Merge with uncropped treelist 
    community_input_file = merge_to_uncropped(postfire_cohorts, landis_path, lp.year, lp.cycle)
    
    ## Write new LANDIS community input file CSV
    community_input_file.to_csv(os.path.join(landis_path,"community-input-file-"+str(lp.year)+".csv"), index = False)
    
    ## Create new Initial Communities file (not necessary I think?)
    write_IC(community_input_file,landis_path,lp.year)
    
    ## End main function ##

def postfire_treelist(postfire_path,treelist_path,year,cycle):
    ## Read in treelist that has been altered by the simulated fire
    newtreelist = pd.read_csv(os.path.join(postfire_path,"AfterFireTrees."+str(year)+".csv")) # the filename will change depending on what we get from Adam's code
    # newtreelist = newtreelist.sample(frac=0.75) # let's pretend some trees burned (temporary)
    ## Read in the pre-QF treelist file that contains addtional information about the trees, crucially biomass, species, and mapcode
    treelist_alldata = pd.read_csv(os.path.join(treelist_path,"Treelist_alldata_"+str(cycle)+".csv"))
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

def merge_to_uncropped(postfire,path,year,cycle):
    prefire = pd.read_csv(os.path.join(path,"Treelist_alldata_"+str(cycle-1)+".csv"))
    prefire_mc = prefire["MapCode"].unique()
    postfire_mc = postfire["MapCode"].unique()
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
    prefire_uncropped = pd.read_csv(os.path.join(path,"community-input-file-"+str(year)+".csv"))
    burndomain_mc = postfire_all["MapCode"].unique()
    uncropped_mc = prefire_uncropped["MapCode"].unique()
    if burndomain_mc != uncropped_mc:
        outside_burndomain = prefire_uncropped[~prefire_uncropped["MapCode"].isin(burndomain_mc)]
        postfire_landis = pd.concat([outside_burndomain, postfire_all])
    else:
        postfire_landis = postfire_all
    return postfire_landis

def write_IC(IC,path,year):
    with open(os.path.join(path,'postfireIC-'+str(year)+".txt"), 'w') as file:
        file.write('LandisData "Initial Communities"\n')
        file.write("\n")
        for i in IC["MapCode"].unique():
            file.write("MapCode {}\n".format(i))
            IC_mc = IC[IC["MapCode"]==i]
            for j in IC_mc["SpeciesName"].unique():
                file.write("{} ".format(j))
                IC_mc_sp = IC_mc[IC_mc["SpeciesName"]==j].reset_index()
                for k in range(0,IC_mc_sp.shape[0]):
                    file.write("{} ({}) ".format(IC_mc_sp.loc[k,"CohortAge"], IC_mc_sp.loc[k,"CohortBiomass"]))
                file.write("\n")
            file.write("\n")
        file.write("\n")

# if __name__=="__main__":
#     toLandis()

