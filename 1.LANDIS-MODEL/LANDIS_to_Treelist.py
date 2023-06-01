# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 14:03:17 2022

@author: Niko Tutland
"""

import numpy as np
from scipy import interpolate
import pandas as pd
import pickle
from random import random
import os
import rasterio as rio
from rasterio.warp import calculate_default_transform, reproject, Resampling
import time

def toTreelist(lp):
    start_time = time.time()
    
    ## Filepaths
    OG_PATH = lp.OG_PATH
    
    ## Create species renaming dict
    spec_rename = dict(zip(lp.landis_spec,lp.fia_spec))
    proc_path = os.path.join(OG_PATH, "1.LANDIS-MODEL/FIA_proc")
    
    if os.path.exists(proc_path):
        print("Importing processed FIA data...")
        FIA_cohorts = pd.read_csv(os.path.join(proc_path,"FIA_cohorts.csv"))
        FIA_all = pd.read_csv(os.path.join(proc_path,"FIA_all.csv"))
    else:
        ## Read in and process FIA data
        print("Processing FIA data...")
        raw_path = os.path.join(OG_PATH, "9.FIA/FIA_raw")
        RF_path = os.path.join(OG_PATH, "9.FIA/RF_models")
        
        FIA_all = process_fia(path=raw_path, states=lp.states)
        
        ## Only include tree species in study area
        FIA_all = aoi_spec(FIA_all,raw_path,lp.fia_spec)
        
        ## Use random forest models to predict age
        print("Predicting tree age...")
        FIA_all = age_calc(FIA_all,RF_path,lp.age_bin)
    
        ## Create cohorts from FIA
        FIA_cohorts = fia_cohorts(FIA_all)
        
        ## Write to file
        os.makedirs(proc_path, exist_ok=True)
        FIA_cohorts.to_csv(os.path.join(proc_path,"FIA_cohorts.csv"), index = False)
        FIA_all.to_csv(os.path.join(proc_path,"FIA_all.csv"), index = False)
        
    
    ## Read in and process LANDIS outputs
    LANDIS_cohorts = process_landis(lp.landis_path,"community-input-file-cycle"+str(lp.cycle)+"_cropped.csv",lp.age_bin,lp.landis_spec,spec_rename)
    
    ## Match LANDIS cohorts to most similar FIA cohort
    print("Matching to LANDIS cohorts...")
    Plotlist = match(LANDIS_cohorts, FIA_cohorts)
    
    ## Pull individual trees from matched cohorts
    Treelist = get_trees(Plotlist,FIA_all)
    
    print("Calculating tree attributes...")
    ## Calculate fields for QUIC-Fire
    Treelist = qf_calcs(Treelist,lp.bulk_density,lp.cl_factor,lp.moisture,lp.sizescale,lp.aoi_elev,lp.region_flag,lp.IC_cropped,lp.landis_path)
    
    ## Replicate trees to fill LANDIS grid cell
    Treelist = replicate_trees(Treelist, lp.L2_res)
    
    ## Assign trees a random location within a LANDIS cell
    print("Assigning tree locations...")
    Treelist = tree_locations(Treelist, lp.landis_path, "output-community-cycle"+str(lp.cycle)+"_cropped.tif", lp.year, lp.L2_res)
    
    ## Clean up for QUIC-Fire
    Treelist_alldata = Treelist.copy()
    Treelist = Treelist[["SPID","X","Y","HT_m","HTLC_m","CD_m","HTMCD_m","CBD","MOIST","SS"]]
    
    ## Import, interpolate, and export LANDIS fuels
    print("Interpolating surface fuels...")
    fuels = get_fuels("SurfaceLitterBiomass-cycle"+str(lp.cycle)+"_cropped.tif", 
                      "ConiferNeedleBiomass-cycle"+str(lp.cycle)+"_cropped.tif", 
                      lp.landis_path, lp.L2_res)
    
    ## Write intermediate files
    # print("Writing files...")
    # file_list = [FIA_all, FIA_cohorts, LANDIS_cohorts, Plotlist, Treelist_alldata]
    # file_names = file_list.copy()
    # i = 0
    # for file in file_list :
    #     obj = [k for k,v in locals().items() if v is file][0]
    #     name = obj+"_"+run+".csv"
    #     file_names[i] = name
    #     i = i+1
    # file_dict = dict(zip(file_names,file_list))
    # write_files(file_dict, path=run_path)
    
    # Write Treelist_alldata for Treelist_to_LANDIS
    Treelist_alldata.to_csv(os.path.join(lp.landis_path,"Treelist_alldata_cycle"+str(lp.cycle)+".csv"), index = False)
    
    # Write final treelist and surface fuels
    os.makedirs("VDM2FM", exist_ok=True)
    treelist_name = "treelist_VDM.dat"
    Treelist.to_csv(path_or_buf = os.path.join("VDM2FM",treelist_name), header=False, index=False, sep=" ")
    fuels_name = "VDM_litter_trees.dat"
    np.savetxt(os.path.join("VDM2FM",fuels_name),fuels,delimiter=" ")
    
    #Calculate domain parameterss to build treelist text file (used in trees script)
    dom_params = get_params(Treelist,lp.landis_path,"output-community-cycle"+str(lp.cycle)+"_cropped.tif",
                            lp.L2_res,lp.QF_res,treelist_name)
    print_fuellist(dom_params, os.path.join(OG_PATH, "5.TREES-QUICFIRE"))
    
    print("Treelist created successfully")
    print(" %s seconds elapsed" % round(time.time() - start_time))
    return dom_params

##### End of main function #####

class Treelist_params:
    """
    Class contains domain parameters
    """
    def __init__(self, dx, dy, dz, nx, ny, nz, csv_name, num_trees, num_spp):    
        self.dx = int(dx)            #x dimensions [m]
        self.dy = int(dy)            #y dimensions [m]
        self.dz = int(dz)            #z dimensions [m]
        self.nx = int(nx)       #number of x cells
        self.ny = int(ny)       #number of y cells
        self.nz = int(nz)            #number of z cells
        self.csv_name = csv_name
        self.num_trees = num_trees
        self.num_spp = num_spp

def process_fia(path, states):
    lbac_to_gm2 = 0.1121
    FIA_all = pd.DataFrame()
    #i = "NC" #for testing loop
    for i in states:
        ## From Trees Table: Pull in attributes we'll need for matching and for building a treelist
        tree = pd.read_csv(os.path.join(path,i+"_TREE.csv"))
        tree.columns = tree.columns.str.lower()
        tree = tree[tree["statuscd"]==1]
        tree = tree[["cn","plt_cn","condid","dia","ht","cclcd","cr","drybio_ag","spcd","tpa_unadj","invyr","bhage","totage"]]
        tree = tree[tree["cn"].notna()]
        tree = tree[tree["dia"].notna()]
        tree = tree.assign(AG_Biomass_gm2 = lambda row: row.drybio_ag*row.tpa_unadj*lbac_to_gm2)
        tree["cr"] = tree["cr"].apply(lambda x: 0.5 if x == 0 else x) #replace any cr==0 with 0.5
        tree["state"] = i
        ## From Condition Table: make sure it's forested with live trees, then get site condition
        cond = pd.read_csv(os.path.join(path,i+"_COND.csv"))
        cond.columns = cond.columns.str.lower()
        cond = cond[(cond["cond_status_cd"]==1) & (cond["live_canopy_cvr_pct"] > 0)]
        cond = cond[["plt_cn","condid","physclcd","sicond","stdage","balive"]]
        ## Join cond and tree tables together
        all_i = tree.merge(cond, how = "left", on=["plt_cn","condid"])
        ## Remove zeros and NAs
        all_i.loc[all_i["sicond"]==0,"sicond"] = all_i.loc[all_i["sicond"]!=0,"sicond"].mean()
        all_i.loc[all_i["stdage"]==0,"stdage"] = all_i.loc[all_i["stdage"]!=0,"stdage"].mean()
        all_i.loc[all_i["sicond"].isna(),"sicond"] = all_i.loc[all_i["sicond"].notna(),"sicond"].mean()
        all_i.loc[all_i["stdage"].isna(),"stdage"] = all_i.loc[all_i["stdage"].notna(),"stdage"].mean()
        all_i = all_i[all_i["balive"].notna()] #live ba should not be na
        all_i = all_i[all_i["balive"] > 0] #or zero
        all_i = all_i[all_i["physclcd"].notna()] #physiognomic class should not be na
        all_i = all_i[all_i["physclcd"] > 0] #or zero
        ## Append each state
        FIA_all = pd.concat([FIA_all,all_i])
    return FIA_all

def aoi_spec(x,path,species):
    REF_SPECIES = pd.read_csv(os.path.join(path,"REF_SPECIES.csv"))
    REF_SPECIES = REF_SPECIES[["SPCD","SPECIES_SYMBOL","MAJOR_SPGRPCD"]]
    REF_SPECIES = REF_SPECIES.rename(columns = {"SPCD" : "spcd"})
    x = x.merge(REF_SPECIES, how = "left", on = "spcd")
    x = x[x["SPECIES_SYMBOL"].isin(species)]
    return x

def age_calc(x,path,age_bin):
    # split df to with/without age
    mask = (x["totage"]>0) | (x["bhage"]>0)
    withage = x[mask]
    noage = x[~mask]
    if len(withage.index > 0):
        # use totage or bhage+10 if age is populated
        withage["AGE"] = withage.apply(lambda x: x.totage if x.totage > 0 else x.bhage + 10, axis = 1)
    # now predict age for everything else
    y = noage[["cn","plt_cn","condid","spcd","drybio_ag","tpa_unadj","invyr","AG_Biomass_gm2","cr","state","SPECIES_SYMBOL","MAJOR_SPGRPCD"]]
    x_vars = noage[["dia","ht","cclcd","physclcd","sicond","stdage","balive","MAJOR_SPGRPCD"]]
    x_vars.loc[:,"cclcd"] = x_vars.loc[:,"cclcd"].astype("category") # make sure canopy class is a categorical variable
    x_vars.loc[:,"physclcd"] = x_vars.loc[:,"physclcd"].astype("category") # make sure physiognomic class is a categorical variable
    FIA_age = pd.DataFrame()
    for spgrp in x["MAJOR_SPGRPCD"].unique():
        RFobj = open(os.path.join(path,'RF_model_spgrp'+str(spgrp)+'.obj'), 'rb')
        RF = pickle.load(RFobj)
        RFobj.close()
        x_spgrp = x_vars[x_vars["MAJOR_SPGRPCD"]==spgrp]
        y_spgrp = y[y["MAJOR_SPGRPCD"]==spgrp]
        x_spgrp_vars = x_spgrp.drop(["MAJOR_SPGRPCD"], axis = "columns")
        x_spgrp_vars["AGE"] = RF.predict(x_spgrp_vars)
        spgrp_i = pd.concat([x_spgrp_vars,y_spgrp], axis = 1)
        FIA_age = pd.concat([FIA_age,spgrp_i])
    FIA_age = pd.concat([FIA_age,withage])
    FIA_age["AGE"] = FIA_age["AGE"].apply(lambda g: roundUp(g, to=age_bin))
    return FIA_age

def fia_cohorts(x):
    x = (x.groupby(["plt_cn","invyr","SPECIES_SYMBOL","spcd","AGE"], as_index = False)
        .sum("AG_Biomass_gm2"))
    x = x.rename(columns = {"plt_cn" : "PLT_CN", "invyr" : "INVYR", "AG_Biomass_gm2" : "AGB_FIA"})
    return x

def process_landis(path,IC_name,age_bin,species,spec_rename):
    x = pd.read_csv(os.path.join(path,IC_name))
    x = x[x["SpeciesName"].isin(species)]
    x["CohortAge"] = x["CohortAge"].apply(lambda g: roundUp(g, to=age_bin))
    x = x.rename(columns = {"SpeciesName" : "SPECIES_SYMBOL", "CohortAge" : "AGE", "CohortBiomass" : "AGB_LANDIS"})
    x = x.replace({"SPECIES_SYMBOL" : spec_rename})
    return x

def match(x,y):
    # Is there a way to do this in parallel?
    #x=LANDIS_cohorts.copy()
    #y=FIA_cohorts.copy()
    mapcodes = pd.unique(x["MapCode"])
    Plotlist = pd.DataFrame()
    #i = 3 #for testing loop
    for i in mapcodes:
        matched = x[x["MapCode"]==i]
        matched = matched.merge(y, how = "inner", on = ["SPECIES_SYMBOL","AGE"])
        matched = matched.assign(AGB_diff = lambda g: abs(g.AGB_LANDIS - g.AGB_FIA))
        min_loc = list(matched.groupby(["SPECIES_SYMBOL","AGE"])["AGB_diff"].idxmin())
        matched = matched.loc[min_loc]
        Plotlist = pd.concat([Plotlist,matched])  
    Plotlist = Plotlist[["PLT_CN","SPECIES_SYMBOL","AGE","MapCode"]]
    return Plotlist

def get_trees(x,y):
    matched_plots = pd.unique(x["PLT_CN"])
    y = y[["cn","plt_cn","dia","ht","cr","tpa_unadj","drybio_ag","balive","SPECIES_SYMBOL","spcd","AGE"]]
    y = y[y["plt_cn"].isin(matched_plots)]
    y = y.rename(columns = {"plt_cn":"PLT_CN","cn":"CN","ht":"HT","dia":"DBH","cr":"C_RATIO","tpa_unadj":"TPA_UNADJ","drybio_ag":"AGB","balive":"BA"})
    Treelist = y.merge(x, how = "inner", on = ["PLT_CN","SPECIES_SYMBOL","AGE"])
    return Treelist

def replicate_trees(x,res):
    ac_per_ha = 2.47105
    x = x.assign(TPH_UNADJ = lambda y: y.TPA_UNADJ*ac_per_ha) #convert to tph
    x = x.assign(tree_multiplier = lambda y: round(y.TPH_UNADJ*res/100)) # find the number to trees per landis cell and round to nearest integer
    x = pd.DataFrame(x.values.repeat(x['tree_multiplier'], axis=0), columns=x.columns) # repeat rows based on the value in tree_multiplier
    return x

# def tree_locations(df,path,file,year,res):
#     ## Import mapcode raster
#     MapCodes = raster_import(os.path.join(path,file))
#     df['X'] = 0.0 #setup columns
#     df['Y'] = 0.0 
#     for y in range(MapCodes.shape[0]):
#         for x in range(MapCodes.shape[1]):
#             mc = MapCodes[y,x]
#             #Calc loc of bottom left corner of MapCode pixel in qf grid
#             qf_x_min = x*res
#             qf_y_min = y*res
#             temp_df = df[df['MapCode']==mc] #make new df with mapcode, and x y of bottom corner of each pixel
#             for i,row in temp_df.iterrows():
#                 #Randomly place trees within larger pixel
#                 x_pixel_loc = random() * res
#                 y_pixel_loc = random() * res
#                 df.X.loc[i] = qf_x_min + x_pixel_loc
#                 df.Y.loc[i] = qf_y_min + y_pixel_loc
#     return df

def tree_locations(df,path,file,year,res):
    ## Import mapcode raster
    MapCodes = raster_import(os.path.join(path,file))
    MapCode = []
    X = []
    Y = []
    for y in range(MapCodes.shape[0]):
        for x in range(MapCodes.shape[1]):
            MapCode.append(MapCodes[y,x])
            X.append(x*res)
            Y.append(y*res)
    df_mc = pd.DataFrame({"MapCode":MapCode,"X_corner":X, "Y_corner":Y})
    df = df.merge(df_mc, how="left", on="MapCode")
    df["X"] = df["X_corner"].apply(lambda x: x + (random()*res))
    df["Y"] = df["Y_corner"].apply(lambda y: y + (random()*res))
    df = df.drop(["X_corner","Y_corner"], axis=1)
    return df

def reproject_raster(IC,path,dst_crs = 'EPSG:4326'):
    with rio.open(os.path.join(path, IC), 'r+') as IC_map :
        src_transform = IC_map.transform
        transform, width, height = calculate_default_transform(
        IC_map.crs, dst_crs, IC_map.width, IC_map.height, *IC_map.bounds)
        kwargs = IC_map.meta.copy()
        kwargs.update({
            'crs': dst_crs,
            'transform': transform,
            'width': width,
            'height': height
        })
        print("Source Transform:\n",src_transform,'\n')
        print("Destination Transform:\n", transform)
        with rio.open(os.path.join(path,"IC_trans.tif"), 'w', **kwargs) as IC_ll:
            for i in range(1, IC_ll.count + 1):
                reproject(
                    source=rio.band(IC_map, i),
                    destination=rio.band(IC_ll, i),
                    src_transform=IC_map.transform,
                    src_crs=IC_map.crs,
                    dst_transform=transform,
                    dst_crs=dst_crs,
                    resampling=Resampling.nearest)

def get_AOI_latlon(IC,path):
    reproject_raster(IC=IC, path=path)
    IC_ll = rio.open(os.path.join(path, "IC_trans.tif"), 'r+')
    band1 = IC_ll.read(1)
    height = band1.shape[0]
    width = band1.shape[1]
    cols, rows = np.meshgrid(np.arange(width), np.arange(height))
    xs, ys = rio.transform.xy(IC_ll.transform, rows, cols)
    lons= np.array(xs)
    lats = np.array(ys)
    IC_long, IC_lat = lons.mean(), lats.mean()
    return IC_long, IC_lat

def roundUp(x,to):
    x = to*(x//to + bool(x%to))
    return x

def raster_import(filepath):    
    # Open the file:
    with rio.open(filepath,"r+") as raster:
        rasterArray = raster.read(1)

    rasterArray = np.flip(rasterArray, axis = 0)
    
    return rasterArray

def qf_calcs(x,cbd,cl,moist,ss,AOI_elev,region_flag,IC,path):
    ft_per_m = 3.28084
    lb_per_g = 0.00220462
    cr_coeff = pd.read_csv("CROWN_WIDTH_TABLE.csv")
    if region_flag==1:
        cr_coeff = cr_coeff[cr_coeff["REGION"]=="SW"]
    elif region_flag==2:
        cr_coeff = cr_coeff[cr_coeff["REGION"]=="NW"]
    else:
        cr_coeff = cr_coeff[cr_coeff["REGION"]=="EC"]
    #^NOTE: It would be more sustainable to do this in a different way for the case when a species occurs outside the region
    #where it was assigned a crown width equation. This would avoid patchwork fixes like for grand fir in California.
    #Would be better to have the code search the species within a region first, and if it's not there, use equations/coefficients
    #from another region.
    x = x.merge(cr_coeff, how = "left", on = "spcd")
    x = x.assign(HT_m = lambda x: x.HT/ft_per_m) # convert ht to meters
    x = x.assign(HTLC_m = lambda x: x.HT_m*(1-(x.C_RATIO/100))) # calculate HTLC
    x = x.assign(AGB_g = lambda x: x.AGB/lb_per_g)
    #x = x.assign(HTLC_m = lambda x: x.HT_m*(1-(x.C_RATIO/100))-(6/ft_per_m)) 
    #^NOTE: Crown Ratio is "ocularly compacted", meaning it can't be used to calculate the lowest
    #part of the crown. There are species-specific coefficients to convert to uncompacted crown ratio,
    #but for now we'll just subract the median difference between uncompacted and compacted from the
    #literature: https://watermark.silverchair.com/sjaf0118.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAscwggLDBgkqhkiG9w0BBwagggK0MIICsAIBADCCAqkGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMK0txij4pSwx-7AuFAgEQgIICenlL-MH-XNIi0H12Ex_FF8lzkeFyNGgriiPFYvCe9zy9ZNpkVtzI_jCiYd1gbWI2NKINIId4MVBbLkvEnCEwCiXklaEDXUUM38Ig8XzUoGU8Hckyc8_RnqPrUXPjAIhPnQz7wnR1MSEQgZnz3zVN8TMCvn5q5REFieC0pY8gv-5SSriEr8T08LVC0Im8B-rK5uE9Ts_r-aO7fTtHI4RVFJlKcfapd4f_sOBLOpwLCF4i_8_Cn2ywQcGQqJmkedhHocj_7RiNrDYkbihqIbKhvDNgpXMlyFHA79BbwsiMkCGxBGUXmGz-mPvyuNr6P8ZCdVSRMj6QQLIcD5mh8-NPv-j4RIktASphM9zNn-YjGdzL2axyRGGzAiGm4Tb_F0bRNBTGbY_b1f2JVljztodAdPzZxXfQXpt5dBOPY5TUgAthR3Ij43AHzAypZ_eYvR5l-kQ0a1thIv0MQv-uaj2E5cafhxAIl0Uq602tblFA_ctWNLM1r-9CZ-Q_bFNkjIX8g5vKGk1v-CsiUvTnjfBSR4j8MOLtJd1kYb5WnCLa1BowHo0CZWiQ_NkTh8yURwP0uyspolrldsALH18iHiDYLkQ2RzB7d1hn22JenqFSeKnHIimVhDLxBInjrO7lsmoTuDpeGZCpwbdfGE3uQjQ-wWHx5HsQK3mqfWDrXsoMaTKCb9jD0pJx3UWzXGNn82hjTPqTLKGjqXMt6Od4W4dxO9oxiyBJ-GeP3ewV6PhFq4Q7yf_1WnKhyPM3NxLob2DOcGwETzJ-EDgQvJ8ZxANoo9Xkr5hS9RY44X3QIoP_OCTd79k593bY9ooKdSXMLSJAIrAFiQtEFRXOUF8
    AOI_long, AOI_lat = get_AOI_latlon(IC,path)
    x["CW"] = x.apply(lambda x: crown_width(x.spcd,x.REGION,x.EQUATION,x.DBH,
                                            x.a1,x.a2,x.a3,x.a4,x.a5,x.a6,x.d1,x.d2,x.s1,
                                            x.HT,x.MinD,x.C_RATIO,x.BA,
                                            AOI_elev,AOI_lat,AOI_long,
                                            x.HI_lat,x.HI_long,x.HI_elev), 
                      axis = 1)
    x = x.assign(CD_m = lambda x: x.CW/ft_per_m)
    x["CBD"] = cbd
    x["MOIST"] = moist
    x["SS"] = ss
    #Create unique species IDs
    x["SPID"] = x["SPECIES_SYMBOL"]
    spids = list(range(1,len(x["SPID"].unique())+1,1))
    sps = list(x["SPID"].unique())
    spid_dict = dict(zip(sps,spids))
    x = x.replace({"SPID" : spid_dict})
    x = x.assign(HTMCD_m = lambda x: x.HT_m - cl*(x.HT_m-x.HTLC_m))
    x = x[["CN","PLT_CN","SPECIES_SYMBOL","SPID","AGE","MapCode","HT_m","HTLC_m","CD_m","HTMCD_m","CBD","MOIST","SS","TPA_UNADJ","DBH","AGB_g"]]
    return x

def crown_width(SPCD,REGION,EQUATION,DBH,a1,a2,a3,a4,a5,a6,d1,d2,s1,HT,MinD,CR,BA,AOI_elev,AOI_lat,AOI_long,HI_lat,HI_long,HI_elev):
    HI = ((AOI_elev - HI_elev)/100)*1 + (AOI_lat - HI_lat)*4 + (HI_long-AOI_long)*1.25
    if REGION == "NW":
        if EQUATION == 1:
            CW = CW_NW_1(MinD,DBH,a1,a2,a3)
        elif EQUATION == 2:
            CW = CW_NW_2(MinD,DBH,CR,BA,HI,a1,a2,a3,a4,a5,a6)
        elif EQUATION == 3:
            CW = CW_NW_3(SPCD,MinD,DBH,HT,CR,BA,a1,a2,a3,a4,a5,a6)
        elif EQUATION == 4:
            CW = CW_NW_4(MinD,DBH,a1,a2)
        elif EQUATION == 5:
            CW = CW_NW_5(MinD,DBH,CR,HT,BA,AOI_elev,a1,a2,a3,a4,a5,a6)
        else:
            CW = CW_NW_6(MinD,DBH,a1,a2)
    elif REGION == "SW":
        if DBH < MinD:
            if SPCD == 17:
                CW = CW_NW_3(SPCD,MinD,DBH,HT,CR,BA,a1,a2,a3,a4,a5,a6)
            else:
                CW = CW_SW_4_5(DBH,HT,d1,d2,s1)
        else:
            if EQUATION == 1:
                CW = CW_SW_1(DBH,a1,a2)
            elif EQUATION == 2:
                CW = CW_SW_1(DBH,a1,a2)
            elif EQUATION == 3:
                CW = CW_SW_3(DBH, a1,a2,a3)
            else:
                CW = CW_NW_3(SPCD,MinD,DBH,HT,CR,BA,a1,a2,a3,a4,a5,a6)
    else:
        if EQUATION == 1:
            CW = CW_EC_1(MinD,DBH,CR,HI,a1,a2,a3,a4,a5)
        else:
            CW = CW_EC_2(MinD,DBH,a1,a2,a3)
    return CW

def CW_EC_1(MinD,DBH,CR,HI,a1,a2,a3,a4,a5):
    if DBH >= MinD:
        CW = a1 + (a2*DBH) + (a3*DBH**2) + (a4*CR) + (a5*HI)
    else:
        CW = (a1 + (a2*MinD) + (a3*MinD**2) + (a4*CR) + (a5*HI))*(DBH/MinD)
    return CW

def CW_EC_2(MinD,DBH,a1,a2,a3):
    if DBH >= MinD:
        CW = a1 + (a2*DBH**a3)
    else:
        CW = (a1 + (a2*MinD**a3))*(DBH/MinD)
    return CW

def CW_SW_1(DBH,a1,a2):
    CW = a1 + a2*DBH
    return CW

def CW_SW_2(DBH,a1,a2):
    CW = a1+DBH**a2
    return CW

def CW_SW_3(DBH,a1,a2,a3):
    CW = a1 + a2*DBH + a3*DBH**2
    return CW

def CW_SW_4_5(DBH,HT,d1,d2,s1):
    if HT < 4.5:
        CW = HT*s1
    else:
        CW = d1 + d2*DBH
    return CW

def CW_NW_1(MinD,DBH,a1,a2,a3):
    if DBH >= MinD:
        CW = a1 + (a2*DBH) + a3*DBH**2
    else:
        CW = (a1 + (a2*MinD) + (a3*MinD**2))+(DBH/MinD)
    return CW

def CW_NW_2(MinD,DBH,CR,BA,HI,a1,a2,a3,a4,a5,a6):
    if DBH >= MinD:
        CW = a1 + (a2*DBH) + (a3*MinD**2) + (a4*CR) + (a5*BA) + (a6*HI)
    else:
        CW = (a1 + (a2*MinD) + (a3*MinD**2) + (a4*CR) + (a5*BA) + (a6*HI))*(DBH/MinD)
    return CW

def CW_NW_3(SPCD,MinD,DBH,HT,CR,BA,a1,a2,a3,a4,a5,a6):
    if SPCD == 264:
        if HT < 5.0:
            CW = (0.8*HT*max(0.5,CR/100))*(1-(HT-5)*0.1)*a1*DBH**a2*HT**a3*((CR/100)*HT)**a4*(HT-5)*0.1
        elif HT > 15:
            CW = a1*(DBH**2)*(HT**3)*(((CR/100)*HT)**a4)
        else:
            CW = 0.8*HT*max(0.5,CR/100)
    else:
        if DBH >= MinD:
            CW = a1*np.exp(a2 + (a3*np.log((CR/100)*HT)) + (a4*np.log(DBH)) + (a5*np.log(HT)) + (a6*np.log(BA)))
        else:
            CW = (a1*np.exp(a2 + (a3*np.log((CR/100)*HT)) + (a4*np.log(MinD)) + (a5*np.log(HT)) + (a6*np.log(BA))))*(DBH/MinD)
    return CW

def CW_NW_4(MinD,DBH,a1,a2):
    if DBH >= MinD:
        CW = a1*DBH**a2
    else:
        CW = (a1*MinD**2)*(DBH/MinD)
    return CW

def CW_NW_5(MinD,DBH,CR,HT,BA,AOI_elev,a1,a2,a3,a4,a5,a6):
    if DBH >= MinD:
        CW = (a1*1)*DBH**a2*HT**a3*((CR/100)*HT)**a4*(BA+1)**a5*np.exp(AOI_elev)**a6
    else:
        CW = ((a1*1)*MinD**a2*HT**a3*((CR/100)*HT)**a4*(BA+1)**a5*np.exp(AOI_elev)**a6)*(DBH/MinD)
    return CW

def CW_NW_6(MinD,DBH,a1,a2):
    if DBH >= MinD:
        CW = a1*DBH**a2
    else:
        CW = (a1*MinD**a2)*(DBH/MinD)
    return CW

def get_fuels(litter_name, needles_name, in_path, L2_res):
    litter = raster_import(os.path.join(in_path,litter_name))
    needles = raster_import(os.path.join(in_path,needles_name))
    fuels = litter + needles
    fuels = fuels/1000 #g to kg
    QF_fact = int(L2_res/2)
    fuels_interp = spline_interp(fuels,QF_fact)
    return fuels_interp
    
def spline_interp(fuels, QF_fact):
    x_L2 = np.linspace(0, fuels.shape[1], fuels.shape[1])
    y_L2 = np.linspace(0, fuels.shape[0], fuels.shape[0])
    x_QF = np.linspace(0, fuels.shape[1], fuels.shape[1]*QF_fact)
    y_QF = np.linspace(0, fuels.shape[0], fuels.shape[0]*QF_fact)
    splinterp = interpolate.RectBivariateSpline(y_L2, x_L2, fuels)
    fuels_interp = splinterp(y_QF, x_QF)
    
    OldMin = fuels_interp.min()
    OldMax = fuels_interp.max()
    NewMin = fuels.min()
    NewMax = fuels.max()
    
    OldRange = (OldMax - OldMin)  
    NewRange = (NewMax - NewMin)  
    squeezed = (((fuels_interp - OldMin) * NewRange) / OldRange) + NewMin
    return squeezed

def write_files(filedict, path):
    for i in range(len(filedict)):
        list(filedict.values())[i].to_csv(path_or_buf = os.path.join(path, list(filedict.keys())[i]),
                                          index = False)

def get_params(df,path,file,L2_res,QF_res,csv_name):
    qf_per_raster = int(L2_res/QF_res)
    MapCodes = raster_import(os.path.join(path,file))
    nx = int(MapCodes.shape[1]*qf_per_raster)
    ny = int(MapCodes.shape[0]*qf_per_raster)
    nz = int(np.ceil(df.HT_m.max()))
    num_spp = len(df["SPID"].unique())
    params = Treelist_params(QF_res, QF_res, 1.0, nx, ny, nz, csv_name, len(df), num_spp)
    return params

def print_fuellist(tp,path):
    with open(os.path.join(path,'fuellist'), 'w') as input_file:
        input_file.write("&fuellist\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("! FIRETEC domain info\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("      nx={} ny={} nz={}           ! Size of HIGRAD/FIRETEC grid [cells]\n".format(tp.nx, tp.ny, tp.nz))
        input_file.write("      dx={} dy={} dz=15           ! Grid Resolution [m]\n".format(tp.dx, tp.dy))
        input_file.write("      aa1=0.1                     ! Vertical stretching component [default=0.1]\n")
        input_file.write("      singlefuel=0                ! Flag forcing single fuel type instead of multiple fuels\n")
        input_file.write("      topofile='flat'  	        ! 'flat' -> no topo, 'name.dat' of topo file for topo\n")
        input_file.write("\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("! Data import from existing files\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("      ifuelin=0                     ! Fuel read-in flag (0 no existing fuel files, 1 read-in exisiting fuel files)\n")
        input_file.write("      inx=1500 iny=1500 inz=128     ! Size of HIGRAD/FIRETEC grid [cells]\n")
        input_file.write("      idx=2.0 idy=2.0 idz=1.0       ! Grid Resolution [m]\n")
        input_file.write("      iaa1=0.0                      ! Vertical stretching component [default=0.1]\n")
        input_file.write("      infuel=1                      ! Number of Fuel Types\n")
        input_file.write("      rhoffile='rhof.dat.orig'      ! Existing rhof fuel file\n")
        input_file.write("      ssfile='sav.dat.orig'         ! Existing sizescale fuel file\n")
        input_file.write("      moistfile='moisture.dat.orig' ! Existing moisture fuel file\n")
        input_file.write("      afdfile='depth.dat.orig'      ! Existing depth fuel file\n")
        input_file.write("\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("! Input dataset info\n")
        input_file.write("! ----------------------------------\n")
        input_file.write("      igrass=0                        ! Grass flag (1 is generalized grass data, 2 is if already have a file with 2d arrays for groundfuels)\n")
        input_file.write("      ngrass=1                        ! Number of Grass Species\n")
        input_file.write("      grassconstant=5                 ! Exponential constant used to determine the decay of grass mass with tree shading\n")
        input_file.write("      grassfile='grass.txt'           ! Grass text file with generalized-data\n")
        input_file.write("\n")
        input_file.write("      itrees=2                        ! Trees flag (1 is generalized tree data, 2 is specific tree data with locations, 3 is specific tree data with randomized locations)\n")
        input_file.write("      ntspecies={}                    ! Number of Tree Species\n".format(tp.num_spp))
        input_file.write("      tfuelbins=1                     ! Number of size bins to distribute branches\n")
        input_file.write("      treefile='{}'                   ! Trees text file with data\n".format(tp.csv_name))
        input_file.write("\n")
        input_file.write("      ilitter=0                       ! Litter flag\n")
        input_file.write("      litterconstant=5                ! Exponential constant to determine increase of litter mass under trees\n")
        input_file.write("      litterfile='litter.txt'         ! Litter text file with generalized-data\n")
        input_file.write("\n")
        input_file.write("      itreatment=0                    ! Treatment flag (0 is no treatment, 1 slashlayer, 2 slashpile, 3 clearing)\n")
        input_file.write("      sdnx=25,75                      ! Range of x cells undergoing treatment\n")
        input_file.write("      sdny=25,75                      ! Range of y cells undergoing treatment\n")
        input_file.write("      sdiameter=25                    ! Diameter of slashpiles\n")
        input_file.write("      sheight=15                      ! Height of slashpiles\n")
        input_file.write("      sprho=30                        ! Bulk density of slashpiles\n")
        input_file.write("      smoist=0.15                     ! Moisture content of slash\n")
        input_file.write("      sdepth=0.2                      ! Depth of slash layer\n")
        input_file.write("/\n")

if __name__=="__main__":
    toTreelist()


