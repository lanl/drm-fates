# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:39:10 2022

@author: Niko Tutland
"""

"""
This takes place after landis to treelist, before trees script

Actually it probably should happen before some steps in landis to treelist
Because the we need to preserve the landis resolution, we can't split landis grid cells
So the fire domain will have to line up with the landis grid
Strategy for this:
    - create burn domain (bbox)
    - convert raster to points (shapefile)
    - find all points that overlap with bbox
    - find farthest point in each cardinal direction (will be a few, pick one)
    - add or subtract resolution/2 to the proper coordinate
    - create a new bounding box with that extent

import community-input-file.csv
import IC.tif
import output-community.img
assign output-community the same crs (etc) as IC. maybe this goes in LandisParams?
crop output-community.tif based on user-supplied shapefile for burn unit
    - use ttrs_quicfire functions to create burn domain
    - use burn domain to crop output-community.tif
filter community-input-file to only include cropped output-community values
use fastfuels to get topo

DO ALL THE FIRE STUFF

updated community-input-file (from postfire treelist) should have all the same mapcodes
"""

import ttrs_quicfire as qf
import LANDIS_to_Treelist as lantotree
import geopandas as gpd
import numpy as np
import os
import shutil
import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
from shapely.geometry import Point, Polygon
# import TTRS_QUICFire_Support as ttrs
# from matplotlib import pyplot as plt
# import numpy as np

def Landis(lp):
    IC_path = os.path.join(lp.landis_path,lp.IC_map)
    OC_path = os.path.join(lp.landis_path,lp.OC_file)
    if lp.spinup:    
        #### Create burn domain that lines up with landis grid cells
        OG_PATH = os.getcwd()
        ## Build domain class from shape:
        shape_paths = qf.Shapefile_Paths()
        dom = qf.dom_from_burn_plot(shape_paths, buffer=200)
        cell_nums = [dom.nx,dom.ny,dom.nz]
        
        qf_arrs = qf.build_ff_domain(dom, FF_request=True) #we need this to use the topo
        filelist = ["Run","Runs","FilesToCopy","Ignitions","CopyToRuns"]
        for i in filelist:
            shutil.rmtree(os.path.join(OG_PATH,i), ignore_errors = True)
        
        ## Import spatial data
        # Burn domain
        burn_domain = gpd.read_file(dom.shape_paths.bbox).to_crs(epsg=5070)
        # Initial Communities raster
        ## Convert landis raster to points
        with rio.open(IC_path, 'r+') as landis_rast:
            val = landis_rast.read(1) # band 1
            no_data=landis_rast.nodata
            geometry = [Point(landis_rast.xy(x,y)[0],landis_rast.xy(x,y)[1]) for x,y in np.ndindex(val.shape) if val[x,y] != no_data]
            v = [val[x,y] for x,y in np.ndindex(val.shape) if val[x,y] != no_data]
            landis_pts = gpd.GeoDataFrame({'geometry':geometry,'data':v})
            landis_pts.crs = landis_rast.crs
        
        ## Find intersection
        domain_mask = landis_pts.within(burn_domain.loc[0,'geometry'])
        domain_pts = landis_pts.loc[domain_mask]
        
        ## Find max and min coordinates to create new burn domain
        N,S,E,W = [domain_pts.geometry.y.max() + lp.L2_res/2,
                   domain_pts.geometry.y.min() - lp.L2_res/2,
                   domain_pts.geometry.x.max() + lp.L2_res/2,
                   domain_pts.geometry.x.min() - lp.L2_res/2]
        
        domain_poly = Polygon([(W,S), (W,N), (E,N), (E,S), (W,S)])
        new_domain = gpd.GeoDataFrame({'col1':['name']}, geometry=[domain_poly], crs = domain_pts.crs)
        
        ## Write to file
        new_domain.to_file(os.path.join(os.path.join(OG_PATH,"Shapefiles","new_bbox.shp")))

    ### Clip landis to new burn domain
    ## Crop the map (raster)
    with fiona.open(os.path.join(os.path.join(OG_PATH,"Shapefiles","new_bbox.shp"))) as shapefile:
        new_domain = [feature["geometry"] for feature in shapefile]
    with rio.open(IC_path, "r+") as IC:
        with rio.open(OC_path, "r+") as OC:
            arr = OC.read(1)
            OC.transform = IC.transform
            OC.crs = IC.crs
            with rio.open(os.path.join(lp.landis_path,lp.OC_cropped), 
                          mode="w",
                          height=OC.height,
                          width=OC.width,
                          count=1,
                          dtype=arr.dtype,
                          crs="EPSG:5070",
                          transform=OC.transform) as new_IC:
                new_IC.write(arr,1)
                IC_path = os.path.join(lp.landis_path,lp.OC_cropped)
    with rio.open(IC_path,"r+") as IC:
        out_image, out_transform = rasterio.mask.mask(IC,new_domain,crop=True)
        out_meta = IC.meta
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
        with rio.open(os.path.join(lp.landis_path,lp.IC_cropped), "w", **out_meta) as IC_cropped:
            IC_cropped.write(out_image)
    
    ## Crop the community input file (csv)
    with rio.open(os.path.join(lp.landis_path,lp.IC_cropped),"r+") as IC:
        cropped_mc = IC.read(1).flatten().tolist()
    cif = pd.read_csv(os.path.join(lp.landis_path,lp.CIF_file))
    cif_cropped = cif[cif["MapCode"].isin(cropped_mc)]
    cif_cropped.to_csv(os.path.join(lp.landis_path,lp.CIF_cropped), index = False)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
