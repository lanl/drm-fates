# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:39:10 2022

@author: Niko Tutland
"""

"""
This takes place after landis to treelist, before trees script

Actually it probably should happen before some steps in landis to treelist
Because the we need to preserve the landis resolution, we can't split landis grid cells
So the fire domain will have to line up with 

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
import os
import shutil
# import TTRS_QUICFire_Support as ttrs
# from matplotlib import pyplot as plt
# import numpy as np

def Landis(lp):
    OG_PATH = os.getcwd()
    ## Build domain class from shape:
    shape_paths = qf.Shapefile_Paths()
    dom = qf.dom_from_burn_plot(shape_paths, buffer=200)
    cell_nums = [dom.nx,dom.ny,dom.nz]
    
    qf_arrs = qf.build_ff_domain(dom, FF_request=True) #we need this to use the topo
    filelist = ["Run","Runs","FilesToCopy","Ignitions","CopyToRuns"]
    for i in filelist:
        shutil.rmtree(os.path.join(OG_PATH,i), ignore_errors = True)
    burn_plot = gpd.read_file(dom.shape_paths.burn_plot).to_crs(epsg=5070)
    burn_domain = gpd.read_file(dom.shape_paths.bbox).to_crs(epsg=5070)
    
    filepath = os.path.join(lp.landis_path,lp.IC_file)
    IC_map, nodata = lantotree.raster_import(filepath)
