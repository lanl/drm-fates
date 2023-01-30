# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:39:10 2022

@author: Niko Tutland
"""

import geopandas as gpd
import numpy as np
import os
import rasterio as rio
import rasterio.mask
import fiona
import pandas as pd
import shutil
from shapely.geometry import Point, Polygon

def Landis(lp):
    if lp.crop_domain:   
        OG_PATH = os.getcwd()
        OC_path = os.path.join(lp.landis_path,"output-community-"+str(lp.year)+".img")
        litter_path = os.path.join(lp.landis_path,"NECN","SurfaceLitterBiomass-"+str(lp.year)+".img")
        needles_path = os.path.join(lp.landis_path,"NECN","ConiferNeedleBiomass-"+str(lp.year)+".img")
        
        if lp.spinup: 
            #### Create burn domain that lines up with landis grid cells
            IC_path = os.path.join(lp.landis_path,lp.IC_map)
            
            if os.path.exists(os.path.join(OG_PATH,"Shapefiles","burn_domain.shp")):
                with fiona.open(os.path.join(OG_PATH,"Shapefiles","burn_domain.shp")) as shapefile:
                    new_domain = [feature["geometry"] for feature in shapefile]
                
                ## Crop the initial communities raster (to get mean lat lon)
                crop_raster(IC_path, new_domain, lp.landis_path, lp.IC_cropped)
            else:
                ## Create initial domain from bounds of burn plot shapefile
                burn_plot = gpd.read_file(os.path.join(OG_PATH,"Shapefiles","burn_plot.shp")).to_crs(epsg=5070)
                tb = burn_plot.total_bounds
                buff = 200
                W = tb[0] - buff
                S = tb[1] - buff
                E = tb[2] + buff
                N = tb[3] + buff            
                domain_poly = Polygon([(W,S), (W,N), (E,N), (E,S), (W,S)])
                burn_domain = gpd.GeoDataFrame({'col1':['name']}, geometry=[domain_poly], crs = 5070)
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
                new_domain.to_file(os.path.join(os.path.join(OG_PATH,"Shapefiles","burn_domain.shp")))
                with fiona.open(os.path.join(os.path.join(OG_PATH,"Shapefiles","burn_domain.shp"))) as shapefile:
                    new_domain = [feature["geometry"] for feature in shapefile]
                
                ## Crop the initial communities raster (to get mean lat lon)
                crop_raster(IC_path, new_domain, lp.landis_path, lp.IC_cropped)

        ### Clip landis to new burn domain
        
        ## Import burn domain polygon
        with fiona.open(os.path.join(os.path.join(OG_PATH,"Shapefiles","new_bbox.shp"))) as shapefile:
            new_domain = [feature["geometry"] for feature in shapefile]
            
        if lp.spinup == False:
            IC_path = os.path.join(lp.landis_path,"IC_original_cropped.tif")
        
        ## Crop the output community raster
        OC_tif = georeference(IC_path,OC_path,"output-community-"+str(lp.year)+".tif",lp.landis_path)
        crop_raster(OC_tif,new_domain,lp.landis_path,"output-community-cycle"+str(lp.cycle)+"_cropped.tif")
        
        ## Crop fuels rasters
        litter_tif = georeference(IC_path,litter_path,"SurfaceLitterBiomass-"+str(lp.year)+".tif",lp.landis_path)
        needles_tif = georeference(IC_path,needles_path,"ConiferNeedleBiomass-"+str(lp.year)+".tif",lp.landis_path)
        crop_raster(litter_tif,new_domain,lp.landis_path,"SurfaceLitterBiomass-cycle"+str(lp.cycle)+"_cropped.tif")
        crop_raster(needles_tif,new_domain,lp.landis_path,"ConiferNeedleBiomass-cycle"+str(lp.cycle)+"_cropped.tif")
        
        ## Crop the community input file (csv)
        with rio.open(os.path.join(lp.landis_path,"output-community-cycle"+str(lp.cycle)+"_cropped.tif"),"r+") as OC:
            cropped_mc = OC.read(1).flatten().tolist()
        cif = pd.read_csv(os.path.join(lp.landis_path,"community-input-file-"+str(lp.year)+".csv"))
        cif_cropped = cif[cif["MapCode"].isin(cropped_mc)]
        cif_cropped.to_csv(os.path.join(lp.landis_path,"community-input-file-cycle"+str(lp.cycle)+"_cropped.csv"), index = False)
        
    else:
        print("LANDIS RUN NOT CROPPED")
        print('   - files called "*_cropped" are same size as LANDIS domain')
        ## Other scripts in the framework are hard-coded to look for cropped files, so if no cropping is 
        ## occurring, we will have to rename them to "_cropped" and copy them to the right location. 
        ## I will need to rewrite things to be more intuitive.
        src = [os.path.join(lp.landis_path,"output-community-"+str(lp.year)+".img"),
               os.path.join(lp.landis_path,"community-input-file-"+str(lp.year)+".csv")]
        dst = [os.path.join(lp.landis_path,"output-community-cycle"+str(lp.cycle)+"_cropped.tif"),
               os.path.join(lp.landis_path,"community-input-file-cycle"+str(lp.cycle)+"_cropped.csv")]
        for i in range(len(src)):
            shutil.copyfile(src[i], dst[i])
        ## georeference .img files
        if lp.spinup:
            shutil.copyfile(os.path.join(lp.landis_path,lp.IC_map),
                            os.path.join(lp.landis_path,lp.IC_cropped))
        IC_path = os.path.join(lp.landis_path,lp.IC_cropped)
        litter_img =  os.path.join(lp.landis_path,"NECN","SurfaceLitterBiomass-"+str(lp.year)+".img")
        litter_tif = "SurfaceLitterBiomass-cycle"+str(lp.cycle)+"_cropped.tif"
        needles_img = os.path.join(lp.landis_path,"NECN","ConiferNeedleBiomass-"+str(lp.year)+".img")
        needles_tif = "ConiferNeedleBiomass-cycle"+str(lp.cycle)+"_cropped.tif"
        georeference(IC_path,litter_img,litter_tif,lp.landis_path)
        georeference(IC_path,needles_img,needles_tif,lp.landis_path)
        
    
def crop_raster(raster_path, bbox, landis_path, out_name):
    with rio.open(raster_path,"r+") as rst:
        out_image, out_transform = rasterio.mask.mask(rst,bbox,crop=True)
        out_meta = rst.meta
        out_meta.update({"driver": "GTiff",
                         "height": out_image.shape[1],
                         "width": out_image.shape[2],
                         "transform": out_transform})
        with rio.open(os.path.join(landis_path,out_name), "w", **out_meta) as cropped:
            cropped.write(out_image)
    
def georeference(IC_path,img_path,tif_name,landis_path):
    with rio.open(IC_path, "r+") as IC:
        with rio.open(img_path, "r+") as img:
            arr = img.read(1)
            img.transform = IC.transform
            img.crs = IC.crs
            with rio.open(os.path.join(landis_path,tif_name), 
                          mode="w",
                          height=img.height,
                          width=img.width,
                          count=1,
                          dtype=arr.dtype,
                          crs="EPSG:5070",
                          transform=img.transform) as tif:
                tif.write(arr,1)
                out_path = os.path.join(landis_path,tif_name)
    return out_path

if __name__=="__main__":
    Landis()
    
    
    
    
    
    
    
    
    
    
    
