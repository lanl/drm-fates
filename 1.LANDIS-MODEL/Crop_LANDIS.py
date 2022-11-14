# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 13:39:10 2022

@author: Niko Tutland
"""

"""
This takes place after a landis run, before landis to treelist

import community-input-file.csv
import IC.tif
import output-community.img
assign output-community the same crs (etc) as IC. maybe this goes in LandisParams?
crop output-community.tif based on user inputs
filter community-input-file to only include cropped output-community values

DO ALL THE FIRE STUFF

updated community-input-file (from postfire treelist) should have all the same mapcodes
"""