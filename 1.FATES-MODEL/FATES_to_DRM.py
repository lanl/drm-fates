# -*- coding: utf-8 -*-
"""
Created on Mon Mar 11 17:13:45 2024

@author: 345578
"""
#import xarray as xr
import netCDF4
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.ndimage
import argparse 

ymin=0
ymax=400
xmin=0
xmax=400
fates_res=100
Fuel_sigma=1.2





################ Parsed arguments###################################
parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter)
### Name of the case you are running 
parser.add_argument('-c', '--casebase', dest='casebase',
                        help='base name of a case (casename=casebase${samplenumber+start_item})', metavar='STRING', required=True)
### First time step
parser.add_argument('-s', '--starttime', dest='startime',
                        help='', metavar='INT', required=True)
###Current iteration
parser.add_argument('-IT', '--ITT', dest='IT',
                        help='base name of a case (casename=casebase${samplenumber+start_item})', metavar='INT', required=True)
###Amount to add per iteration
parser.add_argument('-n', '--number', dest='number',
                        help='base name of a case (casename=casebase${samplenumber+start_item})', metavar='INT', required=True)
### The files you are looking for in casebase
parser.add_argument('-p', '--prefix', dest='prefix',
                        help='base name of a case (casename=casebase${samplenumber+start_item})', metavar='STR', required=True)
#######################################

####Set up arguments 
args = parser.parse_args()
dir1=args.casebase
### Script rules
time=int(args.startime)+(int(args.IT)*int(args.number))

##########################
### run rules
Out='/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/VDM2FM/'


#### Processor for vegetation ####
print(dir1+args.prefix+".elm.r."+str(time)+"-01-01-00000.nc")
Full=netCDF4.Dataset(dir1+args.prefix+".elm.r."+str(time)+"-01-01-00000.nc","r+")

########## Hard coded values the same
# FATES parameters
fates_CWD_frac_twig = 0.045 #CSXM: CWD may indicate Coarse Woody Debris
fates_c2b = 2 #CSXM: carbon to biomass multiplier of bulk structural tissues (unit: ratio)
# fates_leaf_slatop = 0.00662, 0.0189200006425381 # for long-leaf pine & Turkey oak #CSXM: specific Leaf Area (SLA) at top of canopy, projected area basis (unit: m^2/gC)
# denleaf  = -2.3231_r8*sla/prt_params%c2b(ft) + 781.899_r8 L863 kg/m3 #CSXM: leaf dry mass per unit fresh leaf volume (unit: kg/m3) and this formula is from FatesPlantHydraulicsMod.F90
leafdensity = -2.3231*(0.00662+0.0189200006425381)/2/2 + 781.899 # kg/m3
# fates_wood_density: 0.58, 0.65, # g/cm3 #CSXM: mean density of woody tissue in plant (unit: g/cm3)
# 0.615 g/cm3 = 615 kg/m3
wooddensity = (0.58 + 0.65)/2 # kg/m3 #CSxM: the value seems for g/cm3
# fire parameter; a version of fuel surface area to volume ratio
sizescale_pd_df = pd.DataFrame({'fates_pft': [1,2,3], #CSXM: 1=long-leaf pine (needleleaf_evergreen_extratrop_tree), 2=turkey oak (broadleaf_colddecid_extratrop_tree), 3=grass (c4_grass)
                                'sizescale': [0.2,0.6,1]})
grass_pft_index = 3

#CSXM (BGN)
# the following are from the FATES restarting file
# fates_pft:         plant functional type (index)
# fates_dbh:         diameter at breast height (cm)
# fates_height:      plant height (m)
# fates_crown_depth: plant crown depth fraction (fraction) 
# fates_nplant:      number of plants in the cohort (/patch)
# fates_cohort_area: area of the fates cohort (m2)
# leaf_c_val_001:    leaf         carbon, state var, position:001 (kg)
# sapw_c_val_001:    sapwood      carbon, state var, position:001 (kg)
# store_c_val_001:   storage      carbon, state var, position:001 (kg)
# repro_c_val_001:   reproductive carbon, state var, position:001 (kg)
# struct_c_val_001:  structural   carbon, state var, position:001 (kg)
#CSXM (END)
################## Get relevant things from FATES outputs. 

df=pd.DataFrame({"fates_pft":np.array(Full['fates_pft'][:]),
                 "fates_cohort_area":np.array(Full['fates_area'][:]),
                 "fates_height":np.array(Full["fates_height"][:]),
                 'fates_crown_depth' :np.array(Full['fates_canopy_trim'][:]),
                 "fates_nplant":np.array(Full["fates_nplant"][:]),
                 "fates_dbh":np.array(Full["fates_dbh"][:]),
                   "leaf_c_val_001":np.array(Full["leaf_c_val_001"][:]),
                     "sapw_c_val_001":np.array(Full["sapw_c_val_001"][:]),
                      "store_c_val_001":np.array(Full["store_c_val_001"][:]),
                       "repro_c_val_001":np.array(Full["repro_c_val_001"][:]),
                         "struct_c_val_001":np.array(Full[ "struct_c_val_001"][:])
              })




###### Calculating

df["cohort"]=list(range(0,len(df))) ### This saves cohorts for reaggregation 
df["fates_crown_rad"]= df['fates_height']*0.25
df["fates_crown_rad"].loc[df["fates_crown_rad"]<1.0]=.5
df['fates_nplant'] = df['fates_nplant']


df['fates_bagw_twig']=fates_CWD_frac_twig*fates_c2b*sum([df['sapw_c_val_001'], 
                                                         df['store_c_val_001'],
                                                       df['repro_c_val_001'], 
                                       df['struct_c_val_001']])

df['fates_bleaf'] =  df["leaf_c_val_001"] *fates_c2b       
#### Seperate allometry.
pftone=df[df["fates_pft"]==1]
pft2=df[df["fates_pft"]==2]
pftone['fates_height_to_crown_base']  = pftone['fates_height']*.8
pft2['fates_height_to_crown_base']  = pft2['fates_height']*.2
df=pd.concat([pftone,pft2])


df["SizeScale"]= 0.000347222
df["FuelMoisture"]=.3
leaf_twig_bulkd = (leafdensity*df['fates_bleaf'] + wooddensity*df['fates_bagw_twig'])/(df['fates_bleaf']  + df['fates_bagw_twig'])


#For longleaf: 0.2, 1.0, 0.000347222
#For Turkey Oak: 0.6, 1.0, 0.000347222
#### Fuel densities

df["FuelDensity"]=0
df["FuelDensity"]=np.where(df["fates_pft"] == 1, 0.2,df["FuelDensity"])
df["FuelDensity"]=np.where(df["fates_pft"] == 2, 0.6,df["FuelDensity"])


 # biomass weighted moisture content #CSXM: HYDRO is never TRUE in DRM, such that its fuel_moisture_content = 0.5 always
#if (HYDRO == TRUE) {
#    moisture <- read.table(file = file.path(VDM2FM, "livefuel.moisture.txt"), header = TRUE)
#    all.sam.var <- all.sam.var %>%
#      left_join(moisture, by = c("nsam", "fates_pft")) %>%
#      mutate(fuel_moisture_content = c(FATES_LTH_SCPF*fates_bleaf + FATES_STH_SCPF*fates_bagw_twig)/c(fates_bleaf + fates_bagw_twig)) %>%
#      select(-FATES_LTH_SCPF, -FATES_STH_SCPF)
#  } else {
#    all.sam.var <- all.sam.var %>%
#      mutate(fuel_moisture_content = 0.5)
# }
#!

############## Random assortment of plants 
#df["fates_nplant"]=df["fates_nplant"]
wholetrees=df.loc[df["fates_nplant"].round()>1]## Why not just not equal to zero

##make tree density into a list of the repeated value
new_df = wholetrees.loc[ wholetrees.index.repeat(wholetrees["fates_nplant"])].assign(fifo_qty=1).reset_index(drop=True)

### Create a random matrix 
location=pd.DataFrame({})
location['x_i']=np.random.randint(xmin,xmax,size=len(new_df))
location['y_i']=np.random.randint(xmin,xmax,size=len(new_df))
    ### create a unique ID 
location["key"]="x_"+location["x_i"].astype(str)+"y_"+location["y_i"].astype(str)
### Check for repeats. 
location["c"]=1
check9=location.groupby(["key"],as_index=False).count()
check9=check9[check9["y_i"]>1]["key"]
### if repeated Fuzz the location   
sub=location.loc[location["key"].isin(check9)]
for i in list(sub.index):
    location.iloc[i,0]=location.iloc[i,0]+np.random.rand(1,1).flatten().round(2)
    location.iloc[i,1]=location.iloc[i,1]+np.random.rand(1,1).flatten().round(2)
location["key"]="x_"+location["x_i"].astype(str)+"y_"+location["y_i"].astype(str)
### append random locations to data 
location[location["x_i"]>399]['x_i']=399
location[location["y_i"]>399]['y_i']=399

Outputs=pd.concat([location,new_df],axis=1)

#Adam says trees needs. 


#speciesNumber Xlocation, Ylocation, TreeHt, CrownBaseHt, CrownRadius, 
#Ht_toLargestCrown_Radius, FuelDensity, FuelMoisture, SizeScale.

Outputs["fates_height_to_crown_base"].loc[Outputs["fates_height_to_crown_base"]<=0.1]=0.1
Outputs["fates_crown_rad"].loc[Outputs["fates_crown_rad"]<=1.0]=1.0
Outputs["SizeScale"]= 0.000347222
Outputs["fates_cbh2"]=Outputs["fates_height_to_crown_base"]+0.3*(Outputs["fates_height"]-Outputs["fates_height_to_crown_base"])


Outclean=Outputs[["fates_pft","x_i","y_i","fates_height","fates_height_to_crown_base",
            "fates_crown_rad","fates_cbh2","FuelDensity","FuelMoisture","SizeScale"]]
Outclean=Outclean.round(6)

######### write to csv 
Outclean.to_csv(Out+"Pre_fire_trees_"+str(time)+".csv")
np.savetxt(Out+'treelist_VDM.dat',Outclean,delimiter=" ",fmt='%f')




### Save a LUT of starting set 
LUT=Outputs[["x_i","y_i","fates_pft","fates_dbh", 'cohort']]
LUT.to_csv(Out+"LUT.csv")


##############Fuels

xmax=200
ymax=200
sigma_y=1.5
sigma_x=1.5

Fuels=Full['fates_leaf_fines_vec_001'][:]+0.25*Full['fates_ag_cwd_vec_001'][:]
Fuels=Fuels[Fuels>0.00]
if len(Fuels)>3:
    Mu_1=np.mean(Fuels)
    sigma_1=np.std(Fuels)


default=np.zeros([xmax,ymax])
print(default.shape)
for x_i in list(range(0,xmax)):
    for y_i in list(range(0,ymax)):
        default[x_i,y_i]=np.abs(np.random.normal(Mu_1,sigma_1,1))

sigma = [Fuel_sigma , Fuel_sigma ]
y = sp.ndimage.filters.gaussian_filter(default, sigma, mode='constant')

y=np.round(np.abs(y),5)

if int(args.number) > 1:
	memory=pd.read_csv(Out+'MemoryFuels.csv',index_col=0)
	memory=memory.to_numpy()
	y=np.multiply(y, memory)


np.savetxt(Out+"/fuels"+str(time)+".txt",y,delimiter=" ",fmt='%f')
np.savetxt(Out+"/VDM_litter_WG.dat",y,delimiter=" ",fmt='%f')
np.savetxt(Out+"/VDM_litter_trees.dat",y,delimiter=" ",fmt='%f')




