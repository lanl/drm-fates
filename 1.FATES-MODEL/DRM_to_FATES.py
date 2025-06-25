#import xarray as xr
import netCDF4
import pandas as pd
import numpy as np
#import matplotlib.pyplot as plt
import scipy as sp
import scipy.ndimage
import argparse





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

args = parser.parse_args()


dir1=args.casebase
time=int(args.startime)+(int(args.IT)*int(args.number))


###Load full, Load LUT
Full=netCDF4.Dataset(dir1+args.prefix+".elm.r."+str(time)+"-01-01-00000.nc","r+")
#print(Full.keys())
Survived=pd.read_csv('/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/FM2VDM/AfterFireTrees.txt',delimiter=" ",header=None)
LUT=pd.read_csv("/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/VDM2FM/LUT.csv",index_col=0)
print(Survived)
Survived.columns=["fates_pft","x_i","y_i","fates_height","fates_crown_depth",
            "fates_crown_dia","fates_crown_depth","FuelDensity","FuelMoisture",
            "SizeScale"]



Cleaned=Survived.merge(LUT,on=["x_i","y_i","fates_pft"],how="left")

Cleaned=Cleaned.fillna(0)
Cleaned.to_csv("/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/FM2VDM/Afterfire_trees_"+str(time)+".csv")
meme=Cleaned.groupby("cohort",as_index=False).count()[['cohort','fates_dbh']].merge(LUT.groupby("cohort",as_index=False).count()[['cohort','fates_dbh']],on="cohort",how="right")
meme=meme.fillna(0)
meme["Survived"]=meme["fates_dbh_x"]/meme["fates_dbh_y"]
#meme.to_csv("/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/FM2VDM/Afterfire_"+str(time)+".csv")
template=pd.DataFrame({'cohort':list(range(0,int(meme['cohort'].max())))})
output=template.merge(meme,on="cohort",how="left")
Survivor=output.fillna(0)
Mult=Survivor["Survived"].values

Full["fates_nplant"][0:len(Mult)]=np.array(Full["fates_nplant"][0:len(Mult)])*Mult
Fuel=pd.read_csv('/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/FM2VDM/AfterFireLitter.txt',delimiter="    ",header=None)
Fuels_u=Fuel.mean().mean()
#print(Fuel)

print(len(Full['fates_leaf_fines_vec_001']))

Full['fates_leaf_fines_vec_001'][:]=np.repeat(Fuels_u,len(Full['fates_leaf_fines_vec_001']))
Full['fates_ag_cwd_vec_001'][:]=np.repeat(Fuels_u,len(Full['fates_leaf_fines_vec_001']))

minfuels=Fuel+0.1 ###Adj spatial weighting factor.

Fuels_u=minfuels.mean().mean()

Memoryfuels=minfuels/Fuels_u



Memoryfuels.to_csv("/lustre/scratch5/.mdt1/zjrobbins/drm-fates/1.FATES-MODEL/VDM2FM/MemoryFuels.csv")

#print(Full['fates_ag_cwd_vec_001'][:])
Full.close()
