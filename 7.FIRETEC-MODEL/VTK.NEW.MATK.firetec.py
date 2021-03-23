#!/usr/bin/python
# --------------------------------------------------------------
#
#  This script reads binary FIRETEC output and converts it to 
#  VTK format. It can also be used to compute custom metrics and 
#  make simple plots by calling alex_modules. 
#
#  Using alex_modules:
# ------------------- 
#  - a.zheight(ZI)                        returns array of cell heights for each cell
#  - a.volume(Nx,Ny,Nz,dx,dy,zh)          returns array of cell volumes for each cell
#  - a.area_burned(rhof,Nx,Ny,dx)         returns domain-total area burned for one time step 
#  - a.consumption(rhofinit,rhof,volume)  returns domain-total value for fuel consumption for one time step
#  - a.smoke_vol(O2,Nx,Ny,Nz,volume)      returns domain-total volume of "smoke" for one time step
#  - a.vorticity(u,v,w,Nx,Ny,Nz,dx,dy,dz) returns arrays of the three components of vorticity: 
#                                         stream-wise, cross-stream and vertical 
#  - a.plot_reg(x,y,N,xlab,ylab,title,figname)  performs a linear regression between variables x and y, makes a 
#                                         simple scatterplot with the regression line and returns the r2 value
#  - a.TKE(u,v,w,Nx,Ny,Nz)                returns domain-average TKE
#  - a.TKEvert(u,v,w,Nx,Ny,Nz)            returns vector of TKE computed separately for each vertical layer
#
# --------------------------------------------------------------

from evtk.hl import gridToVTK
import numpy as np
from scipy import stats
import struct
import os.path
import sys
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

# --------------------------------------------------------------
# Input fields to read in order they are listed in binary files:
# --------------------------------------------------------------
field_names = ["u","v","w","theta","ka","kb","density","rhof", "rhof0", "rhow", "afd", "ss", "temps", "O2", "sies", "psiwmax", "firad", "fsrad", "convht", "tempg"]

# --------------------------------------------------------------
# These fields need to be divided by density:     
# --------------------------------------------------------------
div_by_dens = ["u", "v", "w", "theta", "ka", "kb", "O2"]

# --------------------------------------------------------------
# Output fields to write to VTK:                                
# --------------------------------------------------------------                                                   
fields_to_write = ["u", "v", "w", "theta", "ka", "kb", "O2", "rhof", "rhof0", "rhos", "temps", "tempg", "convht", "Erhof", "VerEsum", "VerEsum", "StartTime", "ResTime", "Deltrhof", "rhow"]

col      =['red','darkorange','goldenrod','greenyellow','green','mediumturquoise','dodgerblue','magenta','purple','gray']

# --------------------------------------------------------------
# FIRETEC file to extract initial rhowater from (to compute rhosinit)
# --------------------------------------------------------------
initial = 50
incr    = 50
offset  = 0

fname   = './comp/comp.out.'
outname = './comp/water.out.'
metname = '../metric.200x200x41.dat'

strt =  50 + offset
stp = 78500 + offset

# --------------------------------------------------------------
# Read "metric" file that defines grid size, etc
# --------------------------------------------------------------
metfile = open(metname, 'rb')

Nx = struct.unpack("i",metfile.read(4))[0]
Ny = struct.unpack("i",metfile.read(4))[0]
Nz = struct.unpack("i",metfile.read(4))[0]
zt = struct.unpack("d",metfile.read(8))[0]
grdGenScheme = struct.unpack("i",metfile.read(4))[0]
if (grdGenScheme ==2):
    coef1 = struct.unpack("d",metfile.read(8))[0]	
    coef2 = struct.unpack("d",metfile.read(8))[0]	
    
dx = struct.unpack("d",metfile.read(8))[0]
dy = struct.unpack("d",metfile.read(8))[0]
dz = struct.unpack("d",metfile.read(8))[0]
cnt = Nx*Ny*Nz

XI=np.fromstring(metfile.read(Nx*Ny*Nz*8),'d').reshape( (Nx,Ny,Nz), order = 'C')
YI=np.fromstring(metfile.read(Nx*Ny*Nz*8),'d').reshape( (Nx,Ny,Nz), order = 'C')
ZI=np.fromstring(metfile.read(Nx*Ny*Nz*8),'d').reshape( (Nx,Ny,Nz), order = 'C')
topo=np.reshape(np.fromstring(metfile.read(Nx*Ny*8),'d'), (Nx,Ny))

VZF=np.insert(ZI,0,0,axis=2)
VZ=np.diff(VZF)
VOL=np.multiply(VZ,4)
s = (200,200)
nn= (200,200,41)
VerEsum = np.zeros(s)
PercRhofChange = 0.0

metfile.close()

icnt = strt
lrng = (stp-strt)//incr+1

rhosum = np.zeros(lrng)

# ---------------------------------------------------------------
# Open TimeSeries Files
# --------------------------------------------------------------
fname4="AverageIntensity.txt"
fname5="TotalFuelConsumed.txt"
fname7="WindSpeedsCanopy.txt"
fname8="FireSize.txt"
with open(fname4, 'wb') as ff:
 ff.write(b' Time[s],   Ave_EnerygIntensity TreeIntensity [J/s m^2]\n')
with open(fname5, 'wb') as fc:
 fc.write(b' Time[s],  FeulConsumbed [kg\n')
with open(fname7, 'wb') as fca:
 fca.write(b' Time[s], Dom2   Tree2   Patch2   Dom6   Tree6   Patch6   Dom8   Tree8   Patch8 [m/s]\n')
with open(fname8, 'wb') as fap:
 fap.write(b' Time[s] FireArea[m2] ActiveArea FireParimeter[m] ActiveParimeter\n')



# --------------------------------------------------------------
# Read initial FIRETEC file to extract rhofinit
# --------------------------------------------------------------
infile = open(fname+str(initial), 'rb')
outputs = {}

for ifld in range(len(field_names)): #20):                                                                                
    tmpCNT = struct.unpack("i",infile.read(4))[0]
    outputs[field_names[ifld]] = np.fromstring(infile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')
    tmpCNT = struct.unpack("i",infile.read(4))[0]

for ifld in range(len(div_by_dens)):                                                                                      
    the_field = div_by_dens[ifld]                                                                                         
    outputs[the_field] = np.divide(outputs[the_field], outputs["density"])                                                 

outputs["rhos"] = np.add(outputs["rhof"], outputs["rhow"])
outputs["rhofinit"] = outputs["rhof"]
rhofn=(outputs["rhof"])
infile.close()

# --------------------------------------------------------------
# Get Initial Data Once
# --------------------------------------------------------------
CanopyArray=np.zeros(s)
for locy in range(Ny):
    for locx in range(Nx):
        vertfuel = (outputs["rhof0"][locx,locy,2:35])
        CanopyArray[locy,locx] = np.sum(vertfuel, axis = 0)
        #try useing rhof
        #if canopyfuel>0.1:
# Initialize StartTime
StartTime = np.zeros(nn)
EndTime = np.zeros(nn)
ResTime = np.zeros(nn)

# --------------------------------------------------------------
# Loop to read FIRETEC files and write a corresponding .vts file
# --------------------------------------------------------------
time=0.0
nncount=0
for it in range(lrng):
    vtsfile = outname+'.'+str(icnt)+'.vts'
    if os.path.isfile(vtsfile):
        print(vtsfile,' already exists')
    else:
        infile = open(fname+str(icnt), 'rb')
        outputs = {}
        #print("Input File", infile)
        for ifld in range(len(field_names)): #20):
            tmpCNT = struct.unpack("i",infile.read(4))[0]
            outputs[field_names[ifld]] = np.fromstring(infile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')
            tmpCNT = struct.unpack("i",infile.read(4))[0]

        for ifld in range(len(div_by_dens)):
            the_field = div_by_dens[ifld]
            outputs[the_field] = np.divide(outputs[the_field], outputs["density"])

        outputs["rhos"] = np.add(outputs["rhof"], outputs["rhow"])
        Erhof=np.subtract(rhofn,(outputs["rhof"]))
        Erhof=np.multiply(Erhof,VOL)
        Erhof=np.dot(Erhof,19500)
        # Energy Threshold
        threshold_indices = Erhof < 2
        Erhof[threshold_indices] = 0
        outputs["Erhof"] = Erhof  
        # Fire Start and Residance Time
        EndTime[Erhof > 0] = time
        mask = StartTime == 0
        StartTime[mask] = EndTime[mask]
        #StartTime[StartTime < 0.01] = EndTime 
        ResTime = np.subtract(EndTime,StartTime)
        #print ("StartTime:", StartTime[StartTime > 0])
        #print ("ResTime:", ResTime[ResTime > 0])
        outputs["StartTime"] = StartTime 
        outputs["ResTime"] = ResTime

        Deltrhof=np.subtract(outputs["rhof0"],outputs["rhof"])
        Drhofp=np.divide(outputs["rhof"],outputs["rhof0"])
        if it > 2:
           VerEnergy=np.sum(Erhof[:,:,0:20],axis=2)
        else:
           VerEnergy=np.zeros(s)
        FireArea=np.count_nonzero(VerEnergy)
        VerEsum=np.add(VerEsum,VerEnergy)
        VerEsumVis=np.zeros(nn)
        VerEVis=np.zeros(nn)
        for x in range(0, 40):
            VerEsumVis[:,:,x]=VerEsum
            VerEVis[:,:,x]=VerEnergy     
        outputs["VerEsum"] = VerEsumVis
        outputs["VerEVis"] = VerEVis 
        outputs["Deltrhof"] = Deltrhof
        outputs["PercenDR"] = Drhofp
        infile.close()
        rhofn=(outputs["rhof"]) 
# -----------------------------------------
# Put your analysis code here
# -----------------------------------------
        FireArea=FireArea*4
        TotalFireArea= np.count_nonzero(VerEsum) * 4
        Intensity=np.sum(VerEnergy)
        RhoConsumed=np.sum(Deltrhof)
        if FireArea > 0:
           Intensity=Intensity/FireArea
        #print ("InTensity:",Intensity)
        FNtracker = 0
        ActTrack = 0
        P = 0
        AP = 0
        #Find Fire Parimeter
        for locy in range(Ny):
          for locx in range(Nx):
             if FNtracker == 0 and VerEsum[locy,locx]>0:
                P = P + 2
                if VerEsum[locy-1,locx]==0:
                   P = P + 2
                if VerEsum[locy+1,locx]==0:
                   P = P + 2
                FNtracker = 1
             elif FNtracker == 1 and VerEsum[locy,locx]>0:
                if VerEsum[locy-1,locx]==0:
                   P = P + 2
                if VerEsum[locy+1,locx]==0:
                   P = P + 2
             elif FNtracker == 1 and VerEsum[locy,locx]==0:
                P = P + 2
                FNtracker = 0
       
             if ActTrack == 0 and VerEnergy[locy,locx]>0:
                AP = AP + 2
                if VerEnergy[locy-1,locx]==0:
                   AP = AP + 2
                if VerEnergy[locy+1,locx]==0:
                   AP = AP + 2
                ActTrack = 1
             elif ActTrack == 1 and VerEnergy[locy,locx]>0:
                if VerEnergy[locy-1,locx]==0:
                   AP = AP + 2
                if VerEnergy[locy+1,locx]==0:
                   AP = AP + 2
             elif ActTrack == 1 and VerEnergy[locy,locx]==0:
                AP = AP + 2
                ActTrack = 0
        #print ('Perimeter', P, AP)
        with open(fname5, 'a') as fc:
          fc.write('{} {}\n'.format(time,RhoConsumed))

        # Calculate Trees vs patches
        # Average windspeed 2nd cell here
        u_ave2 = np.mean((outputs["u"][:,:,2]))
        u_ave6 = np.mean((outputs["u"][:,:,6]))
        u_ave8 = np.mean((outputs["u"][:,:,8]))
        Tree2sum=0.0
        Patch2sum=0.0
        Tree6sum=0.0
        Patch6sum=0.0
        Tree8sum=0.0
        Patch8sum=0.0
        cfuel2sum=0.0
        cpatch2sum=0.0
        cfuel6sum=0.0
        cpatch6sum=0.0
        cfuel8sum=0.0
        cpatch8sum=0.0
        t2sumcount=0
        p2sumcount=0
        tcount=0
        pcount=0
        TreeIntensity=0
        TreeEsum=0
        TreeArea=0
        fname10="PercentFuelChange" + str(icnt) + ".txt"
        PercRhofChange = 0.0
        nncount = 0
        Drhofp[np.isnan(Drhofp)]=1 
        with open(fname10, 'wb') as pfc: 
         for locz in range (16):
             for locy in range(Ny):
               for locx in range(Nx):
                    #print (nncount, locz, locy, locx, outputs["PercenDR"][locx,locy,locz])
                    PercRhofChange = Drhofp[locx,locy,locz]
                    pfc.write('{}\n'.format(PercRhofChange)) 
                    nncount=nncount + 1
        pfc.close()
        for locy in range(Ny):
            for locx in range(Nx):
                vertfuel = (outputs["rhof0"][locx,locy,2:35])
                #try useing rhof
                if CanopyArray[locy,locx]>0.1:
                   cfuel2sum = cfuel2sum + (outputs["u"][locx,locy,2])
                   cfuel6sum = cfuel6sum + (outputs["u"][locx,locy,6])
                   cfuel8sum = cfuel8sum + (outputs["u"][locx,locy,8])
                   if np.sum(Erhof[locx,locy,:]) > 0.0:
                      TreeEsum = TreeEsum + np.sum(Erhof[locx,locy,:])
                      TreeArea = TreeArea + 1
                      #print ('Are we Ever here?')
                   tcount=tcount + 1
                   pfuelvar = 'fuel_canopy'
        
                else:
                   cpatch2sum = cpatch2sum + (outputs["u"][locx,locy,2])
                   cpatch6sum = cpatch6sum + (outputs["u"][locx,locy,6])
                   cpatch8sum = cpatch8sum + (outputs["u"][locx,locy,8])
                   pcount = pcount + 1
                   pfuelvar = "fuel_patch"
                #print (pstemvar, pfuelvar, locy, locx, u[locx,locy,0])   
        cfuel2mean=cfuel2sum/tcount
        cpatch2mean=cpatch2sum/pcount
        cfuel6mean=cfuel6sum/tcount
        cpatch6mean=cpatch6sum/pcount
        cfuel8mean=cfuel8sum/tcount
        cpatch8mean=cpatch8sum/pcount
        TreeArea=TreeArea*4
        if TreeArea>0:
           TreeIntensity=TreeEsum/(TreeArea)
        print("InTensity:",Intensity, tcount, pcount, FireArea, TreeArea, time)
        with open(fname4, 'a') as ff:
          ff.write('{} {} {}\n'.format(time,Intensity,TreeIntensity))
        with open(fname7, 'a') as fca:
          fca.write('{} {} {} {} {} {} {} {} {} {}\n'.format(time,u_ave2,cfuel2mean,cpatch2mean,u_ave6,cfuel6mean,cpatch6mean,u_ave8,cfuel8mean,cpatch8mean))
        with open(fname8, 'a') as fap:
           fap.write('{} {} {} {} {}\n'.format(time,TotalFireArea,FireArea,P,AP))

# -----------------------------------------
# Write output to VTK (you probably won't need this...
# -----------------------------------------
        output_vars = {}
        print(len(fields_to_write))  
        print(len(outputs))         
        for field in fields_to_write:
            output_vars[field] = outputs[field]
        gridToVTK(outname+str(icnt), XI, YI, ZI, pointData = output_vars)
        time=time+0.5
        icnt += incr  # close loop over multiple files in one simulation
                      # (for trees.temporal.variabiltiy there is an offset of the ignition time, 
                      #  but you can comment this out if you look at trees.spatial.variability)
