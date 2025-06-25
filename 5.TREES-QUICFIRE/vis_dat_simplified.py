# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 16:03:54 2023

@author: FireScience
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0,"../1.LANDIS-MODEL")
import TTRS_QUICFire_Support as ttrs
import shutil

#====================INPUTS====================
#pathfiles------------------
#pf = pathfile where .dat file lives
pf = './'
#of = out pathfile where you save image (as .png)
of = './'
#b = the name of the fuels case you are visualizing
b = 'Adams_method'
#topofile = pathfile and name of topography (if used), if no topography input: ''
topofile = ''#'/Users/joliveto/Desktop/trees_fix_topo/Inputs/rampHill.dat'

#trees/viewing parameters------------------
datfile  = 'treesrhof.dat' #which .dat to make png
nfuel    = 14 #see bottom of trees output, number of output fuels
sum_flag = 1 #if = 0, the topdown figures will be z-slices specified by "plane" ; if = 1 plots topdown figures as a sum in z
plane    = 0 #z-index slice of plotting (0=ground/bottom layer), will *not* be used if sum_flag = 1

#grid parameters------------------
Nx      = 1050 
Ny      = 975 
Nz      = 68 
dx      = 2.0
dy      = 2.0
dz      = 15.0
aa1     = 0.1
f1      = 0.0
stretch = 2 

#======================= DEFINE FUNCTIONS =========================               
def readfield(fuelfile, Nx, Ny, Nz):
    np.frombuffer(fuelfile.read(4),'f')
    return np.frombuffer(fuelfile.read(Nx*Ny*Nz*4), 'f').reshape((Nx,Ny,Nz),order='F')

def plotTopdown(fig,axs,arr,title,X,Y,sum_flag,plane):
    if sum_flag == 1:
        arr = np.sum(arr,axis=2)
        sp1 = axs.pcolormesh(X[:,:,0],Y[:,:,0],arr,cmap='jet',shading='auto', )
        #sp1 = axs.imshow(arr,cmap='Greens')
    else:
        arr = arr[:,:,plane]
        sp1 = axs.pcolormesh(X[:,:,0],Y[:,:,0],arr,cmap='jet',shading='auto')
        #sp1 = axs.imshow(arr,cmap='Greens')
    cbar = fig.colorbar(sp1, ax=axs)
    #cbar.ax.set_ylabel(ylabel, rotation=270)
    axs.set_title(title)
    axs.set_xlabel('X [m]')
    axs.set_ylabel('Y [m]')
    return 0


#=============================MAIN===================================
rhof= np.zeros(nfuel*Nx*Ny*Nz).reshape(nfuel,Nx,Ny,Nz)
rhoffile = open(pf+datfile,'rb')
for ift in range(nfuel):
    print('Reading fuel type:',ift)
    rhof[ift,:,:,:] = readfield(rhoffile,Nx,Ny,Nz)
    trhof = rhof[ift,:,:,:]
    print( 'SPECIES ',ift+1,' MIN = ',np.min(trhof) ,' ; MAX = ',np.max(trhof))
rhoffile.close()
print(rhof.shape)

x = 12
spx = rhof[x,:,:,:]
spx_sum = np.sum(spx, axis = 2)
spx_sum = np.flip(spx_sum, axis=0)
spx_sum = np.rot90(spx_sum, 3)


sp_all = np.sum(rhof, axis=0)
# sp_all = np.sum(sp_all, axis=2)
sp_all = np.flip(sp_all, axis=1)
sp_all = np.rot90(sp_all,3)
sp_all = np.swapaxes(sp_all,0,2)
sp_all = np.swapaxes(sp_all, 1, 2)

plt.figure(1)
plt.set_cmap('viridis')
plt.imshow(sp_all[0,:,:],origin='lower')
plt.colorbar()

def sp_sum(nsp, nx, ny, nz, ii):
    dat_list = ["rhof","moist","ss","fueldepth"]
    for i in dat_list:
        print("Importing trees",str(i)," 4D dat file")
        shutil.copyfile("trees"+str(i)+".dat","trees"+str(i)+"4D_cycle"+str(ii)+".dat")
        rhof= np.zeros(nsp*nx*ny*nz).reshape(nsp,nx,ny,nz)
        datfile = "trees"+str(i)+".dat"
        rhoffile = open("./"+datfile,'rb')
        for ift in range(nsp):
            print('Reading species ',ift)
            rhof[ift,:,:,:] = readfield(rhoffile,nx,ny,nz)
            trhof = rhof[ift,:,:,:]
            print( 'SPECIES ',ift+1,' MIN = ',np.min(trhof) ,' ; MAX = ',np.max(trhof))
        rhoffile.close()
        print(rhof.shape)
        # flip, rotate, and swap axes
        sp_all = np.sum(rhof, axis=0)
        sp_all = np.flip(sp_all, axis=1)
        sp_all = np.rot90(sp_all,3)
        sp_all = np.swapaxes(sp_all,0,2)
        sp_all = np.swapaxes(sp_all, 1, 2)
        
        ttrs.export_fortran_dat_file(sp_all, datfile)
        print("3D array created for trees",str(i),".dat")
    return  

sp_sum(nfuel,Nx,Ny,Nz,1)

cell_nums = [Nx,Ny,Nz]
arr = ttrs.import_fortran_dat_file("treesrhof_summed.dat", cell_nums)

plt.figure(1)
plt.set_cmap('viridis')
plt.imshow(arr[3,:,:],origin='lower')
plt.colorbar()

