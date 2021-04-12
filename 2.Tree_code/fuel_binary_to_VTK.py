#!/usr/bin/python

from evtk.hl import gridToVTK
import numpy as np
import struct
import os.path
import sys

metname       = './metric'
readfilename  = 'fuels.dat'
outname       = 'output'

nfuel = 5
lfuel = 12

# careful - listed in order of how they appear in binary output                                                                                        
fuel_field_names= ["rhof", "moist", "ss", "afd"]

#======================= DEFINE FUNCTIONS =========================                                                                                            

def readfield(infile, Nx, Ny, Nz, lfuel):
  output = np.zeros(Nx*Ny*Nz).reshape(Nx,Ny,Nz)
  temp = np.fromstring(infile.read(Nx*Ny*lfuel*4), 'f').reshape((Nx,Ny,lfuel),order='F')
  output[:,:,0:lfuel]+= temp[:,:,:]
  return output

def zheight(ZI):
# generates array of cell heights from z-index array                                                                                
        Z = np.copy(ZI, order='K')
        ZItemp = Z[0,0,:]
        ZItemp[0] = ZItemp[0] * 2
        for ii in range(1,len(ZItemp)):
                ZItemp[ii] = (ZItemp[ii] - sum(ZItemp[:ii]))*2
        for ii in range(len(ZItemp)):
                Z[:,:,ii] = ZItemp[ii]

        print(Z[1,1,:])
        return Z

def metrics(metname):
# returns Nx, Ny, Nz, volume                                                                                                                                                 
        f = open(metname, 'rb')
        Nx = struct.unpack("i",f.read(4))[0]
        Ny = struct.unpack("i",f.read(4))[0]
        Nz = struct.unpack("i",f.read(4))[0]
        zt = struct.unpack("d",f.read(8))[0]
        grdGenScheme = struct.unpack("i",f.read(4))[0]
        if (grdGenScheme ==2):
                coef1 = struct.unpack("d",f.read(8))[0]
                coef2 = struct.unpack("d",f.read(8))[0]
        dx = struct.unpack("d",f.read(8))[0]
        dy = struct.unpack("d",f.read(8))[0]
        dz = struct.unpack("d",f.read(8))[0]
        cnt = Nx*Ny*Nz
        XI = np.fromstring(f.read(Nx*Ny*Nz*8),'d').reshape((Nx,Ny,Nz), order='C')
        YI = np.fromstring(f.read(Nx*Ny*Nz*8),'d').reshape((Nx,Ny,Nz), order='C')
        ZI = np.fromstring(f.read(Nx*Ny*Nz*8),'d').reshape((Nx,Ny,Nz), order='C')
        Z = zheight(ZI)
        volume = np.multiply(dx,dy,Z)
        f.close()
        return Nx, Ny, Nz, XI, YI, ZI, volume

def read_fields(fname, Nx, Ny, Nz, lfuel):
  outputs = {}
  infile = open(fname, 'rb')

  for ii in range(len(fuel_field_names)):
    for iii in range(nfuel):
      tmpCNT = struct.unpack("i",infile.read(4))[0]
      if fuel_field_names[ii] == "afd":
        outputs[fuel_field_names[ii]+'_'+str(iii+1)] = readfield(infile, Nx, Ny, Nz, 1)
      else:
        outputs[fuel_field_names[ii]+'_'+str(iii+1)] = readfield(infile, Nx, Ny, Nz, lfuel)
      tmpCNT = struct.unpack("i",infile.read(4))[0]

  return outputs

def select_data(pD):
  output = {}
  for field in fuel_field_names:
      for i in range(nfuel):
          output[field+'_'+str(i+1)] = pD[field+'_'+str(i+1)]
  return output

def main():
        offset = 0
        point_data = {}

        Nx, Ny, Nz, XI, YI, ZI, volume = metrics(metname)
         
        f = open(readfilename, 'rb')
        point_data.update(read_fields(readfilename, Nx, Ny, Nz, lfuel))
        if not os.path.isfile(outname):
            gridToVTK(outname, XI, YI, ZI, pointData = select_data(point_data))
main()

