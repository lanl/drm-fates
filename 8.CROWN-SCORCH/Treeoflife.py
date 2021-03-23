import numpy as np
from scipy import stats
import struct
import os.path
import sys
import itertools

fnameIn1   = 'PercentFuelChange.txt'
fnameIn2   = 'TreeTracker.txt'
fnameIn3   = 'treelist_LLM-400x400.txt'
fnameIn4   = 'LLM_litter_WG.txt'
fnameIn5   = 'LLM_litter_trees.txt'
fnameOut   = 'AfterFireTrees.txt'
fnameOut2  = 'AfterFireWG.txt'
fnameOut3  = 'AfterFireLitter.txt'

# DEFINE DOMAIN
Nx = 200
Ny = 200
Nz = 16
s = (200,200)
cellnum = 0
cellptr = 0
conccell = 0
array1d = 0
array1d = Nx*Ny*Nz
groundf = Nx*Ny
percentFuelChang1d = np.zeros(array1d)
grassfuel = np.zeros(groundf)
litterfuel = np.zeros(groundf) 
newgrassfuel = np.zeros(s)
newlitterfuel = np.zeros(s) 

# READ IN PERCENT FUEL AFTER Burn
pfc = open(fnameIn1, 'r')
count = 0
for line in pfc:
 percentFuelChang1d[count] = line
 count = count + 1

pfc.close()
print(percentFuelChang1d[8], percentFuelChang1d[47777])
print (type(percentFuelChang1d[8]))

## JUST CHECKING HERE TO SEE HOW MANY CELL'S REDUCED FUEL
ccc = 0
cc = 0
c = 0
for locz in range (Nz):
    for locy in range(Ny):
      for locx in range(Nx):
        ccc = locx + (locy*Nx) + (Nx*Ny*locz)
        #print locx, locy, locz, ccc
        if locz == 0:
           if percentFuelChang1d[ccc] < 1:
              c = c + 1
        if percentFuelChang1d[ccc] < 1:
           cc = cc + 1 
           #print locx, locy, locz, cc
print (cc, c)


## READ IN GROUND FUEL
gf = open(fnameIn4, 'r')
grassfuel = gf.read().split()
gf.close()
lf = open(fnameIn5, 'r')
litterfuel = lf.read().split()
lf.close()
gloc = 0
planarloc = 0
print (grassfuel[0],grassfuel[1],grassfuel[2],grassfuel[200],type(grassfuel[200]))
### ASSUMED A PLANAR VIEW OF LITTER AND GRASS ARRAYS FROM LLM **** BUT NEED TO CHECK !!!!!!!
for locy in range(Ny):
  for locx in range(Nx):
    gloc = locx + (locy * Nx)
    ### THIS IS TO CONVERT CARTISIAN TO PLANAR **** MAY NOT BE NECESSARY **** 
    planarloc = (Nx*(Ny-1)-(locy*Nx)) + locx 
    newgrassfuel[locy,locx] = float(grassfuel[gloc]) * percentFuelChang1d[planarloc] 
    newlitterfuel[locy,locx] = float(litterfuel[gloc]) * percentFuelChang1d[planarloc]
    #print (gloc, locx, locy, planarloc, percentFuelChang1d[planarloc], grassfuel[gloc], newgrassfuel[locx,locy])
#print (newgrassfuel[:,:])
#afw = open(fnameOut2, 'wb')
#afw.write(newgrassfuel[:,:])
#afw.close() 
np.savetxt(fnameOut2, newgrassfuel,fmt='%10.2f')
np.savetxt(fnameOut3, newlitterfuel,fmt='%10.2f')

cc = 0
c6 = 0
c7 = 0
## READ IN TREE TRACKER FILE
newfuel = 0.0
sppflag = 0
smallturk = 1
smallllp = 1 
tf = open(fnameIn2, 'r')
td = open(fnameIn3, 'r')
aft = open(fnameOut, 'wb')
for e_tf, e_td in itertools.izip(tf, td):
  line_tf = e_tf.split()
  line_td = e_td.split()
  #print (line_td)
  #print (line_tf)
  cellnum = int(line_tf[1])
  #print (type(line_tf[1]), type(cellnum), cellnum)
  totfuel = 0.0
  sppflag = int(line_td[0])
  for cid in range(cellnum):
    cellptr = int(line_tf[cid+2])
    conccell = float(line_tf[cid+2+cellnum])
    newfuel = percentFuelChang1d[cellptr] * conccell 
    totfuel = totfuel + newfuel
    cc = cc + 1
    #print (line_tf[0], cid, cellptr, percentFuelChang1d[cellptr], conccell, newfuel)
  #print (line_tf[0], line_td[0], totfuel)
  if sppflag == 1:
     if smallllp > totfuel:
        smallllp = totfuel
     if totfuel > 0.30:
        aft.write(e_td)
        #print ("LLP", totfuel)
     else:
        c6 = c6 + 1 
  elif sppflag == 2:
     if smallturk > totfuel:
        smallturk = totfuel
     if totfuel > 0.75:
        aft.write(e_td)
        #print ("TurkeyOak", totfuel)
     else:
       c7 = c7 + 1
  #print ("")

print ("")
print ("")
print ('Small LLP',smallllp,c6,'Small Turk',smallturk,c7, cc)

#with open(fnameIn2, 'r') as treetkr:
#  with open(fnameIn3, 'r') as treedta:
#       
#tf = open(fnameIn2, 'r')
#for line in tf:
#  line = line.strip()  
#  print (line)

#for locz in range (Nz):
#    for locy in range(Ny):
#      for locx in range(Nx):
        







