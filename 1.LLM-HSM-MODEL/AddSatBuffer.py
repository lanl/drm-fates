import numpy as np
import os.path
import subprocess
from shutil import copyfile
import sys

meters = 10
num_cell = meters
sidecells = num_cell / 2

nx = 200
ny = 200

file_in='../5.TREES-QUICFIRE/Saturation.txt'
#file_out='../5.TREES-QUICFIRE/BuffSaturation.txt'
file_out='../5.TREES-QUICFIRE/Saturation.txt'

file = open(file_in, 'r')
Lines = file.readlines()
count = len(Lines)
a = Lines[0].split(' ')
a = [int(i) for i in a]

nx = a[0]
ny = a[1]
nz = a[2]
cells = nx*ny*nz
nnx = int(nx + sidecells + sidecells)
nny = int(ny + sidecells + sidecells)
satarray = []
with open(file_out, 'wb') as ss:
  ss.write('{} {} {}\n'.format(nnx,nny,nz))
c = 0
cc = 0
for z in range(0, nz):
    for y in range(0, nny):
        for x in range(0, nnx):
            if (y>sidecells-1)and(x>sidecells-1)and(y<(nny-sidecells))and(x<(nnx-sidecells)):
                satarray.append(float(Lines[cc+1]))
                cc = cc + 1
            else:
                satarray.append(0.252525)
            c = c + 1
with open(file_out,'ab') as ss:
  np.savetxt(ss,satarray,fmt="%s")
