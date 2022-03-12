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

file_in='LLM2FT/treelist_LLM.dat'
file_in2='LLM2FT/LLM_litter_trees.dat'
file_in3='LLM2FT/LLM_litter_WG.dat'

file_out='LLM2FT/treelist_LLM.dat'
file_out2='LLM2FT/LLM_litter_trees.dat'
file_out3='LLM2FT/LLM_litter_WG.dat'

file1 = open(file_in, 'r')
Lines = file1.readlines()
count = len(Lines)
c = 0
data = []
for line in Lines:
    a = Lines[c].split(' ')
    a = [float(i) for i in a]
    a[1] = a[1] + meters
    a[2] = a[2] + meters
    data.append(a)
    c = c + 1
BuffData = np.array(data)
np.savetxt(file_out,BuffData,fmt='%s %s %s %s %s %s %s %s %s %s')

file2 = open(file_in2, 'r')
substrings = file2.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,(200,200))
ddd = data_array
nxx = np.zeros(nx)
nxx[:] = 0.1
nyy = np.zeros(ny + num_cell)
nyy[:] = 0.1
for x in range(0, sidecells):
   ddd = np.insert(ddd, (ny + x), nxx, 1)
   ddd = np.insert(ddd, 0, nxx, 1)
for y in range(0, sidecells):
    ddd = np.insert(ddd, (nx + y), nyy, 0)
    ddd = np.insert(ddd, 0, nyy, 0)
np.savetxt(file_out2,ddd,fmt='%s')

file3 = open(file_in3, 'r')
substrings = file3.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,(200,200))
ddd = data_array
nxx = np.zeros(nx)
nxx[:] = 0.1
nyy = np.zeros(ny + num_cell)
nyy[:] = 0.1
for x in range(0, sidecells):
   ddd = np.insert(ddd, (ny + x), nxx, 1)
   ddd = np.insert(ddd, 0, nxx, 1)
for y in range(0, sidecells):
    ddd = np.insert(ddd, (nx + y), nyy, 0)
    ddd = np.insert(ddd, 0, nyy, 0)
np.savetxt(file_out3,ddd,fmt='%s')


