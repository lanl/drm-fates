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

file_in=  '../1.LLM-HSM-MODEL/FT2LLM/AfterFireTrees.txt'
file_in2= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireLitter.txt'
file_in3= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireWG.txt'

file_out= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireTrees.txt'
file_out2='../1.LLM-HSM-MODEL/FT2LLM/AfterFireLitter.txt'
file_out3='../1.LLM-HSM-MODEL/FT2LLM/AfterFireWG.txt'

file1 = open(file_in, 'r')
Lines = file1.readlines()
count = len(Lines)
c = 0
data = []
for line in Lines:
    a = Lines[c].split(' ')
    a = [float(i) for i in a]
    a[1] = a[1] - meters
    a[2] = a[2] - meters
    data.append(a)
    c = c + 1
Re_BuffData = np.array(data)
np.savetxt(file_out,Re_BuffData,fmt='%s %s %s %s %s %s %s %s %s %s')

file2 = open(file_in2, 'r')
substrings = file2.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,((nx + num_cell),(ny + num_cell)))
ddd = data_array
for x in range(0, sidecells):
    ddd = np.delete(ddd, 0, 0)
    ddd = np.delete(ddd, (ny + sidecells - x -1), 0)
for y in range(0, sidecells):
    ddd = np.delete(ddd, 0, 1)
    ddd = np.delete(ddd, (nx + sidecells - y -1), 1)
np.savetxt(file_out2,ddd,fmt='%s')

file3 = open(file_in3, 'r')
substrings = file3.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,((nx + num_cell),(ny + num_cell)))
ddd = data_array
for x in range(0, sidecells):
    ddd = np.delete(ddd, 0, 0)
    ddd = np.delete(ddd, (ny + sidecells - x -1), 0)
for y in range(0, sidecells):
    ddd = np.delete(ddd, 0, 1)
    ddd = np.delete(ddd, (nx + sidecells - y -1), 1)
np.savetxt(file_out3,ddd,fmt='%s')

