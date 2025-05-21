import numpy as np
import os.path
import subprocess
from shutil import copyfile
import sys


#def NoZeros():
nx = 200
ny = 200
file_in2='AfterFireLitter.0.txt'
file_in3='AfterFireWG.0.txt'
file_out2='LLM_litter_trees-NoZeros.dat'
file_out3='LLM_litter_WG-NoZeros.dat'
file2 = open(file_in2, 'r')
substrings = file2.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,(200,200))
data_array = data_array + 0.03
ddd = data_array
print(file_out2)
np.savetxt(file_out2,ddd,fmt='%s')

file3 = open(file_in3, 'r')
substrings = file3.read().split()
data = [float(i) for i in substrings]
data_array = np.reshape(data,(200,200))
data_array = data_array + 0.03
ddd = data_array
print(file_out3)
np.savetxt(file_out3,ddd,fmt='%s')


#    return

