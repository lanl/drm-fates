import numpy as np
import os.path
import subprocess
from shutil import copyfile
import sys

def add_sat():
    meters = 10
    num_cell = meters
    sidecells = int(num_cell / 2)
    nx = 200
    ny = 200
    file_in='../5.TREES-QUICFIRE/Saturation.txt'
    #file_out='../5.TREES-QUICFIRE/BuffSaturation.txt'
    file_out='../5.TREES-QUICFIRE/Saturation.txt'
    print('ADDING SATURATION BUFFTER')    
    
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
    with open(file_out, 'w') as ss:
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
    with open(file_out,'a') as ss:
      np.savetxt(ss,satarray,fmt="%s")

    print('SATURATION BUFFTERI ADDED')
    return


def add_tree_buff():
    meters = 10
    num_cell = meters
    sidecells = int(num_cell / 2)
    nx = 200
    ny = 200
    print('ADDING TREE BUFFTER')
    file_in='LLM2FT/treelist_LLM.dat'
    file_out='LLM2FT/treelist_LLM.dat'
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
    print('TREE BUFFER ADDED')
    return


def add_surf_buff():
    meters = 10
    num_cell = meters
    sidecells = int(num_cell / 2)
    nx = 200
    ny = 200
    file_in2='LLM2FT/LLM_litter_trees.dat'
    file_in3='LLM2FT/LLM_litter_WG.dat'
    file_out2='LLM2FT/LLM_litter_trees.dat'
    file_out3='LLM2FT/LLM_litter_WG.dat'
    copyfile('LLM2FT/LLM_litter_trees.dat','LLM2FT/LLM_litter_trees.OG.dat')
    copyfile('LLM2FT/LLM_litter_WG.dat','LLM2FT/LLM_litter_WG.OG.dat')
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

    return


def remove_tree_buff():
    meters = 10
    num_cell = meters
    sidecells = int(num_cell / 2)
    nx = 200
    ny = 200
    file_in=  '../1.LLM-HSM-MODEL/FT2LLM/AfterFireTrees.txt'
    file_out= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireTrees.txt'
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

    return


def remove_surf_buff():
    meters = 10
    num_cell = meters
    sidecells = int(num_cell / 2)
    nx = 200
    ny = 200
    file_in2= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireLitter.txt'
    file_in3= '../1.LLM-HSM-MODEL/FT2LLM/AfterFireWG.txt'
    file_out2='../1.LLM-HSM-MODEL/FT2LLM/AfterFireLitter.txt'
    file_out3='../1.LLM-HSM-MODEL/FT2LLM/AfterFireWG.txt'
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

    return

