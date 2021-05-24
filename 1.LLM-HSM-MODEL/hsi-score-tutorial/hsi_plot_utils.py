#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import pandas as pd
import sys

def plot_area_matrix(xx,title):
    fig, ax = plt.subplots(figsize=(18, 12))
    im = ax.imshow(xx, interpolation = 'nearest', cmap='jet')#, vmin=0, vmax=200) 
    cbar=fig.colorbar(im, ax=ax, extend='both')
    cbar.ax.tick_params(labelsize=20) 
    #plt.pcolor(v, r, z, cmap='jet')
    plt.title(title,fontsize=20)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    ## Loop over data dimensions and create text annotations.
    #for i in range(len(xx[:,0])):
    #    for j in range(len(xx[0,:])):
    #        text = ax.text(j, i, round(xx[i, j],1),
    #                       ha="center", va="center", color="w")
    return ax

def gridon(ax):
    ax.plot(range(26),np.ones(26)*4.5,'w--',linewidth=3)
    ax.plot(range(26),np.ones(26)*9.5,'w--',linewidth=3)
    ax.plot(range(26),np.ones(26)*14.5,'w--',linewidth=3)
    ax.plot(range(26),np.ones(26)*19.5,'w--',linewidth=3)
    ax.plot(np.ones(27)*4.5, range(-1,26),'w--',linewidth=3)
    ax.plot(np.ones(27)*9.5, range(-1,26),'w--',linewidth=3)
    ax.plot(np.ones(27)*14.5, range(-1,26),'w--',linewidth=3)
    ax.plot(np.ones(27)*19.5, range(-1,26),'w--',linewidth=3)
    plt.xlim(0.5,24.5)
    plt.ylim(-0.5,24.5)
    
def plot_species_scores(pp):
    sc_rcw=np.asarray(pp.age_sc)+np.asarray(pp.hw_sc)+np.asarray(pp.ageHW_sc)+np.asarray(pp.hwHW_sc)
    plt.figure(figsize=(8, 6))
    plt.plot(sc_rcw*0.25,'r',linewidth=2)
    plt.plot(pp.gt_sc,'g',linewidth=2)
    plt.plot(pp.sq_sc,'b',linewidth=2)
    plt.xlabel('Time [years]',fontsize=16)
    plt.ylabel('Scores [-]',fontsize=16);
    plt.legend(['RCW','GT','SQ',],loc=0);

def plot_rcw_scores(pp):
    plt.figure(figsize=(8, 6))
    plt.plot(pp.age_sc,'r')
    plt.plot(pp.hw_sc,'b')
    plt.plot(pp.ageHW_sc,'k')
    plt.plot(pp.hwHW_sc,'k--')
    plt.xlabel('Time [years]',fontsize=16)
    plt.ylabel('Score [-]',fontsize=16);
    plt.legend(['LLP age','LLP height','HW age', 'HW height'],loc=0);

    
def plot_tree_count(pp):
    plt.plot(pp.lp_gt_age,'b--')
    plt.plot(pp.lp_tot_age,'b')
    plt.plot(pp.hw_gt_age,'r--')
    plt.plot(pp.hw_tot_age,'r')
    constr=str(pp.tree_mature_age)
    plt.legend(['llp > '+constr,'llp tot','hw > '+constr,'hw tot'])
    plt.ylabel('Number of trees [-]')
    plt.xlabel('Time [years]')

def plot_llp2hw_ratio(pp):
    constr=str(pp.tree_mature_age)
    print 'HW trees > '+constr+':',pp.hw_gt_age[-10:-1]
    print 'HW total:',pp.hw_tot_age[-1]
    print 'LLP trees > '+constr+':',pp.lp_gt_age[-10:-1]
    print 'LLP total:',pp.lp_tot_age[-1]
    print 'hw/llp ratio:',round((pp.lp_gt_age[-1]/pp.hw_gt_age[-1]),2)
    ratio=np.asarray(pp.lp_gt_age)/np.asarray(pp.hw_gt_age)
    plt.title('Ratio of mature LLP to HW for the entire domain')
    plt.ylabel('Ratio [-]')
    plt.xlabel('Time [years]')
    plt.plot(ratio,'b-');

def plot_multi_output(sc_gt,sc_sq,sc_rcw):
    plt.figure(figsize=(8, 6))
    plt.plot(sc_rcw[:,0],'r',linewidth=2);
    plt.plot(sc_gt[:,0],'b',linewidth=2);
    plt.plot(sc_sq[:,0],'g',linewidth=2);
    plt.legend(['RCW','GT','SQ',],loc=0);
    plt.plot(sc_rcw,'r',linewidth=2);
    plt.plot(sc_gt,'b',linewidth=2);
    plt.plot(sc_sq,'g',linewidth=2);
    plt.xlabel('Time [years]',fontsize=20)
    plt.ylabel('Score [-]',fontsize=20)

def plot_area_numerics(xx,title,icheck):
    ynames=[0.0,0.1,0.2,0.3]
    fig, ax = plt.subplots(figsize=(10, 8))
    im = ax.imshow(xx, interpolation = 'nearest', cmap='seismic', vmin=0, vmax=1) 
    if icheck: # equal 1 for readslopes case
        im = ax.imshow(xx, interpolation = 'nearest', cmap='seismic', vmin=-0.01, vmax=0.01) 
    cb=fig.colorbar(im, ax=ax, extend='both')
    
    cb.ax.set_yticklabels(cb.ax.get_yticklabels(), fontsize=20)
    #plt.pcolor(v, r, z, cmap='jet')
    plt.title(title, fontsize=20)
    [nx, ny]=np.shape(xx)    
    #print nx,ny
    ## Loop over data dimensions and create text annotations.
    #print xx[1, 1]
    for i in range(nx):
        for j in range(ny):
            text = ax.text(j, i, round(xx[i, j],2),
                           ha="center", va="center", color="y",fontsize=18)
    
    xnames = [0.0,0.05,0.1,0.15]
    xticks = np.arange(0,nx,1)
    yticks = np.arange(0,ny,1)
    ax.set_xticks(yticks)
    ax.set_yticks(xticks)
    ax.set_xticklabels(xnames, fontsize=22)
    ax.set_yticklabels(ynames, fontsize=22)
    plt.xlabel('mast prob', fontsize=22)
    plt.ylabel('fire prob', fontsize=22)

def readslopes(filename):

    with open(filename) as f:
        content = f.readlines()

    str1 = ''.join(content)
    str2 = str1.replace("{", "\n")
    str3 = str2.replace("'slp1':", "")
    str4 = str3.replace("'slp2':", "")
    str5 = str4.replace("'slp3':", "")
    str6 = str5.replace("}", "")
    str7 = str6.replace(" ", "")
    text_file = open("slope_out.csv", "w")
    text_file.write(str7)
    text_file.close()

    data = pd.read_csv("slope_out.csv",header=None) 
    vis=0
    if vis:
        Slrcw=fillscore(data[0])       
        Slsq=fillscore(data[1])
        Slgt=fillscore(data[2]) 
        plot_area_numerics(Slgt,'GT score',1);
        #plt.savefig('gt_score', bbox_inches="tight")
        plot_area_numerics(Slrcw,'RCW score',1);
        #plt.savefig('rcw_score', bbox_inches="tight")
        plot_area_numerics(Slsq,'SQ score',1);
    
    print data.head()
    return data#[Slrcw, Slsq, Slgt]

def fillscore(x):
    k=0
    sc=np.zeros([4,4])
    for i in range(4):
        for j in range(4):
            sc[i,j]=x[k]
            k+=1
    return sc 

def readdats(datfilename):

    m=np.loadtxt(datfilename, skiprows=3, usecols=[1,2,4,5,6]) 
    Srcw=fillscore(m[:,2])       
    Ssq=fillscore(m[:,3])
    Sgt=fillscore(m[:,4]) 

    plot_area_numerics(Sgt,'GT score',0)
    #plt.savefig('gt_score', bbox_inches="tight")
    plot_area_numerics(Srcw,'RCW score',0)
    #plt.savefig('rcw_score', bbox_inches="tight")
    plot_area_numerics(Ssq,'SQ score',0)
    

def barchart(datfilename,slopefile,title):
    
    f = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
    ax = f.add_subplot(gs[0])
    ax2 = f.add_subplot(gs[1])
    
    m=readslopes(slopefile)
    N = len(m[2])
    rcw = m[0]
    gt = m[2]
    sq = m[1]

    ind = np.arange(N)*4 
    width = 0.75
    xwidth = 1.55   
    xlocs=ind + 2*width / 2

    ax.bar(ind, rcw, width, label='RWC',align='center', alpha=0.8)
    ax.bar(ind + width, sq, width,  label='SHS',align='center', alpha=0.8)
    ax.bar(ind + 2*width, gt, width,  label='GT',align='center', alpha=0.8)
    ax.set_ylabel('Slope',fontsize=16)
    ax.tick_params(labelsize=14)
    plt.setp(ax, xticks=xlocs, xticklabels=('', '', '', '', '', '', '', '', \
                                   '', '', '', '', '', '', '', '' ))
    #---------------------------------------------

    m=np.loadtxt(datfilename, skiprows=3, usecols=[4,5,6])
    N = len(m[:,2])
    rcw = m[:,2]
    gt = m[:,0]
    sq = m[:,1]
    ax2.bar(ind, rcw, width, label='RWC',align='center', alpha=0.8)
    ax2.bar(ind + width, sq, width,  label='SHS',align='center', alpha=0.8)
    ax2.bar(ind + 2*width, gt, width,  label='GT',align='center', alpha=0.8)

    plt.setp(ax2, xticks=xlocs, xticklabels=('f:0 \n m:0', '0 \n 0.05', '0 \n 0.1', '0 \n 0.15', \
                                   '0.1 \n 0', '0.1 \n 0.05', '0.1 \n 0.1', '0.1 \n 0.15', \
                                   '0.2 \n 0', '0.2 \n 0.05', '0.2 \n 0.1', '0.2 \n 0.15', \
                                   '0.3 \n 0', '0.3 \n 0.05', '0.3 \n 0.1', '0.3 \n 0.15' ))

    ax2.set_ylabel('Score (mean)',fontsize=16)
    ax2.tick_params(labelsize=14)
    ax2.set_ylim(0, 1)
    ax2.legend(loc='best',fontsize=14)
    plt.text(15,1.42,title,fontsize=16)
    plt.tight_layout()

