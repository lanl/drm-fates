#!/usr/bin/env python
import numpy as np
import sys
import random

def midstory(pl_height,hw_height):
    a=pl_height.flatten()
    a[a < 2]=0
    a[a > 20]=0
    a[a != 0]=1

    b=hw_height.flatten()
    b[b < 2]=0
    b[b > 15]=0
    b[b != 0]=1

    return (np.sum(a)+np.sum(b))

class hsi_score:
    def norm_gaussian(self):
        # mean is 0
        # sigma defines the spread of Gaussian distribution
        x = self.x0
        s = self.sigma
        #print 'sigma:', s
        y0=np.exp(- x*x / (2.0 * s*s))   
        return y0
    
    def gscore(self):
        # num: is a current number
        # bn: benchmarking number 
        # sm: smoothing
        # sigma: is a distribution (bigger sigma, wider the distribution)
        num = self.number
        bmn = self.bench
        sm = self.smooth
        r=float(abs(bmn-num)) # absolute distance from 0 mean
        #print 'num',num,' bn',bmn,' r',r,' sm',self.smooth,' sigma',self.sigma
        if r > sm:
            self.x0=r-sm
            y0=self.norm_gaussian()
        else:
            y0=1
        #print 'score:', y0
        return y0
    
    def gscore_range(self, x):
        # x: a 2D array (matrix)
        # lb,ub: prescribed range
        # bn: benchmark number corresponds to the ideal number of elements 
        # smooth : relaxes the benchmark constrain 
        # (e.g. benchmark number of pines for RCW is 200 trees per ha)
        x1=np.asarray(x).reshape(-1)
        xx=x1[self.lb<=x1]
        self.number=len(xx[xx<=self.ub])
        self.sigma=1
        sc=self.gscore()
        #print 'num:',self.number,'bn:',self.bench,'sc:',sc
        return sc
    
    def subset_matrix(self, a, i, j):
        [n,m]=a.shape
        b=a[i:i+self.nsize,j:j+self.msize]
        #print b.shape
        return b
    
    def score_subset(self, x):#, a, b, bn, nsize=5, msize=5, smooth=1):
        # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
        # only works when matrix can be equally divided to subsets
        # nsize : vertical size of the window
        # msize : horizonal size of the window 

        #print  nsize, msize   
        [n,m]=x.shape
        v=range(0,int(n),self.nsize)
        #h=range(0,m-msize+1) # if steping size is 1
        h=range(0,int(m),self.msize)

        if np.mod(int(m),self.msize): # NOTE, mod(4,2) = 0, i.e. false 
            print 'ERROR: matrix horizonal dim:',m,' cannot be div by',self.msize
        if np.mod(int(n),self.nsize):
            print 'ERROR: matrix vertical dim:',n,' cannot be div by',self.nsize
        s=len(v)*len(h)
        subset_sc = np.zeros(s, 'float') 
        #print subset_sc
        k=0
        for i in v:
            for j in h:
                #print i,j
                sub=self.subset_matrix(x, i, j)
                #print self.lb,self.ub,self.bench,self.smooth
                subset_sc[k]=self.gscore_range(sub)#, a, b, bn, smooth)
                #print sub
                #print 'score: ',subset_sc[k],'bn: ',self.bench
                #print '-------------'
                k+=1

        sc=np.mean(subset_sc)
        #print 'mean score:',sc
        #print '-------------'
        return sc

    def prob_of_use_sq(self, lp_age, lp_count,hw_age, hw_count, mature_tree_age=5):
        # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
        # only works when matrix can be equally divided to subsets
        # nsize : vertical size of the window
        # msize : horizonal size of the window 

        # get the dimensions of the matrix   
        [n,m]=np.shape(lp_age)
        #NOTE: size of the subset is fixed
        self.nsize = 5
        self.msize = 5

        # check if matrix can be evenly divided
        if np.mod(int(m),self.msize):
            raise ValueError('matrix horizonal dim:',m,' cannot be divided by',self.msize)
        if np.mod(int(n),self.nsize):
            raise ValueError('ERROR: matrix vertical dim:',n,' cannot be divided by',self.nsize)  

        # divide matrix to equal subsets
        v=range(0,int(n),self.nsize)
        #h=range(0,m-msize+1) # if steping size is 1
        h=range(0,int(m),self.msize)

        # loop through subsets and calculate the #lp/#hw ratio
        ratio=[]
        prob_of_use=0.0
        for i in v:
            for j in h:
                # subset each matrix
                sub_pl_age=self.subset_matrix(lp_age, i, j)
                sub_lp_count=self.subset_matrix(lp_count, i, j)
                sub_hw_age=self.subset_matrix(hw_age, i, j)
                sub_hw_count=self.subset_matrix(hw_count, i, j)
                # find indexes of the mature LP and HW trees
                idl=np.nonzero(sub_pl_age > mature_tree_age) 
                idh=np.nonzero(sub_hw_age > mature_tree_age)
                # sum number of trees per cell
                num_of_lp=np.sum(sub_lp_count[idl[0][:],idl[1][:]])
                num_of_hw=np.sum(sub_hw_count[idh[0][:],idh[1][:]])
                # check division by zero 
                if num_of_hw!=0:
                    ratio.append(num_of_lp/num_of_hw)
                    #print num_of_lp,num_of_hw,num_of_lp/num_of_hw

        if ratio==[]:
            print 'WARNING: in prob_of_use_sq, no hardwood trees older than ',mature_tree_age,\
                    'years.'
            if idl!=[]: 
                prob_of_use=0.2
            else:
                print 'ERROR: in prob_of_use_sq, all mature trees killed by fire'
                sys.exit(1)
        else:
            # Probability of use curve calculated according to the Fig.2 (Perkins et al., 2008)
            x0=np.mean(np.asarray(ratio))
            #print 'x0:',x0
            self.number = x0
            self.bench = 8 #check this!!!
            self.smooth = 0
            self.sigma = 6
            #prob_of_use=self.gscore(8-x0, 8, 0, sigma=6)-.2
            prob_of_use=self.gscore()-.2
            #print 'sfs score:',prob_of_use
            if prob_of_use<0:
                prob_of_use=0

        return prob_of_use



    def prob_of_use_gt(self, lp_height, lp_count, hw_height, hw_count):
        # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
        # only works when matrix can be equally divided to subsets
        # nsize : vertical size of the window
        # msize : horizonal size of the window

        # get the dimensions of the matrix
        [n,m]=np.shape(lp_height)
        self.nsize=3
        self.msize=3

        ntot=self.nsize*self.msize

        # divide matrix to equal subsets
        v=range(0,int(n),self.nsize)
        #h=range(0,m-msize+1) # if steping size is 1
        h=range(0,int(m),self.msize)

        # loop through subsets and calculate the #lp/#hw ratio
        percent_canopy=[]
        prob_of_use=[] # new version
        for i in v[:-1]:
            for j in h[:-1]:
                # subset each matrix
                sub_pl_height=self.subset_matrix(lp_height, i, j)
                sub_lp_count=self.subset_matrix(lp_count, i, j)
                sub_hw_height=self.subset_matrix(hw_height, i, j)
                sub_hw_count=self.subset_matrix(hw_count, i, j)
                # find indexes of the mature LP and HW trees
                # overstory
                # max height the LLP ~ 30-35m, HW ~ 20-25m 
                idl=np.nonzero(sub_pl_height > 20)
                idh=np.nonzero(sub_hw_height > 15)
                # sum number of trees per cell
                num_of_lp=np.sum(sub_lp_count[idl[0][:],idl[1][:]])
                num_of_hw=np.sum(sub_hw_count[idh[0][:],idh[1][:]])
                
                # adding midstory assuming that it each cell that has young trees
                # contributing to 3% to the total canopy cover
                mid=3*midstory(sub_pl_height,sub_hw_height)

                # percent canopy within 9 cells (10 llp, 10 hw)
                # 4 border cell only half way in the 15m circle
                # number of trees produce stable HSI for canopy cover
                #percent_canopy.append((num_of_lp+num_of_hw)/140)
                # here I changed it to percent, assuming that mature trees take about 4.4%
                # which 3.3x3.3m2 of the total area of 9 cells 15x15m2=225m2
                # eliminating corner cells, allowing 175m2, leads to 5.7%~6% area covered by 1 tree
                over=6*(num_of_lp+num_of_hw)
                percent_canopy.append((mid+over)/100)

        # Mean probability burrow status = 'Abandoned' according to the Fig.4 (Catano et al., 2014)
        # Here flipped that formula to the probability of use for consistency with SFS

        self.bench = 0
        self.smooth = 0
        self.sigma = 0.9
        for i in percent_canopy:
            self.number = float(i)
            prob_of_use.append(self.gscore ()-.3)
            #print 'per_can:',round(i,3)
            #print 'prob of use:',self.gscore()-.3
        
        prob=np.asarray(prob_of_use)
        prob[prob<0]=0

        return np.mean(prob)#prob_of_use

    def prob_of_use_gt_new(self, lp_height, lp_count, hw_height, hw_count):
        # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
        # only works when matrix can be equally divided to subsets
        # nsize : vertical size of the window
        # msize : horizonal size of the window

        # get the dimensions of the matrix
        [n,m]=np.shape(lp_height)
        self.nsize=3
        self.msize=3

        ntot=self.nsize*self.msize

        # divide matrix to equal subsets
        v=range(0,int(n),self.nsize)
        #h=range(0,m-msize+1) # if steping size is 1
        h=range(0,int(m),self.msize)

        # loop through subsets and calculate the #lp/#hw ratio
        percent_canopy=[]
        prob_of_use=[] # new version
        for i in v[:-1]:
            for j in h[:-1]:
                # subset each matrix
                sub_pl_height=self.subset_matrix(lp_height, i, j)
                sub_lp_count=self.subset_matrix(lp_count, i, j)
                sub_hw_height=self.subset_matrix(hw_height, i, j)
                sub_hw_count=self.subset_matrix(hw_count, i, j)
                a=np.zeros(ntot)
                a[sub_pl_height.flatten()>1]=1
                b=np.zeros(ntot)
                b[sub_hw_height.flatten()>1]=1
                c=a+b
                c[c==2]=1
                nrand=0.5+abs(0.5-np.random.rand(ntot))
                percent_canopy.append(np.sum(c)/ntot)
                #percent_canopy.append((np.count_nonzero(a+b))/ntot)
																

        # Mean probability burrow status = 'Abandoned' according to the Fig.4 (Catano et al., 2014)
        # Here flipped that formula to the probability of use for consistency with SFS

        #print percent_canopy[-1]
        self.bench = 0
        self.smooth = 0
        self.sigma = 0.9
        for i in percent_canopy:
            self.number = float(i)
            prob_of_use.append(self.gscore ()-.3)
            #print 'per_can:',round(i,3)
            #print 'prob of use:',self.gscore()-.3
        
        prob=np.asarray(prob_of_use)
        prob[prob<0]=0

        return np.mean(prob)#prob_of_use
    
    def __init__(self):
        """HSI score variables"""
        self.x0 = 0         # mean value in GSD
        self.sigma = 2      # standard deviation 
        self.bench = 0      # benchmarking number 
        self.nsize = 5      # subset grid cell area vertical
        self.msize = 5      # subset grid cell area horizontal
        self.smooth = 1     # smoothness factor
        self.number = 0     # number compared against benchmarking number
        self.lb = 2         # lower boundary constraint
        self.ub = 6         # upper boundary constraint



