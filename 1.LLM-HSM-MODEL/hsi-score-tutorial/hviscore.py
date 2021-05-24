#!/usr/bin/env python
import numpy as np
import sys

def norm_gaussian(x0,sigma):
    # mean is 0
    # sigma defines the spread of Gaussian distribution
    #print 'sigma:', sigma
    y0=np.exp(- x0*x0 / (2.0 * sigma*sigma))   
    return y0


def gscore(num, bmn, sm, sigma=2):
    # num: is a current number
    # bn: benchmarking number 
    # sm: smoothing
    # sigma: is a distribution (bigger sigma, wider the distribution)

    r=float(abs(bmn-num)) # absolute distance from 0 mean
    #print 'num',num,' bn',bmn,' r',r,' sm',sm,' sigma',sigma
    if r > sm:
        y0=norm_gaussian(r-sm,sigma)
    else:
        y0=1
    #print 'score:', y0
        
    return y0

def age2dbh(age):
	# Note, need a better formula
	# age [yrs]
	# dbh [cm], the final result is in meters
	
	#dbh=0
	#if (age > 0):
		#Greenberg, Cathryn H.; Simons, Robert W. 1999. Age, composition, and stand structure 
		#of old-growth oak sites in the Florida high pine landscape: implications for ecosystem 
		#management and restoration. Natural Areas Journal. 19(1): 30-40. (Figure 7)
	dbh=2.77+0.31*age+0.004*age*age
	
	# Another equation by James et al., (2004) [see HSI webpage]
	# The problem with this formulation is that at some point older trees 
	# will have smaller dbh
	# dbh=(0.5-(0.817-0.004*age)**2)/0.008
	
	return dbh
	
def dbh2ba(dbh):
	# Ref: https://en.wikipedia.org/wiki/Basal_area
	# Note dbh is in [cm], ba in [m]
	ba = 0.00007854*dbh*dbh
	return ba
	
def subset_matrix(a, i, j, nsize, msize):
    [n,m]=np.shape(a)
    b=a[i:i+nsize,j:j+msize]
    #print b.shape
    return b
    
def score_subset(x, a, b, bn, nsize=5, msize=5, smooth=1):
	# divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
    # only works when matrix can be equally divided to subsets
    # nsize : vertical size of the window
    # msize : horizonal size of the window 
    
    #print  nsize, msize   
    [n,m]=x.shape
    v=range(0,int(n),nsize)
    #h=range(0,m-msize+1) # if steping size is 1
    h=range(0,int(m),msize)
#    try:
#    	raise(np.mod(int(m),msize))
#    except ValueError:
#    	print 'ERROR: matrix horizonal dim:',m,' cannot be div by',msize
#    	raise
    if np.mod(int(m),msize): # NOTE, mod(4,2) = 0, i.e. false 
        print 'ERROR: matrix horizonal dim:',m,' cannot be div by',msize
    if np.mod(int(n),nsize):
        print 'ERROR: matrix vertical dim:',n,' cannot be div by',nsize
    s=len(v)*len(h)
    subset_sc = np.zeros(s, 'float') 
    #print subset_sc
    k=0
    for i in v:
        for j in h:
            #print i,j
            sub=subset_matrix(x, i, j, nsize, msize)
            #print a, b, bn, smooth
            subset_sc[k]=gscore_range(sub, a, b, bn, smooth)
            #print sub
            #print 'score: ',subset_sc[k],'bn: ',bn
            #print '-------------'
            k+=1
    
    sc=np.mean(subset_sc)
    #print 'mean score:',sc
    #print '-------------'
    return sc
           
def gscore_range(x, a, b, bench, smooth):
    # x: a one dimensional array
    # a,b: prescribed range
    # bn: benchmark number corresponds to the ideal number of elements 
    # smooth : relaxes the benchmark constrain 
    # (e.g. benchmark number of pines for RCW is 200 trees per ha)
    x1=np.asarray(x).reshape(-1)
    xx=x1[a<=x1]
    num=len(xx[xx<=b])
    sc=gscore(num, bench, smooth, sigma=1)
    #print 'num:',num,'bn:',bench,'sc:',sc
    return sc
    
def prob_of_use_sq(lp_age, lp_count,hw_age, hw_count, mature_tree_age=5,nsize=5, msize=5):
    # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
    # only works when matrix can be equally divided to subsets
    # nsize : vertical size of the window
    # msize : horizonal size of the window 

    # get the dimensions of the matrix   
    [n,m]=np.shape(lp_age)
    
    # check if matrix can be evenly divided
    if np.mod(int(m),msize):
        raise ValueError('matrix horizonal dim:',m,' cannot be divided by',msize)
    if np.mod(int(n),nsize):
        raise ValueError('ERROR: matrix vertical dim:',n,' cannot be divided by',nsize)  
    
    # divide matrix to equal subsets
    v=range(0,int(n),nsize)
    #h=range(0,m-msize+1) # if steping size is 1
    h=range(0,int(m),msize)

    # loop through subsets and calculate the #lp/#hw ratio
    ratio=[]
    prob_of_use=0.0
    for i in v:
        for j in h:
            # subset each matrix
            sub_pl_age=subset_matrix(lp_age, i, j, nsize, msize)
            sub_lp_count=subset_matrix(lp_count, i, j, nsize, msize)
            sub_hw_age=subset_matrix(hw_age, i, j, nsize, msize)
            sub_hw_count=subset_matrix(hw_count, i, j, nsize, msize)
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
		prob_of_use=gscore(x0, 8, 0, sigma=6)-.2
		#print 'sfs score:',prob_of_use
		if prob_of_use<0:
			prob_of_use=0	

    return prob_of_use
    
#    try:
#        1/len(ratio)
#    except ZeroDivisionError as err:
#        print 'WARNING: number of hardwood trees is equal to ',len(ratio),\
#              '. Consider increasing resolution or HW age constrain.'
#        #raise
    
    
def midstory(pl_age,hw_age,lp_count,hw_count):
        a=pl_age.flatten()
        a[a < 6]=0
        a[a > 50]=0
        a[a != 0]=1

        b=hw_age.flatten()
        b[b < 6]=0
        b[b > 50]=0
        b[b != 0]=1

        return (np.sum(a)+np.sum(b))*2

def prob_of_use_gt(lp_age, lp_count,hw_age, hw_count, mature_tree_age=5,nsize=3, msize=3):
    # divides matrix x(n,m) to subsets of matrices of size (nsize, msize)
    # only works when matrix can be equally divided to subsets
    # nsize : vertical size of the window
    # msize : horizonal size of the window

    # get the dimensions of the matrix
    [n,m]=np.shape(lp_age)

    ntot=nsize*msize

    # divide matrix to equal subsets
    v=range(0,int(n),nsize)
    #h=range(0,m-msize+1) # if steping size is 1
    h=range(0,int(m),msize)

    # loop through subsets and calculate the #lp/#hw ratio
    percent_canopy=[]
    prob_of_use=[] # new version
    for i in v[:-1]:
        for j in h[:-1]:
            # subset each matrix
            sub_pl_age=subset_matrix(lp_age, i, j, nsize, msize)
            sub_lp_count=subset_matrix(lp_count, i, j, nsize, msize)
            sub_hw_age=subset_matrix(hw_age, i, j, nsize, msize)
            sub_hw_count=subset_matrix(hw_count, i, j, nsize, msize)
            # find indexes of the mature LP and HW trees
            # overstory
            idl=np.nonzero(sub_pl_age > mature_tree_age)
            idh=np.nonzero(sub_hw_age > mature_tree_age)
            # sum number of trees per cell
            num_of_lp=np.sum(sub_lp_count[idl[0][:],idl[1][:]])
            num_of_hw=np.sum(sub_hw_count[idh[0][:],idh[1][:]])
            over=5.7*(num_of_lp+num_of_hw)

            # adding midstory assuming that it each cell that has young trees
            # contributing to 2% to the total canopy cover
            mid=midstory(sub_pl_age,sub_hw_age,sub_lp_count,sub_hw_count)

            percent_canopy.append((mid+over)/100)
            # percent canopy within 9 cells (10 llp, 10 hw)
            # 4 border cell only half way in the 15m circle
            # number of trees produce stable HSI for canopy cover
            #percent_canopy.append((num_of_lp+num_of_hw)/140)
            # here I changed it to percent, assuming that mature trees take about 4.4%
            # which 3.3x3.3m2 of the total area of 9 cells 15x15m2=225m2
            # eliminating corner cells, allowing 175m2, leads to 5.7% areas covered by 1 tree
            #percent_canopy.append(5.7*(num_of_lp+num_of_hw)/100)

    # Mean probability burrow status = 'Abandoned' according to the Fig.4 (Catano et al., 2014)
    # Here flipped to the probability of use
    # old version
    # x0=np.mean(np.asarray(percent_canopy))
    # prob_of_use=gscore(x0, 0, 0, sigma=0.9)-.3
    # new version
    [prob_of_use.append(gscore(float(i), 0, 0, sigma=0.9)-.3)for i in percent_canopy]
    prob=np.asarray(prob_of_use)
    prob[prob<0]=0

    return np.mean(prob)#prob_of_use