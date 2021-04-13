#!/usr/bin/env python
import numpy as np
import random
import pandas as pd
import scipy.interpolate.ndgriddata as ndgriddata
import matplotlib.pyplot as plt

def regrid_LLM2FT(nx,ny,nx1,ny1,imethod="linear",data=[]):
    # methods: nearest, linear cubic 
    
    if data==[]:
        data=np.random.rand(nx*ny)
        data_reshape=np.reshape(data,(nx,ny))
    
    x_old = np.linspace(1.0, 5*nx, nx) # every cell is 25 m2
    y_old = np.linspace(1.0, 5*ny, ny) # total 25 cells in each direction
    X_old, Y_old = np.meshgrid(x_old, y_old)

    x = np.linspace(1.0, 2*nx1, nx1) # every cell is 4 m2
    y = np.linspace(1.0, 2*ny1, ny1) # total 25 cells in each direction
    X, Y = np.meshgrid(x, y)

    data_new=ndgriddata.griddata((X_old.flatten(), Y_old.flatten()), data.flatten(), (X, Y), method=imethod)
     
    return data_new

def regrid_FT2LLM(nx,ny,nx1,ny1,imethod="linear",data=[]):
    # methods: nearest, linear cubic 
    
    if data==[]:
        data=np.random.rand(nx*ny)
        data_reshape=np.reshape(data,(nx,ny))
    
    x_old = np.linspace(1.0, 2*nx, nx) # every cell is 25 m2
    y_old = np.linspace(1.0, 2*ny, ny) # total 25 cells in each direction
    X_old, Y_old = np.meshgrid(x_old, y_old)

    x = np.linspace(1.0, 5*nx1, nx1) # every cell is 4 m2
    y = np.linspace(1.0, 5*ny1, ny1) # total 25 cells in each direction
    X, Y = np.meshgrid(x, y)

    data_new=ndgriddata.griddata((X_old.flatten(), Y_old.flatten()), data.flatten(), (X, Y), method=imethod)
     
    return data_new

def plot_area_matrix(xx,title,pbar='yes'):
    fig, ax = plt.subplots(figsize=(14, 8))
    im = ax.imshow(xx, interpolation = 'nearest', cmap='jet')#, vmin=0, vmax=200) 
    if pbar=='yes':
        cbar=fig.colorbar(im, ax=ax, extend='both')
        cbar.ax.tick_params(labelsize=20) 
    plt.title(title,fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    #plt.gridon()
    return ax
    
def count_tress(tt,llm_model):
    x0=2.5
    y0=2.5
    
    [n,m]=llm_model.shape
    tree_count=np.zeros((n,m))
    [n_tree,m_tree]=np.shape(tt)
    print ('number FT of trees in treelist: ',n_tree)

    for i in range(n+1):#n
        x=x0+i*5
        for j in range(m+1):#m
            y=y0+j*5
            for ii in range(n_tree):#n_tree
                xx=tt[ii][1]
                yy=tt[ii][2]
                if (xx>x-2.5) and (xx<=x+2.5):
                    if (yy>y-2.5) and (yy<=y+2.5):
                        tree_count[i,j]=tree_count[i,j]+1
    print ('total LLM tree count: ',np.sum(tree_count))

    return tree_count


def dbh1_model(H):   #  example: n=500, alpha=2
    # Gonzalez-Benecke et al., (2014), Table 3
    #a1,a2: parameter estimates
    a1=0.4031
    a2=0.97
    
    y = a1+a2*np.log(H-1.37)      
    dbh=np.exp(y)
    
    return dbh

def dbh2_model(H,dbh):   #  example: n=500, alpha=2
    # Gonzalez-Benecke et al., (2014), Table 3
    #a1,a2,a3: parameter estimates
    a1=0.342191
    a2=0.804174
    a3=0.185289
    
    y = np.log(dbh)-a1-a2*np.log(H-1.37)   
    lnCA=y/a3
    CA=np.exp(lnCA)
        
    return CA

def dbh_hw_model(H):   #  example: n=500, alpha=2
    # Lynch et al., (2005), Table 2
    #a1,a2: parameter estimates for HWs
    beta0=3.7997
    beta1=-14.917
    power=-10/9
    
    y = (np.log(H-1.37)-beta0)/beta1      
    dbh=y**power
    
    return dbh

def dbh2cr_hw(dbh):   #  example: n=500, alpha=2
    # Goelz (1996), Table 2
    #lnb1,b2: parameter estimates for Willow Oak
    lnb1=0.342191
    b2=0.804174
    b1=np.exp(lnb1)
    y = b1*dbh**b2   
        
    return y

def FT_LLM_FT_tlist(tid,tlist,new_count,CR,dbh,ht):
    x0=2.5
    y0=2.5
    
    [n,m]=new_count.shape
    tree_count=np.zeros((n,m))
    #[n,m]=[2,2]
    [n_tree,m_tree]=np.shape(tlist)
    #[n,m]=new_count.shape
    #tree_count=np.zeros((n,m))
    
    #print (np.shape(tlist))
    #print (np.shape(new_count))
    
    tlist_array=np.asarray(tlist)
    
    for i in range(n):#n
        x=x0+i*5
        for j in range(m):#m
            y=y0+j*5
            index=[]
            for ii in range(n_tree):#n_tree
                xx=tlist[ii][1]
                yy=tlist[ii][2]
                if (xx>x-2.5) and (xx<=x+2.5):
                    if (yy>y-2.5) and (yy<=y+2.5):
                        tree_count[i,j]=tree_count[i,j]+1
                        index.append(ii)
            #print (x,y,index)
            #print ('old:',tree_count[i,j],'new',new_count[i,j])
            
            for ii in np.asarray(index):
                tlist[ii][3]=CR[i,j]
                tlist[ii][4]=dbh[i,j]
                tlist[ii][5]=ht[i,j]      
            
            if tree_count[i,j]>new_count[i,j]:
                #print ('difference:',tree_count[i,j]-new_count[i,j])
                for ii in np.asarray(index):
                    #print ('ii',ii) 
                    tlist[ii][0]=-999
                    
            elif tree_count[i,j]<new_count[i,j]:
                line=np.round([CR[i,j],dbh[i,j],ht[i,j]],3)
                #print (line)
                xx=random.uniform(-2.5,2.5)+x
                yy=random.uniform(-2.5,2.5)+y
                line0=np.round([tid,xx,yy],3)
                tlist.append(np.hstack((line0,line)) )
    return tlist

def save_litter_LLM_FT(filename,ftitle,litter,fplot):
    [nx,ny]=litter.shape
    print ('shape of the litter matrix:',nx,ny)
    new_litter=regrid_LLM2FT(nx,ny,200,200,"linear",litter);
    if fplot=='plot':
        axx=plot_area_matrix(new_litter,ftitle,'yes')
    print ('sum of the old litter matrix:',litter.sum())
    print ('sum of the old litter matrix:', new_litter.sum())
    np.savetxt(filename,  (1/(4*1.5))*new_litter, fmt='%.2f')
    print ('litter file:',filename,' saved!')
    #checking the size of the saved matrix
    #a=np.loadtxt(filename)
    #print (a.shape)

def read_FT_2_LLM(file_litter,file_wg,file_tlist,pp,graph=0):
    # converting WG and tree litters and tree to LLM
    # updating litters and tree counts in the LLM from FT
    aflitter=np.loadtxt(file_litter)
    aflitter=4*1.5*aflitter
    if graph:
        axx=plot_area_matrix(aflitter,'LLP tree litter kg/4m2','yes') 
    [nx,ny]=aflitter.shape
    print ('FT tree litter size :',nx,ny)
    old_tot_litter=regrid_FT2LLM(nx,ny,80,80,"linear",aflitter);
    if graph:
        axx=plot_area_matrix(old_tot_litter,'LLP tree litter [kg/25m2]','yes')
    print ('FT tree litter sum (2m2 res):',aflitter.sum())
    print ('LLM tree litter sum (5m2 res):',old_tot_litter.sum())
    
    afwg=np.loadtxt(file_wg)
    afwg=4*1.5*afwg
    if graph:
        axx=plot_area_matrix(afwg,'LLP tree litter kg/4m2','yes') 
    [nx,ny]=afwg.shape
    print ('FT WG litter size :',nx,ny)
    old_WG=regrid_FT2LLM(nx,ny,80,80,"linear",afwg);
    if graph:
        axx=plot_area_matrix(old_WG,'LLP tree litter [kg/25m2]','yes')
    print ('FT WG litter sum (2m2 res):',afwg.sum())
    print ('LLM WG litter sum (5m2 res):',old_WG.sum())
    

    # saving to LLM, need to double check, here I assume that HW litter is only 5% from total
    pp.litterWG=old_WG
    #0.993381053918314 0.0066 percent difference
    pp.litterHW=old_tot_litter*0.006
    pp.litter=old_tot_litter*0.993

    #converting treelist to treenumber 
    treelist=np.loadtxt(file_tlist)
    LLP_trees=treelist[treelist[:,0]==1]
    HW_trees=treelist[treelist[:,0]==2]

    tt_llp=LLP_trees.tolist()
    tt_hw=HW_trees.tolist()

    LLPcount=count_tress(tt_llp,pp.old_LPcount);
    HWcount=count_tress(tt_hw,pp.old_LPcount);

    pp.old_LPcount=LLPcount
    pp.old_HWcount=HWcount
    
    return pp

def save_FT_treelist(fin,fout,graph):
    data=np.loadtxt(fin)
    [n,m]=data.shape
    htlc=np.zeros(n)
    hmaxcr=np.zeros(n)
    canopydensity=np.zeros(n)
    newdata=np.zeros((n,8))

    htlc=data[:,-1]*0.7649-4.1628
    htlc[htlc<0]=0
    htlc[data[:,0]==2]=0
    hmaxcr=htlc+.07
    hmaxcr[data[:,0]==2]=0
    canopydensity[data[:,0]==2]=0.6
    canopydensity[data[:,0]==1]=0.2

    newdata[:,0:3]=data[:,0:3] #tid,x,y
    newdata[:,3]=data[:,-1] #height
    newdata[:,4]=htlc
    newdata[:,5]=2*data[:,3] #Cdimater
    newdata[:,6]=hmaxcr 
    newdata[:,7]=canopydensity

    df_new = pd.DataFrame(newdata)
    if graph:
        df_new.plot(subplots=True, layout=(4,2),figsize=(12, 10));
    df_new.to_csv(fout, sep=' ',header=False,index=False)
    print (df_new.shape)

def update_tree_info_per_location(pp,ftreelist,graph):
    # updating CR and dbh for LLPs and HWs
    lp_height=pp.old_ht.copy()
    lp_height[lp_height<1.37]=1.37
    lp_dbh=dbh1_model(lp_height)
    #print dbh
    if graph:
        axx=plot_area_matrix(lp_dbh,'LLP tree DBH [cm]','yes')

    lp_CA=dbh2_model(lp_height,lp_dbh)
    lp_CR=np.sqrt(lp_CA/np.pi) #LLP Crown Area
    all_NaNs = np.isnan(lp_CR)
    lp_CR[all_NaNs] = 0
    if graph:
        axx=plot_area_matrix(lp_CR,'LLP tree CR [m]','yes')

    hw_height=pp.old_htHW.copy()
    hw_height[hw_height<1.37]=1.37
    hw_dbh=dbh1_model(hw_height)
    #print dbh
    if graph:
        axx=plot_area_matrix(hw_dbh,'HW tree DBH [cm]','yes')

    hw_CR=dbh2cr_hw(hw_dbh/2.54) # note dbh is in inch
    hw_CR=hw_CR/3.281            # CR is in feet convert to meters
    all_NaNs = np.isnan(hw_CR)
    hw_CR[all_NaNs] = 0
    if graph:
        axx=plot_area_matrix(hw_CR,'HW tree CR [m]','yes') 
    
    #count only trees higher than 1.37m and dbh>1
    lp_count=pp.old_LPcount.copy()
    lp_count[pp.old_ht<1.37]=0
    lp_count[lp_dbh<1.001]=0
    if graph:
        axx=plot_area_matrix(lp_count,'LLP tree number ','yes') 
    print ('number of LLPs:',np.sum(lp_count))

    hw_count=pp.old_HWcount.copy()
    hw_count[pp.old_htHW<1.37]=0
    hw_count[hw_dbh<1.001]=0
    if graph:
        axx=plot_area_matrix(hw_count,'LLP tree number ','yes') 
    print ('number of HWs:',np.sum(hw_count))
    
    treelist=np.loadtxt(ftreelist)

    LLP_trees=treelist[treelist[:,0]==1]
    HW_trees=treelist[treelist[:,0]==2]
    #id, x,y, lp_CR[i,j],lp_dbh[i,j],p.old_ht[i,j]
    tt_llp=LLP_trees[:,0:6].tolist()
    tt_hw=HW_trees[:,0:6].tolist()

    print ('tt_llp:',np.shape(tt_llp))    
    lplist=FT_LLM_FT_tlist(1,tt_llp,lp_count,lp_CR,lp_dbh,pp.old_ht);
    print ('lplist:',np.shape(lplist))

    print ('tt_hw:',np.shape(tt_hw))    
    hwlist=FT_LLM_FT_tlist(2,tt_hw,hw_count,hw_CR,hw_dbh,pp.old_htHW);
    print ('hwlist:',np.shape(hwlist))

    #remove dead trees
    lp_newlist=[]
    for i in range(len(lplist)):
        if lplist[i][0]!=-999:
            lp_newlist.append(lplist[i])        
    print ('lp_newlist',np.shape(lp_newlist))

    hw_newlist=[]
    for i in range(len(hwlist)):
        if hwlist[i][0]!=-999:
            hw_newlist.append(hwlist[i])        
    print ('hw_newlist',np.shape(hw_newlist))
    
    return lp_newlist,hw_newlist

def create_treelist(p,filename):
    # Creates the treefile for Tree code and provides a quick look at the prduced dataset
    # NOTE: this step is required only when we first time create treelist, then we used the 
    #       same treelist.

    lp_count=p.old_LPcount.copy()
    lp_count[p.old_ht<1.37]=0
    #counting HWs
    hw_count=p.old_HWcount.copy()
    hw_count[p.old_htHW<1.37]=0

    x0=2.5 # start in the [0,0,5,5] cell
    y0=2.5
    llp_xy=[]
    hw_xy=[]

    [n,m]=lp_count.shape

    for i in range(n):
        x=x0+i*5
        for j in range(m):
            y=y0+j*5
            if lp_count[i,j]!=0: #add long-leaf pines
                if p.lp_dbh[i,j]>1:
                    line=np.round([p.lp_CR[i,j],p.lp_dbh[i,j],p.old_ht[i,j]],3)
                    for item in range(int(lp_count[i,j])):
                        xx=random.uniform(-2.5,2.5)+x
                        yy=random.uniform(-2.5,2.5)+y
                        line0=np.round([1,xx,yy],3)
                        llp_xy.append( np.hstack((line0,line)) )

            if hw_count[i,j]!=0: #add hardwoods
                if p.hw_dbh[i,j]>1:
                    line=[p.hw_CR[i,j].round(2),p.hw_dbh[i,j].round(2),p.old_htHW[i,j].round(2)]
                    for item in range(int(hw_count[i,j])):
                        xx=random.uniform(-2.5,2.5)+x
                        yy=random.uniform(-2.5,2.5)+y
                        line0=np.round([2,xx,yy],3)
                        hw_xy.append( np.hstack((line0,line)) )
    print ('shape of llp (x,y):',np.shape(llp_xy))
    print ('shape of hw (x,y):',np.shape(hw_xy))

    df_hw = pd.DataFrame(hw_xy)
    df = pd.DataFrame(llp_xy)
    df=df.append(df_hw)
    df.to_csv(filename, sep=' ',header=False,index=False)

    data=np.loadtxt(filename)
    [n,m]=data.shape
    htlc=np.zeros(n)
    hmaxcr=np.zeros(n)
    canopydensity=np.zeros(n)
    newdata=np.zeros((n,10))

    htlc=data[:,-1]*0.7649-4.1628
    htlc[htlc<0]=0
    htlc[data[:,0]==2]=0
    hmaxcr=htlc+.07
    hmaxcr[data[:,0]==2]=0
    canopydensity[data[:,0]==2]=0.6
    canopydensity[data[:,0]==1]=0.2

    newdata[:,0:3]=data[:,0:3]
    newdata[:,3]=data[:,-1]
    newdata[:,4]=htlc
    newdata[:,5]=2*data[:,3] #Cdimater
    newdata[:,6]=hmaxcr #Cdimater
    newdata[:,7]=canopydensity
    newdata[:,8]=np.ones(n)
    newdata[:,9]=np.ones(n)*0.000347222
    
    df_new = pd.DataFrame(newdata)
    df_new.to_csv(filename, sep=' ',header=False,index=False)
    print (filename,'is created!')
    
    return