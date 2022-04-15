# This a master coupling script, that streamlines all the components of the DRM framework
#
# (c) Elchin Jafarov 03/30/2021

import numpy as np
import matplotlib.pyplot as plt
import time
import pandas as pd
import random
import sys,os
import os.path
import subprocess
from shutil import copyfile
from subprocess import call
import shutil
sys.path.insert(0, '1.LLM-HSM-MODEL/')
import LLM_model_class as llm
import LLM_FT_utils as llmft
import hsiscore_class as HSI
import hsi_plot_utils as hsi_plt
import LLM_display
sys.path.insert(0, '7.QUICFIRE-MODEL/projects/Tester')
import postfuelfire_new as pff
import Buffer as buff

def LLMspinup(nyears):
    # --spinup run ---
    p = llm.LLM()     # assign p to the llm class
    p.dim = 80
    p.instantiate(0)  # 1: reads input data from file, 0: generate inputs internally
    p.readfireprobfromfile=0
    p.readmastprobfromfile=0
    p.verbose=0

    start_time = time.time()
    p.run(nyears) # here 200 is a number of years
    print("--- %s seconds ---" % (time.time() - start_time))
    p.save_pickle() #saves the results
    
    return
   
def LLMtransient(nyears):
    # --transient run ---
    #del p
    p = llm.LLM() 
    p.dim = 80
    p.randfromfile=0
    p.instantiate(1)  # 1: reads input data from file, 0: generate inputs internally
    p.verbose=0       # 0: do not print out scores, 1: print scores on the screen
    p.tree_mature_age=10
    p.readfireprobfromfile=0
    p.readmastprobfromfile=0

    start_time = time.time()
    p.run(nyears)
    print("--- %s seconds ---" % (time.time() - start_time))
    if np.sum(p.litter)==0:
        print('no litter, make an extra run...')
        p.fire_prob=0
        p.run(1)
    #p.save_randfromfile()
    return p

def dbh_cr(p):
    lp_count=p.old_LPcount.copy()
    lp_count[p.old_ht<1.37]=0

    lp_height=p.old_ht.copy()
    lp_height[lp_height<1.37]=1.37
    p.lp_dbh=llmft.dbh1_model(lp_height)
    #print dbh
    #axx=llmft.plot_area_matrix(p.lp_dbh,'LLP tree DBH [cm]','yes')

    lp_CA=llmft.dbh2_model(lp_height,p.lp_dbh)
    p.lp_CR=np.sqrt(lp_CA/np.pi) #LLP Crown Area
    all_NaNs = np.isnan(p.lp_CR)
    p.lp_CR[all_NaNs] = 0
    #axx=llmft.plot_area_matrix(p.lp_CR,'LLP tree CR [m]','yes')

    hw_height=p.old_htHW.copy()
    hw_height[hw_height<1.37]=1.37
    p.hw_dbh=llmft.dbh1_model(hw_height)
    #print dbh
    #axx=llmft.plot_area_matrix(p.hw_dbh,'HW tree DBH [cm]','yes')

    hw_CR=llmft.dbh2cr_hw(p.hw_dbh/2.54) # note dbh is in inch
    p.hw_CR=hw_CR/3.281            # CR is in feet convert to meters
    all_NaNs = np.isnan(p.hw_CR)
    p.hw_CR[all_NaNs] = 0
    #axx=llmft.plot_area_matrix(p.hw_CR,'HW tree CR [m]','yes') 
    
    return p

def savelittersLLMQF(p,i):
    filename='LLM2FT/LLM_litter_WG.dat'
    ftitle='WG litter [kg/4m2]'
    llmft.save_litter_LLM_FT(filename,ftitle,p.litterWG,'plot','grass')
    newname = 'litter_WG.' + str(i) + '.png'
    os.rename('litter.png', newname)

    filename='LLM2FT/LLM_litter_trees.dat'
    ftitle='LLP + HW litter [kg/4m2]'
    tree_litter=p.litterHW+p.litter
    llmft.save_litter_LLM_FT(filename,ftitle,tree_litter,'plot','litter')
    newname = 'litter_Tree.' + str(i) + '.png'
    os.rename('litter.png', newname)

    percent_LP_litter=np.sum(p.litter)/np.sum(p.litterHW+p.litter)
    percent_HW_litter=np.sum(p.litterHW)/np.sum(p.litterHW+p.litter)
    print ('lit_LLP%, lit_HW%:',percent_LP_litter,percent_HW_litter)
    
    return

def runTreeQF():
# Note: Adam has a QF Tree code in '5.TREES-QUICFIRE'
    src='LLM2FT/'
    dst='../5.TREES-QUICFIRE/'
    print(os.getcwd())
    copyfile(src+'LLM_litter_WG.dat',dst+'LLM_litter_WG.dat')
    copyfile(src+'treelist_LLM.dat',dst+'treelist_LLM.dat')
    copyfile(src+'LLM_litter_trees.dat',dst+'LLM_litter_trees.dat')

    os.chdir(dst)
    status=subprocess.call(["./trees"])
    while status != 0:
       print('Tree program failed to execute...')
       status=subprocess.call(["./trees"])

    if status==0:
        print('Tree program run successfully!')
        ### Copying Tree Files to Fire Affects Assessment
        copyfile('TreeTracker.txt','../8.CROWN-SCORCH/TreeTracker.txt')
        copyfile('treelist_LLM.dat','../8.CROWN-SCORCH/treelist_LLM.dat')
        copyfile('LLM_litter_WG.dat','../8.CROWN-SCORCH/LLM_litter_WG.dat')
        copyfile('LLM_litter_trees.dat','../8.CROWN-SCORCH/LLM_litter_trees.dat')

    return

def runQF(i): 
    #copy produced by Tree program files to the QF folder
    #os.chdir("/Users/elchin/Documents/Adams_project/llm-hsm-ft/")
    print(os.getcwd())
    src=''
    dst='../7.QUICFIRE-MODEL/projects/ftFiles/'
    copyfile(src+'treesfueldepth.dat',dst+'treesfueldepth.dat')
    copyfile(src+'treesmoist.dat',dst+'treesmoist.dat')
    copyfile(src+'treesrhof.dat',dst+'treesrhof.dat')
    copyfile(src+'treesss.dat',dst+'treesss.dat')
    
    #Before running the shell script make sure that 
    # file path '/projects/Tester/QUIC_fire.inp' 
    # is pointing to the gridlist file in the projects/ftFiles directory.
    os.chdir("../7.QUICFIRE-MODEL/mac_compile/")
    import subprocess
    status=subprocess.call(["./compile_and_run.sh"])
    if status==0:
        print('QF run successfully!')
    else:
        print('QF failed to execute...')
    #Successful run should produce bunch of binary files in 
    #7.QUICFIRE-MODEL/projects/Tester. Now run the postfire script 
    #that will generate PercentFuelChange.txt file required for the next step.
    os.chdir("../projects/Tester")
    pff.main(200) #ASK ADAM! # 2200 seconds for a normal run, 200 sec for a quick-test
    # MAtch this value at Line 5 of 7.QUICFIRE-MODEL/projects/Tester/QUIC_fire.inp
    direc = "Plots"
    dd = direc + str(i)
    print (os.getcwd())
    if os.path.exists(dd):
       shutil.rmtree(dd)
    os.rename('Plots', dd)
    os.mkdir('Plots')
    
    return

def runCrownScorch():
    os.chdir("../../../8.CROWN-SCORCH")
    copyfile('../7.QUICFIRE-MODEL/projects/Tester/PercentFuelChange.txt','../8.CROWN-SCORCH/PercentFuelChange.txt')
    LiveDead = []
    file_names = ['PercentFuelChange.txt', 
                  'TreeTracker.txt', 
                  'treelist_LLM.dat',
                  'LLM_litter_WG.dat', 
                  'LLM_litter_trees.dat', 
                  '../1.LLM-HSM-MODEL/FT2LLM/AfterFireTrees.txt', 
                  '../1.LLM-HSM-MODEL/FT2LLM/AfterFireWG.txt',
                  '../1.LLM-HSM-MODEL/FT2LLM/AfterFireLitter.txt']

    for i in range(len(file_names)-3):
        # check if all input files exist
        llmft.check_file_exists(file_names[i])

    LiveDead = llmft.Treeoflife(file_names)
    
    return LiveDead
    
def runLLMcyclical(p,nyears):
    os.chdir("../1.LLM-HSM-MODEL")
    flitter='FT2LLM/AfterFireLitter.txt'
    fwg='FT2LLM/AfterFireWG.txt'
    ftlist='FT2LLM/AfterFireTrees.txt'
    p=llmft.read_FT_2_LLM(flitter,fwg,ftlist,p)

    #run the LLM-HSI for nyears years
    p.fire_prob=0
    start_time = time.time()
    p.run(nyears)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    return p

def updateTreelist(p,ii):
    ftlist='FT2LLM/AfterFireTrees.txt'
    [lp_list,hw_list]=llmft.update_tree_info_per_location(p,ftlist,0)

    df_hw = pd.DataFrame(hw_list)
    df = pd.DataFrame(lp_list)
    df=df.append(df_hw)
    df.plot(subplots=True, layout=(4,2),figsize=(12, 10));
    df.to_csv('treelist_LLM.dat', sep=' ',header=False,index=False)
    print(os.getcwd())
    file_in='treelist_LLM.dat'
    file_out='LLM2FT/treelist_LLM.dat'
    llmft.save_FT_treelist(file_in,file_out,0)
    
    df = pd.read_csv('LLM2FT/treelist_LLM.dat',sep=' ',
                names=["Tree id","x coord [m]","y coord [m]","Ht [m]",
                      "htlc [m]","CRDiameter [m]","hmaxcr [m]",
                      "canopydensity  [kg/m3]", "CR fuel moist [frac]",
                      "CR fuel size scale [m]" ])
    df.plot(subplots=True, layout=(5,2),figsize=(12, 10));
    plt.tight_layout()
    plt.savefig('TreeData.png')

    plt.figure(figsize=(8, 6))
    plt.plot(df["x coord [m]"].values,df["y coord [m]"].values,'.')
    plt.title('Tree distribution in the FT domain');   
    print("Total number of trees: ",df["x coord [m]"].size )
    newname = 'TreeMap.' + str(ii) + '.png'
    plt.savefig('TreeMap.png')
    os.rename('TreeMap.png', newname)

    return

#-----main------
#
nyears=3      # number of years for spinup and transient runs
ncycyear=2    # number of cyclical year run
ncycle=2      # number of loops

#Build Trees
os.chdir('5.TREES-QUICFIRE')
ierr = call('make', shell=True)
os.chdir('../1.LLM-HSM-MODEL')

LLMspinup(nyears)          # temporary llm class
llm=LLMtransient(nyears)   # permanent llm class
llm=dbh_cr(llm)            # calculates dbh and crown radius 
savelittersLLMQF(llm,0)
llmft.create_treelist(llm,'LLM2FT/treelist_LLM.dat')

#### MAKE INTO FUNCTION
df = pd.read_csv('LLM2FT/treelist_LLM.dat',sep=' ',
                    names=["Tree id","x coord [m]","y coord [m]","Ht [m]",
                        "htlc [m]","CRDiameter [m]","hmaxcr [m]",
                        "canopydensity  [kg/m3]", "CR fuel moist [frac]",
                        "CR fuel size scale [m]" ])
df.plot(subplots=True, layout=(5,2),figsize=(12, 10));
plt.tight_layout()
plt.savefig('TreeInfo.png')

plt.figure(figsize=(8, 6))
plt.plot(df["x coord [m]"].values,df["y coord [m]"].values,'.')
plt.title('Tree distribution in the FT domain');
print("Total number of trees: ",df["x coord [m]"].size )
plt.savefig('TreePlot.0.png')

hsi_plt.plot_species_scores(llm)
plt.savefig('HVI.0.png')
#### MAKE ABOVE INTO FUNTION

## Change Coordinates for QUICFIRE HERE ###
#buff.add_sat()
#buff.add_tree_buff()
#buff.add_surf_buff()

LiveDead=[]
for i in range(ncycle):
    ii = i + 1
    runTreeQF()                       # runs the tree program to create QF inputs
    runQF(i)                           # runs QuickFire
    L=np.array(runCrownScorch())                  # runs the tree program to create LLM inputs
    L=np.insert(L,0,ii)
    LiveDead.append(L)    
    ## Change Coordinates Back to Eco system model HERE ###
    #buff.remove_tree_buff()
    #buff.remove_surf_buff()
    print('Loop Number: ',i)
    llm=runLLMcyclical(llm,ncycyear)  # runs LLM-HSM with no fire 
    hsi_plt.plot_species_scores(llm)  # Plotting HVI
    plt.savefig('HVI.png')
    sc_rcw=np.asarray(llm.age_sc)+np.asarray(llm.hw_sc)+np.asarray(llm.ageHW_sc)+np.asarray(llm.hwHW_sc)
    savelittersLLMQF(llm,ii)
    updateTreelist(llm,ii)               # this also updates dbh and cr 
    ## Change Coordinates for QUICFIRE HERE ###
    #buff.add_tree_buff()
    #buff.add_surf_buff()    
    dd = 'HVI.' + str(i) + '.png'
    print (dd)
    #plt.savefig('HVI.png') 
    os.rename('HVI.png', dd)

LiveDead=np.array(LiveDead)
np.savetxt('LiveDead.txt',LiveDead,fmt='%i',header='Fire LLP(L/D) Turk(L/D)')
print ('ADAM SQ', llm.sq_sc)
print ('ADAM GT', llm.gt_sc)
print ('ADAM RCW',sc_rcw)
np.savetxt('HVI-score.txt', np.c_[sc_rcw, llm.sq_sc, llm.gt_sc], fmt='%1.4e')        
