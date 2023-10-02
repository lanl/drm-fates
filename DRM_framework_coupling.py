
# This is a master coupling script, that streamlines all the components of the DRM framework

# (c) Elchin Jafarov 03/30/2021

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import time
import pandas as pd
import random
import sys,os
import os.path
import subprocess
from shutil import copyfile
from subprocess import call
import shutil
import yaml
sys.path.insert(0, '1.LLM-HSM-MODEL/')
import LLM_model_class as llm
import LLM_FT_utils as llmft
import hsiscore_class as HSI
import hsi_plot_utils as hsi_plt
import LLM_display
sys.path.insert(0, '7.QUICFIRE-MODEL/projects/Tester')
#TBK import postfuelfire_new as pff
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
    p.run(nyears)   # here 200 is a number of years
    print("--- %s seconds ---" % (time.time() - start_time))
    p.save_pickle() # saves the results
    
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

def dbh_cr(p): #CSXM: LLM only
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

def savelittersLLMQF(p,i): #CSXM: LLM only
    filename='VDM2FM/VDM_litter_WG.dat'
    ftitle='WG litter [kg/4m2]'
    llmft.save_litter_LLM_FT(filename,ftitle,p.litterWG,'plot','grass')
    newname = 'litter_WG.' + str(i) + '.png'
    os.rename('litter.png', newname)

    filename='VDM2FM/VDM_litter_trees.dat'
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
    if VDM == "LLM":
       VDM_folder = "1.LLM-HSM-MODEL"
    elif VDM == "FATES":
       VDM_folder = "1.FATES-MODEL"
    src='../'+VDM_folder+'/VDM2FM/'
    dst='../5.TREES-QUICFIRE/'
    print(os.getcwd())
    copyfile(src+'VDM_litter_WG.dat',dst+'VDM_litter_WG.dat')
    copyfile(src+'treelist_VDM.dat',dst+'treelist_VDM.dat')
    copyfile(src+'VDM_litter_trees.dat',dst+'VDM_litter_trees.dat')

    os.chdir(dst)
    status=subprocess.call(["./trees"])
    while status != 0:
       print('Tree program failed to execute...')
       status=subprocess.call(["./trees"])

    if status==0:
        print('Tree program run successfully!')
        ### Copying Tree Files to Fire Affects Assessment
        copyfile('TreeTracker.txt','../8.CROWN-SCORCH/TreeTracker.txt')
        copyfile('treelist_VDM.dat','../8.CROWN-SCORCH/treelist_VDM.dat')
        copyfile('VDM_litter_WG.dat','../8.CROWN-SCORCH/VDM_litter_WG.dat')
        copyfile('VDM_litter_trees.dat','../8.CROWN-SCORCH/VDM_litter_trees.dat')

    return

def runQF(i): 
    # copy produced by Tree program files to the QF folder
    #os.chdir("/Users/elchin/Documents/Adams_project/llm-hsm-ft/")
    src=''
    dst='../7.QUICFIRE-MODEL/projects/ftFiles/'
    copyfile(src+'treesfueldepth.dat',dst+'treesfueldepth.dat')
    copyfile(src+'treesmoist.dat',dst+'treesmoist.dat')
    copyfile(src+'treesrhof.dat',dst+'treesrhof.dat')
    copyfile(src+'treesss.dat',dst+'treesss.dat')
    
    # Before running the shell script make sure that 
    # file path '/projects/Tester/QUIC_fire.inp' 
    # is pointing to the gridlist file in the projects/ftFiles directory.
    os.chdir("../7.QUICFIRE-MODEL/mac_compile/")
    # ASXM (BGN)
    # https://www.geeksforgeeks.org/how-to-search-and-replace-text-in-a-file-in-python/
    if i<= 1: # only needed for i = 0 and 1
        with open("adv_compile_and_run.sh", "r") as file:
            data_adv_compile_and_run = file.read()
            if i == 0:
                data_adv_compile_and_run = data_adv_compile_and_run.replace("compile=0 run=1", "compile=1 run=1")
            if i == 1:
                data_adv_compile_and_run = data_adv_compile_and_run.replace("compile=1 run=1", "compile=0 run=1")
        with open("adv_compile_and_run.sh", "w") as file:
            file.write(data_adv_compile_and_run)
    # ASXM (END)
    import subprocess 
    status=subprocess.call(["./adv_compile_and_run.sh"])
    if status==0:
        print('QF run successfully!')
    else:
        print('QF failed to execute...')
    # Successful run should produce bunch of binary files in 
    # 7.QUICFIRE-MODEL/projects/Tester. Now run the postfire script 
    # that will generate PercentFuelChange.txt file required for the next step.
    os.chdir("../projects/Tester")
    # Outputs for paraview plotting 
    #python3 drawfire.py /usr/projects/higrad/rutuja/E3SM_cases/proj1/7.QUICFIRE-MODEL/projects/Tester/ 1 0
    # Match this value at Line 5 of 7.QUICFIRE-MODEL/projects/Tester/QUIC_fire.inp
    direc = "Plots"
    dd = direc + str(i)
    if os.path.exists(dd):
        shutil.rmtree(dd)
    os.rename('Plots', dd)
    os.mkdir('Plots')
    dd = "fuels-dens-00000." + str(i) + ".vin" #CSXM: vin is the renamed bin
    os.rename('fuels-dens-00000.bin', dd)
    dd = "fire_indexes." + str(i) + ".vin"     #CSXM: vin is the renamed bin
    os.rename('fire_indexes.bin', dd)

    return

def runCrownScorch(ii):
    os.chdir("../../../8.CROWN-SCORCH")
    copyfile('../7.QUICFIRE-MODEL/projects/Tester/PercentFuelChange.txt','../8.CROWN-SCORCH/PercentFuelChange.txt')
    LiveDead = []
    if VDM == "LLM":
       VDM_folder = "1.LLM-HSM-MODEL"
    elif VDM == "FATES":
       VDM_folder = "1.FATES-MODEL" 
    os.makedirs('../'+VDM_folder+'/FM2VDM', exist_ok=True)
    file_names = ['PercentFuelChange.txt', 
                  'TreeTracker.txt', 
                  'treelist_VDM.dat',
                  'VDM_litter_WG.dat', 
                  'VDM_litter_trees.dat', 
                  '../'+VDM_folder+'/FM2VDM/AfterFireTrees.txt', 
                  '../'+VDM_folder+'/FM2VDM/AfterFireWG.txt',
                  '../'+VDM_folder+'/FM2VDM/AfterFireLitter.txt']

    for i in range(len(file_names)-3):
        # check if all input files exist
        llmft.check_file_exists(file_names[i])

    LiveDead = llmft.Treeoflife(file_names)
    # saving output files with loop index
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireTrees.txt','../'+VDM_folder+'/FM2VDM/AfterFireTrees.'+ str(ii) +'.txt')
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireWG.txt','../'+VDM_folder+'/FM2VDM/AfterFireWG.'+ str(ii) +'.txt')
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireLitter.txt','../'+VDM_folder+'/FM2VDM/AfterFireLitter.'+ str(ii) +'.txt')
 
    return LiveDead
    
def runLLMcyclical(p,nyears): #CSXM: LLM only
    os.chdir("../1.LLM-HSM-MODEL")
    flitter='FM2VDM/AfterFireLitter.txt'
    fwg='FM2VDM/AfterFireWG.txt'
    ftlist='FM2VDM/AfterFireTrees.txt'
    p=llmft.read_FT_2_LLM(flitter,fwg,ftlist,p)

    #run the LLM-HSI for nyears years
    p.fire_prob=0
    start_time = time.time()
    p.run(nyears)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    return p

def updateTreelist(p,ii): #CSXM: LLM only
    ftlist='FM2VDM/AfterFireTrees.txt'
    [lp_list,hw_list]=llmft.update_tree_info_per_location(p,ftlist,0)

    df_hw = pd.DataFrame(hw_list)
    df = pd.DataFrame(lp_list)
    df=df.append(df_hw)
    df.plot(subplots=True, layout=(4,2),figsize=(12, 10));
    df.to_csv('treelist_VDM.dat', sep=' ',header=False,index=False)
    file_in='treelist_VDM.dat'
    file_out='VDM2FM/treelist_VDM.dat'
    llmft.save_FT_treelist(file_in,file_out,0)
    
    df = pd.read_csv('VDM2FM/treelist_VDM.dat',sep=' ',
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

VDM = "FATES" # Vegetation Demography Model: "LLM" or "FATES"

# Adam suggested values
# CSXM: nyears and ncycyear may need to be the same because the workflow looks for restarting files (double check)
# nyears   = 100 # number of years for spinup (turn on SPITFIRE)                   Rutuja: 10 Adam: 100
# ncycyear = 5   # number of years for VDM to run in each loop (turn off SPITFIRE) Rutuja: 5  Adam: 5
# ncycle   = 10  # number of loop                                                  Rutuja: 20 Adam: 10 (Adam decided to go from 10 to 20 on Apr. 11 2023)
with open('config.yaml', 'r') as file:
    config_dict = yaml.safe_load(file)
nmonths = config_dict['NMONTHS']
ncycyear = config_dict['NCYCYEAR']
ncycle = config_dict['NCYCLE']

#Build Trees
os.chdir('5.TREES-QUICFIRE')
ierr = call('make', shell=True)

# SPINUP
if VDM == "LLM":
    os.chdir('../1.LLM-HSM-MODEL')
    LLMspinup(nyears)          # temporary llm class
    llm=LLMtransient(nyears)   # permanent llm class
    llm=dbh_cr(llm)            # calculates dbh and crown radius 
    savelittersLLMQF(llm,0)
    llmft.create_treelist(llm,'VDM2FM/treelist_VDM.dat')
elif VDM == "FATES":
    RESTART="FALSE"
    os.chdir('../1.FATES-MODEL')
    with open('../config.yaml', 'r') as file:
        y = yaml.safe_load(file)
        y['STOP_N'] = nmonths
        y['FINAL_TAG_YEAR'] = int(y['DATM_CLMNCEP_YR_START']) + nmonths//12 - 1
        #y['FINAL_TAG_YEAR'] = confi_dict['FINAL_TAG_YEAR'] #+ nyears - 1 
        # y['REST_N'] = nyears # commented out by SXM
        #y['REST_N'] = 1 #ASXM: hardwired as 1 as otherwise do not generate the restarting file at nyears+ncycyear for the first ncycyear (i.e., 160+5 years) (to fix later by working on nyears, ncycyear and ncycle)
        #y['FINAL_TAG_YEAR'] = confi_dict['FINAL_TAG_YEAR'] #+ nyears - 1 
        y['CYCLE_INDEX'] = 0
    with open('../config.yaml', 'w') as file:
        yaml.dump(y, file, default_flow_style=False, sort_keys=False)
    dir = '../1.FATES-MODEL/VDM2FM'
    shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir)
    # last visit reading: stopped here by SXM
    #subprocess.call(['sh', './src/prep_elm_parallel.sh'])
    subprocess.call(['sh', './src/run_elm_parallel.sh', RESTART])

    #### MAKE INTO FUNCTION
df = pd.read_csv('VDM2FM/treelist_VDM.dat',sep=' ',
                          names=["Tree id","x coord [m]","y coord [m]","Ht [m]",
                              "htlc [m]","CRDiameter [m]","hmaxcr [m]",
                              "canopydensity  [kg/m3]", "CR fuel moist [frac]",
                              "CR fuel size scale [m]", "treeid" ])
df.drop("treeid", inplace=True, axis = 1)
df.plot(subplots=True, layout=(5,2),figsize=(12, 10));
plt.tight_layout()
os.makedirs('figures', exist_ok=True)
plt.savefig('figures/TreeInfo.png')

plt.title('Tree distribution in the FT domain');
print("Total number of trees: ",df["x coord [m]"].size )
plt.savefig('figures/TreePlot.0.png')

if VDM == "LLM":
    hsi_plt.plot_species_scores(llm)
    plt.savefig('figures/HVI.0.png')
#### MAKE ABOVE INTO FUNTION

## Change Coordinates for QUICFIRE HERE ###

#buff.add_surf_buff()

LiveDead=[]
for i in range(ncycle):
    ii = i + 1
    runTreeQF()                    # runs the tree program to create QF inputs
    runQF(i)                       # runs QuickFire
    L=np.array(runCrownScorch(ii)) # runs the tree program to create LLM inputs
    L=np.insert(L,0,ii)
    LiveDead.append(L)    
    ## Change Coordinates Back to Eco system model HERE ###
    #buff.remove_tree_buff()
    #buff.remove_surf_buff()
    print('Loop Number: ',ii)
    if VDM == "LLM":
        llm=runLLMcyclical(llm,ncycyear)  # runs LLM-HSM with no fire 
        hsi_plt.plot_species_scores(llm)  # Plotting HVI
        plt.savefig('figures/HVI.png')
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
        print ('SQ', llm.sq_sc)
        print ('GT', llm.gt_sc)
        print ('RCW',sc_rcw)
        np.savetxt('HVI-score.txt', np.c_[sc_rcw, llm.sq_sc, llm.gt_sc], fmt='%1.4e')
    elif VDM == "FATES":
        os.chdir('../1.FATES-MODEL')
        subprocess.call(['sh', './src/update.restart.treelist.sh']) 
        RESTART="TRUE"
        with open('../config.yaml', 'r') as file:
            y = yaml.safe_load(file)
            y['STOP_N'] = ncycyear #ncycyear*(1 + (ncycle - 1))
            # y['REST_N'] = ncycyear # commented out by SXM
            #y['REST_N'] = 1  #ASXM: hardwired as 1, not nec. though (to fix later by working on nyears, ncycyear and ncycle)
            y['FINAL_TAG_YEAR'] = y['FINAL_TAG_YEAR'] + y['NCYCYEAR']
            y['CYCLE_INDEX'] = ii
        with open('../config.yaml', 'w') as file:
            yaml.dump(y, file, default_flow_style=False, sort_keys=False)
        subprocess.call(['sh', './src/run_elm_parallel.sh', RESTART])

LiveDead=np.array(LiveDead)
os.makedirs('output', exist_ok=True)
np.savetxt('LiveDead.txt', LiveDead, fmt='%i', header='Fire LLP(L/D) Turk(L/D)')
