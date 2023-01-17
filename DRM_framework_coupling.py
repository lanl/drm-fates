# This a master coupling script, that streamlines all the components of the DRM framework
#
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
import shutil
from subprocess import call
from time import sleep
import yaml
sys.path.insert(0, '1.LLM-HSM-MODEL/')
import LLM_model_class as llm
import LLM_FT_utils as llmft
import hsiscore_class as HSI
import hsi_plot_utils as hsi_plt
import LLM_display
sys.path.insert(0, '7.QUICFIRE-MODEL/projects/Tester')
import postfuelfire_new as pff 
import Buffer as buff


#VDM = "LLM" # Vegetation Demography Model: "LLM" or "FATES" or "LANDIS" #why is this here?

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
    elif VDM == "LANDIS":
        VDM_folder = "1.LANDIS-MODEL"
    src='../'+VDM_folder+'/VDM2FM/'
    dst='../5.TREES-QUICFIRE/'

    file_list = ["VDM_litter_WG.dat","treelist_VDM.dat","VDM_litter_trees.dat"]
    for i in file_list:
        if os.path.isfile(os.path.join(os.getcwd(),"VDM2FM",i)):
            copyfile(src+i,dst+i)
    os.chdir(dst)
    with subprocess.Popen(
        ["wsl","./trees"], stdout=subprocess.PIPE
    ) as process:

        def poll_and_read():
            print(f"{process.stdout.read1().decode('utf-8')}")
        
        while process.poll() != 0:
            poll_and_read()
            sleep(1)
        if process.poll()==0:
            print('Tree program run successfully!')
            ### Copying Tree Files to Fire Affects Assessment
            file_list = ["TreeTracker.txt","treelist_VDM.dat","VDM_litter_WG.dat","VDM_litter_trees.dat"]
            for i in file_list:
                if os.path.isfile(os.path.join(os.getcwd(),i)): 
                    copyfile(i,'../8.CROWN-SCORCH/'+i)
    
    # while status != 0:
    #    print('Tree program failed to execute...')
    #    status=subprocess.call(["wsl","./trees"])

    
    return

def runQF(i,VDM):
    #copy produced by Tree program files to the QF folder
    #os.chdir("/Users/elchin/Documents/Adams_project/llm-hsm-ft/")
    src=''
    dst='../7.QUICFIRE-MODEL/projects/LandisTester/'
    copyfile(src+'treesfueldepth.dat',dst+'treesfueldepth.dat')
    copyfile(src+'treesmoist.dat',dst+'treesmoist.dat')
    copyfile(src+'treesrhof.dat',dst+'treesrhof.dat')
    copyfile(src+'treesss.dat',dst+'treesss.dat')
    
    if VDM == "LANDIS":
        os.chdir("../1.LANDIS-MODEL")
        import TTRS_QUICFire_Support as ttrs
        import re
        os.chdir("../5.TREES-QUICFIRE")
        with open('fuellist', 'r', encoding='utf-8') as file:
            filelist = file.readlines()
        line_id = "nx="
        lines = [match for match in filelist if line_id in match]
        line = lines[0]
        cell_nums = list(map(int, re.findall(r'\d+', line)))
        rhof = ttrs.import_fortran_dat_file("treesrhof.dat", cell_nums)
        os.chdir("../1.LANDIS-MODEL/VDM2FM")
        surf = np.loadtxt("VDM_litter_trees.dat")
        rhof[0,:,:] = surf
        os.chdir("../../7.QUICFIRE-MODEL/projects/LandisTester/")
        ttrs.export_fortran_dat_file(rhof,"treesrhof.dat")
        os.chdir("../../../5.TREES-QUICFIRE")
    
    os.chdir("../7.QUICFIRE-MODEL/mac_compile/")
    # HAD TO CHANGE adv_compile_and_run.sh ARGUMENT testcase TO MATCH dst IN LINE 165
    # MUST CHANGE QF INPUTS TO MATCH DOMAIN SIZE
    status=subprocess.call(["wsl","./adv_compile_and_run.sh"])
    if status!=0: #when QF fails it still gives an exit status of 0...
        print('QF failed to execute...') 
        return
    else:
        print('QF run successfully!')
        #Successful run should produce bunch of binary files in 
        #7.QUICFIRE-MODEL/projects/Tester. Now run the postfire script 
        #that will generate PercentFuelChange.txt file required for the next step.
        os.chdir("../projects/LandisTester")
        pff.main(0)
        # MAtch this value at Line 5 of 7.QUICFIRE-MODEL/projects/Tester/QUIC_fire.inp
        direc = "Plots"
        dd = direc + str(i)
        if os.path.exists(dd):
           shutil.rmtree(dd)
        os.rename('Plots', dd)
        os.mkdir('Plots')
        dd = "fuels-dens-00000." + str(i) + ".vin"
        os.rename('fuels-dens-00000.bin', dd)
        dd = "fire_indexes." + str(i) + ".vin"
        os.rename('fire_indexes.bin', dd)    
        return

def runCrownScorch(ii):
    ii = 1
    os.chdir("../../../8.CROWN-SCORCH")
    copyfile('../7.QUICFIRE-MODEL/projects/LandisTester/PercentFuelChange.txt','../8.CROWN-SCORCH/PercentFuelChange.txt')
    LiveDead = []
    if VDM == "LLM":
        VDM_folder = "1.LLM-HSM-MODEL"
    elif VDM == "FATES":
        VDM_folder = "1.FATES-MODEL" 
    elif VDM == "LANDIS":
        VDM_folder = "1.LANDIS-MODEL"
    os.makedirs('../'+VDM_folder+'/FM2VDM', exist_ok=True)
    if VDM == "LANDIS":
        dostuff()
    else:
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
    
def runLLMcyclical(p,nyears):
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

def updateTreelist(p,ii):
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

VDM = "LANDIS" # Vegetation Demography Model: "LLM" or "FATES" or "LANDIS"

nyears=4      # number of years for spinup and transient runs
ncycyear=5    # number of cyclical year run
ncycle=4      # number of loops

#Build Trees
os.chdir('5.TREES-QUICFIRE')
ierr = call('make clean', shell=True)
ierr = call('make', shell=True)

# SPINUP
if VDM == "LLM":
    os.chdir('../1.LLM-HSM-MODEL')
    LLMspinup(nyears)          # temporary llm class
    llm=LLMtransient(nyears)   # permanent llm class
    llm=dbh_cr(llm)            # calculates dbh and crown radius 
    os.makedirs('VDM2FM', exist_ok=True)
    savelittersLLMQF(llm,0)
    llmft.create_treelist(llm,'VDM2FM/treelist_VDM.dat')
elif VDM == "FATES":
    RESTART="FALSE"
    os.chdir('../1.FATES-MODEL')
    with open('../config.yaml', 'r') as file:
        y = yaml.safe_load(file)
        y['STOP_N'] = nyears
        y['REST_N'] = nyears
        y['FINAL_TAG_YEAR'] = y['DATM_CLMNCEP_YR_START'] + nyears - 1
        y['CYCLE_INDEX'] = 0
    with open('../config.yaml', 'w') as file:
        yaml.dump(y, file, default_flow_style=False, sort_keys=False)
    dir = '../1.FATES-MODEL/VDM2FM'
    shutil.rmtree(dir, ignore_errors=True)
    os.makedirs(dir)
    subprocess.call(['sh', './src/prep_elm_parallel.sh'])
    subprocess.call(['sh', './src/run_elm_parallel.sh', RESTART])
elif VDM == "LANDIS":
    os.chdir('../1.LANDIS-MODEL')
    import LANDIS_to_Treelist as Landis
    import Run_LANDIS as Run
    import Crop_LANDIS as Crop
    os.chdir("..")
    OG_PATH = os.getcwd()
    cycle = 0      # current iteration (will be looped through range(0,ncycle))
    # Build Landis Parameters object for spinup
    L2_params = Run.LandisParams(OG_PATH, nyears, ncycyear, ncycle, cycle, spinup=True)
    # Run LANDIS
    Run.Landis(L2_params)
    os.chdir("..")
    # Crop to fire domain
    Crop.Landis(L2_params)
    # Build Treelist
    Treelist_params = Landis.toTreelist(L2_params)  
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
i = 0
for i in range(ncycle):
    ii = i + 1
    runTreeQF()                       # runs the tree program to create QF inputs
    runQF(i,VDM)                           # runs QuickFire
    L=np.array(runCrownScorch(ii))                  # runs the tree program to create LLM inputs
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
        plt.savefig('HVI.png') 
        os.rename('HVI.png', dd)
        print ('ADAM SQ', llm.sq_sc)
        print ('ADAM GT', llm.gt_sc)
        print ('ADAM RCW',sc_rcw)
        np.savetxt('HVI-score.txt', np.c_[sc_rcw, llm.sq_sc, llm.gt_sc], fmt='%1.4e')

    elif VDM == "FATES":
        os.chdir('../1.FATES-MODEL')
        subprocess.call(['sh', './src/update.restart.treelist.sh']) 
        RESTART="TRUE"
        with open('../config.yaml', 'r') as file:
            y = yaml.safe_load(file)
            y['STOP_N'] = ncycyear #ncycyear*(1 + (ncycle - 1))
            y['REST_N'] = ncycyear 
            y['FINAL_TAG_YEAR'] = y['FINAL_TAG_YEAR'] + ncycyear
            y['CYCLE_INDEX'] = ii
        with open('../config.yaml', 'w') as file:
            yaml.dump(y, file, default_flow_style=False, sort_keys=False)
        subprocess.call(['sh', './src/run_elm_parallel.sh', RESTART])
        
    elif VDM == "LANDIS":
        os.chdir("../1.LANDIS-MODEL")
        import Treelist_to_LANDIS as Treelist
        os.chdir("..")
        OG_PATH = os.getcwd()
        cycle = ii      # current iteration (will be looped through range(0,ncycle))
        # Build Landis Parameters object for cycles
        L2_params = Run.LandisParams(OG_PATH, nyears, ncycyear, ncycle, cycle, spinup=False)
        # Update Landis run with new treelist
        Treelist.toLandis(L2_params)
        # Run landis
        Run.Landis(L2_params)
        # Build another treelist
        Landis.toTreelist(L2_params)

LiveDead=np.array(LiveDead)
os.makedirs('output', exist_ok=True)
np.savetxt('LiveDead.txt',LiveDead,fmt='%i',header='Fire LLP(L/D) Turk(L/D)')
