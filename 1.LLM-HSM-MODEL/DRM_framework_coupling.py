# This a master coupling script, that should streamline all the components of the DRM framework
#
#(c) Elchin Jafarov 03/30/2021

import numpy as np
import hsiscore_class as HSI
import LLM_model_class as llm
import hsi_plot_utils as hsi_plt
import time
import pandas as pd
import random
import LLM_FT_utils as llmft
import sys

def first_LLM_run(nyears,ffolder):
    # nyears: number of years for the spinup and transient runs
    # ffolder: the location of the folder where the files will be saved
    # Runs the spinup and transient LLM runs
    # Produces 3 files: litter_WG.dat, litter_tree.dat, treelist_LLM.dat
    # in the ffolder
    hsi=HSI.hsi_score()
    # --spinup run ---
    p = llm.LLM()     # assign p to the llm class
    #NOTE: p.instantiate(1) assumes that runs correspoding input files do exist if not then run
    p.dim = 80
    p.instantiate(0)  # 1: reads input data from file, 0: generate inputs internally
    p.readfireprobfromfile=0
    p.readmastprobfromfile=0
    p.verbose=0

    start_time = time.time()
    print('Running LLM spinup for',nyears,' years....')
    p.run(nyears) # here 200 is a number of years
    print("--- spinup time %s seconds ---" % (time.time() - start_time))
    p.save_pickle() #saves the results

    # --transient run ---
    del p
    p = llm.LLM() 
    p.dim = 80
    p.randfromfile=0
    p.instantiate(1)  # 1: reads input data from file, 0: generate inputs internally
    p.verbose=0       # 0: do not print out scores, 1: print scores on the screen
    p.tree_mature_age=10
    p.readfireprobfromfile=0
    p.readmastprobfromfile=0

    start_time = time.time()
    print('Running LLM transient for',nyears,' years....')
    p.run(nyears)
    print("--- transient time %s seconds ---" % (time.time() - start_time))
    if np.sum(p.litter)==0:
        print('no litter, making an extra run...')
        p.fire_prob=0
        p.run(1)
    # --end of transient run ---
    print ('End of spinup and transient runs')
    print ('--------------------------')
    print ('')

    # --calculate dbh and CR for LLPs and HWs ---
    print ('Calculate dbh and CR for LLPs and HWs')
    lp_count=p.old_LPcount.copy()
    lp_count[p.old_ht<1.37]=0

    lp_height=p.old_ht.copy()
    lp_height[lp_height<1.37]=1.37
    p.lp_dbh=llmft.dbh1_model(lp_height)

    lp_CA=llmft.dbh2_model(lp_height,p.lp_dbh)
    p.lp_CR=np.sqrt(lp_CA/np.pi) #LLP Crown Area
    all_NaNs = np.isnan(p.lp_CR)
    p.lp_CR[all_NaNs] = 0

    hw_height=p.old_htHW.copy()
    hw_height[hw_height<1.37]=1.37
    p.hw_dbh=llmft.dbh1_model(hw_height)

    hw_CR=llmft.dbh2cr_hw(p.hw_dbh/2.54) # note dbh is in inch
    p.hw_CR=hw_CR/3.281            # CR is in feet convert to meters
    all_NaNs = np.isnan(p.hw_CR)
    p.hw_CR[all_NaNs] = 0

    # --Regrid LLM litters (80x80), each cell area 25m2 and save them in FT format (200x200), 4m2 cell area --
    filename=ffolder+'/litter_WG.dat'
    ftitle='WG litter [kg/4m2]'
    llmft.save_litter_LLM_FT(filename,ftitle,p.litterWG,'noplot')

    filename=ffolder+'/litter_tree.dat'
    ftitle='LLP + HW litter [kg/4m2]'
    tree_litter=p.litterHW+p.litter
    llmft.save_litter_LLM_FT(filename,ftitle,tree_litter,'noplot')

    percent_LP_litter=np.sum(p.litter)/np.sum(p.litterHW+p.litter)
    percent_HW_litter=np.sum(p.litterHW)/np.sum(p.litterHW+p.litter)
    print ('lit_LLP%, lit_HW%:',percent_LP_litter,percent_HW_litter)

    # Make a treelist file for the first run only
    llmft.create_treelist(p,ffolder+'/treelist_LLM.dat')
    
    return

if len(sys.argv) == 1:
    n_years = 1
    folder_name = 'LLM2FT'
else:
    n_years = sys.argv[1]
    folder_name = sys.argv[2]
    
print ('Run time [years]:',n_years)
print ('Folder name:',folder_name)
print ('--------------------------')

first_LLM_run(int(n_years),folder_name)
