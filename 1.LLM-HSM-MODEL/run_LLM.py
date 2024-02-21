"""
This script houses the LLM functions previously in DRM_framework_coupling,
streamlining the imports and simplifying organization. No functions have 
been altered.
"""

# Core imports
import os
from time import time

# External imports
import numpy as np

# Internal imports
import LLM_model_class as llm
import LLM_FT_utils as llmft
import hsiscore_class as HSI
import hsi_plot_utils as hsi_plt
import LLM_display


def LLMspinup(nyears):
    # --spinup run ---
    p = llm.LLM()  # assign p to the llm class
    p.dim = 80
    p.instantiate(0)  # 1: reads input data from file, 0: generate inputs internally
    p.readfireprobfromfile = 0
    p.readmastprobfromfile = 0
    p.verbose = 0

    start_time = time()
    p.run(nyears)  # here 200 is a number of years
    print("--- %s seconds ---" % (time() - start_time))
    p.save_pickle()  # saves the results

    return


def LLMtransient(nyears):
    # --transient run ---
    # del p
    p = llm.LLM()
    p.dim = 80
    p.randfromfile = 0
    p.instantiate(1)  # 1: reads input data from file, 0: generate inputs internally
    p.verbose = 0  # 0: do not print out scores, 1: print scores on the screen
    p.tree_mature_age = 10
    p.readfireprobfromfile = 0
    p.readmastprobfromfile = 0

    start_time = time()
    p.run(nyears)
    print("--- %s seconds ---" % (time() - start_time))
    if np.sum(p.litter) == 0:
        print("no litter, make an extra run...")
        p.fire_prob = 0
        p.run(1)
    # p.save_randfromfile()
    return p


def dbh_cr(p):
    lp_count = p.old_LPcount.copy()
    lp_count[p.old_ht < 1.37] = 0

    lp_height = p.old_ht.copy()
    lp_height[lp_height < 1.37] = 1.37
    p.lp_dbh = llmft.dbh1_model(lp_height)
    # print dbh
    # axx=llmft.plot_area_matrix(p.lp_dbh,'LLP tree DBH [cm]','yes')

    lp_CA = llmft.dbh2_model(lp_height, p.lp_dbh)
    p.lp_CR = np.sqrt(lp_CA / np.pi)  # LLP Crown Area
    all_NaNs = np.isnan(p.lp_CR)
    p.lp_CR[all_NaNs] = 0
    # axx=llmft.plot_area_matrix(p.lp_CR,'LLP tree CR [m]','yes')

    hw_height = p.old_htHW.copy()
    hw_height[hw_height < 1.37] = 1.37
    p.hw_dbh = llmft.dbh1_model(hw_height)
    # print dbh
    # axx=llmft.plot_area_matrix(p.hw_dbh,'HW tree DBH [cm]','yes')

    hw_CR = llmft.dbh2cr_hw(p.hw_dbh / 2.54)  # note dbh is in inch
    p.hw_CR = hw_CR / 3.281  # CR is in feet convert to meters
    all_NaNs = np.isnan(p.hw_CR)
    p.hw_CR[all_NaNs] = 0
    # axx=llmft.plot_area_matrix(p.hw_CR,'HW tree CR [m]','yes')

    return p


def savelittersLLMQF(p, i):
    filename = "VDM2FM/VDM_litter_WG.dat"
    ftitle = "WG litter [kg/4m2]"
    llmft.save_litter_LLM_FT(filename, ftitle, p.litterWG, "plot", "grass")
    newname = "litter_WG." + str(i) + ".png"
    os.rename("litter.png", newname)

    filename = "VDM2FM/VDM_litter_trees.dat"
    ftitle = "LLP + HW litter [kg/4m2]"
    tree_litter = p.litterHW + p.litter
    llmft.save_litter_LLM_FT(filename, ftitle, tree_litter, "plot", "litter")
    newname = "litter_Tree." + str(i) + ".png"
    os.rename("litter.png", newname)

    percent_LP_litter = np.sum(p.litter) / np.sum(p.litterHW + p.litter)
    percent_HW_litter = np.sum(p.litterHW) / np.sum(p.litterHW + p.litter)
    print("lit_LLP%, lit_HW%:", percent_LP_litter, percent_HW_litter)

    return
