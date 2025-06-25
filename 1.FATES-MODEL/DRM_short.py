
import os
import numpy as np
import pandas
import sys
import os.path
import subprocess
from shutil import copyfile
from subprocess import call
import shutil
#import yaml
sys.path.insert(0, '1.LLM-HSM-MODEL/')
#import LLM_model_class as llm
import LLM_FT_utils as llmft
#import hsiscore_class as HSI
#import hsi_plot_utils as hsi_plt
#import LLM_display
import matplotlib
#sys.path.insert(0, '8.CROWN-SCORCH/')
#import Treeoflife as tol 


sys.path.insert(0, '7.QUICFIRE-MODEL/projects/Tester')



def runQF(i):
    # copy produced by Tree program files to the QF folder
    #os.chdir("/Users/elchin/Documents/Adams_project/llm-hsm-ft/")
    os.chdir(os.path.dirname(__file__))
    src='5.TREES-QUICFIRE/'
    dst='7.QUICFIRE-MODEL/projects/Tester/'
    copyfile(src+'treesfueldepth.dat',dst+'treesfueldepth.dat')
    copyfile(src+'treesmoist.dat',dst+'treesmoist.dat')
    copyfile(src+'treesrhof.dat',dst+'treesrhof.dat')
    copyfile(src+'treesss.dat',dst+'treesss.dat')

    # Before running the shell script make sure that
    # file path '/projects/Tester/QUIC_fire.inp'
    # is pointing to the gridlist file in the projects/ftFiles directory.
    os.chdir("7.QUICFIRE-MODEL/mac_compile/")
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
    #python drawfire.py $RUN5/drm-fates/7.QUICFIRE-MODEL/projects/Tester/
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
	
	
	
	
def runTreeQF():
	# Note: Adam has a QF Tree code in '5.TREES-QUICFIRE'
    os.chdir(os.path.dirname(__file__))
    VDM_folder = "1.FATES-MODEL"
    src=VDM_folder+'/VDM2FM/'
    print(src)
    dst='5.TREES-QUICFIRE/'
    print(os.getcwd())
    copyfile(src+'VDM_litter_WG.dat',dst+'LLM_litter_WG.dat')
    copyfile(src+'treelist_VDM.dat',dst+'treelist_LLM.dat')
    copyfile(src+'VDM_litter_trees.dat',dst+'LLM_litter_trees.dat')

    os.chdir(dst)
    status=subprocess.call(["./trees"])
    while status != 0:
       print('Tree program failed to execute...')
       status=subprocess.call(["./trees"])

    if status==0:
        print('Tree program run successfully!')
        copyfile('TreeTracker.txt','../8.CROWN-SCORCH/TreeTracker.txt')
        copyfile('treelist_LLM.dat','../8.CROWN-SCORCH/treelist_VDM.dat')
        copyfile('LLM_litter_WG.dat','../8.CROWN-SCORCH/VDM_litter_WG.dat')
        copyfile('LLM_litter_trees.dat','../8.CROWN-SCORCH/VDM_litter_trees.dat')
    return

def runCrownScorch(ii):
    os.chdir(os.path.dirname(__file__))
    os.chdir("8.CROWN-SCORCH")
    copyfile('../7.QUICFIRE-MODEL/projects/Tester/PercentFuelChange.txt','../8.CROWN-SCORCH/PercentFuelChange.txt')
    LiveDead = []
    #if VDM == "LLM":
    #   VDM_folder = "1.LLM-HSM-MODEL"
    #elif VDM == "FATES":
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

    llmft.Treeoflife(file_names)
    # saving output files with loop index
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireTrees.txt','../'+VDM_folder+'/FM2VDM/AfterFireTrees.'+ str(ii) +'.txt')
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireWG.txt','../'+VDM_folder+'/FM2VDM/AfterFireWG.'+ str(ii) +'.txt')
    copyfile('../'+VDM_folder+'/FM2VDM/AfterFireLitter.txt','../'+VDM_folder+'/FM2VDM/AfterFireLitter.'+ str(ii) +'.txt')

    #return LiveDead








i=3
ii = i + 1
print("runTreeOF")
runTreeQF()                    # runs the tree program to create QF inputs
print("RunQF")
runQF(i)                       # runs QuickFire
L=np.array(runCrownScorch(ii)) # runs the tree program to create LLM inputs
print(L)
