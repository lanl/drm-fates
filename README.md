WORKFLOW INSTRUCTIONS:
---------------------

 * Run LLM
 * Postprocess LLM resutls
 * Run VegMap-ParFlow (optional)
 * Run Parflow (optional)
 * Run Tree code
 * Run CanopyEnergyBalance code 
 * Run FIRETEC or QUICFIRE
 * Run Fuel Density code


Run LLM
------------

Go to LLM-HSM-MODEL and follow the instructions in the "Coupling_script_LLM-HSM.ipynb" script. 

Postprocess LLM resutls
------------

Run Elchins' LLM postprocessing to get trees code, and hydrology input.
Produces 3 files:
1) 'treelist.txt' Containg  tree species and location<br/>
   Format in colums:<br/>
   <em> spp	x	y	Ht	HTLC	CanDiameter	HttoMaxCrDi	CanDensity	Canmoisture	sizescale </em>
2) 'LLM_litter_trees.txt'  Array in the domeintions of the domain containg litter from trees density
3) 'LLM_litter_WG.txt' Array in the domeintions of the domain containg grass density  

Run VegMap-ParFlow (optional)
------------

Running this code will map trees for Parflow. Requires running 'MakeTreeMap_LLM.cpp' first to make 'CellTreeMap.txt', a long column of species number <br/>
   <em>1 = LLP, 2 = Turkey Oak, 0 = grass</em><br/>
   \*NOTE: that number FIRETEC and PARFLOW domain cells are hard coded in 'MakeTreeMap_LLM.cpp'\*<br/>
   <em>  L64         int nx = 200, ny = 200;</em><br/>
   Then run 'makeslope.cpp' to produce the 'drv_vegm.dat' file need by ParFlow-CLM<br/> 
   This Also makes slope files for ParFlow in topography as needed<br/>
   \*NOTE that number FIRETEC and PARFLOW domain cells are hard coded 'makeslope.cpp'\*<br/>
   <em>  L62         int nx=200, ny=200,....</em>

Run Parflow (optional)
------------

This run will simulate for subsurface moisture conditons. Produced 'Saturation.txt' that contains all the subsurface moistures. Alex's runs Spin-up runs are in /scratch/wildfire1/ajonko/Eglin-Parflow-OneSoil. Adam's runs are in es38:/lclscratch/aatchley/Eglin-Parflow/CLIMATE-RUNS

Run Tree code
------------

Run Trees_EJ code aka 'trees' using the files 'LLM_litter_trees.txt', 'LLM_litter_WG.txt' and the 'treelist.txt'  and ''Saturation.txt'' as input. Produces 'treesmoist.dat'  'treesss.dat' 'treesfueldepth.dat'  'treesrhof.dat'. \*Note that input file are in parent directory while, fuellest (genral input) is in subdirectory 'Input' -- Weird\*. \* 'Trees_EJ_working-9-14-WORKING' for FIRETECH \* 'TREES-QUICFIRE' For QUICFire

Run CanopyEnergyBalance code 
------------

Use 'treesrhof.dat' to run CanopyEnergyBalance code to calculate canopy & Dead Fuel moisture. Also uses 'treesmoist.dat' and rewrites 'treesmoist.dat' with added dead fuel moisture.

Run FIRETEC or QUICFIRE
------------

A) \* QUICFIRE run notes \* <br/>
QUICFIRE Diretory Tree <br/>
mac_compile  projects  source_code<br/>
                        ``/`` ``\``<br/>
                      ``/``     `` \``<br/>
             ftFiles          Tester<br/>

Put fuel data (all the trees*.dat) to run QUICFIRE in the ftFiles directory along with the ignite.dat 
The input decks are in the Tester directory, the main input deck is 'QUIC_fire.inp'. To run the model exicute './compile_and_run.sh 0' in the mac_compile directory. That .sh file has flag for whether or not you want to recompile source and also lets you point the executable to a project folder.That's the testcase variable at the top, which points quicfire to your .inp files in the project/Tester folder defined by testcase

B) /* FIRETEC NOTES /*
FIRETEC postprocess file is 'VTK.NEW.MATK.firetec.py'.  In addition to the 'PercentFuelChange.txt' file it produces all kinds of other goodies.

 Run Fuel Density code
------------

This code kill trees and convert inout to LLM. See directory CROWN_SCORCH for description of Treeoflife.py. Once complete go back to step 1 (Run LLM). 
