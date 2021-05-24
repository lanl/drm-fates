WORKFLOW FOR TRACKING FIRE AFFECTS TO TREES
USED FOR COMBUSTION AND OR CROWN SCORCH


1) Create a tree map that includes:
    a. Species Identification     [-]
    b. X and Y location in a plot [m,m]
    c. Tree height                [m]
    d. Bottom of crown height     [m]
    e. Crown diameter             [m]
    f. Height to max crown diameter [m]
    g. Fuel bulk density          [kg/m3]
    h. Fuel moisture content      [fraction]
    i. Fuel size scale            [m]


2) Use the fuel map to populate fuel in FIRETEC or QUICFIRE
   This is done by Alex Josephson’s trees code.
   In there I used the code to create a file "TreeTracker.txt"
   that contains (See TreeTracker.txt-LABELS for an example):

Tree#  #ofCellsW/Fuel  ---------------Cell_ID-------------------      ------------- %ofTree_in_each_tree---------

   Each row represents an individual tree and contains the number of cells in the domain that contains fuel,
   the cell_ID of those trees, and the % of the tree in each of those cells.  For example,
   if an individual has 4 cellsW/Fuel in a domain, there will be 4 cell ID's listed, and 4 %fuel
   value corresponding to the cell ID's in order.
   The cell_ID is created from the spatial domain loop by:
   z --> NZ
    Y --> NY
      X --> NX
        CELL_ID
      End X
    End Y
   End Z


3) Simulate the Fire using FIRETEC or QUICFIRE


4) Postprocess the binary output files for parameters on interest.
   For now, I only use change in fuel density to get at direct combustion.
   Air temperature, moisture, ect can also be used for scorch.
   This postprocessing step creates a file called 'PercentFuelChange.txt'
   and contains one column of data that is the percent change in fuel.
   The data is created from the spatial doing loop by:
   z --> NZ
    Y --> NY
      X --> NX
        DATA
      End X
    End Y
   End Z
   **** Note that now CELL_ID and DATA are in the same spatial format ****


5) Evaluate if a tree lives or dies.
   This is done by 'Treeoflife.py'
   I first read in PercentFuelChange.txt to populate a 1D array of the change in fuel.
   Then I read in the TreeTracker.txt file and loop through each tree individually.
   A subloop is then set by the #ofCellsW/Fuel (L109) that sets a cell pointer ('cellptr'[L110])
   and reads the percent fuel of tree for the specific cell ('conccell'[L110]).
   A 'newfuel' is then a direct comparison to the mapped 1D array of change in fuel.
   #### Here is where I would put some crown scorch function based off of temperature ####
   #### This function based off of a threshold needs to provide a change of fuel density ####
   If the total new fuel is higher than some threshold, then the tree survives and a new
   tree map (step 1) can be created with only the surviving trees.

   In addition, I evaluate what is the post fire surface fuel load and use that to populate LLM or FATES ecosystem model.


