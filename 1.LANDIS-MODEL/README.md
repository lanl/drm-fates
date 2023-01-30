## Setting up the DRM with LANDIS and QUIC-Fire

LANDIS-II can be used in the Disturbance Response Model framework to incorporate the impacts of fine-scale fire
behavior into a broad-scale model of landscape change. The framework assumes that the LANDIS model has already
been built and parameterized by the user. There is flexibility with what the LANDIS model can include, but there
are two extension requirements:

- Net Ecosystem Carbon & Nitorgen (NECN) must be used as the succession extension. 
- Biomass Community Output extension must be used for summarizing cohort biomass.

Other extensions that simulate disturbance may be used (e.g. harvest or wind disturbance extensions), but extensions
that incorporate fire disturbance should most likely be avoided since they are redundant with the goals of the fire
simulation models in the DRM.

#### Tips for setting up the LANDIS run:

+ The DRM will use outputs from NECN and Biomass Community Output, so make sure those outputs are printed at a
frequency that aligns with the lengths of the spinup and fire cycles. For example, if the spinup will last 4 years
and each cycle will be 5 years, then both NECN and Biomass Community Output Timestep should be set to 1.
However, if for example the spinup is 3 years and each cycle is 6 years, then outputs could be printed at 
3-year intervals (Timestep = 3). For best results:

	+ Set Timestep in NECN input file to 1
	+ Set Timestep in Community Biomass Output input file to 1

+ For a few key parameters in the NECN input file, make sure there are no comments (or even spaces!) after the 
input value. The framework references these values as the last item in their line. Relevant lines are:

	+ InitialCommunities
	+ InitialCommunitiesMap
	+ InitialDeadWoodSurfaceMapName
	+ InitialDeadCoarseRootsMapName
	+ InitialFineFuels

### FIA Data

The fire simululation models require a 3D fuels array that is created from a list of individual trees with unique
attributes. Because LANDIS does not track individual trees, we have developed a method for assigning these attributes
by matching LANDIS cohorts to trees from USFS Forest Inventory and Analysis (FIA) plots. Thus, before running the DRM
with LANDIS, the user must download the necessary FIA data.

+ FIA data can be accessed at the 
[FIA DataMart](https://experience.arcgis.com/experience/3641cea45d614ab88791aef54f3a1849/)

+ Data should be downloaded in CSV format and placed in 9.FIA/FIA_raw.

+ Download data from the state the AOI is in, as well as all (or some) surrounding states. Data from FIA plots in these
states will be used to approximate the LANDIS cohorts, so download any states with forests that have tree species present
in the AOI. Likewise, data from states that do not have those species (even if they border the AOI state) do not need to
be downloaded.

+ For each state XX, download **XX_TREE.csv** and **XX_COND.csv**. No other datasets are necessary.

### Cropping the LANDIS domain.

The user might wish to run each iteration of the LANDIS model at a landscape extent, but conduct fire simulations on a
smaller section of that domain. Included in the DRM is the option to crop the LANDIS run before each fire simulation to
save time on the more computationally-intensive fire simulation. If cropping is desired, the user should do the following:

1. Create a shapefile of the desired burn plot, called **burn_plot.shp**
2. Put the shapefile in 1.LANDIS-MODEL/Shapefiles
3. Change the "crop_domain" input in LANDIS_options.py to True *(see below)*

### LANDIS_options.py

Once the LANDIS run is set up, parameters must be input into 1.LANDIS-MODEL/LANDIS_options.py so that a treelist 
can be generated. All parameters are input as values in a python dictionary, so lists must be enclosed in brackets 

+ Run Info:

	+ **states**: list of states including and surrounding the area of interest (AOI). FIA data from these states
will be compiled for matching to LANDIS cohorts when building a treelist, so states that do not include trees species
in the LANDIS run do not need to be included (e.g. a simulation of Colorado subalpine forests would not need to 
include Kansas in the list of states, but should include Wyoming)

	+ **fia_spec**: list of the tree species in the LANDIS run as they appear in the SPECIES_SYMBOL field of FIA data.
These are the codes used by the USDA PLANTS (Plant List of Attributes, Names, Taxonomy, and Symbols) database.

	+ **landis_spec**: list of the corresponding codes used in the LANDIS run, whether or not they are different from
the PLANTS codes. Species must be listed in the same order as fia_spec.

	+ **region_flag**: flag indicating where the AOI is located. 

		+ 1 = California
		+ 2 = Other western state (WA, OR, ID, MT, WY, NV, UT, AZ, NM, CO)
		+ 3 = Midwest, eastern, or southern state (any state not listed above)

	+ **age_bin**: integer indicating how age cohorts should be grouped. Default is 10. Tree ages are ceiling rounded
to the nearest age bin to assign their cohorts (e.g. a 9 year-old tree is in the 0-10 year cohort, and an 11 year-old tree
is in the 11-20 year cohort).

	+ **aoi_ele**v: the approximate average elevation of the AOI. We are working on a method to calculate this in the
code.

+ Fuels:

	+ **bulk_density**: value of canopy bulk density for all trees (kg/m<sup>3</sup>). Default is 0.7. Future versions 
of the framework will include the option to set different bulk densities for each tree species.

	+ **cl_factor**: value indicating the fraction of the crown above the maximum crown diameter for all trees. Default 
is 0.8. Future versions of the framework will include the option to set different CL factors for each tree species.

	+ **moisture**: value of canopy moisture content (proportion) for all trees. Default is 1. Future version of the 
framework will include the option to set different moisture values for each tree species.

	+ **sizescale**: value for the size scale of canopy fuels. Default is 0.0005.

+ Spatial Info:

	+ **crop_domain**: boolean indicating whether to crop the LANDIS domain to a smaller burn domain before each fire
simulation.

+ Fire Effects:

	+ **mortality_thresholds**: threshold of percent of the tree canopy remaining after fire, under which a tree will
not survive. Enter one value for each tree species, in the order corresponding to the species in the fia_spec and
landis_spec lists. 

### QUICFire_options.py

The DRM has built-in functions to write QUIC-Fire input files based on user inputs in QUICFire_options.py. Like in
LANDIS_options.py, all values are input to a python dictionary.

+ Simulation Parameters:

	+ **QFVD**: version of QUICFire to use. Currently versions 4 and 5 are supported. 
		+ *NOTE: Differences may exist between versions within a major release. One in particular might be
the omegarelax input in the QU_simparams.inp file. If this is causing QUIC-Fire to crash, comment or uncomment line 320 
in QFVD4/print_outputs.py*

	+ **PROJ_FOLDER**: name of the folder containing the QUIC-Fire input files. 
		+ *NOTE: this **must** match the "testcase" argument in 
7.QUICFIRE-MODEL/mac_compile/adv_comiple_and_run.sh*

	+ **SimTime**: fire simulation time (s).

	+ **print_times**: time interval (s) at which to print output arrays of fuel density, etc. 

+ Domain Settings:

	+ **nx**: horizontal fire grid x size (number of 2m cells)
	+ **ny**: horizontal fire grid y size (number of 2m cells)
	+ **nz**: vertical grid size (number of 1m cells)

		+ If running the DRM with LANDIS, the domain size will be altered in the code based on the size of either
the LANDIS domain or the user-provided burn domain, so these are placeholder values. If using LLM or FATES, domain size 
should be set here.

+ Topo Settings:

	+ **topo_custom**: boolean indicating whether a custom topo.inp file will be used. Currently the DRM only supports
flat topography, but future versions will include the option to use custom topography.

	+ **max_topo**: DO NOT CHANGE. Set to zero for flat topography. Future versions of the DRAM will alter this value 
based on the difference between the highest and lowest elevations in the AOI.

+ Wind Settings:

	+ **windspeed**: windspeed in m/s.
	+ **winddir**: wind direction in degrees.

+ Ignition Settings:

	+ **custom_ig**: boolean indicating whether or not custom ignitions should be used. Currently the
DRM does not produce custom ignition patterns internally, but they may be provided by the user through an ignite.dat file.
If set to true, next four inputs will not be used.

	+ **ig_xmin**: x coordinate (m) of lower-right corner of rectangular ignition
	+ **ig_ymin**: y coordinate (m) of lower-right corner of rectangular ignition
	+ **ig_xlen**: length (m) of rectangular ignition, in positive x direction 
	+ **ig_ylen**: length (m) of rectangular ignition, in positive y direction 

		+ If using LANDIS, the above values will be calculated based on the size of the burn domain. Ignitions
will be 100m from the west edge of the domain, will be 1/2 the length of the y dimension of the domain, and start 1/4 the
y length of the domain from the origin. If using LLM or FATES, ignition will be determined from the above values.

	+ **ig_method**, **ig_pattern**, **ig_spacing**: inputs for potential future implementation of driptorch.


