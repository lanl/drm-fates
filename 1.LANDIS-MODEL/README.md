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

### LANDIS_options.py

Once the LANDIS run is set up, parameters must be input into LANDIS_options.py so that a treelist can be generated.
All parameters are input as values in a python dictionary, so lists must be enclosed in brackets 

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

	+ **moisture**: value of canopy moisture content (%/100) for all trees. Default is 1. Future version of the 
framework will include the option to set different moisture values for each tree species.

	+ **sizescale**: value for the size scale of canopy fuels. Default is 0.0005.

+ Spatial Info:

	+ **QF_res**: DO NOT CHANGE. The horizontal QUIC-Fire resolution of 2m.

+ Fire Effects:

	+ **mortality_thresholds**: threshold of percent of the tree canopy remaining after fire, under which a tree will
not survive. Enter one value for each tree species, in the order corresponding to the species in the fia_spec and
landis_spec lists. 


