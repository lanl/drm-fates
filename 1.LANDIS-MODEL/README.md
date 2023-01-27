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
and each cycle will be 5 years, then both NECN and Biomass Community Output should be set to print output files every
year. However, if for example the spinup is 3 years and each cycle is 6 years, then outputs could be printed at 
3-year intervals.

+ 