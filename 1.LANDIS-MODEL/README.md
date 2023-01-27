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

