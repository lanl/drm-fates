# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin
package require parflow
namespace import Parflow::*


#set tcl_precision 16

#-----------------------------------------------------------------------------
# File input version number
#-----------------------------------------------------------------------------
pfset FileVersion 4

set runname Spin
# set number of time steps
set tt 30

set porosity         [pfload $runname.out.porosity.silo]
set mask             [pfload $runname.out.mask.silo]
set top              [pfcomputetop $mask
set top              [pfcomputetop $mask]]
set specific_storage [pfload $runname.out.specific_storage.silo]

for {set i 0} {$i <= $tt} {incr i} {
  set total_water_in_domain 0.0
  set filename [format "%s.out.press.%05d.pfb" $runname $i]
  set pressure [pfload $filename]
  pfsave $pressure -sa "Pressure.$i.txt"
  set surface_storage [pfsurfacestorage $top $pressure]
  pfsave $surface_storage -silo "surface_storage.$i.silo"
  set total_surface_storage [pfsum $surface_storage]
  set total_water_in_domain [expr $total_water_in_domain + $total_surface_storage]
  pfsave $surface_storage -sa "Surface_Storage.$i.txt"
  set filename [format "%s.out.satur.%05d.pfb" $runname $i]
  set saturation [pfload $filename]
  pfsave $saturation SaturationNEW.$i.silo
  pfsave $saturation -sa "Saturation.$i.txt"
}





