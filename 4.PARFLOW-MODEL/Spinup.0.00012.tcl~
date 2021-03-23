
#This runs the tilted-v catchment problem
#  similar to that in Kollet and Maxwell (2006) AWR

set tcl_precision 17

set runname Spin 
            
#
# Import the ParFlow TCL package
#
lappend auto_path $env(PARFLOW_DIR)/bin 
package require parflow
namespace import Parflow::*

pfset FileVersion 4

pfset Process.Topology.P 8
pfset Process.Topology.Q 4
pfset Process.Topology.R 1

#---------------------------------------------------------
# Computational Grid
#---------------------------------------------------------
pfset ComputationalGrid.Lower.X           0.0
pfset ComputationalGrid.Lower.Y           0.0
pfset ComputationalGrid.Lower.Z           0.0

pfset ComputationalGrid.NX                200
pfset ComputationalGrid.NY                200
pfset ComputationalGrid.NZ                10

pfset ComputationalGrid.DX	          2 
pfset ComputationalGrid.DY                2 
pfset ComputationalGrid.DZ	          0.01

#---------------------------------------------------------
# The Names of the GeomInputs
#---------------------------------------------------------
pfset GeomInput.Names                  "domaininput litterinput"
pfset GeomInput.domaininput.GeomName   "domain"
pfset GeomInput.domaininput.InputType  Box

pfset GeomInput.litterinput.GeomName   "litter"
pfset GeomInput.litterinput.InputType  Box


#---------------------------------------------------------
# Domain Geometry 
#---------------------------------------------------------
pfset Geom.domain.Lower.X                        0.0
pfset Geom.domain.Lower.Y                        0.0
pfset Geom.domain.Lower.Z                        0.0
 
pfset Geom.domain.Upper.X                        400.0
pfset Geom.domain.Upper.Y                        400.0
pfset Geom.domain.Upper.Z                        0.10
pfset Geom.domain.Patches             "x-lower x-upper y-lower y-upper z-lower z-upper"


#---------------------------------------------------------
# Litter Geometry 
#---------------------------------------------------------
pfset Geom.litter.Lower.X                        0.0
pfset Geom.litter.Lower.Y                        0.0
pfset Geom.litter.Lower.Z                        0.9

pfset Geom.litter.Upper.X                        400.0
pfset Geom.litter.Upper.Y                        400.0
pfset Geom.litter.Upper.Z                        0.10

#---------------------------------------------------------
# Variable dz Assignments
#---------------------------------------------------------
pfset Solver.Nonlinear.VariableDz   True
pfset dzScale.GeomNames             domain
pfset dzScale.Type                  nzList
pfset dzScale.nzListNumber          10
## 9 layers, starts at 0 for the bottom to 8 for the top
## Soil layers (0-6) are 0.8, 0.4, 0.4, 0.2, 0.1, 0.05, 0.05 meters
## Litter/duff (7 and 8) 0.025, 0.025 meters  This will vary from run to run
pfset Cell.0.dzScale.Value          100.0
pfset Cell.1.dzScale.Value          40.0
pfset Cell.2.dzScale.Value          20.0
pfset Cell.3.dzScale.Value          20.0
pfset Cell.4.dzScale.Value          10.0
pfset Cell.5.dzScale.Value          3.0
pfset Cell.6.dzScale.Value          3.0
pfset Cell.7.dzScale.Value          2.0
pfset Cell.8.dzScale.Value          1.0
pfset Cell.9.dzScale.Value          1.0

#-----------------------------------------------------------------------------
# Perm
#-----------------------------------------------------------------------------

pfset Geom.Perm.Names                 "domain  litter"




# Values in m/hour
#
pfset Geom.domain.Perm.Type            TurnBands
pfset Geom.domain.Perm.Value           0.239
pfset Geom.domain.Perm.LambdaX            20.0
pfset Geom.domain.Perm.LambdaY            20.0
pfset Geom.domain.Perm.LambdaZ            3.0
pfset Geom.domain.Perm.GeomMean            0.91
pfset Geom.domain.Perm.Sigma               2.5
pfset Geom.domain.Perm.Seed                1
pfset Geom.domain.Perm.NumberLines         150
pfset Geom.domain.Perm.RZeta               5.0
pfset Geom.domain.Perm.KMax                100.0000001
pfset Geom.domain.Perm.DelK                0.2
pfset Geom.domain.Perm.MaxNPts             5
pfset Geom.domain.Perm.MaxCpts             200
pfset Geom.domain.Perm.LogNormal           Log
pfset Geom.domain.Perm.StratType           Bottom


pfset Perm.TensorType               TensorByGeom

pfset Geom.Perm.TensorByGeom.Names  "domain  litter"
pfset Geom.domain.Perm.TensorValX  20.0
pfset Geom.domain.Perm.TensorValY  20.0
pfset Geom.domain.Perm.TensorValZ  1.0d0

# Values in m/hour
#
pfset Geom.litter.Perm.Type            Constant
pfset Geom.litter.Perm.Value           0.325
#pfset Geom.litter.Perm.Value           0.005410

pfset Geom.litter.Perm.TensorValX      1.0
pfset Geom.litter.Perm.TensorValY      1.0
pfset Geom.litter.Perm.TensorValZ      1.0

#-----------------------------------------------------------------------------
# Specific Storage
#-----------------------------------------------------------------------------

pfset SpecificStorage.Type            Constant
pfset SpecificStorage.GeomNames       "domain  litter"
pfset Geom.domain.SpecificStorage.Value 1.0e-4
pfset Geom.litter.SpecificStorage.Value 1.0e-4

#-----------------------------------------------------------------------------
# Phases
#-----------------------------------------------------------------------------

pfset Phase.Names "water"

pfset Phase.water.Density.Type	        Constant
pfset Phase.water.Density.Value	        1.0

pfset Phase.water.Viscosity.Type	Constant
pfset Phase.water.Viscosity.Value	1.0

#-----------------------------------------------------------------------------
# Contaminants
#-----------------------------------------------------------------------------

pfset Contaminants.Names			""

#-----------------------------------------------------------------------------
# Retardation
#-----------------------------------------------------------------------------

pfset Geom.Retardation.GeomNames           ""

#-----------------------------------------------------------------------------
# Gravity
#-----------------------------------------------------------------------------

pfset Gravity				1.0

#-----------------------------------------------------------------------------
# Setup timing info
#-----------------------------------------------------------------------------
# 
# Base unit is hour.  Met input is Daily i.e 24hours/per met input = timestep.Value = 24
pfset TimingInfo.BaseUnit        1
pfset TimingInfo.StartCount      0
pfset TimingInfo.StartTime       0.0
# 8760 hours in a year.           
pfset TimingInfo.StopTime        262800 
#pfset TimingInfo.StopTime         1488
pfset TimingInfo.DumpInterval    8760
pfset TimeStep.Type              Constant
#pfset TimeStep.Value             0.1666666666667
pfset TimeStep.Value             100.0
 
#-----------------------------------------------------------------------------
# Porosity
#-----------------------------------------------------------------------------

pfset Geom.Porosity.GeomNames          "domain litter"

pfset Geom.domain.Porosity.Type          Constant
pfset Geom.domain.Porosity.Value         0.45
pfset Geom.litter.Porosity.Type          Constant
pfset Geom.litter.Porosity.Value         0.65

#-----------------------------------------------------------------------------
# Domain
#-----------------------------------------------------------------------------

pfset Domain.GeomName domain

#-----------------------------------------------------------------------------
# Relative Permeability
#-----------------------------------------------------------------------------

pfset Phase.RelPerm.Type               VanGenuchten
pfset Phase.RelPerm.GeomNames          "domain litter"

pfset Geom.domain.RelPerm.Alpha         4.0
pfset Geom.domain.RelPerm.N             2.0 
pfset Geom.litter.RelPerm.Alpha         10.0
pfset Geom.litter.RelPerm.N             1.70

#---------------------------------------------------------
# Saturation
#---------------------------------------------------------

pfset Phase.Saturation.Type              VanGenuchten
pfset Phase.Saturation.GeomNames         "domain litter"

pfset Geom.domain.Saturation.Alpha        4.0
pfset Geom.domain.Saturation.N            2.0
pfset Geom.domain.Saturation.SRes         0.05
pfset Geom.domain.Saturation.SSat         1.0
pfset Geom.litter.Saturation.Alpha        10.0
pfset Geom.litter.Saturation.N            1.70
pfset Geom.litter.Saturation.SRes         0.035
pfset Geom.litter.Saturation.SSat         1.0



#-----------------------------------------------------------------------------
# Wells
#-----------------------------------------------------------------------------
pfset Wells.Names                           ""

#-----------------------------------------------------------------------------
# Time Cycles
#-----------------------------------------------------------------------------
pfset Cycle.Names "constant rainrec"
pfset Cycle.constant.Names              "alltime"
pfset Cycle.constant.alltime.Length      1
pfset Cycle.constant.Repeat             -1

# rainfall and recession time periods are defined here
# rain for 1 hour, recession for 2 hours

pfset Cycle.rainrec.Names                 "rain rec"
pfset Cycle.rainrec.rain.Length           2
pfset Cycle.rainrec.rec.Length            3
pfset Cycle.rainrec.Repeat                -1
 
#-----------------------------------------------------------------------------
# Boundary Conditions: Pressure
#-----------------------------------------------------------------------------
pfset BCPressure.PatchNames                   [pfget Geom.domain.Patches]

pfset Patch.x-lower.BCPressure.Type		      FluxConst
pfset Patch.x-lower.BCPressure.Cycle		      "constant"
pfset Patch.x-lower.BCPressure.alltime.Value	      0.0

pfset Patch.y-lower.BCPressure.Type		      FluxConst
pfset Patch.y-lower.BCPressure.Cycle		      "constant"
pfset Patch.y-lower.BCPressure.alltime.Value	      0.0

pfset Patch.z-lower.BCPressure.Type		      FluxConst
pfset Patch.z-lower.BCPressure.Cycle		      "constant"
pfset Patch.z-lower.BCPressure.alltime.Value	      0.0

pfset Patch.x-upper.BCPressure.Type		      FluxConst
pfset Patch.x-upper.BCPressure.Cycle		      "constant"
pfset Patch.x-upper.BCPressure.alltime.Value	      0.0

pfset Patch.y-upper.BCPressure.Type		      FluxConst
pfset Patch.y-upper.BCPressure.Cycle		      "constant"
pfset Patch.y-upper.BCPressure.alltime.Value	      0.0

## overland flow boundary condition with very heavy rainfall then slight ET
pfset Patch.z-upper.BCPressure.Type		      OverlandFlow

pfset Patch.z-upper.BCPressure.Cycle		      "constant"
#pfset Patch.z-upper.BCPressure.alltime.Value	       0.0
#pfset Patch.z-upper.BCPressure.alltime.Value          -0.00005
pfset Patch.z-upper.BCPressure.alltime.Value           -0.00010

#pfset Patch.z-upper.BCPressure.Cycle		      "rainrec"
#pfset Patch.z-upper.BCPressure.rain.Value	      -0.05
#pfset Patch.z-upper.BCPressure.rec.Value	      0.000001

#---------------------------------------------------------
# Topo slopes in x-direction
#---------------------------------------------------------
## Converting heterogeneous slopes for Parflow
set xslope [pfload -sa Xgrid.sa.txt]
pfsave $xslope -pfb "Xgrid.pfb"
pfsave $xslope -silo "Xgrid.silo"
#
pfset TopoSlopesX.Type "PFBFile"
#pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.GeomNames "domain"
pfset TopoSlopesX.FileName   Xgrid.pfb

#pfset TopoSlopesX.GeomNames          "litter"
#pfset TopoSlopesX.Type               "Constant"
#pfset TopoSlopesX.Geom.litter.Value  0.000



#---------------------------------------------------------
# Topo slopes in y-direction
#---------------------------------------------------------
## Converting heterogeneous slopes for Parflow
set yslope [pfload -sa Ygrid.sa.txt]
pfsave $yslope -pfb "Ygrid.pfb"
pfsave $yslope -silo "Ygrid.silo"
#
pfset TopoSlopesY.Type "PFBFile"
#pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.GeomNames "domain"
pfset TopoSlopesY.FileName Ygrid.pfb

#pfset TopoSlopesY.GeomNames          "litter"
#pfset TopoSlopesY.Type               "Constant"
#pfset TopoSlopesY.Geom.litter.Value   0.20

#-------------------------------
# Save DEM.sa.txt as .pfb For post-process imagining
#-------------------------------
set dem [pfload -sa DEM.sa.txt]
pfsave $dem -pfb "Slope.pfb"

#---------
##  Distribute slopes
#---------

pfset ComputationalGrid.NZ                1


pfdist Xgrid.pfb
pfdist Ygrid.pfb

pfset ComputationalGrid.NZ                10


#---------------------------------------------------------
# Mannings coefficient 
#---------------------------------------------------------

pfset Mannings.Type "Constant"
pfset Mannings.GeomNames "litter"
pfset Mannings.Geom.litter.Value 5.8e-5
#pfset Mannings.Geom.litter.Value 0.018
#-----------------------------------------------------------------------------
# Phase sources:
#-----------------------------------------------------------------------------

pfset PhaseSources.water.Type                         Constant
pfset PhaseSources.water.GeomNames                    domain
pfset PhaseSources.water.Geom.domain.Value        0.0


#-----------------------------------------------------------------------------
# Exact solution specification for error calculations
#-----------------------------------------------------------------------------

pfset KnownSolution                                    NoKnownSolution


#-----------------------------------------------------------------------------
# Set solver parameters
#-----------------------------------------------------------------------------

pfset Solver                                             Richards
pfset Solver.MaxIter                                    500000
pfset Solver.TerrainFollowingGrid                        True
pfset Solver.Nonlinear.MaxIter                           130
pfset Solver.Nonlinear.ResidualTol                       1e-7
pfset Solver.Nonlinear.EtaChoice                         Walker1 
pfset Solver.Nonlinear.EtaChoice                         EtaConstant
pfset Solver.Nonlinear.EtaValue                          0.001
pfset Solver.Nonlinear.UseJacobian                       False
pfset Solver.Nonlinear.DerivativeEpsilon                 1e-17
pfset Solver.Nonlinear.StepTol				 1e-5
pfset Solver.Nonlinear.Globalization                     LineSearch
pfset Solver.Linear.KrylovDimension                      20
pfset Solver.Linear.MaxRestart                           12
pfset Solver.MaxConvergencFailures                       6

pfset Solver.Linear.Preconditioner                       PFMG
pfset Solver.Linear.Preconditioner.SymmetricMat          Nonsymmetric
pfset Solver.Linear.Preconditioner.MGSemi.MaxIter        1
pfset Solver.Linear.Preconditioner.MGSemi.MaxLevels      10
pfset Solver.PrintSubsurf				False
pfset  Solver.Drop                                      1E-20
pfset Solver.AbsTol                                     1E-6
 
#pfset Solver.LSM                                         CLM
set useclm false   
#set useclm true 

#pfset Solver.CLM.MetForcing                             1D
#pfset Solver.CLM.MetFileName                            Narr-Parflow-July-Aug-2011.txt 
#pfset Solver.WriteSiloCLM                               True	
#pfset Solver.PrintCLM                                   True
##pfset Solver.CLM.ReuseCount                             600
#pfset Solver.CLM.CLMDumpInterval                        1 

pfset Solver.WriteSiloSubsurfData                       True
pfset Solver.WriteSiloPressure                          True
pfset Solver.WriteSiloSaturation                        True
pfset Solver.WriteSiloConcentration                     True

pfset Solver.WriteSiloSlopes                            True
pfset Solver.WriteSiloMask                              True
pfset Solver.WriteSiloEvapTrans                         True
pfset Solver.WriteSiloEvapTransSum                      True
pfset Solver.WriteSiloOverlandSum                       True
pfset Solver.WriteSiloMannings                          True
pfset Solver.WriteSiloSpecificStorage                   True


#---------------------------------------------------------
# Initial conditions: water pressure
#---------------------------------------------------------
# set water table to be at the bottom of the domain, the top layer is initially dry
#set PressDomain [pfload -sa SpinPressVdZ.txt]
#pfsave $PressDomain -pfb "Spin.Press.pfb"
#
#pfset ICPressure.Type                                   PFBFile
#pfset ICPressure.GeomNames                              domain
#pfset Geom.domain.ICPressure.FileName                   Spin.Press.pfb
#puts "Before press dist "
#pfdist Spin.Press.pfb
#puts "After press dist"

pfset ICPressure.Type                                   HydroStaticPatch
pfset ICPressure.GeomNames                              domain
pfset Geom.domain.ICPressure.Value                      -3.95

pfset Geom.domain.ICPressure.RefGeom                    domain
pfset Geom.domain.ICPressure.RefPatch                   z-upper

#-----------------------------------------------------------------------------
# Run and Unload the ParFlow output files
#----------------------------------------------------------------------------

# This prepares the output directory.
set outputdir Spinup_PET_0.00010
file mkdir "$outputdir"

puts "before force copy"
#file copy -force Input_1_slope_short.pfsol $outputdir/Input_1_slope_short.pfsol
file copy -force Xgrid.pfb $outputdir/Xgrid.pfb 
file copy -force Ygrid.pfb.dist $outputdir/Xgrid.pfb.dist
file copy -force Ygrid.pfb $outputdir/Ygrid.pfb
file copy -force Ygrid.pfb.dist $outputdir/Ygrid.pfb.dist
#file copy -force drv_clmin.dat $outputdir/drv_clmin.dat
#file copy -force drv_vegm.dat $outputdir/drv_vegm.dat
file copy -force Spin.Press.pfb $outputdir/Spin.Press.pfb
file copy -force Spin.Press.pfb.dist $outputdir/Spin.Press.pfb.dist
puts "After force copy"

## With CLM:
## This enables CLM with its own output directory and divides outputs into subdirectories.
#if { $useclm == "true" } {
#file mkdir "$outputdir/$runname.clm_output"
#file mkdir "$outputdir/$runname.clm_output/eflx_lh_tot"
#file mkdir "$outputdir/$runname.clm_output/eflx_lwrad_out"
#file mkdir "$outputdir/$runname.clm_output/eflx_sh_tot"
#file mkdir "$outputdir/$runname.clm_output/eflx_soil_grnd"
#file mkdir "$outputdir/$runname.clm_output/qflx_evap_grnd"
#file mkdir "$outputdir/$runname.clm_output/qflx_evap_soi"
#file mkdir "$outputdir/$runname.clm_output/qflx_evap_tot"
#file mkdir "$outputdir/$runname.clm_output/qflx_evap_veg"
#file mkdir "$outputdir/$runname.clm_output/qflx_infl"
#file mkdir "$outputdir/$runname.clm_output/qflx_top_soil"
#file mkdir "$outputdir/$runname.clm_output/qflx_tran_veg"
#file mkdir "$outputdir/$runname.clm_output/swe_out"
#file mkdir "$outputdir/$runname.clm_output/t_grnd"
#file mkdir "$outputdir/$runname.clm_output/diag_out"
#
#pfset Solver.CLM.CLMFileDir   "$runname.clm_output/"
#
#set np [expr [pfget Process.Topology.P] * [pfget Process.Topology.Q] * [pfget Process.Topology.R]]
#for {set fileno 0} {$fileno <= [expr $np - 1]} {incr fileno} {
#	file copy -force drv_clmin.dat $outputdir/drv_clmin.dat.$fileno
#        file copy -force drv_vegm.dat $outputdir/drv_vegm.dat.$fileno
#        file copy -force Narr-Parflow-July-Aug-2011.txt $outputdir/Narr-Parflow-July-Aug-2011.txt.$fileno
#        file copy -force Narr-Parflow-July-Aug-2011.txt $outputdir/Narr-Parflow-July-Aug-2011.txt
#}
#file copy -force drv_vegp.dat $outputdir/drv_vegp.dat
#}

cd "$outputdir"

#pfdist Xgrid.pfb
#pfdist Ygrid.pfb
pfwritedb $runname
#pfdist litter20sl.out.press.03651.pfb
pfrun $runname
pfundist $runname
pfundist Xgrid.pfb
pfundist Ygrid.pfb

if { $useclm == "true" } {
file delete "drv_vegp.dat"
for {set fileno 0} {$fileno <= [expr $np - 1]} {incr fileno} {   
	file delete "drv_clmin.dat.$fileno"
	file delete "drv_vegm.dat.$fileno"
	file delete "Narr-Parflow-July-Aug-2011.txt.$fileno"
}}
#file delete Input_1_slope_short.pfsol

##===============================================
######  WATER BALANCE----------------------------
##===============================================

# SWITCH - uncomment "set runwaterbalance true" to run the water balance post-processor.
set runwaterbalance false
#set runwaterbalance true


if { $runwaterbalance == "true" } {
set slope_x          [pfload $runname.out.slope_x.silo]
set slope_y          [pfload $runname.out.slope_y.silo]
set mannings         [pfload $runname.out.mannings.silo]
set specific_storage [pfload $runname.out.specific_storage.silo]
set porosity         [pfload $runname.out.porosity.silo]
set mask             [pfload $runname.out.mask.silo]
set top              [pfcomputetop $mask]
set nx [expr [pfget ComputationalGrid.NX]]
set ny [expr [pfget ComputationalGrid.NY]]
set nz [expr [pfget ComputationalGrid.NZ]]
set tt [expr ([pfget TimingInfo.StopTime] - [pfget TimingInfo.StartTime]) / [pfget TimeStep.Value]]
set surface_area_of_domain [expr [pfget ComputationalGrid.DX] * [pfget ComputationalGrid.DY] * [pfget ComputationalGrid.NX] * [pfget ComputationalGrid.NY]]

set CLMinterval 6
set sum_total_surface_runoff 0

set prev_total_water_balance 0.0
puts "======================================================"
puts "Timestep  SurfStor  SubsurfStor  Runoff       BCFlux     Total_Water  Expected_waterBalance Percent Diff"

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
    pfsave $saturation -sa "Saturation.$i.txt"    

    set subsurface_storage [pfsubsurfacestorage $mask $porosity $pressure $saturation $specific_storage]
    pfsave $subsurface_storage -silo "subsurface_storage.$i.silo"
    set total_subsurface_storage [pfsum $subsurface_storage]
    set total_water_in_domain [expr $total_water_in_domain + $total_subsurface_storage]


    set surface_runoff [pfsurfacerunoff $top $slope_x $slope_y $mannings $pressure]
    pfsave $surface_runoff -silo "surface_runoff.$i.silo"
    pfsave $surface_runoff -sa "Surface_runoff.$i.txt"
    set total_surface_runoff [expr [pfsum $surface_runoff] * [pfget TimeStep.Value]]
    set sum_total_surface_runoff [expr $sum_total_surface_runoff + $total_surface_runoff]
    #puts [format "Surface runoff\t\t\t\t\t : %.16e" $total_surface_runoff]

##|||Define BC Flux based on simulation inputs and outputs

## Use constant flux from defined boundary conditions:
#    set bc_flux [pfget Patch.top.BCPressure.alltime.Value]
#    set boundary_flux [expr $bc_flux * $surface_area_of_domain * [pfget TimeStep.Value]]

#Use ET flux array with CLM:
if { $useclm == "true" } {
    set bc_flux 0.0
    if { $i > 0 } {
	set filename [format "%s.out.evaptranssum.%05d.silo" $runname $i]
	set evaptrans [pfload $filename]
        pfsave $evaptrans -sa "EvapTransSum.$i.txt"
	set bc_flux [expr $bc_flux - [pfsum $evaptrans]]

        set filename [format "%s.out.evaptrans.%05d.silo" $runname $i]
        set evaptransNotsum [pfload $filename]
        pfsave $evaptransNotsum -sa "EvapTrans.$i.txt"

       set filename [format "%s.out.qflx_evap_grnd.%05d.silo" $runname $i]
       set fexist [file exist $filename]

       if {$fexist==1} {
         set filename [format "%s.out.qflx_evap_grnd.%05d.silo" $runname $i]
         set EvapGrnd [pfload $filename]
         pfsave $EvapGrnd -sa "EvapGrnd.$i.txt"
         set filename [format "%s.out.qflx_evap_soi.%05d.silo" $runname $i]
         set EvapSoil [pfload $filename]
         pfsave $EvapSoil -sa "EvapSoil.$i.txt"
         set filename [format "%s.out.qflx_evap_veg.%05d.silo" $runname $i]
         set EvapCanopy [pfload $filename]
         pfsave $EvapCanopy -sa "EvapCanopy.$i.txt"
         set filename [format "%s.out.qflx_evap_tot.%05d.silo" $runname $i]
         set EvapTot [pfload $filename]
         pfsave $EvapTot -sa "EvapTot.$i.txt"
         set filename [format "%s.out.qflx_tran_veg.%05d.silo" $runname $i]
         set Tran [pfload $filename]
         pfsave $Tran -sa "VegTrans.$i.txt"
       }
    }
    
    set boundary_flux [expr $bc_flux * [pfget TimeStep.Value]]
} else {
# Use fluxes defined for timesteps in defined boundary conditions:
    if [expr $i < 7] {
	set bc_flux [pfget Patch.z-upper.BCPressure.rain.Value]
    } {
	set bc_flux [pfget Patch.z-upper.BCPressure.rec.Value]
    }
    set boundary_flux [expr $bc_flux * [pfget TimeStep.Value]]
}


    # Note flow into domain is negative
    set expected_difference [expr $boundary_flux + $total_surface_runoff]

    if { $i > 0 } {


	if [expr $expected_difference != 0.0] {
	    set percent_diff [expr (abs(($prev_total_water_balance - $total_water_in_domain) - $expected_difference)) / abs($expected_difference) * 100]
	}

	set expected_water_balance [expr $prev_total_water_balance - $expected_difference]
	set percent_diff [expr abs(($total_water_in_domain - $expected_water_balance)) / $expected_water_balance * 100]

#	if [expr $percent_diff > 0.005] {
#	    puts [ format "Error: Water balance is not correct %05d" $i]
#	    set passed 0                          
#	} 

    }

    set prev_total_water_balance [expr $total_water_in_domain]
    if { $i > 0 } {
	    puts [format "%05d  %.05e  %.05e  %.05e  %.05e  %.05e %.05e %.05e %.05e" $i $total_surface_storage $total_subsurface_storage $total_surface_runoff $boundary_flux $total_water_in_domain $expected_water_balance $percent_diff $sum_total_surface_runoff]
    } else {
	    puts [format "%05d  %.05e  %.05e  %.05e  %.05e  %.05e" $i $total_surface_storage $total_subsurface_storage $total_surface_runoff $boundary_flux $total_water_in_domain]
    }
}
puts "\n\n"
puts [format "%.05e" $sum_total_surface_runoff]
set outfile [open "Sum_Total_Surface_runoff.txt" w]
puts $outfile "$sum_total_surface_runoff"
close $outfile
puts "\n\n"
}



# To finish, we return to the original directoryE.
cd ".."


