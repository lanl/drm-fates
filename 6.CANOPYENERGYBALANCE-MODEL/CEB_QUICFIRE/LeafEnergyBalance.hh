#ifndef SNOW_ENERGY_BALANCE_
#define SNOW_ENERGY_BALANCE_
#include <sstream>
#include "LeafEnergyBalance.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFRACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vapor pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
***************************** */

namespace SurfaceEnergyBalance {

struct VaporPressure {
  double temp;
  double relative_humidity;
  double saturated_vaporpressure;
  double actual_vaporpressure;
  double dewpoint_temp;
};

struct EnergyBalance {
  double air_temp;
  double relative_rumidity;
  double Us;
  double fQswIn;
  double Pr;

  double Dt;
  double WaterStorage;
  double max_water;

  double dewpoint_temp;
  double albedo_value;

  double stephB;
  double Apa;
  double SEs;
  double Dhe;
  double Sqig;
  double gZr;
  double rowaCp;
  double rowaLe;
  double density_w;
  //    double density_air;
  double Hf;
  double Ls;
  double Le;
  double VKc;
  //    double Cp;
  double Zr;
  double Zo;

  double fQlwIn;
  double QswIn;
  double fQlwOut;
  double fQh;
  double fQe;
  double Ts;

  double EvR;
  double EvL;
  double lai;

  double varvar;
  std::string funcall;

};

//vapor_pressure * vp_pointer;
struct LocalData {
  VaporPressure vp_air;
  VaporPressure vp_surface;
  EnergyBalance st_energy;

};

void CalcEFluxTempIndependent (LocalData& seb);
void CalcEFluxTempDependent (LocalData& seb, double T);

// #### FUNCTIONS TO CALCULATE SATURATION VAPOR PRESSURE & ACTUALL VAPOR PRESSURE
void VaporCalc (VaporPressure& vp);


// #### FUNCTIONS TO CALCULATE ALBEDO
void AlbedoCalc (EnergyBalance& eb);


/*  FUNCTION TO SOLVE SNOW SURFACE ENERGY BALANCE FOR SURFACE TEMP (Ts)
    ###############################################################################################*
    Bisection Method
    (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) + Qc(Ts) +  Qe(Ts) = 0
    Substitute Xx for all Ts
    ############################################################################################## */
double BisectionZeroFunction(LocalData& dat, double Xx);
void BisectionEnergyCalc (LocalData& dat);

//  FUNCTION TO CALCULATE EVAPORATION & CONDENSATION RATE 
void EvapCalc (EnergyBalance& eb);

// FUNCTION TO CHANGE WATER STORAGE BASED ON EVAPORATION OR CONDENSATION AND PRECIPITATION.
void StorageCalc (EnergyBalance& eb);

// FUNCTION TO MAKE SURE PROPER MASS OF WATER IS DELEVERED TO ATS WHEN SNOWPACK DISAPEARS
void MaxStor (EnergyBalance& eb);

// MAIN LEAF ENERGY BALANCE FUNCTION
void LeafEnergyBalance (LocalData& seb);

}// Namespace

#endif
