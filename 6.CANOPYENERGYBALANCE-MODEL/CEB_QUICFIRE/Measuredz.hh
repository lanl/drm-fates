#ifndef SNOW_ENERGY_BALANCE_
#define SNOW_ENERGY_BALANCE_
#include <sstream>
#include "Measuredz.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFRACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vapor pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
***************************** */

namespace Measuredz {


struct EnergyBalance {

  double Dt;

};

// MAIN LEAF ENERGY BALANCE FUNCTION
void Measuredz (LocalData& mdz);

}// Namespace

#endif
