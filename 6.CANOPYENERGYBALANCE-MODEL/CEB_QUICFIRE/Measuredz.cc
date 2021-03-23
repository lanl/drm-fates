#include <iostream>
#include <cmath>

//#include "dbc.hh"
#include "Measuredz.hh"

/* 
**********************   */

void Measuredz::Measuredz (LocalData& mdz) {
  // 
  double EmissivitySky = 0.787 + (0.7641*log(seb.vp_air.dewpoint_temp/273)); // From http://www.ibpsa.org/proceedings/BS2017/BS2017_569.pdf  Equation 19.
  double Ca = 1;  // 
  double EmissivitySkyClouds = EmissivitySky * Ca; 
  seb.st_energy.fQlwIn = EmissivitySkyClouds*seb.st_energy.stephB*(std::pow(seb.st_energy.air_temp,4)); 
  seb.st_energy.Dhe=(((std::pow(seb.st_energy.VKc,2)*seb.st_energy.Us))/(std::pow(log(seb.st_energy.Zr/seb.st_energy.Zo),2)));
}


