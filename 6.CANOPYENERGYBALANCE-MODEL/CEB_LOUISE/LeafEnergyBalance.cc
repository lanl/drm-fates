#include <iostream>
#include <cmath>

//#include "dbc.hh"
#include "LeafEnergyBalance.hh"

/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFACE TEMPERUATURE.

**** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****

+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
100 pa++++

*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vapor pressure over snow is taken from Buck, 1996 'Buck Research Manual'
See: http://cires.colorado.edu/~voemel/vp.html
**********************   */

void SurfaceEnergyBalance::CalcEFluxTempIndependent (LocalData& seb) {
  // Calculate incoming long-wave radiation
  double EmissivitySky = 0.787 + (0.7641*log(seb.vp_air.dewpoint_temp/273)); // From http://www.ibpsa.org/proceedings/BS2017/BS2017_569.pdf  Equation 19.
  double Ca = 1;  // Find Cloud factor. See example http://www.ibpsa.org/proceedings/BS2017/BS2017_569.pdf  Equation 21.
  double EmissivitySkyClouds = EmissivitySky * Ca; 
  seb.st_energy.fQlwIn = EmissivitySkyClouds*seb.st_energy.stephB*(std::pow(seb.st_energy.air_temp,4)); 
  //seb.st_energy.fQlwIn = 0; 
  // Calculate D_h, D_e
  seb.st_energy.Dhe=(((std::pow(seb.st_energy.VKc,2)*seb.st_energy.Us))/(std::pow(log(seb.st_energy.Zr/seb.st_energy.Zo),2)));
}

void SurfaceEnergyBalance::CalcEFluxTempDependent (LocalData& seb, double T) {
  double Sqig;
  if (seb.st_energy.Us == 0.) {
    Sqig = 0.;
  } else {
    double Ri  = ((seb.st_energy.gZr*(seb.st_energy.air_temp-T))/(seb.st_energy.air_temp*std::pow(seb.st_energy.Us,2)));
    Sqig = (1/(1+10*Ri));
  }
   seb.st_energy.Sqig = Sqig;
  //std::cout<<Sqig<<"  "<<seb.st_energy.air_temp<<"  "<<T<<"  "<<seb.st_energy.air_temp-T<<std::endl; 

  // Calculate outgoing long-wave radiation
  seb.st_energy.fQlwOut = -seb.st_energy.SEs*seb.st_energy.stephB*(std::pow(T,4));
  //std::cout<<"StephB*SEs, SurfTemp, longwaveIn LongwaveOut  "<<seb.st_energy.SEs*seb.st_energy.stephB<<"  "<<T<<"  "<<seb.st_energy.air_temp-T<<"  "<<seb.st_energy.fQlwIn<<"  "<<seb.st_energy.fQlwOut<<std::endl;
  //seb.st_energy.fQlwOut = 0;
  // Calculate sensible heat flux
  seb.st_energy.fQh = seb.st_energy.rowaCp*seb.st_energy.Dhe*Sqig*(seb.st_energy.air_temp-T);

  // Update vapore pressure of leaf
  seb.vp_surface.temp=T;
  VaporCalc (seb.vp_surface);

// CHECK TO SEE IF INTERCEPTED WATER IS PRESENT ON FOLIAGE TO GET CORRECT SURFACE VAPOR PRESSURE & CALCULATE LATENT HEAT EXCHANGE
  if (seb.st_energy.WaterStorage>0.0) {// Checking for standing water
    seb.st_energy.fQe = seb.st_energy.rowaLe*seb.st_energy.Dhe*Sqig*0.622*((seb.vp_air.actual_vaporpressure - seb.vp_surface.saturated_vaporpressure)/seb.st_energy.Apa);
    // ***** Think about limiting fQe here.......*******
    EvapCalc(seb.st_energy);   
    } else {// no standing water
    seb.st_energy.fQe = seb.st_energy.rowaLe*seb.st_energy.Dhe*Sqig*0.622*((seb.vp_air.actual_vaporpressure - seb.vp_surface.actual_vaporpressure)/seb.st_energy.Apa);
    if (seb.st_energy.fQe < 0){
        seb.st_energy.fQe = 0;
    }
    EvapCalc(seb.st_energy);
  }
}

// #### FUNCTIONS TO CALCULATE SATURATION VAPOR PRESSURE & ACTUALL VAPOR PRESSURE #########################
void SurfaceEnergyBalance::VaporCalc (VaporPressure& vp) {
  double temp;
  //Convert from Kelvin to Celsius
  temp = vp.temp-273.15;
  // Sat vap. press o/water Dingman D-7 (Bolton, 1980)
  vp.saturated_vaporpressure = 0.611*std::exp((17.67*temp)/(temp+243.5));
  // (Bolton, 1980)
  vp.actual_vaporpressure=vp.saturated_vaporpressure*vp.relative_humidity;
  // Find dewpoint Temp Dingman D-11
  vp.dewpoint_temp=(std::log(vp.actual_vaporpressure)+0.4926)/(0.0708-0.00421*std::log(vp.actual_vaporpressure));
  // Convert Tdp from Celsius to Kelvin
  vp.dewpoint_temp=vp.dewpoint_temp + 273.15;
}
// #### FUNCTIONS TO CALCULATE ALBEDO
void SurfaceEnergyBalance::AlbedoCalc (EnergyBalance& eb) {
    double perSnow = 0.0, perTundra=0.0, perWater=0.0;
    double AlNeedle=0.14;
    //double TransitionVal=eb.AlbedoTrans;  // Set to 2 cm
    // weighted vaverage function for albedo
    //eb.albedo_value = ((AlSnow*perSnow)+(AlTundra*perTundra)+(AlWater*perWater))/(perSnow + perTundra + perWater);
    eb.albedo_value = AlNeedle; 
}

/*  FUNCTION TO SOLVE SNOW SURFACE ENERGY BALANCE FOR SURFACE TEMP (Ts)
    ###############################################################################################*
    Bisection Method
    (1-albedo)QswIn + Qlwin + QlwOut(Ts) + Qh(Ts) +  Qe(Ts) = 0
    Substitute Xx for all Ts
    ####################################################################################################
*/
double SurfaceEnergyBalance::BisectionZeroFunction(LocalData& seb, double Xx) {
  CalcEFluxTempDependent(seb, Xx);
  // BALANCE EQUATION
  double ZERO = seb.st_energy.fQswIn + seb.st_energy.fQlwIn + seb.st_energy.fQlwOut + seb.st_energy.fQh
      + seb.st_energy.fQe;
  //std::cout<<"Energy Balance  "<<seb.st_energy.fQswIn<<"  "<<seb.st_energy.fQlwIn<<"  "<<seb.st_energy.fQlwOut<<"  "<<seb.st_energy.fQh<<"  "<<seb.st_energy.fQe<<"  "<<Xx<<"  "<<ZERO<<std::endl;
  return ZERO;
}
// BISECTION METHOD
void SurfaceEnergyBalance::BisectionEnergyCalc (LocalData& seb) {
  double tol = 1.e-6;
  double deltaX = 0.05;

  double Xx = seb.st_energy.air_temp;
  double FXx = BisectionZeroFunction(seb, Xx);
  // NOTE: decreasing function
  double a,b,Fa,Fb;
  if (FXx > 0) {
    b = Xx;
    Fb = FXx;
    a = Xx;
    Fa = FXx;
    while (Fa > 0) {
      b = a;
      Fb = Fa;
      a += deltaX;
      Fa = BisectionZeroFunction(seb,a);
    }
  } else {
    a = Xx;
    Fa = FXx;
    b = Xx;
    Fb = FXx;
    while (Fb < 0) {
      a = b;
      Fa = Fb;
      b -= deltaX;
      Fb = BisectionZeroFunction(seb,b);
    }
  }
  int maxIterations = 200;
  double res;
  int iter;
  for (int i=0; i<maxIterations; ++i) { //Besection Iterations Loop Solve for Ts using Energy balance equation #######
    Xx = (a+b)/2;
    res = BisectionZeroFunction(seb, Xx);
    if (res>0) {
      b=Xx;
    } else {
      a=Xx;
    }
    if (std::abs(res)<tol) {
      break;
    }
    iter=i;
  }//End Bisection Interatiion Loop  Solve for Ts using Energy balance equation ##################################
  if (std::abs(res) >= tol) {
  //  ASSERT(0);
  }
  seb.st_energy.Ts=Xx;
}
// FUNCTION TO CALCULATE EVAPORATION & CONDENSATION RATE 
void SurfaceEnergyBalance::EvapCalc (EnergyBalance& eb) {
  eb.EvR=eb.fQe/(eb.density_w*eb.Le); // EvR=Qe/(ROWw*Le); // Evaporation or condensation rate [m/s]
  eb.EvL=eb.EvR*eb.Dt;      // EvL=EvR*Dt;           // One dimmential Evap or condensation lenght [m] for a timestep
  //std::cout<<"IN THE ENERGY BALANCE:  "<<eb.fQe<<"  "<<eb.EvL;
  if ((eb.EvL * -1) > eb.WaterStorage) {// Evaporates all the water, limits Lantent Heat to avialable water 
    eb.EvR= - eb.WaterStorage/eb.Dt; // Limiting rate based on water avialalbe
    eb.fQe = eb.EvR * eb.density_w * eb.Le;  // Setting Lantent Heat
    //std::cout<<"  Water Limited:  "<<eb.fQe<<"  "<<eb.EvL<<std::endl;
  }
}
// FUNCTION TO UPDATE WATER STORAGE BASED ON LATENT HEAT
void SurfaceEnergyBalance::StorageCalc (EnergyBalance& eb) {
  eb.WaterStorage = eb.WaterStorage + eb.EvL + eb.Pr; // Adjusts water storage to evaporation or condensation
  if (eb.WaterStorage<0) {eb.WaterStorage = 0;} 
}
// FUNCTION TO LIMIT WATER STORAGE BASED ON LAI
void SurfaceEnergyBalance::MaxStor (EnergyBalance& eb) {
   eb.max_water = eb.lai * 0.001; // Assumes all surfaces can store 1mm of water.
   if (eb.WaterStorage > eb.max_water){
      eb.WaterStorage = eb.max_water;
   }
}


// MAIN SNOW ENERGY BALANCE FUNCTION
void SurfaceEnergyBalance::LeafEnergyBalance (LocalData& seb) {
// Caculate Vapor pressure and dewpoint temperature from Air
  VaporCalc(seb.vp_air);

// Find Albedo
  AlbedoCalc(seb.st_energy); // Calculating Surface Albedo
  seb.st_energy.fQswIn=(1-seb.st_energy.albedo_value)*seb.st_energy.QswIn;

// Update temperature-independent fluxes
  CalcEFluxTempIndependent(seb);

//Solving for surface temperature and Energy balance equation
  BisectionEnergyCalc(seb);

// Update water Storage
  StorageCalc(seb.st_energy);   
 
// Update Check water storage
  MaxStor(seb.st_energy);
}
