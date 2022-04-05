#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdio.h>

#include "LeafEnergyBalance.hh"
#include "spa.h"  //include the SPA header file
/* THIS CODE WILL CALCULATE THE SNOWSURFACE ENERGY BALANCE THE SNOW SURFRACE TEMPERUATURE.

 **** Incomming Longwave radation is cacualted in this version, but if data is avaialbe we could incoroerate it with the avialable met data****
 
+++Atmospheric pressure is often used in snow models, If data is available we could incorperate it but for now Pa is held constant at
 100 pa++++
  
*** Equation for saturated vapor pressure over water is taken from Bolton, 1980 'Monthly Weather Review'
*** Equation for saturated vapor pressure over snow is taken from Buck, 1996 'Buck Research Manual'
 See: http://cires.colorado.edu/~voemel/vp.html
 *****************************
 */

using namespace std;

double MeasureDZ (int Zcell){ // Function to determine the vertical size (dz) of a cell.
  double a = 0.1, l = 41.0, dz = 15.0;
  double aaa = (1 - a)/pow((l*dz),2); 
  double bb = (Zcell)*dz;
  double Height = (aaa * pow(bb,3)) + a*bb;
  return Height;
}

struct Direction
{
  double sin_X, sin_Y, azimuth, XBottomCorner, YBottomCorner, XtopCorner, YtopCorner, DX, DY;
  string corner;
};

void XYagnle (Direction& ang_){
  double PI=3.14159265;
  if (ang_.azimuth<90) {
     ang_.sin_X = -sin(ang_.azimuth*PI/180), ang_.sin_Y = -sin((90-ang_.azimuth)*PI/180);
     ang_.corner = "NE", ang_.XtopCorner = -ang_.DX, ang_.YtopCorner = -ang_.DY;
     ang_.XBottomCorner = ang_.DX, ang_.YBottomCorner = ang_.DY;
  }
  if ((ang_.azimuth>=90)&&(ang_.azimuth<180)) {
     ang_.sin_X = -sin((90-(ang_.azimuth-90))*PI/180), ang_.sin_Y = sin((ang_.azimuth-90)*PI/180);
     ang_.corner = "SE", ang_.XtopCorner = -ang_.DX, ang_.YtopCorner = ang_.DY;
     ang_.XBottomCorner = ang_.DX, ang_.YBottomCorner = 0;
  }
  if ((ang_.azimuth>=180)&&(ang_.azimuth<270)) {
     ang_.sin_X = sin((ang_.azimuth-180)*PI/180), ang_.sin_Y = sin((90-(ang_.azimuth-180))*PI/180);
     ang_.corner = "SW", ang_.XtopCorner = ang_.DX, ang_.YtopCorner = ang_.DY;
     ang_.XBottomCorner = 0, ang_.YBottomCorner = 0;
  }
  if (ang_.azimuth>=270) {
     ang_.sin_X = sin((90-(ang_.azimuth-270))*PI/180), ang_.sin_Y = -sin((ang_.azimuth-270)*PI/180);
     ang_.corner = "NW", ang_.XtopCorner = ang_.DX, ang_.YtopCorner = -ang_.DY;
     ang_.XBottomCorner = 0, ang_.YBottomCorner = ang_.DY;
  }
}

int main(int argc, char* argv[])
{
    //Required input values into SPA structure
    spa_data spa;  //declare the Sun Position Analysis (SPA) structure
    // Time Parameters
    spa.year          = 2019;
    spa.month         = 0;
    spa.day           = 0;
    spa.hour          = 0;
    spa.minute        = 0;
    spa.second        = 0;
    // Location Parameters
    spa.longitude     = -86.500379;
    spa.latitude      = 30.606949;
    spa.elevation     = 45.7;
    spa.timezone      = -7.0;
    //Other Stuff
    spa.delta_ut1     = 0;
    spa.delta_t       = 67;
    spa.pressure      = 820;
    spa.temperature   = 11;
    spa.slope         = 0;
    spa.azm_rotation  = -10;
    spa.atmos_refract = 0.5667;
    spa.function      = SPA_ALL;

    Direction ang_; 

    //Constant Parameters for Energy Equaions
    double ROWa = 1.275;                  // Density of Air ------------------------- [kg/m^3]
    double ROWw = 1000;                   // Density of Water ----------------------- [kg/m^3]
    double Ls = 2834000.000000;           // Latent heat of sublimation ------------- [J/kg]  *** Check this kJ to J conversion In equations***
    double Hf = 333500.0;                 // Heat of fusion for melting snow -------- [J/kg] *** Check this kJ to J conversion In equations***
    double Le = 2497848.;                 // Latent heat of vaporization ------------ [J/kg]  *** Check this kJ to J conversion In equations***
    double StephBolt = 0.000000056703730; // Stephan-boltzmann Constant ------------- [W/m^2 K^4]
    double VKc = 0.41;                    // Von Karman Constant -------------------- [-]
    double Cp = 1004.0;                   // Specific heat of air ------------------- [J/K kg]
    double g = 9.807;                     // Acceleration due to gravity ------------ [m/s^2]
    double Zr = 2.0;                      // Referance ht of wind speed ------------- [m]
    double Zo = 0.005;                    // Roughness length  ---------------------- [m]  *** I don't know what this number should be ****
    double Apa = 100;                     // Atmospheric Pressure ------------------- [KPa]
    double SEs = 0.98;                    // Emissivity for grass  ----------- [-] ** https://kopernio.com/viewer?doi=10.1016/S0034-4257(96)00123-X&route=6
    double SEtun = 0.92;                  // Surface Emissivity for tundtra  -------- [-] ** Taken from P. ReVelle (Thesis)
    
    double Ta = 0.0;
    double RH = 0.0;
    double WaterStorage = 0.0;
    double QswIn = 0.0;
    double QlwIn = 0.0;
    double Us = 0.0;
    double Dt = 600;
    double Pr = 0.0;
    double LAI = 0.0, ShadeF = 0.0, FuelDensity = 0.0;
    int fcellnum, numRel=0, NZ=0, NY=0, NX=0; 
    double TIME = 0.0;
    double djunk=0.0;    

    // Domain Size and Dimensions **** Need to figure out how to do changing DZ 
    double DZ = 2.0, DY = 2.0, DX = 2.0;  
    double rhomicro = 500.0, sizescale = 0.0005; 
 
string trash, AMout, LD, MD, SD, FD, TTS1, TTS2, TTS3, FS;

SurfaceEnergyBalance::LocalData st_data;
    
//FILL IN DATA STRUCTURE NEEDED FOR PARAMETERS IN ENERGY BALANCE
st_data.st_energy.density_w=ROWw;
st_data.st_energy.Hf=Hf;
st_data.st_energy.Apa=Apa;
st_data.st_energy.stephB=StephBolt;
st_data.st_energy.gZr=g*Zr;
st_data.st_energy.rowaCp=ROWa*Cp;
st_data.st_energy.rowaLe=ROWa*Le;
st_data.st_energy.Ls=Ls;
st_data.st_energy.Le=Le;
st_data.st_energy.SEs=SEs;
st_data.st_energy.VKc=VKc;
st_data.st_energy.Zr=Zr;
st_data.st_energy.Zo=Zo;

st_data.vp_surface.relative_humidity = 1;

//MAKE MOCK DOMAIN of 4x3x3 cells. [z,x,y]
ifstream FuelData;
FD = "FuelData.txt";
FuelData.open(FD.c_str());
if (!FuelData) {
    cout <<"ERROR opening " <<FD<<" input file"<<endl;
} getline(FuelData, trash);
FuelData>>NZ>>NY>>NX;

/* SECTION 1. USE FIRETEC DOMAIN TO DETERMIN ARRAY STRUCTUR
   #########################################################################
   Loops through FireTEC domain to find all cells with fuel.
   Sets array sizes and locations for Shade, LAI, and WaterStored.
   Each array element will have an x,y,z location
   #########################################################################
*/
std::cout<<"FuelDensity  dn  X-loc  Y-loc  Z-loc"<<std::endl;
int ArrayCount = 0;
for (int dn; dn < (NZ*NY*NX); dn++) { // LOOP TO COUNT NUMBER OF CELLS WITH FUEL
    FuelData>>FuelDensity;  
    if (FuelDensity>0){
    ArrayCount++;
    }
}
FuelData.close();
double LaiArray [ArrayCount];
double Xloc [ArrayCount];
double Yloc [ArrayCount];
double Zloc [ArrayCount];
double DZarray [ArrayCount];
double UnderZ [ArrayCount];
double FuelMoisture [ArrayCount];
double QswFactor [ArrayCount];
FuelData.open(FD.c_str());
if (!FuelData) {
    cout <<"ERROR opening " <<FD<<" input file"<<endl;
} getline(FuelData, trash);
getline(FuelData, trash);
int ArrayIndex = 0;
for (int dn; dn < (NZ*NY*NX); dn++) { // LOOP TO ASSING DATA, LOCATION & LAI OF EACH CELL WITH FUEL
    FuelData>>FuelDensity; // ADD FUNCTION HERE *#! *#! *#! *#! *#! *#! *#! *#! *#!
    if (FuelDensity>0){
       Xloc[ArrayIndex] = (dn%NX)*DX;
       Yloc[ArrayIndex] = ((dn/NX)%NY)*DY;
       Zloc[ArrayIndex] = ((dn/(NX*NY))%NZ); // This is not the streched Z location Needs to measure DZ.
       Zloc[ArrayIndex] = MeasureDZ(Zloc[ArrayIndex]); // Bottom of cell.
       UnderZ[ArrayIndex] = ((dn/(NX*NY))%NZ);
       FuelMoisture[ArrayIndex] = 0.0, QswFactor[ArrayIndex] = 0.0;
       // Call Function to determine DZ
       DZarray[ArrayIndex] = MeasureDZ(((dn/(NX*NY))%NZ)+1) - Zloc[ArrayIndex];  
       // Determine LAI   
       LaiArray[ArrayIndex] = FuelDensity*DZarray[ArrayIndex] / (rhomicro * sizescale); // Consider Reading Sizescale form file. 
       //std::cout<<ArrayIndex<<"   "<<(dn%NX)*DX<<"     "<<((dn/NX)%NY)*DY<<"     "<<((dn/(NX*NY))%NZ)*DZ<<"  "<<Zloc[ArrayIndex]<<"  "<<DZarray[ArrayIndex]<<"  "<<LaiArray[ArrayIndex]<<std::endl;
       //std::cout<<Zloc[ArrayIndex]<<"  "<<DZarray[ArrayIndex]<<std::endl;
       ArrayIndex++;
    }
}
FuelData.close();

/* SECTION 2.  CAST SHADE 
   #########################################################################
   Loops through fuel arrays to cast shade based on LAI and the angle of the sun.
   Locates affected cells.
   Sums the affect (in Shade array) of all shade cast from other cells.
   -- Establish the track of the sun
   -- Analize the angle of the sun and fuel is casts shade on
      ++ Sum Shade Factor.
https://rredc.nrel.gov/solar/codesandalgorithms/spa/
   #########################################################################
*/
int EndTime = 14;
    spa.month         = 6;
    spa.day           = 10;
    spa.hour          = 0;
double min=0.0, sec=0.0, daylight=0.0, zenith=0.0, azimuth=0.0, sunrise=0.0, sunset=0.0;
double ZsdB=0.0, YsdB=0.0, XsdB=0.0;
double ZsdT=0.0, YsdT=0.0, XsdT=0.0; 

int sunresult;
int f = 0;
ang_.DX = DX, ang_.DY = DY; 
for (int tt; tt < EndTime; tt++) {
    // FINDING POSITION OF SUN BASED ON TIME
    spa.hour = tt%24;
    spa.day = 10 + (tt/24)%24;
    //std::cout<<tt<<"  "<<tt%24<<"  "<<1 + (tt/24)%24<<std::endl; 
    sunresult = spa_calculate(&spa);
     if (sunresult == 0) {
     zenith = spa.zenith;
     azimuth = spa.azimuth;
     sunrise = spa.sunrise;
     sunset = spa.sunset;
     daylight = spa.sunset - spa.sunrise;
     min = 60.0*(daylight - (int)(daylight));
     sec = 60.0*(min - (int)min);
     ang_.azimuth=spa.azimuth;
      //std::cout<<tt<<"  "<<spa.day<<"  "<<tt%24<<"  Azimuth: "<<azimuth<<"  Zenith: "<<zenith;
     //printf("  Daylight:   %02d:%02d:%02d \n", (int)(daylight), (int)min, (int)sec); 
     }else printf("SPA Error Code: %d\n", sunresult);
     if (zenith<91) { // Sun is above the horizon Shade will happen

        XYagnle(ang_); // Function to find the sin for X and Y based on angle of the sun. --> Also finds corner (side sun is shining) 

        std::cout<<tt<<"  "<<spa.day<<"  "<<tt%24<<"  Azimuth: "<<azimuth<<"  Zenith: "<<zenith<<"  "<<ang_.sin_X<<"  "<<ang_.sin_Y<<"  ********************************"<<std::endl;
        //printf("  Daylight:   %02d:%02d:%02d \n", (int)(daylight), (int)min, (int)sec);

        // LOOP THROUGH FUEL HERE & CAST SHADDOW
	for (f=0; f < ArrayCount; f++){// Fuel Loop  --> Calculate Shade by shooting lines from ground to sun.  Intercepted fuel causes shade.
            // PUT FUNCTION TO CAST SHADOW HERE *#! *#! *#! *#! *#! *#! *#! *#! *#! *#! *#! *#! *#! *#!
            /*      --- Still Need To Do ---
              1) Looping Through Fuel cells.
              2) Make a sub Z loop -- to find shade at each level below shadding cell.
            */ 
            //std::cout<<"Z cell, hieght:  "<<Zloc[f]<<"  "<<DZarray[f]<<"  Xbot to Xtop   Ybot to Y top      X  Y  Z"<<std::endl; 
            // Make it so that it finds the X, Y shade length at Z cell center.
            ZsdB = Zloc[f]*tan((zenith)*3.14159265/180);
            XsdB = ZsdB*ang_.sin_X, YsdB = ZsdB*ang_.sin_Y; // Find length of cell bottom X and Y shadow side, ZsdB is hypotenuse
            ZsdT = (Zloc[f]+DZarray[f])*tan((zenith)*3.14159265/180);
            XsdT = ZsdT*ang_.sin_X, YsdT = ZsdT*ang_.sin_Y; // Find length of cell top X and Y shadow side, ZsdT is hypotenuse
            XsdB = XsdB + Xloc[f] + ang_.XBottomCorner, XsdT = XsdT + Xloc[f] + ang_.XtopCorner; 
            YsdB = YsdB + Yloc[f] + ang_.YBottomCorner, YsdT = YsdT + Yloc[f] + ang_.YtopCorner;
            int uz = 0;

std::cout<<f<<"  "<<Zloc[f]<<"  "<<DZarray[f]<<"  "<<ZsdB<<"  "<<Xloc[f]<<" "<<XsdB<<" "<<XsdT<<"  "<<Yloc[f]<<" "<<YsdB<<" "<<YsdT<<std::endl;
            //for (uz = 0; uz < ArrayCount; uz++) {
//std::cout<<f<<" "<<uz<<"  "<<Zloc[f]<<"  "<<DZarray[f]<<"  "<<ZsdB<<"  "<<Xloc[f]<<" "<<XsdB<<" "<<XsdT<<"  "<<Yloc[f]<<" "<<YsdB<<" "<<YsdT<<"      "<<Xloc[uz]<<"  "<<Yloc[uz]<<"  "<<Zloc[uz]<<std::endl;
                //if (Xloc[uz]>) {
                //
                //}

//            } 
            /*     ****   NEED TO DO   *****
             1) CALCUALTE SHADE FACTOR BASED ON LAI
             2) SUM SHADE FACTOR
            */ 
        } // End Fuel Loop
     } 
}// End Time Lopp

/* SECTION 3.  SOLVE SURFACE ENERGY BALANCE 
   #########################################################################
   Solves the surface energy balance for each cell in domain.  
   Tracks water stored as interception or dew and leaf temperature.
   ######################################################################### 
*/
for (int f; f < ArrayIndex; f++) {
  stringstream fcellconverter;
  fcellconverter << f;
  string strFcell = fcellconverter.str(); // Creates naming structure for relization files used to open input files.

// METEROLOGICAL INPUT FILE
  ifstream MetData;
  MD = "MetData.3.txt";
  MetData.open(MD.c_str());
  if (!MetData){
      cout <<"ERROR opening "<<MD<<" input file"<<endl;
  } getline (MetData,trash);
  MetData>>numRel;
  
  //SHADE FACTOR FILE
  ifstream ShadeData;
  SD = "ShadeData." + strFcell + ".txt";
  ShadeData.open(SD.c_str());
  if (!ShadeData){
      cout <<"ERROR opening "<<SD<<" input file"<<endl; 
  } getline (ShadeData,trash);  

  // OUTPUT FILES
  //Creating output file
  ofstream Amanzistream;
  AMout = "Amanzi." + strFcell + ".out.txt";
  Amanzistream.open(AMout.c_str());
  Amanzistream<<"TIME   Ts     WaterStorage"<<endl;
  
  ofstream TimeSeires1stream;
  TTS1="TempSeires." + strFcell + ".out.txt";
  TimeSeires1stream.open(TTS1.c_str());
  TimeSeires1stream<<"Tstep   Ta(C)    Ts(C)  WaterStored[m]    Qe    TIME[D] "<<endl;
  TimeSeires1stream.setf(ios::fixed, ios::floatfield);
  TimeSeires1stream.precision(3);
  TimeSeires1stream.width(10);
  
  ofstream TimeSeires2stream;
  TTS2="EvapSeires." + strFcell + ".out.txt";
  TimeSeires2stream.open(TTS2.c_str());
  TimeSeires2stream<<"Tstep   RH   Ea    ESS   Us   TIME[D] "<<endl;
  TimeSeires2stream.setf(ios::fixed, ios::floatfield);
  TimeSeires2stream.precision(6);
  TimeSeires2stream.width(10);
  
  ofstream TimeSeires3stream;
  TTS3="EnergySeires." + strFcell + ".out.txt";
  TimeSeires3stream.open(TTS3.c_str());
  TimeSeires3stream<<"Tstep   QswIn     QlwIn  QlwOut  Qh  Qe  Dh"<<endl;
  TimeSeires3stream.setf(ios::fixed, ios::floatfield);
  TimeSeires3stream.precision(3);
  TimeSeires3stream.width(10);

  // Set Parameters For TIME LOOP iterations 
  //Initial Conditions = Initial WaterStored[m] ~ Assume Fully Saturated
  //Fuel Data
  st_data.st_energy.lai = LaiArray[f];
  //Initial Conditions = Initial WaterStored[m] ~ Assume Fully Saturated
  st_data.st_energy.WaterStorage = st_data.st_energy.lai * 0.001;
  for (int rel=1;rel<=numRel;rel++) { //TIME LOOP LOOP LOOP #####################
      MetData  >>djunk>>QswIn>>QlwIn>>Ta>>RH>>Us>>Pr;
      // Read Shade Factor here
      ShadeData>>ShadeF;  
      // Shade Factor Calculation ********************************** !!!!!!!!!!!!!!!!!!
      QswIn = QswIn*(1-ShadeF); 
      // Assing all meterological data -- Converting Celcius to Kelvin for all temp Inputs
      st_data.vp_air.temp = Ta + 273.15;                 // MET DATA 13
      st_data.vp_air.relative_humidity=RH;      // MET DATA 12
      st_data.st_energy.air_temp= Ta + 273.15; // MET DATA         5
      st_data.st_energy.Us=Us;  // MET DATA              6
      st_data.st_energy.Pr=Pr;  // MET DATA              8
      st_data.st_energy.Dt=Dt;  // ATS CACLC             10
      st_data.st_energy.QswIn=QswIn; // MET DATA 15
      st_data.st_energy.fQlwIn=QlwIn;

      // Run SEB        
      SurfaceEnergyBalance::LeafEnergyBalance(st_data);
          
      // Saving Snow data for next timestep
          
      // Calculating Data for ATS
      
      TIME = TIME + (Dt/600);  // Time in Days
   
      Amanzistream<<setw(4)<<rel<<"  "<<TIME<<"  "<<setw(3)<<st_data.st_energy.Ts-273.15<<"   "<<st_data.st_energy.WaterStorage<<endl;
      TimeSeires1stream<<setw(4)<<rel<<"    "<<setw(3)<<Ta<<"  "<<st_data.st_energy.Ts-273.15<<"   "<<st_data.st_energy.WaterStorage<<"  "<<st_data.st_energy.fQe<<"  "<<TIME<<endl;
      TimeSeires2stream<<setw(4)<<rel<<"    "<<setw(4)<<st_data.vp_air.relative_humidity<<"   "<<st_data.vp_air.actual_vaporpressure<<"  "<<st_data.st_energy.EvL<<"  "<<Us<<"  "<<TIME<<endl;
      TimeSeires3stream<<setw(4)<<rel<<"    "<<setw(3)<<st_data.st_energy.fQswIn<<"  "<<st_data.st_energy.fQlwIn<<"  "<<st_data.st_energy.fQlwOut<<"  "<<st_data.st_energy.fQh<<"  "<<st_data.st_energy.fQe<<"  "<<st_data.st_energy.Dhe<<endl;
      
      // Resetting variables for next iteration
      st_data.vp_air.actual_vaporpressure=0.0;
      st_data.st_energy.fQswIn=0.0, st_data.st_energy.fQlwIn=0.0, st_data.st_energy.fQlwOut=0.0, st_data.st_energy.fQe=0.0, st_data.st_energy.fQh=0.0;
      st_data.st_energy.fQe=0.0;
      st_data.st_energy.varvar=0.0;
      
      //Updateing loop variables
    
  }//  END TIME LOOP LOOP LOOP LOOP #####################
  //std::cout<<f<<"  "<<st_data.st_energy.Ts-273.15<<"   "<<st_data.st_energy.WaterStorage<<endl;
  Amanzistream.close();
  ShadeData.close();
  TimeSeires1stream.close();
  TimeSeires2stream.close();
  TimeSeires3stream.close();
  TIME = 0;
}  // End Fuel Cell Loop
 
	return 0;
}
