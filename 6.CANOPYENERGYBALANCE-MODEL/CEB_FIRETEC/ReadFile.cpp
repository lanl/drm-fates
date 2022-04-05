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


struct Direction
{
  double sin_X, sin_Y, azimuth, XBottomCorner, YBottomCorner, XtopCorner, YtopCorner, DX, DY;
  string corner;
};


int main(int argc, char* argv[])
{
    //Required input values into SPA structure


    //Constant Parameters for Energy Equaions
    double LAI = 0.0, ShadeF = 0.0, FuelDensity = 0.0;
    int fcellnum, numRel=0, NZ=0, NY=0, NX=0; 
    double TIME = 0.0;
    double djunk=0.0;    

    // Domain Size and Dimensions **** Need to figure out how to do changing DZ 
    double DZ = 2.0, DY = 2.0, DX = 2.0;  
    double rhomicro = 500.0, sizescale = 0.0005; 
 
string trash, AMout, LD, MD, SD, FD, TTS1, TTS2, TTS3, FS;

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


/* SECTION 3.  SOLVE SURFACE ENERGY BALANCE 
   #########################################################################
   Solves the surface energy balance for each cell in domain.  
   Tracks water stored as interception or dew and leaf temperature.
   ######################################################################### 
*/
 
	return 0;
}
