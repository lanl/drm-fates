#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <stdio.h>
#include "allocate.h"

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
  double sin_X, sin_Y, tan_X, tan_Y, azimuth, XBottomCorner, YBottomCorner, XtopCorner, YtopCorner, DX, DY, DZ;
  double x_raySt, y_raySt, z_raySt, x_rayEnd, y_rayEnd, z_rayEnd, GndAngle, T_length;
  int Xdir, Ydir, Zdir;
  string corner;
};

void XYangle (Direction& ang_){
  double PI=3.14159265;
  if (ang_.azimuth<90) {
     ang_.sin_X = sin(ang_.azimuth*PI/180), ang_.sin_Y = sin((90-ang_.azimuth)*PI/180);
     ang_.tan_X = tan(ang_.azimuth*PI/180), ang_.tan_Y = tan((90-ang_.azimuth)*PI/180);
     ang_.corner = "NE", ang_.XtopCorner = -ang_.DX, ang_.YtopCorner = -ang_.DY;
     ang_.XBottomCorner = ang_.DX, ang_.YBottomCorner = ang_.DY;
  }
  if ((ang_.azimuth>=90)&&(ang_.azimuth<180)) {
     ang_.sin_X = sin((90-(ang_.azimuth-90))*PI/180), ang_.sin_Y = -sin((ang_.azimuth-90)*PI/180);
     ang_.tan_X = tan((90-(ang_.azimuth-90))*PI/180), ang_.tan_Y = -tan((ang_.azimuth-90)*PI/180);
     ang_.corner = "SE", ang_.XtopCorner = -ang_.DX, ang_.YtopCorner = ang_.DY;
     ang_.XBottomCorner = ang_.DX, ang_.YBottomCorner = 0;
  }
  if ((ang_.azimuth>=180)&&(ang_.azimuth<270)) {
     ang_.sin_X = -sin((ang_.azimuth-180)*PI/180), ang_.sin_Y = -sin((90-(ang_.azimuth-180))*PI/180);
     ang_.tan_X = -tan((ang_.azimuth-180)*PI/180), ang_.tan_Y = -tan((90-(ang_.azimuth-180))*PI/180);
     ang_.corner = "SW", ang_.XtopCorner = ang_.DX, ang_.YtopCorner = ang_.DY;
     ang_.XBottomCorner = 0, ang_.YBottomCorner = 0;
  }
  if (ang_.azimuth>=270) {
     ang_.sin_X = -sin((90-(ang_.azimuth-270))*PI/180), ang_.sin_Y = sin((ang_.azimuth-270)*PI/180);
     ang_.tan_X = -tan((90-(ang_.azimuth-270))*PI/180), ang_.tan_Y = tan((ang_.azimuth-270)*PI/180);
     ang_.corner = "NW", ang_.XtopCorner = ang_.DX, ang_.YtopCorner = -ang_.DY;
     ang_.XBottomCorner = 0, ang_.YBottomCorner = ang_.DY;
  }
}

void CELLexit (Direction& ang_){
  double X_ = 0.0,  Y_ = 0.0,  Z_ = 0.0;   
  double C_1 = 0.0, X1 = 0.0, Y1 = 0.0, Z1 = 0.0;
  ang_.Xdir = 0, ang_.Ydir = 0, ang_.Zdir = 0;
  Z_ = ang_.DZ - ang_.z_raySt;
  if (ang_.azimuth<90) { // Early Morning Sun
    X_ = ang_.DX - ang_.x_raySt, Y_ = ang_.DY - ang_.y_raySt;
    C_1 = X_  / ang_.sin_X;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    Y1 = C_1 * ang_.sin_Y;
    X1 = X_;  
    if ((abs(Y1)<Y_)&&(Z1<Z_)){  // Search X-Face Exit
       ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
       ang_.Xdir = 1;
       //std::cout<<"X-FACE (<90):  "<<C_1<<"   "<<X_<<"  "<<Y_<<"  "<<Z_; 
    }else{
      C_1 = Z_ / tan((ang_.GndAngle)*3.14159265/180);
      Y1 = C_1 * ang_.sin_Y;
      X1 = C_1 * ang_.sin_X;
      Z1 = Z_;
      if ((abs(Y1)<Y_)&&(X1<X_)){ // Search Z top face Exit
        ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
        ang_.Zdir = 1;
       // std::cout<<"Z-FACE (<90):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
      }else{
         C_1 = Y_  / ang_.sin_Y;
         Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
         X1 = C_1 * ang_.sin_X;
         Y1 = Y_;
         if ((abs(X1)<X_)&&(Z1<Z_)){ // Search Y-Face Exit
           ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
           ang_.Ydir = 1;
           //std::cout<<"Y-FACE (<90):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
         }}}}

  if ((ang_.azimuth>=90)&&(ang_.azimuth<180)) { // Late Morning Sun
    X_ = ang_.DX - ang_.x_raySt, Y_ = abs(ang_.y_raySt);
    if (Y_==0){Y_ = ang_.DY;}
    C_1 = X_  / ang_.sin_X;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    Y1 = C_1 * ang_.sin_Y;
    X1 = X_;
    if ((abs(Y1)<Y_)&&(Z1<Z_)) { // Search X-Face Exit
       ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
       ang_.Xdir = 1;
       //std::cout<<"X-FACE (>90):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
    }else{
       C_1 = Z_ / tan((ang_.GndAngle)*3.14159265/180);
       Y1 = C_1 * ang_.sin_Y;
       X1 = C_1 * ang_.sin_X; 
       Z1 = Z_;      
       if ((abs(Y1)<Y_)&&(X1<X_)){ // Search Z top face Exit
         ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
         ang_.Zdir = 1;
         //std::cout<<"Z-FACE (>90):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
       }else{
         C_1 = Y_  / - ang_.sin_Y;
         Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
         X1 = C_1 * ang_.sin_X;
         Y1 = -Y_;
         if ((abs(X1)<X_)&&(Z1<Z_)){ // Search Y-Face Exit
           ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
           ang_.Ydir = -1;
           //std::cout<<"Y-FACE (>90):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
         }}}}

  if ((ang_.azimuth>=180)&&(ang_.azimuth<270)) { // Early Afternoon Sun
    X_ = abs(ang_.x_raySt), Y_ = abs(ang_.y_raySt);
    if (Y_==0){Y_ = ang_.DY;}
    if (X_==0){X_ = ang_.DX;}
    C_1 = Y_  / - ang_.sin_Y;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    X1 = C_1 * ang_.sin_X;
    Y1 = - Y_;
    if ((abs(X1)<X_)&&(Z1<Z_)){ // Search Y-Face Exit
      ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
      ang_.Ydir = -1;
      //std::cout<<"Y-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
    }else{
      C_1 = Z_ / tan((ang_.GndAngle)*3.14159265/180);
      Y1 = C_1 * ang_.sin_Y;
      X1 = C_1 * ang_.sin_X;
      Z1 = Z_;
      if ((abs(Y1)<Y_)&&(abs(X1)<X_)){ // Search Z top face Exit
        ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
        ang_.Zdir = 1;
        //std::cout<<"Z-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
      }else{
    C_1 = X_  / - ang_.sin_X;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    Y1 = C_1 * ang_.sin_Y;
    X1 = - X_;
    if ((abs(Y1)<Y_)&&(Z1<Z_)) { // Search X-Face Exit
       ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
       ang_.Xdir = -1;
       //std::cout<<"X-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
    }}}}

  if (ang_.azimuth>=270) { // Late Afternoon Sun
    X_ = abs(ang_.x_raySt), Y_ = ang_.DY - ang_.y_raySt;
    if (X_==0){X_ = ang_.DX;}
    C_1 = X_  / - ang_.sin_X;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    Y1 = C_1 * ang_.sin_Y;
    X1 = - X_;
    if ((abs(Y1)<Y_)&&(Z1<Z_)) { // Search X-Face Exit
       ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
       ang_.Xdir = -1;
       //std::cout<<"X-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
    }else{
      C_1 = Z_ / tan((ang_.GndAngle)*3.14159265/180);
      Y1 = C_1 * ang_.sin_Y;
      X1 = C_1 * ang_.sin_X;
      Z1 = Z_;
      if ((abs(Y1)<Y_)&&(abs(X1)<X_)){ // Search Z top face Exit
        ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
        ang_.Zdir = 1;
        //std::cout<<"Z-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
      }else{
    C_1 = Y_  / ang_.sin_Y;
    Z1 = (C_1 * tan((ang_.GndAngle)*3.14159265/180));
    X1 = C_1 * ang_.sin_X;
    Y1 =  Y_;
    if ((abs(X1)<X_)&&(Z1<Z_)){ // Search Y-Face Exit
      ang_.x_rayEnd =X1, ang_.y_rayEnd = Y1, ang_.z_rayEnd = Z1, ang_.T_length = C_1;
      ang_.Ydir = 1;
      //std::cout<<"Y-FACE (>180):  "<<C_1<<"  "<<X_<<"  "<<Y_<<"  "<<Z_;
    }}}}
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
    double WaterMass = 0.0;
    double QswIn = 0.0;
    double QlwIn = 0.0;
    double Us = 0.0;
    double Dt = 60;
    double Pr = 0.0;
    double LAI = 0.0, ShadeF = 0.0, FuelDensity = 0.0;
    int fcellnum, NZ=41, NY=200, NX=200; 
    double TIME = 0.0;
    double djunk=0.0;    

    // Domain Size and Dimensions **** Need to figure out how to do changing DZ 
    double DZ = 2.0, DY = 2.0, DX = 2.0;  
    double rhomicro = 500.0, sizescale = 0.0005; 
 
string trash, AMout, AVS_out, Mois_out, LD, MD, SD, FD, TTS1, FS;

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


/* SECTION 1. USE FIRETEC DOMAIN TO DETERMIN ARRAY STRUCTUR
   #########################################################################
   Loops through FireTEC domain to find all cells with fuel.
   Sets array sizes and locations for Shade, LAI, and WaterStored.
   Each array element will have an x,y,z location
   #########################################################################
*/
//////MAKE MOCK DOMAIN of 4x3x3 cells. [z,x,y]
//ifstream FuelData;
//FD = "FuelData.txt";
//FuelData.open(FD.c_str());
//if (!FuelData) {
//    cout <<"ERROR opening " <<FD<<" input file"<<endl;
//} getline(FuelData, trash);
//FuelData>>NZ>>NY>>NX;
////// MOCK DOMAIN ABOVE
int ArrayCount = 0;
int count = 0;
const int NData = (NZ * NY * NX);
const int Header_Length = 4;
const int Trailer_Length = 4;
char header[Header_Length];
float FuelPointer[NData];
char trailer[Trailer_Length];

int EndTime = 288;
double fueltol = 0.000005;


int x_count=0, y_count=0;
//READ BINARY UNFORMATED FORTRAN FILE --> treesrhof.dat file
ifstream ifs("treesrhof.dat");             // open file
ifs.read(header, sizeof(header));      // read header
ifs.read((char*)FuelPointer, sizeof(FuelPointer));  // read the data
ifs.read(trailer, sizeof(trailer));   // read the trailer 

double * ZZarray;
double * DDZZarray;
double * DomainXloc;
double * DomainYloc;
double * DomainZloc;
double * moistcontent;
int * FuelMAP;
int * DomainMAP;
double ** FuelMoisture;
double ** LeafTemp;
ZZarray = new double [NZ]; 
DDZZarray = new double [NZ]; 
DomainXloc = new double [NX*NY*NZ];
DomainYloc = new double [NX*NY*NZ];
DomainZloc = new double [NX*NY*NZ];
moistcontent = new double [NX*NY*NZ];
FuelMAP = new int [NX*NY*NZ]; 
DomainMAP = new int [NX*NY*NZ];
FuelMoisture = allocate_2d_double(EndTime, NData);  
LeafTemp = allocate_2d_double(EndTime, NData);

for (int zz = 0; zz < NZ; zz++) {
  ZZarray[zz] = MeasureDZ(zz);
  DDZZarray[zz] = MeasureDZ(zz+1) - MeasureDZ(zz);
  for (int yy = 0; yy < NY; yy++) {
    for (int xx = 0; xx < NX; xx++){
      FuelMAP[count] = 0;
      DomainMAP[count] = -99;
      DomainXloc[count] = (count%NX)*DX;
      DomainYloc[count] = ((count/NX)%NY)*DY;
      DomainZloc[count] = ((count/(NX*NY))%NZ);
      DomainZloc[count] = MeasureDZ(DomainZloc[count]);
      //FuelData>>FuelDensity; // MOCK DOMAIN READING
      //if (FuelDensity>0){  // MOCK DOMAIN READING 
      if(FuelPointer[count] > fueltol){  // BIG DOMAIN
         //cout<<zz<<"  "<<yy<<"  "<<xx<<"  "<<FuelPointer[count]<<endl;
         FuelMAP[count] = 1; 
         ArrayCount++;
         y_count++;
      }count++;
    } 
  } 
cout<<zz<<"  Number of Cells per level: "<<y_count<<endl;
y_count=0;
}

//FuelData.close(); // MOCK DOMAIN

std::cout<<count<<"   Array Length: "<<ArrayCount<<std::endl;
int * Address;
double * LaiArray;
double * Xloc;
double * Yloc;
double * Zloc;
double * DZarray;
double * OverZ;
double * DZoverZcells;
double * CellsPerZlevel;
double ** QswFactor;
Address = new int [ArrayCount];
LaiArray = new double [ArrayCount];
Xloc = new double [ArrayCount];
Yloc = new double [ArrayCount];
Zloc = new double [ArrayCount];  //= {-99};
DZarray = new double [ArrayCount]; // = {0};
OverZ = new double [ArrayCount]; // Array size large than needed.
DZoverZcells = new double [ArrayCount]; // Array size large than needed.
CellsPerZlevel = new double [ArrayCount]; // Array size large than needed.
QswFactor = allocate_2d_double(EndTime, ArrayCount); // = {0};

int Zlevelcount=0, index=0;
cout <<"Here 4"<<endl;

////// MOCK DOMAIN *****
//FuelData.open(FD.c_str());
//if (!FuelData) {
//    cout <<"ERROR opening " <<FD<<" input file"<<endl;
//} getline(FuelData, trash);
//getline(FuelData, trash);
////// MOCK DOMAIN

int ArrayIndex = 0;
for (int dn=0; dn < (NZ*NY*NX); dn++) { // LOOP TO ASSING DATA, LOCATION & LAI OF EACH CELL WITH FUEL
    //FuelData>>FuelDensity; // ADD FUNCTION HERE *#! *#! *#! *#! *#! *#! *#! *#! *#! MOCK DOMAIN
    //if (FuelDensity>0){  // MOCK DOMAIN *********
    if (FuelPointer[dn] > fueltol){ // BIG DOMAIN
       Address[ArrayIndex] = dn;
       DomainMAP[dn] = ArrayIndex;
       Xloc[ArrayIndex] = (dn%NX)*DX;
       Yloc[ArrayIndex] = ((dn/NX)%NY)*DY;
       Zloc[ArrayIndex] = ((dn/(NX*NY))%NZ); // This is not the streched Z location Needs to measure DZ.
       Zloc[ArrayIndex] = MeasureDZ(Zloc[ArrayIndex]); // Bottom of cell.
       //UnderZ[ArrayIndex] = ((dn/(NX*NY))%NZ);
       // Call Function to determine DZ
       DZarray[ArrayIndex] = MeasureDZ(((dn/(NX*NY))%NZ)+1) - Zloc[ArrayIndex];  
       // Determine LAI  
       //LaiArray[ArrayIndex] = FuelDensity*DZarray[ArrayIndex] / (rhomicro * sizescale); // MOCK DOMAIN READING   
       LaiArray[ArrayIndex] = FuelPointer[dn] * DZarray[ArrayIndex] / (rhomicro * sizescale); // Consider Reading Sizescale form file. // BIG DOMAIN 
       if ( DZarray[ArrayIndex]>DZarray[ArrayIndex-1]){
          OverZ[Zlevelcount] = Zloc[ArrayIndex], DZoverZcells[Zlevelcount] = DZarray[ArrayIndex];
          std::cout<<"Z-level counts  "<<Zlevelcount<<"  "<<Zloc[ArrayIndex]<<"  "<<DZarray[ArrayIndex]<<"  "<<ArrayIndex<<"  "<<CellsPerZlevel[Zlevelcount-1]<<"  "<<OverZ[Zlevelcount]<<std::endl;
          //cout<<"LAI: "<<LaiArray[ArrayIndex]<<"   FuelDensity: "<<FuelPointer[dn]<<"    DZ: "<<DZarray[ArrayIndex]<<endl;
          Zlevelcount++; 
       }
       ArrayIndex++;
    }
}
//FuelData.close(); // MOCK DOMAIN
// Count the number of cells in each level
double cell_data = DZarray[0];
Zlevelcount=0;
int counter=0;
for (int i=0; i < ArrayCount; i++) {
    if (cell_data<DZarray[i]){
        cell_data = DZarray[i];
        CellsPerZlevel[Zlevelcount] = counter;
        Zlevelcount++;
        counter=0;
    }
    counter++;
}
CellsPerZlevel[Zlevelcount] = counter;
Zlevelcount++;

for (int zl = 0; zl < Zlevelcount; zl++) {
    cout<<zl<<"   Cells Per Level: "<<CellsPerZlevel[zl]<<endl;
}

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
    spa.month         = 3;
    spa.day           = 2;
    spa.hour          = 0;
    spa.minute        = 0;
double min=0.0, sec=0.0, daylight=0.0, zenith=0.0, azimuth=0.0, sunrise=0.0, sunset=0.0;
double Z_ShaddowLength=0.0, Y_Shaddow=0.0, X_Shaddow=0.0;
double Z_topShaddow=0.0, Y_topShaddow=0.0, X_topShaddow=0.0;
double Xpp=0.0, Ypp=0.0, Zpp=0.0, Zsp=0.0, Ysp=0.0, Xsp=0.0;
double ZspTop=0.0, YspTop=0.0, XspTop=0.0; 
double X1=0.0, X2=0.0, Y1=0.0, Y2=0.0, Z1=0.0, Z2=0.0;
double GroundAngle=0.0, C_1 = 0.0, C_2=0.0, T_length=0.0;
double TransectLAI = 0.0;
double DeepShadeLAI = 12.0;  //An LAI of 10 means that only ambiant daylight gets in 0.5% of the SWradation/.
int sunresult;
int tt = 0, f = 0, z = 0, in = 0, ID=0, sumloop = 0, levelsum = 0, StartLevelSearch = 0, EndLevelSearch = 0;
ang_.DX = DX, ang_.DY = DY;
ang_.x_raySt = 0.0, ang_.y_raySt = 0.0, ang_.z_raySt = 0.0;
ang_.x_rayEnd = 0.0, ang_.y_rayEnd = 0.0, ang_.z_rayEnd = 0.0; 

for (tt=0; tt < EndTime; tt++) {
    // FINDING POSITION OF SUN BASED ON TIME
    spa.minute = (tt%6)*10;
    spa.hour = (tt/6)%24;
    spa.day = 2 + (tt/144)%144;

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
     }else printf("SPA Error Code: %d\n", sunresult);
     //std::cout<<tt<<"  "<<spa.day<<"  "<<spa.hour<<"  "<<spa.minute<<"  Azimuth: "<<azimuth<<"  Zenith: "<<zenith<<"  "<<ang_.sin_X<<"  "<<ang_.sin_Y<<"  ***********************  "<<Zlevelcount<<std::endl;
     if (zenith<85) { // Sun is above the horizon Shade will happen
      XYangle(ang_); // Function to find the sin for X and Y based on angle of the sun. --> Also finds corner (side sun is shining) 
      GroundAngle = 90 - zenith;
      ang_.GndAngle = GroundAngle;
      std::cout<<tt<<"  "<<spa.day<<"  "<<spa.hour<<"  "<<spa.minute<<"  Azimuth: "<<azimuth<<"  Zenith: "<<zenith<<"  "<<ang_.sin_X<<"  "<<ang_.sin_Y<<"  ************************************  "<<Zlevelcount<<std::endl;

      // LOOP THROUGH FUEL HERE & CAST SHADDOW
      for (int l=0; l < Zlevelcount; l++) { 
      for (f=0; f < CellsPerZlevel[l]; f++){// Fuel Loop  --> Calculate Shade by shooting lines from ground to sun.  Intercepted fuel causes shade.
          index = f+sumloop; // Index sets the cell being evaluated for shade.
        z = l; // Should limit the number of internal loops so that we are only evaluating above the referanced cell.
        if (Zloc[index] < OverZ[Zlevelcount-1]) {  // The highest level of fuel never gets shade

// Find point that ray exits cell
          ang_.x_raySt = 0.75*DX, ang_.y_raySt = 0.75*DY, ang_.z_raySt = 0, ang_.DZ = DZarray[index];
          in = 0;
          // MAP TO LARGE DOMAIN
          ID = Address[index];
          //ID = DomainMAP[ID];
          //cout<<index<<"    "; 
// Make   Function  
          CELLexit(ang_); 
          Xpp = ang_.x_raySt+ang_.x_rayEnd, Ypp = ang_.y_raySt+ang_.y_rayEnd, Zpp = ang_.z_rayEnd;
          Xsp = Xpp + Xloc[index], Ysp = Ypp + Yloc[index] , Zsp = Zpp + Zloc[index];          
          //cout<<"Cell ID: "<<ID<<"    "<<Xsp<<"  "<<Ysp<<"  "<<Zsp<<"     "<<ang_.Xdir<<"  "<<ang_.Ydir<<"  "<<ang_.Zdir;
//New Cell
          if ((ang_.Ydir==0)&&(ang_.Zdir==0)) { // X - FACE
             ang_.x_raySt = 0;
             ang_.y_raySt = Ypp, ang_.z_raySt = Zpp;
             ID = ID + ang_.Xdir;
          }else if ((ang_.Xdir==0)&&(ang_.Ydir==0)){ // Z - FACE
            ang_.z_raySt = 0;  
            ang_.y_raySt = Ypp, ang_.x_raySt = Xpp;
            z++;
            ang_.DZ = DDZZarray[z];
            ID = ID + (NX*NY); // Going Up a level  ************************************************
          } else if ((ang_.Xdir==0)&&(ang_.Zdir==0)) { // Y - FACE
             ang_.y_raySt = 0;
             ang_.x_raySt = Xpp, ang_.z_raySt = Zpp;
             ID = ID + (ang_.Ydir * NX);
          }
          if ((ID<0)||(ID >= (NX*NY*NZ))) {
           break;}
          //cout<<"   Cell ID: "<<ID<<"            Transect Length: "<<ang_.T_length<<endl;
          // Cell Can't shade itself !!!!!!!!!!
  
          if ((Xsp <= 0)||(Ypp <= 0)||(Xsp >= (NX*DY))||(Ysp >= (NY*DY))) {in = 2;}
          while ((in < 1)&&(z<NZ)) {
              //cout<<index<<"        "; 
              CELLexit(ang_);
              Xpp = ang_.x_raySt + ang_.x_rayEnd, Ypp = ang_.y_raySt + ang_.y_rayEnd, Zpp = ang_.z_raySt + ang_.z_rayEnd;
              Xsp = Xpp + DomainXloc[ID], Ysp = Ypp + DomainYloc[ID] , Zsp = Zpp + DomainZloc[ID];  // Index here is wrong
              //cout<<index<<"   Cell ID: "<<ID<<"    "<<Xsp<<"  "<<Ysp<<"  "<<Zsp<<"     "<<ang_.Xdir<<"  "<<ang_.Ydir<<"  "<<ang_.Zdir;
              if (FuelMAP[ID]>0){ // FUEL EXISTS **** CAST SHADE 
                // Normalize transect lenght by DZ to wieght LAI of shadding cell
                T_length = ang_.T_length;
                if (T_length > DZarray[DomainMAP[ID]]){T_length = DZarray[DomainMAP[ID]];}  // Change to array address
                TransectLAI = T_length/DZarray[DomainMAP[ID]];
                TransectLAI = TransectLAI * LaiArray[DomainMAP[ID]];
                  // Sum Shade Factor in Affected Cell
                QswFactor[tt][index] = QswFactor[tt][index] + (TransectLAI/DeepShadeLAI);
//if ((tt>75)&&(tt<79)) {
//   if ((index==20050)||(index==20055)||(index==58358)||(index==58362)) {
// cout<<index<<"   SHADE FACTOR: "<<QswFactor[tt][index]<<"  TransectLAI: "<<TransectLAI<<"  LaiArray: "<<LaiArray[DomainMAP[ID]]<<"  DZarray: "<<DZarray[DomainMAP[ID]]<<"  T_length: "<<T_length<<"      "<<ang_.Xdir<<"  "<<ang_.Ydir<<"  "<<ang_.Zdir<<endl;
//}}
              }
              if ((ang_.Ydir==0)&&(ang_.Zdir==0)) { // X - FACE
                 ang_.x_raySt = 0;
                 ang_.y_raySt = Ypp, ang_.z_raySt = Zpp;
                 ID = ID + ang_.Xdir;
              }else if ((ang_.Xdir==0)&&(ang_.Ydir==0)){ // Z - FACE
                ang_.z_raySt = 0;  
                ang_.y_raySt = Ypp, ang_.x_raySt = Xpp;
                z++;
                ang_.DZ = DDZZarray[z];
                ID = ID + (NX*NY); // Going Up a level ************************************************
              } else if ((ang_.Xdir==0)&&(ang_.Zdir==0)) { // Y - FACE
                 ang_.y_raySt = 0;
                 ang_.x_raySt = Xpp, ang_.z_raySt = Zpp;
                 ID = ID + (ang_.Ydir * NX);
              }
              if ((ID<0)||(ID >= (NX*NY*NZ))) { 
           break;}
              //cout<<"   Cell ID: "<<ID<<"            Transect Length: "<<ang_.T_length<<endl;
              if ((Xsp <= 0)||(Ysp <= 0)||(Xsp >= (NX*DY))||(Ysp >= (NY*DY))) {in = 2;}
          }
          if (QswFactor[tt][index]>0.95){QswFactor[tt][index]=0.95;} // 0.05% of SW will always get through for ambiant light.
        } // If top level statement
        levelsum = 0;
      }sumloop = sumloop + CellsPerZlevel[l]; } // End Fuel Loop
     sumloop = 0;
     } 
}// End Time Loop
cout<<endl<<endl;

//for (int ii=0; ii<ArrayCount; ii++) {
//cout<<ii<<"  Address[ii]: "<<Address[ii]<<"    "<<DomainMAP[Address[ii]]<<"   FuelMAP[Address[ii]]: "<<DomainMAP[Address[ii]]<<endl;
//}

cout<<"tt: "<<tt<<"  index: "<<index<<"    "<<ArrayIndex<<endl<<endl;

/* SECTION 3.  SOLVE SURFACE ENERGY BALANCE 
   #########################################################################
   Solves the surface energy balance for each cell in domain.  
   Tracks water stored as interception or dew and leaf temperature.
   ######################################################################### 
*/

// OUTPUT FILES
//Creating output file
ofstream Amanzistream;
AMout = "Amanzi.ALL.out.txt";
Amanzistream.open(AMout.c_str());
Amanzistream<<"TIME   Ts     WaterStorage"<<endl;
double CellVol = 0.0;
for (int f; f < ArrayIndex; f++) {
  stringstream fcellconverter;
  fcellconverter << f;
  string strFcell = fcellconverter.str(); // Creates naming structure for relization files used to open input files.

// METEROLOGICAL INPUT FILE
  ifstream MetData;
  MD = "EglinMetData.txt";
  MetData.open(MD.c_str());
  if (!MetData){
      cout <<"ERROR opening "<<MD<<" input file"<<endl;
  } getline (MetData,trash);

  // Set Parameters For TIME LOOP iterations 
  //Initial Conditions = Initial WaterStored[m] ~ Assume Fully Saturated
  //Fuel Data
  st_data.st_energy.lai = LaiArray[f];
  //Initial Conditions = Initial WaterStored[m] ~ Assume Fully Saturated
  st_data.st_energy.WaterStorage = st_data.st_energy.lai * 0.00025;

  Amanzistream<<"Fuel Cell: "<<f<<"    NEW CELL  *#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#*#"<<endl;

  for (int rel=1; rel<=EndTime; rel++) { //TIME LOOP LOOP LOOP #####################
      MetData  >>djunk>>QswIn>>QlwIn>>Ta>>RH>>Us>>Pr;
      // Read Shade Factor here
      // Shade Factor Calculation ********************************** !!!!!!!!!!!!!!!!!!
      QswIn = QswIn*(1-QswFactor[rel-1][f]);
      //if (QswFactor[rel-1][f]>0) {
         //cout<<f<<"  "<<rel<<"  "<<QswFactor[rel-1][f]<<"  "<<QswIn<<endl; 
      //} 
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

      TIME = TIME + (Dt/60);  // Time in Days

      Amanzistream<<setw(4)<<rel<<"  "<<TIME<<"  "<<setw(3)<<st_data.st_energy.Ts-273.15<<"   "<<st_data.st_energy.WaterStorage;

      // Convert Depth of water[m] to mass of water[kg] in cell 
      WaterMass = st_data.st_energy.WaterStorage * (DX * DY) * 1000; 
      // Convert mass of water to water content (mass of water / Dry mass of fuel)
      CellVol = DX * DY * DZarray[f];
       
      FuelMoisture[rel-1][Address[f]] = WaterMass / (FuelPointer[Address[f]] * CellVol);
      LeafTemp[rel-1][Address[f]] = st_data.st_energy.Ts-273.15;

      Amanzistream<<"     "<<WaterMass<<"  "<<CellVol<<"  "<<FuelPointer[Address[f]]<<"  "<<FuelMoisture[rel-1][Address[f]]<<endl;

      // Resetting variables for next iteration
      st_data.vp_air.actual_vaporpressure=0.0;
      st_data.st_energy.fQswIn=0.0, st_data.st_energy.fQlwIn=0.0, st_data.st_energy.fQlwOut=0.0, st_data.st_energy.fQe=0.0, st_data.st_energy.fQh=0.0;
      st_data.st_energy.fQe=0.0;
      st_data.st_energy.varvar=0.0;

  }//  END TIME LOOP LOOP LOOP LOOP #####################

  TIME = 0;
}  // End Fuel Cell Loop

Amanzistream.close();

/* SECTION 4.  MAKE VISULAIZATION AND FILES FOR FIRETEC
   #########################################################################
   Makes AVS files and binary fuel moisture files  
   Output variables are water stored in canopy and leaf temperature.
   ######################################################################### 
*/

double z_position = 0.0;
int nodeZ=41;
int nodeNum = (nodeZ + 1) * (NY + 1) * (NX + 1);

//Create AVS file
ofstream AVSstream;
ofstream Moiststream;

AVS_out = "AVS_Domain.MAIN.inp";
AVSstream.open(AVS_out.c_str());
AVSstream<<"       "<<nodeNum<<"     "<<(nodeZ * NY * NX)<<"     0      2      0"<<endl;
count=0;
for (int z = 0; z < (nodeZ+1); z++) {
  for (int y = 0; y < (NY + 1); y++) {
    for (int x = 0; x < (NX + 1); x++){
     z_position = (z%(nodeZ+1)); // This is not the streched Z location Needs to measure DZ.
     z_position = MeasureDZ(z_position); // Bottom of cell.
     count++;
     AVSstream<<(count)<<"  "<<x<<"  "<<y<<"  "<<z_position<<endl;
     //cout<<(count)<<"  "<<x<<"  "<<y<<"  "<<z_position<<endl;
}}}

count=0;
int node1=0, node2=0, node3=0, node4=0,  node5=0, node6=0, node7=0, node8=0, node_=0;
for (int z = 0; z < (nodeZ); z++) {
  for (int y = 0; y < NY; y++) {
    for (int x = 0; x < NX; x++){
    count++;
    node_ = (((count-1)%NX) + 1) + ((y%NY) * NX) + (y%NY);
    node1 = node_ + (((z%(nodeZ+1))) * ((NX + 1) * (NY + 1)));
    node2 = node1 + 1;
    node4 = node1 + (NX+1); //node_ + (z%(nodeZ+1)) * ((NY+1)*(NX+1));
    node3 = node4 + 1;  // --> Needs to be be x=0, y=1, z=0
    node5 = node_ + ((NX + 1) * (NY + 1)) + (((z%(nodeZ+1))) * ((NX + 1) * (NY + 1)));
    node6 = node5 + 1;
    node8 = node5 + (NX+1);
    node7 = node8 + 1;   // --> Needs to be be x=0, y=1, z=1
    AVSstream<<(count)<<"    1  hex  "<<node5<<"  "<<node6<<"  "<<node7<<"  "<<node8<<"  "<<node1<<"  "<<node2<<"  "<<node3<<"  "<<node4<<endl;
    //cout<<(count)<<"    1  hex  "<<node1<<"  "<<node2<<"  "<<node3<<"  "<<node4<<"  "<<node5<<"  "<<node6<<"  "<<node7<<"  "<<node8<<endl;
}}}
AVSstream<<"002  1  1"<<endl;
AVSstream<<"WaterStored, real"<<endl;
AVSstream<<"LeafTemp, real"<<endl;

AVSstream.close();

ifstream Orig;
ofstream NewFile ;

for (tt=0; tt < EndTime; tt++) { // Time Step Loop
  stringstream ttconverter;
  ttconverter << tt;
  string strTT = ttconverter.str(); // Creates naming structure for relization files used to open input files. 

  //if (tt%3 == 0) {
  AVS_out = "AVS_Domain." + strTT + ".inp";
  Mois_out = "MoistureContent." + strTT + ".txt";
  cout<<tt<<"  "<<strTT<<"  "<<AVS_out<<endl;
  Moiststream.open(Mois_out.c_str());  

  //Open the file
  Orig.open("AVS_Domain.MAIN.inp");
  NewFile.open(AVS_out.c_str());
  NewFile<<Orig.rdbuf();
  //NewFile.close();
  Orig.close();
   
  //AVSstream.open(AVS_out.c_str());

  z_position = 0.0;
  count=0;
  for (int z = 0; z < (nodeZ); z++) {
    for (int y = 0; y < NY; y++) {
      for (int x = 0; x < NX; x++){
      //cout<<(count+1)<<"  "<<tt<<"  ";
      //if (tt>43){
      //cout<<FuelMoisture[tt][count]<<"  "<<LeafTemp[tt][count]<<endl;
      //}
      NewFile<<(count+1)<<"  "<<FuelMoisture[tt][count]<<"  "<<LeafTemp[tt][count]<<endl;
      Moiststream<<FuelMoisture[tt][count]<<endl;
      count++;
  }}}
  NewFile.close();
  Moiststream.close(); 
  //}
  std::cout<<"Array Length: "<<ArrayCount<<std::endl;
//ofstream outfile;
//outfile.open("moisture.dat", ios::binary | ios::out);
//outfile.write(header, sizeof(header));
//outfile.write((char*) FuelMoisture[1], sizeof(FuelMoisture[1]));
//outfile.write(trailer, sizeof(trailer));
//outfile.close();

}


delete [] ZZarray; 
delete [] DDZZarray;
delete [] DomainXloc;
delete [] DomainYloc;
delete [] DomainZloc;
delete [] FuelMAP;
delete [] DomainMAP;
delete [] Address;
delete [] LaiArray;
delete [] Xloc;
delete [] Yloc;
delete [] Zloc;
delete [] DZarray;
delete [] OverZ;
delete [] DZoverZcells;
delete [] CellsPerZlevel;
delete [] moistcontent;
deallocate_2d_double(QswFactor, EndTime, ArrayCount);
deallocate_2d_double(FuelMoisture, EndTime, NData);
deallocate_2d_double(LeafTemp, EndTime, NData);
 
	return 0;
}
