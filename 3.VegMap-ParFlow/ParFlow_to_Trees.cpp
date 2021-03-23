#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <sstream>
/* #############################################################################################################################################
 T#E S#########E.
 BECOMES THE CELL LENGTH IN THE CRUNCFLOW 1D GEOCHEMICAL SIMULATION).*******
 
 NO OUTPUT IS CREATED REGARDING STREAMLINES THAT DO NOT PASS THROUGH THE SOURCE ZONE.
 
 
 HARDWARED INFORMATION NEEDED IN CODE:
                             VAlUSE PASSED IN:
                                         ####################################################
                             HARDWIRED IN
 
                            FILES NEEDED:
 
                            OUTPUT FILES CREATED:

                            OUTPUT DIRECTORY:

 ********************************************************************
  */
using namespace std;

int linecounter (string fname) {
	// This function simply counts the number of lines in an input file.
	
	ifstream input;
	
	int n;
	string line;
	
	input.open(fname.c_str());
	if (!input) {
		cout << "ERROR: Unable to open " << fname << "." << endl;
		exit (1);
	}
	n = 0;
	while (getline(input,line)) {
		n++;
	}
	
	input.close();
	
       return n;
}

int main(int argc, char* argv[])
{
	double SLnumber = 0.0;
	double FI=0.0, SR=0.0;
	double djnk=0.0, Time=0.0;
	double Xlow=0.0, Xhigh=0.0;
        double Ylow=0.0, Yhigh=0.0;
			
	int l=0, i=0;
	int lcount=0;
	int sl=0, cell=0;
	int x=0, y=0, z=0;
        int nx = 200, ny = 200;
        int dx = 2, dy = 2;

	char cjunk; 
	double djunk;
	
	string trash, file, InIenF, InSpRateF, OutIenF, CHECK, IenF, Xloc, Yloc;
	
//Allocate 2D data arrays
        double *Xarray;
        double *Yarray;
        double *spparray;

//Identify streamline files from slim to read in & Output files created
	ofstream Intensefile;
	OutIenF = "CellTreeMap.txt";
        Intensefile.open(OutIenF.c_str());
        Intensefile.setf(ios::fixed, ios::floatfield);
        Intensefile.precision(3);
        Intensefile.width(10);
        Intensefile<<nx<<"  "<<ny<<endl;

        ofstream checkfile;
        CHECK = "CellTreeMapCHECK.txt";
        checkfile.open(CHECK.c_str());
        checkfile.setf(ios::fixed, ios::floatfield);
        checkfile.precision(3);
        checkfile.width(10);
        checkfile<<nx<<"  "<<ny<<endl;        
	
  // Tree locations input
        ifstream locstream; 
        IenF = "CellTreeMap.txt";
        lcount = linecounter (IenF);
       
        Xarray = new double[lcount];
        Yarray = new double[lcount];
        spparray = new double[lcount];
    
  //Opening Files
        locstream.open(IenF.c_str());
        if(!locstream)
        {cout << "Error opening '"<<IenF<<" input file"<<endl;
                return -1;  }
        cout<<"Line Number of trees "<<lcount<<endl;
 
        for (int sl=0; sl < lcount; sl++) {// Loop through each tree location
           // getline(locstream,Xloc,',');
           // getline(locstream,Yloc,'\n');
           // Xarray[sl] = atof(Xloc.c_str());
           // Yarray[sl] = atof(Yloc.c_str());
            locstream>>spparray[sl]>>Xarray[sl]>>Yarray[sl]>>djnk>>djnk; 
            //cout<<Xloc<<endl;
            
        }//End tree loop

        cout<<"Between Loops"<<endl;

       for (int y=0; y < ny; y++) {//Domain cell loop
           Ylow = y*dy, Yhigh = (y*dy)+dy;
           cout<<"Ylow, Yhigh: "<<Ylow<<"  "<<Yhigh<<endl;
           for (int x=0; x < nx; x++) {
               Xlow = x*dx, Xhigh = (x*dx)+dx;
               cell = 0;
               for (int sl=0; sl < lcount; sl++) {// Checing Tree location
                   if ((Xarray[sl]>Xlow)&&(Xarray[sl]<Xhigh)&&(Yarray[sl]>Ylow)&&(Yarray[sl]<Yhigh)) {
                       cell = spparray[sl];
                       cout<<"X and Y fits: "<<Ylow<<"  "<<Yarray[sl]<<"  "<<Xarray[sl]<<"  "<<Yhigh<<"     "<<cell<<endl;
                   }
               } 
               checkfile<<x<<"  "<<y<<"  "<<cell<<endl;
               Intensefile<<cell<<endl; 
           }
       }
       cout<<"After Domain Loop"<<endl;
       Intensefile.close();
       locstream.close();
       checkfile.close();                
       cout<<"About to delete arrays"<<endl;
       delete[] Xarray;
       delete[] Yarray;	
       delete[] spparray; 
          	
	return 0;
}
	
	
	
