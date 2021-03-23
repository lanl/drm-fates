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
	double tFA=0.0, aFA=0.0, tP=0.0, aP=0.0;
	double djnk=0.0, time=0.0;
	double tree=0.0, shrub=0.0, grass=0.0;
			
	int l=0, i=0;
	int lcount=0, lcount_SpRt=0;
	int sl=0, rel=0;
	int nx=200, ny=200, nz=1, hibin=0;

	char cjunk; 
	double djunk;
 //Set number of relizations to Test. = number of PARFLOW relizations	
	
	string trash, file, VegOut, IenF, OutF, OutFF, TIME, DISTANCE;

//Identify streamline files from slim to read in & Output files created
        ifstream locstream;
        IenF = "CellTreeMap.txt";
        locstream.open(IenF.c_str());
        if(!locstream)
        {cout << "Error opening '"<<IenF<< " input file"<<endl;
                return -1;  }
        lcount = linecounter (IenF);
        getline(locstream, trash); 

	ofstream outfile;
	OutF = "Xslope.txt";
        outfile.open(OutF.c_str());
        outfile<<nx<<"  "<<ny<<"  "<<nz<<endl;
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(3);
        outfile.width(10);

        ofstream youtfile;
        OutFF = "Yslope.txt";
        youtfile.open(OutFF.c_str());
        youtfile<<nx<<"  "<<ny<<"  "<<nz<<endl;
        youtfile.setf(ios::fixed, ios::floatfield);
        youtfile.precision(3);
        youtfile.width(10);

        ofstream vegfile;
        VegOut = "NEW_drv_vegm.dat";
        vegfile.open(VegOut.c_str());
        vegfile<<" x  y  lat    lon    sand clay color  fractional coverage of grid by vegetation class (Must/Should Add to 1.0)"<<endl;
        vegfile<<"       (Deg)     (Deg)  (%/100)   index  1    2    3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18"<<endl;
        vegfile.setf(ios::fixed, ios::floatfield);
        vegfile.precision(1);
        vegfile.width(2);


        for (int y=0; y<ny; y++){
          for (int x=0; x<nx; x++){
            outfile<<0.0<<endl;
            youtfile<<0.02<<endl;
            locstream>>tree;
            if (tree==0.0){shrub=0.0, grass=1.0;}
            if (tree==1.0){shrub=0.0, grass=0.0;}
            if (tree==2.0){shrub=1.0, grass=0.0, tree=0.0;}
            vegfile<<y+1<<"  "<<x+1<<"  35.837 -106.417  0.16 0.265   2  "<<tree<<"  0.0  0.0  0.0  0.0  "<<shrub<<"  0.0  0.0  0.0  "<<grass<<"  0.0  0.0  0.0  0.0  0.0  0.0  0.0  0.0"<<endl;
          }
        }
        locstream.close();
        vegfile.close();		
        youtfile.close();   	
       	outfile.close();
	return 0;
}
	

