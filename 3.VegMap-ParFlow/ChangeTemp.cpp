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
	double rads=0.0, radl=0.0, rain=0.0, wind_u=0.0, wind_v=0.0, temp=0.0, press=0.0, rh=0.0;
	double djnk=0.0, time=0.0;
			
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
        IenF = "Narr-Parflow-July-Aug-2014.txt";
        locstream.open(IenF.c_str());
        if(!locstream)
        {cout << "Error opening '"<<IenF<< " input file"<<endl;
                return -1;  }
        lcount = linecounter (IenF);

	ofstream outfile;
	OutF = "Hot-Narr.txt";
        outfile.open(OutF.c_str());
        outfile.setf(ios::fixed, ios::floatfield);
        outfile.precision(4);
        outfile.width(10);



        for (int y=0; y<(lcount-1); y++){
          locstream>>rads>>radl>>rain>>temp>>wind_u>>wind_v>>press>>rh;
          temp=temp+4;
          outfile<<rads<<"  "<<radl<<"  "<<rain<<"  "<<temp<<"  "<<wind_u<<"  "<<wind_v<<"  "<<press<<"  "<<rh<<endl;
        }
        locstream.close();
       	outfile.close();
	return 0;
}
	

