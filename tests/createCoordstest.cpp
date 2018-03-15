#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 createCoordstest.cpp -l:libgsm.so.2 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {0, 360, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(6.5,38.5,0.25,nlats);
	vector <float> lons = createCoord(11.45,18.55,0.1,nlats);
	vector <double> times(1); 
	times[0]=ymd2gday("2003-6-1")-ymd2gday("2000-1-1");

	cout << lats.size()  << ": "; printArray(lats);
	cout << lons.size()  << ": "; printArray(lons);

//	// indexC test
//	vector <float> gr(10);
//	for (int i=0; i<10; ++i) gr[i] = -89.75+0.5*i;
//	
//	vector <float> gr100(100);
//	for (int i=0; i<100; ++i) gr100[i] = -90+0.05*i;

//	printArray(gr);
//	for (int i=0; i<100; ++i){
//		cout << gr100[i] << " --> " << gr[indexC(gr, gr100[i])] << "\n";
//	}

//	cout << "----------------------------------------\n";


	return 0;

}


