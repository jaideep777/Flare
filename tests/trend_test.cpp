#include <iostream>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
#include "../include/gsm.h"
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v3/include -L/home/jaideep/codes/libgsm_v3/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 trend_test.cpp -l:libgsm.so.3 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
//	ofstream gsml("gsm_log.txt");
//	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {0, 360, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
//	int nlons, nlats, nlevs, ntimes;
//	vector <float> lons = createCoord(-180+.5/2,180-0.5/2,0.5,nlons);
//	vector <float> lats = createCoord(-90+.5/2,90-.5/2,0.5,nlats);
//	vector <float> levs = createCoord(1,1,1,nlevs);
//	vector <double> times(16*24); 
//	for (int i=0; i<times.size(); ++i) times[i]= ymd2gday("2001-1-1")+ i*365.2524/24 + 6 - ymd2gday("2000-1-1");

	string files[] = 
	{
		"/media/jaideep/Totoro/Data/GPP_modis/gpp.2000-2015.nc",
	};

	vector <string> infiles(files, files+1); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.createNcInputStream(infiles, glim);
	hires.printGrid();

	gVar slope = hires.trend(ymd2gday("2000-1-1"), ymd2gday("2015-12-31"));
	
	slope.writeOneShot("/tmp/npp.slope.nc");
	
	return 0;

}


