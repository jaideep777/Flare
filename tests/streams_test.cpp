#include <iostream>
#include "../include/gsm.h"
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 inputvar_test.cpp -lgsm -lnetcdf_c++ 

int main(){
	
	// ~~~~~~ Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// speficy log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {0, 150, -60, 60};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	string files[] = {
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2000.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2001.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2002.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2003.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2004.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2005.nc"
		};
			
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(-89.75,89.75,0.5,nlats);
	vector <float> lons = createCoord(0.25,359.75,0.5,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=ymd2gday("2003-6-1")-ymd2gday("2000-1-1");
//	times[1]=ymd2gday("2005-11-2")-ymd2gday("2000-1-1");

	gVar v("ts", "deg C", "days since 2000-1-1 6:0:0");
	v.setCoords(times, levs, lats, lons);
//	v.printGrid();
	
	vector <string> filenames(files, files+6);

//	v.createNcInputStream(filenames, glim);	
//	v.readVar_gt(ymd2gday("2003-06-01")+hms2xhrs("0:0:0"),1);
//	v.closeNcInputStream();

	
	
//	gVar w;
//	w.createOneShot(files[1]);
//	w.printGrid();
//	w.printValues();
	
	v.readOneShot(files[1]);
	v.printGrid();

	v.createNcOutputStream("testnc1.nc");
	v.writeVar(0);	
	v.closeNcOutputStream();
	
	return 0;

}


