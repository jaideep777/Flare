#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v2/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o exec lesson.cpp -l:libgsm.so.2 -lnetcdf_c++ 


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


	// read data from NC files
	vector <string> files;
	files.push_back("/media/jaideep/WorkData/Fire_G/precip_imd/rf_imd.2000.nc");
	files.push_back("/media/jaideep/WorkData/Fire_G/precip_imd/rf_imd.2003.nc");

	gVar pr;
	pr.initMetaFromFile("/media/jaideep/WorkData/Fire_G/precip_imd/rf_imd.2000.nc");
	pr.createNcInputStream(files, glim);
	pr.readVar_reduce_sd( ymd2gday("2000-1-1"), ymd2gday("2000-2-1") );	

	
	vector <double> times(1, ymd2gday("2000-1-15"));
	pr.setCoords(times, pr.levs, pr.lats, pr.lons);
	pr.printGrid();
	
	pr.createNcOutputStream("sample.nc");
	pr.writeVar(0);	
	
	
	pr.closeNcOutputStream();
	pr.closeNcInputStream();
	
	return 0;
}



