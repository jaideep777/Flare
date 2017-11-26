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

	// define coordinates
	int nlons, nlats, nlevs, ntimes;
	vector <float> lons = createCoord(66.5, 100.5, 0.5, nlons);
	vector <float> lats = createCoord(6.5, 38.5, 0.5, nlats);
	vector <float> levs = createCoord(1, 1, 1, nlevs);
	vector <float> times = createCoord( ymd2gday("2000-1-1") - ymd2gday("1950-1-1"),  ymd2gday("2000-1-1") - ymd2gday("1950-1-1"), 1, ntimes );

	cout << "lons = "; 
	printArray(lons);

	cout << "times = "; 
	printArray(times);

	// create georef variable
	gVar T("ts", "deg C", "days since 1950-1-1 0:0:0");
	T.printGrid();
	
	vector <double> t(times.size()); for (int i=0; i< t.size(); ++i) t[i] = times[i];
	T.setCoords(t, levs, lats, lons);
	T.printGrid();

	// read data from NC files
	vector <string> files;
	files.push_back("/media/jaideep/WorkData/Fire_G/precip_imd/rf_imd.2000.nc");
	files.push_back("/media/jaideep/WorkData/Fire_G/precip_imd/rf_imd.2003.nc");
	
	T.createNcInputStream(files, glim);
	T.printGrid();

	// read data
	T.readVar_gt( ymd2gday("2000-1-1") ,0);
	T.readVar_gt( ymd2gday("2003-6-23") ,0);
	
	gVar U;
	U.copyMeta(T);
	
	gVar W = U+T;
	
	return 0;
}




