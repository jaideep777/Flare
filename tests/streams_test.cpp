#include <iostream>
#include "../include/gsm.h"
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 inputvar_test.cpp -lgsm -lnetcdf_c++ 

int main(){
	
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);

	ofstream gsml("gsm_log1.txt");
	gsm_log = &gsml;

	float glimits[] = {0, 150, -60, 60};
	vector <float> glim(glimits, glimits+4);

	string files[] = {
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2000.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2001.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2002.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2003.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2004.nc",
		"/media/jaideep/WorkData/Fire_G/ncep_20cen/temp_sfc/air.sfc.2005.nc"
		};
			
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(-90,90,0.5,nlats);
	vector <float> lons = createCoord(0,360,0.5,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=ymd2gday("2003-6-1")-ymd2gday("2000-1-1");
//	times[1]=ymd2gday("2005-11-2")-ymd2gday("2000-1-1");

	gVar v("ts", "deg C", "days since 2000-1-1 6:0:0");
	v.setCoords(times, levs, lats, lons);
	v.values.resize(nlons*nlats*nlevs);
//	v.printGrid();
	
	vector <string> filenames(files, files+6);

	v.createNcInputStream(filenames, glim);	
	v.readVar_gt(ymd2gday("2003-06-01")+hms2xhrs("0:0:0"),1);
	v.closeNcInputStream();

	v.createNcOutputStream("testnc1.nc");
	v.writeVar(0);	
	v.closeNcOutputStream();
	
		
	return 0;

}


