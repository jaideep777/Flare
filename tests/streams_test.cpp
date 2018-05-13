#include <iostream>
#include "../include/gsm.h"
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v3/include -L/home/jaideep/codes/libgsm_v3/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 streams_test.cpp -l:libgsm.so.3 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {-180, 360, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// set your data files
	string files[] = {
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2001.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2002.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2003.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2004.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2005.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2006.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2007.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2008.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2009.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2010.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2011.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2012.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2013.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2014.nc",
		"/media/jaideep/WorkData/Fire_G/fire_BA/burned_area.2015.nc"
	};

//	string files[] = {
//		"/media/jaideep/Totoro/Data/precip_trmm/pr.trmm.2003.nc"
//			};

	vector <string> filenames(files, files+15);

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(-89.75,89.75,0.5,nlats);
	vector <float> lons = createCoord(-180+0.25,180-0.25,0.5,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=ymd2gday("2003-6-1")-ymd2gday("2000-1-1");


	// test bilinear-interpolating stream
	gVar v("ba", "%", "days since 2000-1-1 6:0:0");
	v.setCoords(times, levs, lats, lons);
	v.createNcInputStream(filenames, glim);	 // default stream is bilinear-interpolating stream
	v.printGrid();
	v.readVar_gt(ymd2gday("2008-6-25"), 0);
	v.writeOneShot("ba20080625_0.5.nc");

	// test non-interpolating stream
	gVar w;
	w.initMetaFromFile(filenames[0]);
	w.createNcInputStream(filenames, glim, "none");	// specify no interpolation
	w.printGrid();
	w.readVar_gt(ymd2gday("2008-6-25"), 0);
	w.writeOneShot("ba20080625.nc");

	v.readVar_reduce_mean(ymd2gday("2001-1-1"), ymd2gday("2001-12-31"));

	w.readVar_reduce_mean(ymd2gday("2002-1-1"), ymd2gday("2002-12-31"));


//	// create output stream and write variable to output	
//	v.createNcOutputStream("testnc1.nc");
//	v.writeVar(0);	

//	// close streams when done
//	v.closeNcOutputStream();
//	v.closeNcInputStream();
//	v.writeOneShot("trmm20030101.nc");



//	// demonstration that operators copy NcStream data and duplicate pointers
//	gVar w("haha", "-", "days since 2000-1-1 6:0:0");
//	w.setCoords(times, levs, lats, lons);
	
//	gVar z = v+w;

//	v.printGrid();
//	w.printGrid();
//	z.printGrid();
//	
//	*gsm_log << sizeof(gVar) << "=======================================================" << sizeof(NcFile_handle) << endl;
//	
//	z.readVar_gt(ymd2gday("2003-06-01")+hms2xhrs("0:0:0"), 0);
//	v.readVar_gt(ymd2gday("2003-06-01")+hms2xhrs("0:0:0"), 0);
	

	// demonstration of readVar_reduce()
	
	


//	v.closeNcInputStream();



//	v.closeNcInputStream();
//	z.closeNcInputStream();
	
//	w.createOneShot(files[1]);
//	w.printGrid();
//	w.printValues();
	
//	v.readOneShot(files[1]);
//	v.printGrid();

//	v.createNcOutputStream("testnc1.nc");
//	v.writeVar(0);	
//	v.closeNcOutputStream();
	
	return 0;

}


