#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test4.cpp -l:libgsm.so.2 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {-180, 180, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lons = createCoord(-180+.5/2,180-0.5/2,0.5,nlons);
	vector <float> lats = createCoord(-90+.5/2,90-.5/2,0.5,nlats);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(16*24); 
	for (int i=0; i<times.size(); ++i) times[i]= ymd2gday("2001-1-1")+ i*365.2524/24 + 6 - ymd2gday("2000-1-1");

	string fireFiles[] = 
	{
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2001.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2002.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2003.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2004.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2005.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2006.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2007.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2008.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2009.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2010.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2011.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2012.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2013.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2014.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2015.nc",
		"/media/jaideep/Totoro/Data/Fire_BA/burned_area.2016.nc"
	};

	vector <string> infiles(fireFiles, fireFiles+16); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.createNcInputStream(infiles, glim);
	hires.printGrid();

	gVar lores("ba", "m2", "days since 2000-1-1");
	lores.setCoords(times, levs, lats, lons);
	lores.printGrid();
	lores.createNcOutputStream("/media/jaideep/Totoro/Data/Fire_BA/burned_area.2001-2016_0.5deg.nc");

	for (int t=0; t<lores.ntimes; ++t){
		hires.readVar_gt(lores.ix2gt(t)+4, 0);
		gVar tmp = coarseGrain_sum(hires, lores.lons, lores.lats);
		lores.copyValues(tmp);
		lores.writeVar(t);
		cout << "t = " << gt2string(lores.ix2gt(t)) << endl;
	}

	lores.closeNcOutputStream();
	hires.closeNcInputStream();	
	
//	lores.writeOneShot("GHS_pop_GPW42000_0.5deg_world.nc");


//	// test gsm_upper_bound and indexC 
//	int nlons;
//	vector <float> a = createCoord(ilim[0]+0.125, ilim[1]-0.125, 0.25, nlons);
//	
//	int wrong_count = 0;
//	for (int i=0; i<7200; ++i){
//		int my_upperbound = gsm_upper_bound(a, hires.lons[i]);
//		vector <float>::iterator stl_upperbound = upper_bound(a.begin(), a.end(), hires.lons[i]);
////		cout << hires.lons[i] << " " <<  *stl_upperbound << " " <<  a[my_upperbound] << "\n";
////		cout << hires.lons[i] << " " <<  distance(a.begin(), stl_upperbound) << " " <<  my_upperbound << "\n";
//		if (distance(a.begin(), stl_upperbound) != my_upperbound) ++wrong_count;
//	}
//	cout << "Wrong upper bounds = " << wrong_count << endl;
//
//
//	// indexC time test
//	for (int i=0; i<7200*3600; ++i){
//		indexC(a, -177);
//		indexC(a, -177);
//	}
	
	return 0;

}


