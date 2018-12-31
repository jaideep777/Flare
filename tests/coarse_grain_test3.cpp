#include <iostream>
#include <gsm.h>
#include <netcdf>
#include <vector>
#include <algorithm>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test3.cpp -l:libgsm.so.2 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
//	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
//	ofstream gsml("gsm_log.txt");
//	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {-180, 180, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lons = createCoord(-180+1./2,180-1./2,1,nlons);
	vector <float> lats = createCoord(-90+1./2,90-1./2,1,nlats);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=0;

	vector <string> infiles(1, "/home/jaideep/Data/Fire_BA_GFED4.1s/nc/GFED_4.1s_0.5deg.1997-2016.nc"); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.createNcInputStream(infiles, glim);
	hires.printGrid();

	gVar lores;
	lores.copyMeta(hires);
	lores.setCoords(lores.times, levs, lats, lons);
	lores.printGrid();
	lores.createNcOutputStream("/home/jaideep/Data/Fire_BA_GFED4.1s/nc/GFED_4.1s_1deg.1997-2016.nc");

	for (int t=0; t<lores.ntimes; ++t){
		hires.readVar_it(t);
		gVar tmp = coarseGrain_mean(hires, lores.lons, lores.lats);
		lores.copyValues(tmp);
		lores.writeVar(t);
		cout << "t = " << t << endl;
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


