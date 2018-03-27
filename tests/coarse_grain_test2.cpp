#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test2.cpp -l:libgsm.so.2 -lnetcdf_c++ 


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
	vector <float> lats = createCoord(5.25,49.75,0.5,nlats);
	vector <float> lons = createCoord(60.25,99.75,0.5,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=0;

	gVar hires;
	hires.createOneShot("/media/jaideep/Totoro/Data/World_population_density/GHS_POP_GPW42000_GLOBE_R2015A_54009_1k_v1_0/wpd_2000_reprojected.nc", glim);
	hires.printGrid();

	gVar lores = coarseGrain_mean(hires, lons, lats);
	lores.printGrid();
	
	lores.writeOneShot("GHS_pop_GPW42000_0.5deg_world	.nc");


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


