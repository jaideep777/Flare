#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;

// g++ -O3 -I/usr/local/netcdf-c/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v3/include -L/home/jaideep/codes/libgsm_v3/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 lage_file_test.cpp -l:libgsm.so.3 -lnetcdf_c++ 


int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
//	ofstream gsml("gsm_log.txt");
//	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {-180, 180, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(-90+0.25,90-0.25,0.5,nlats);
	vector <float> lons = createCoord(-180+0.25,180-0.25,0.5,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=0;

	gVar hires;
	hires.createOneShot("data/GHS_pop_GPW42000_reprojected_Globe.nc", glim);
	hires.printGrid();

	clock_t start, end;
	start = clock();
	hires.fill(0);
	end = clock();
	cout << "Execution time for coarsegraining is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

//	clock_t start, end;
	start = clock();
	fill(hires.values.begin(), hires.values.end(), 0);
	end = clock();
	cout << "Execution time for coarsegraining is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

//	clock_t start, end;
	start = clock();
	float mv = hires.missing_value;
	for_each(hires.values.begin(), hires.values.end(), [mv](float &x){ x = (x== mv)? mv:0;});
	end = clock();
	cout << "Execution time for coarsegraining is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

	vector <float> v(5, 0);
	v[2] = hires.missing_value;
	for_each(v.begin(), v.end(), [mv](float &x){ x = (x== mv)? mv:5;});
	for_each(v.begin(), v.end(), [](float x){ cout << x << " ";});
	cout << endl;

//	clock_t start, end;
//	start = clock();
//	gVar lores = coarseGrain_mean(hires, lons, lats);
//	end = clock();
//	cout << "Execution time for coarsegraining is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
//	lores.printGrid();
//	
//	lores.writeOneShot("GHS_pop_GPW42000_0.5deg.nc");


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


