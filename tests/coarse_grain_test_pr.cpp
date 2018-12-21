#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;



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
	int nlons, nlats, nlevs, ntimes;
	vector <float> lats = createCoord(11.0,12.0,0.25,nlats);
	vector <float> lons = createCoord(77.0,78.0,0.25,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1); 
	times[0]=ymd2gday("2003-6-1")-ymd2gday("2000-1-1");


	// indexC test
	vector <float> gr(10);
	for (int i=0; i<10; ++i) gr[i] = -89.75+0.5*i;
	
	vector <float> gr100(100);
	for (int i=0; i<100; ++i) gr100[i] = -90+0.05*i;

	printArray(gr);
	for (int i=0; i<100; ++i){
		cout << gr100[i] << " --> " << gr[indexC(gr, gr100[i])] << "\n";
	}

	cout << "----------------------------------------\n";

	// test coarsegraining on a simple dummy variable
	gVar v("ts", "deg C", "days since 2000-1-1 6:0:0");
	v.setCoords(times, levs, lats, lons);
	float a[] = {-1,0,0,1,1,
				 -2,1,1,2,2,
				 -2,1,1,2,2,
				 -3,3,3,4,4,
				 -3,3,3,4,4};
	v.values = vector <float> (a,a+25);
	v.printGrid();
	v.printValues();
	
	vector <float> y(4); y[0] = 11; y[1] = 11.5; y[2] = 12.0; y[3] = 12.5;
	vector <float> x(4); x[0] = 77; x[1] = 77.5; x[2] = 78.0; x[3] = 78.5;
	gVar vc = coarseGrain_sum(v, x, y);
	vc.printGrid();
	vc.printValues();

	
	// test coarsegraining on a real large variable
	vector <float> ilim(4);
	ilim[0] = 60.25;
	ilim[1] = 99.75;
	ilim[2] = 5.25;
	ilim[3] = 49.75;

	vector <float> lons1 = createCoord(60.25,99.75,0.5,nlons);
	vector <float> lats1 = createCoord(5.25,49.75,0.5,nlats);
	vector <float> levs1 = vector<float>(1,0);
	vector <double> t1 = vector<double>(1,0);
	
//	gVar hires;
//	hires.createOneShot("/home/jaideep/Data/precip_trmm/combined/reordered_dims/pr.2002.monthly.nc", ilim);

	gVar hires("pr0.5", "mm", "days since 1950-1-1");
	hires.setCoords(t1, levs1, lats1, lons1);
	hires.createNcInputStream(vector<string>(1,"/home/jaideep/Data/precip_trmm/combined/reordered_dims/pr.trmm-perm.2002.nc"), ilim, "bilinear");
	hires.readVar_reduce_mean(ymd2gday("2002-1-1"), ymd2gday("2002-1-31"));
	hires.printGrid();


//	vector<int> indices = cgIndices(hires.lons, hires.lats, lons1, lats1);
//	cout << "indices: " << indices.size() << endl;

//	gVar lores = coarseGrain_mean(hires, lons1, lats1, indices);
//	gVar lores = lterp(hires, lons1, lats1);
//	lores.printGrid();
	

	hires.writeOneShot("out_pr.nc");



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


