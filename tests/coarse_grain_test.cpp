#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v2/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test.cpp -l:libgsm.so.2 -lnetcdf_c++ 


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

//	vector <float> v1(1e9);
//	int a1; cin >> a1;
	
	// test coarsegraining on a real large variable
	vector <float> ilim(4);
	ilim[0] = 66.5;
	ilim[1] = 100.5;
	ilim[2] = 6.5;
	ilim[3] = 38.5;
	
	gVar hires;
	hires.createOneShot("/home/jaideep/codes/Rajiv_carbon_project/MOD17A2_GPP.2000.M01.nc", ilim);
//	hires.ntimes=1; hires.times=vector <double> (1,0);	// This is somehow still required else program segfaults
	hires.printGrid();
//	hires.printValues();

//	hires.createNcOutputStream("out.nc");
//	hires.writeVar(0);
//	hires.closeNcOutputStream();

	vector <float> lats1 = createCoord(ilim[2],ilim[3],0.25,nlats);
	vector <float> lons1 = createCoord(ilim[0],ilim[1],0.25,nlons);
	gVar lores = coarseGrain_mean(hires, lons1, lats1);
	lores.printGrid();
	
	lores.createNcOutputStream("out.nc");
	lores.writeVar(0);
	lores.closeNcOutputStream();
	
	return 0;

}


