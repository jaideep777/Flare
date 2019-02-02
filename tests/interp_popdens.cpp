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
	vector <float> lons = createCoord(-180+0.25/2,180-0.25/2,0.25,nlons);
	vector <float> lats = createCoord(-90+0.25/2,90-0.25/2,0.25,nlats);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(16*12); 
	for (int i=0; i<times.size(); ++i) times[i] = ymd2gday("2000-1-1") - ymd2gday("1950-1-1") + i*365.2524/12 + 15;
	
	double t2000 = ymd2gday("2000-6-15") - ymd2gday("1950-1-1");
	double t2015 = ymd2gday("2015-6-15") - ymd2gday("1950-1-1");
	
	gVar pop2000;
	pop2000.createOneShot("pop2000_0.25deg.nc");
//	pop2000 = coarseGrain_mean(pop2000, lons, lats);
	pop2000.printGrid();
//	pop2000.writeOneShot("pop2000_0.25deg.nc");

	gVar pop2015;
	pop2015.createOneShot("pop2015_0.25deg.nc");
//	pop2015 = coarseGrain_mean(pop2015, lons, lats);
	pop2015.printGrid();
//	pop2015.writeOneShot("pop2015_0.25deg.nc");

	gVar lores("pop", "/km2", "days since 1950-1-1");
	lores.setCoords(times, levs, lats, lons);
	lores.printGrid();
	lores.createNcOutputStream("/media/jaideep/Totoro/Data/World_population_density/GHS_POP_GPW4.2000-2015.nc");

	for (int t=0; t<lores.ntimes; ++t){
		
		for (int i=0; i<nlons*nlats*nlevs; ++i){
			if (pop2000[i] == pop2000.missing_value) lores[i] = lores.missing_value;
			else if (pop2000[i] < 1e-10 || pop2015[i] < 1e-10) lores[i] = 0;
			else lores[i] = pow(pop2015[i]/(pop2000[i]+1e-10), (times[t]-t2000)/(t2015-t2000))*pop2000[i];
		}
		lores.writeVar(t);
		cout << "t = " << t << endl;
	}

	lores.closeNcOutputStream();
//	hires.closeNcInputStream();	
//	
////	lores.writeOneShot("GHS_pop_GPW42000_0.5deg_world.nc");


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


