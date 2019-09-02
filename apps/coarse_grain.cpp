#include <iostream>
#include <gsm.h>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <cmath>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test3.cpp -l:libgsm.so.2 -lnetcdf_c++ 


int main(int argc, char** argv){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
//	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
//	ofstream gsml("gsm_log.txt");
//	gsm_log = &gsml;

	float res = str2float(argv[1]);

	// create a grid limits vector for convenience
	float glimits[] = {-180, 180, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	vector <string> infiles(1, argv[2]); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.printGrid();
	
	float lat0 = floor(hires.lats[0]/res)*res;
	float latf = ceil(hires.lats[hires.lats.size()-1]/res)*res;
	float lon0 = floor(hires.lons[0]/res)*res;
	float lonf = ceil(hires.lons[hires.lons.size()-1]/res)*res;
	
	cout << lon0 << " " << lonf << " " << lat0 << " " << latf << endl;

	// create the coordinates for our georeferenced variable
	int nlons, nlats, nlevs, ntimes;
	vector <float> lons = createCoord_from_edges(lon0,lonf,res,nlons);
	vector <float> lats = createCoord_from_edges(lat0,latf,res,nlats);
	
	gVar lores;
	lores.copyMeta(hires);
	lores.setCoords(hires.times, hires.levs, lats, lons);
	lores.printGrid();
	lores.createNcInputStream(infiles, glim, "coarsegrain");

	if (hires.times.empty()){
		lores.readVar_it(0);
		lores.writeOneShot(argv[3]);
	}
	else{
		lores.createNcOutputStream(argv[3]);
		for (int t=0; t<lores.ntimes; ++t){
			lores.readVar_it(t);
			lores.writeVar(t);
			cout << "t = " << t << endl;
		}
		lores.closeNcOutputStream();
	}
	
	lores.closeNcInputStream();	
		
	return 0;

}




