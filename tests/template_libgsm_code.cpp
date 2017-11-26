#include <iostream>
#include "../include/gsm.h"
#include <netcdfcpp.h>
#include <vector>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v2/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 streams_test.cpp -l:libgsm.so.2 -lnetcdf_c++ 


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

	return 0;
}
