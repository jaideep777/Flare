#include <iostream>
#include <gsm.h>
#include <fstream>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v2/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 precip_generator.cpp -l:libgsm.so.2 -lnetcdf_c++ 

int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {66.5, 100.5, 6.5, 38.5};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	gVar pr;
	pr.createOneShot("/media/jaideep/Totoro/Data/precip_IMD/precip_imd/slope_1901-2015.nc", glim);
	pr.printGrid();

	ifstream fin("/media/jaideep/WorkData/Rajiv_carbon_project/Export_Output_india_rp.txt");	
	string s;
	getline(fin, s);

	ofstream fout("/media/jaideep/WorkData/Rajiv_carbon_project/india_rp_rf_1901-2015.txt");
	fout << "X\tY\tNETCARBON1\tRF\n";
	float lat, lon, carb, rf;

	int count = 0;
	while (!fin.eof()){	// keep reading from txt file until file ends
		fin >> 	lon >> lat >> carb;		// read lon, lat, carbon from txt file
		rf = pr.getValue(lon, lat, 0);	// read corresponding precip from nc file (getValue will interpolate and give value for lon,lat)
		fout << lon << "\t" << lat << "\t" << carb << "\t" << rf << endl;
		++count;
		if (count % 1000 == 0) cout << count << endl;
	}


	return 0;

}




