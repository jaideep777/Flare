#include <iostream>
#include <gsm.h>
#include <netcdf>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <cmath>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 coarse_grain_test3.cpp -l:libgsm.so.2 -lnetcdf_c++ 
vector <string> split(string str, char at = '.'){
	std::replace(str.begin(), str.end(), at, ' ');  // replace ':' by ' '

	vector<string> array;
	stringstream ss(str);
	string temp;
	while (ss >> temp)	array.push_back(temp); 
	return array;
}

vector <string> readlines(string ifile){

	vector<string> array;
	ifstream fin(ifile);
	string temp;
	while (fin >> temp)	array.push_back(temp); 
	return array;
}


string extension(string filename){
	vector <string> array = split(filename, '.');
	return *array.rbegin();
}


int main(int argc, char** argv){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
//	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
//	ofstream gsml("gsm_log.txt");
	gsm_log = &cout;

	string ifile = argv[1];
	string ofile = argv[2];
	
	string ext = extension(ifile);
	cout << "Extension of input file is: " << ext << endl;	
	
	
	// create a grid limits vector for convenience
//	float glimits[] = {-180, 180, -90, 90};
	vector <float> glim = {-180, 180, -90, 90}; //(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	vector <string> infiles;
	
	if (ext == "nc") 		infiles = {ifile};
	else if (ext == "txt")  infiles = readlines(ifile);
	
	gVar in;
	in.initMetaFromFile(infiles[0]);
	in.printGrid();
	in.createNcInputStream(infiles, glim, "none");
	
	if (in.times.empty()){
		CERR << "Input file has no time dimension" << endl;
	}
	else{
		gVar P;
		gVar trend = in.yearlytrend(gt2year(in.ix2gt(0)), gt2year(in.ix2gt(in.ntimes-1)), &P);
		trend.writeOneShot(ofile);
		P.writeOneShot("P_"+ofile);
	}
	
	in.closeNcInputStream();	
		
	return 0;

}




