#include <iostream>
#include <gsm.h>
#include <fstream>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v2/include -L/home/jaideep/codes/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 correl_nc_and_nc.cpp -l:libgsm.so.2 -lnetcdf_c++ 

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

	
	// to correlate 2 nc files, they must be at the same resolution.
	// let us define coordinate system in which to read both files (interpolated)
	int nlons, nlats, nlevs, ntimes;
	float res = 0.25;
	vector <float> lats = createCoord(glim[2]+res/2, glim[3]-res/2, res,nlats);
	vector <float> lons = createCoord(glim[0]+res/2, glim[1]-res/2, res,nlons);
	vector <float> levs = createCoord(1,1,1,nlevs);
	vector <double> times(1, ymd2gday("2003-6-1")-ymd2gday("2000-1-1")); 

	// since trends are only 1 time slice, we can use the readOneShot() function instead of creating a stream.

	gVar pr("pr", "mm/day", "days since 2000-1-1");	// trend file #1
	pr.setCoords(times, levs, lats, lons);	// set coordinate system for the gVar
	pr.readOneShot("/home/jaideep/codes/Rajiv_carbon_project/gpp_2000-01_0.25.nc", glim);	// readOneShot reads the slice from NC file and interpolates to specifed coordinates
	pr.printGrid();

	// you can write out an nc file in one command if you want to check that it was read correctly
//	pr.writeOneShot("out3.nc");

	gVar lai("lai", "NA", "days since 2000-1-1");	// trend file #2
	lai.setCoords(times, levs, lats, lons);
	lai.readOneShot("/home/jaideep/codes/Rajiv_carbon_project/gpp_2000-01_0.25.nc", glim);
	lai.printGrid();


	// open text file for writing pr-lai pairs
	ofstream fout("ncnc.txt");
	for (int ilev=0; ilev< lai.nlevs; ++ilev){
		for (int ilat=0; ilat< lai.nlats; ++ilat){
			for (int ilon=0; ilon< lai.nlons; ++ilon){
				fout << levs[ilev] << "\t" << lats[ilat] << "\t" << lons[ilon] << "\t";	// write coord
				fout << pr(ilon,ilat,ilev) << "\t" << lai(ilon, ilat, ilev) << endl; // write pr, lai
			}
		}
	}
	fout.close();
	

	return 0;

}




