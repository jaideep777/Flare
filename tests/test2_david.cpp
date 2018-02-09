#include <iostream>
#include <gsm.h>
#include <fstream>
using namespace std;

//g++ -I/usr/local/netcdf-c-4.4.1/include -I/usr/local/netcdf-cxx-legacy/include -I/home/david/Documents/libgsm/include -L/home/david/Documents/libgsm/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 test2.cpp -l:libgsm.so.2 -lnetcdf_c++ 

int main(){
	
	NcError err(NcError::silent_nonfatal);
	
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	float latlons[] = {74.8657202830595,77.4094184394495,8.29401732590651, 12.7914921513535};
	vector <float> latlonlimit(latlons, latlons+4);

	
	int nlons, nlats, nlevs, ntimes;
	int ilon, ilat, ilev;
	int lon, lat, lev;

	
	
	gVar pr;	
	pr.createOneShot("/home/david/slope_03_14.nc", latlonlimit);	
	pr.printGrid();
	
	gVar lai;	
	lai.createOneShot("/home/david/LAI_kerala_time.nc", latlonlimit);
	lai.printGrid();
	lai.writeOneShot("/home/david/manasi.nc");
	

   ofstream fout("trendoutput.txt");
	for (int ilev=0; ilev< lai.nlevs; ++ilev){
		for (int ilat=0; ilat< lai.nlats; ++ilat){
			for ( int  ilon=0; ilon< lai.nlons; ++ilon){
			   
			   fout<<lai.lons[ilon]<<"\t"<<lai.lats[ilat]<<"\t"<<lai(ilon,ilat,ilev)<<"\t";
			   fout<<pr.getValue(lai.lons[ilon],lai.lats[ilat])<<endl;
			}
		}
	}
	fout.close();
 
	
	    
	

	return 0;

}




