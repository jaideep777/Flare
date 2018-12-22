#include <iostream>
#include <math.h>
#include <netcdfcpp.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <map>
using namespace std;

#include <gsm.h>

/* --------------------------------------------------
compile command:
g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o sl standardize_lons.cpp -l:libgsm.so.2 -lnetcdf_c++  
-----------------------------------------------------*/


//int str2int(string s){
//	istringstream sin;
//	sin.str(s);
//	int val;
//	sin >> val;
//	return val;
//}


void read(int ilon0, int ilat0, int ilev0, int nlons, int nlats, int nlevs, float* source, float* dst){
//	cout << "reading: \n";
	for (int ilev=ilev0; ilev < ilev0+nlevs; ++ilev){
		for (int ilat=ilat0; ilat < ilat0+nlats; ++ilat){
			for (int ilon=ilon0; ilon < ilon0+nlons; ++ilon){
//				cout << IX3(ilon-ilon0,ilat-ilat0,ilev-ilev0, nlons, nlats) << " <-- " << IX3(ilon,ilat,ilev,12,2) << endl;
				dst[IX3(ilon-ilon0,ilat-ilat0,ilev-ilev0, nlons, nlats)] = source[IX3(ilon,ilat,ilev,12,2)];			
			}
		}
	}
}



int main(int argc, char ** argv){

	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);


//	string filename = argv[1]; 
	
//	gsm_info_on = true;
//	gsm_debug_on = true;

	float l[] = {0,30,60,90,120,150,180,210,240,270,300,330};

	float arr[] = {0,30,60,90,120,150,180,210,240,270,300,330,
				   0,3,6,9,12,15,18,21,24,27,30,33, 
				   
				   0,.30,.60,.90,.120,.150,.180,.210,.240,.270,.300,.330,
				   0,3.0,6.0,9.0,12.0,15.0,18.0,21.0,24.0,27.0,30.0,33.0};
	
	
	vector <float> lons(l,l+12);
	vector <float> lats(2); lats[0] = -60; lats[1] = 60;
	vector <float> levs(2); levs[0] = 0; levs[1] = 1;
	int nlons = 12, nlats = 2, nlevs = 2;
	
	vector <float> lons_orig = lons;
	
	cout << "   ";
	for (int i=0; i<nlons; ++i) cout << i << "  "; 
	cout << "\n   ";
	printArray(lons);	
	printCube(arr, nlons,nlats,nlevs);

	float glimits[] = {-180, -40, -40, 40};
	vector <float> glim(glimits, glimits+4);

	if (glim[0] > 180) glim[0] -= 360;
	if (glim[1] > 180) glim[1] -= 360;

	if (glim[1] < glim[0]){
		cout << "Incorrect Lon bounds: [" << glim[0] << ", " << glim[1] << "]\n";
		exit (1);
	}

	auto it = find_if(lons.begin(), lons.end(), [](float x){return x > 180;});	// find first lon > 180
	for_each(it, lons.end(), [](float &x){x -= 360;});							// bring all lons > 180 t0 principle range
	int shift = lons.size() - distance(lons.begin(), it);						// array to be shifted right by (n-it) elements

	cout << "Shift right by: " << shift << endl;

	if (it != lons.begin() && it != lons.end()){								// if shift is > 0 and < N, 
		shiftRight(lons.data(), lons.size(), shift);							//    shift lons 
	}

	cout << "   ";
	for (int i=0; i<nlons; ++i) cout << i << "  "; 
	cout << "\n   ";
	printArray(lons);	
	printCube(arr, nlons,nlats,nlevs);
	
	int wlonix1 = ncIndexLo(lons, glim[0]);
	int elonix1 = ncIndexHi(lons, glim[1]);


	vector <float> lons_trim(&lons[wlonix1], &lons[elonix1]+1);

	int wlonix = (wlonix1-shift+lons.size()) % lons.size();
	int elonix = (elonix1-shift+lons.size()) % lons.size();

	cout << "Calculated bound indices: " << wlonix << " " << elonix << endl;
	
	bool splitRead = (wlonix > elonix);
	
	vector <float> arr_trim;
	if (!splitRead){
		int nlons_trim = elonix-wlonix+1;
		arr_trim.resize(nlons_trim*nlats*nlevs);
		read(wlonix, 0, 0, nlons_trim, nlats, nlevs, arr, arr_trim.data());

		cout << "\n   ";
		printArray(lons_trim);	
		printCube(arr_trim.data(), nlons_trim, nlats, nlevs);
	}
	else{
		cout << "Reading in splits: [0, " << elonix << "], [" << wlonix << ", " << nlons-1 << "]\n";   
		int nlons1 = elonix-0+1;
		vector <float> seg1(nlons1*nlats*nlevs); 
		read(0, 0, 0, nlons1, nlats, nlevs, arr, seg1.data());

		int nlons2 = nlons-1-wlonix+1;
		vector <float> seg2(nlons2*nlats*nlevs); 
		read(wlonix, 0, 0, nlons2, nlats, nlevs, arr, seg2.data());
		
		printCube(seg1.data(), nlons1, nlats, nlevs);
		printCube(seg2.data(), nlons2, nlats, nlevs);
		
		int nlons_trim = nlons1+nlons2;
		arr_trim.resize(nlons_trim*nlats*nlevs);
		cout << "Resized array lons: " <<  nlons_trim << endl;
		
	
		for (int ilev=0; ilev < nlevs; ++ilev){
			for (int ilat=0; ilat < nlats; ++ilat){
				for (int ilon=0; ilon < nlons2; ++ilon){
					arr_trim[IX3(ilon,ilat,ilev, nlons_trim, nlats)] = seg2[IX3(ilon,ilat,ilev, nlons2, nlats)];
				}
			}
		}
		for (int ilev=0; ilev < nlevs; ++ilev){
			for (int ilat=0; ilat < nlats; ++ilat){
				for (int ilon=0; ilon < nlons1; ++ilon){
					arr_trim[IX3(ilon+nlons2,ilat,ilev, nlons_trim, nlats)] = seg1[IX3(ilon,ilat,ilev, nlons1, nlats)];
				}
			}
		}
		
		cout << "\n   ";
		printArray(lons_trim);	
		printCube(arr_trim.data(), nlons_trim, nlats, nlevs);

	}



//	// create a grid limits vector for convenience

	gVar vin;
	vin.createOneShot("data/gpp.intercept.2001-2010.nc", glim);
	vin.writeOneShot("gpp_out.nc");


}






