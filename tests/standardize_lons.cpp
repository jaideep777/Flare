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




int main(int argc, char ** argv){

	float lon0 = 0.25, lonf =  359.75, lat0 = -89.75, latf = 89.75, dx = 0.5, dy = 0.5, dt0 = 1, dlev = 1;
	int nlons, nlats, nlevs = 1;
	//const float glimits_custom[4] = {lon0, lonf, lat0, latf};

	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);


//	string filename = argv[1]; 
	
//	gsm_info_on = true;
//	gsm_debug_on = true;

	float arr[] = {1,2,3,4,5,6,7,8,9,
				   11, 22, 33, 44, 55, 66, 77, 88, 99, 
				   
				   1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9,
				   11.1, 22.2, 33.3, 44.4, 55.5, 66.6, 77.7, 88.8, 99.9};
	
	shiftCubeRight(arr, 9, 2, 2, 4);
	
	for (int k=0; k<2; ++k){
		for (int j=0; j<2; ++j){
			for (int i=0; i<9; ++i) cout << arr[IX3(i,j,k,9,2)] << " ";
			cout << "\n";	
		}
		cout << "\n";	
	}
	cout << endl;

	// create a grid limits vector for convenience
	float glimits[] = {0, 50, -90, 90};
	vector <float> glim(glimits, glimits+4);

	gVar vin;
	vin.createOneShot("data/gpp.intercept.2001-2010.nc", glim);
	vector <float> &lons = vin.lons;
//	printArray(lons);

	float wlon = 30, elon = 180;

	float l[] = {0,30,60,90,120,150,180,210,240,270,300,330};
//	vector <float> lons(l, l+12);

	auto it = find_if(lons.begin(), lons.end(), [](float x){return x > 180;});	// find first lon > 180
	int shift = lons.size() - distance(lons.begin(), it);						// array to be shifted right by that many elements
	for_each(it, lons.end(), [](float &x){x -= 360;});							// bring all lons > 180 t0 principle range
	
	printArray(lons);
	cout << "shift by: " << shift << endl;
	
	if (it != lons.begin() && it != lons.end()){								// if shift is > 0 and < N, shift both
		shiftRight(lons.data(), lons.size(), shift);							//    lons and data
//		shiftCubeRight(vin.values.data(), vin.nlons, vin.nlats, vin.nlevs, shift);
	}

	int wlonix1 = ncIndexLo(lons, wlon);
	int elonix1 = ncIndexHi(lons, elon);

	cout << "lon bnds = " << lons[wlonix1] << ", " << lons[elonix1] << endl;
	
	bool partialread = (lons[wlonix1]*lons[elonix1] < 0); // range includes 0, then we will have to read in 2 passes

	int wlonix = (wlonix1-shift+lons.size())%lons.size();
	int elonix = (elonix1-shift+lons.size())%lons.size();
	cout << "lon bnds = " << l[wlonix] << ", " << l[elonix] << endl;

	if (it != lons.begin() && it != lons.end()){								// if shift is > 0 and < N,
//		shiftRight(lons.data(), lons.size(), shift);							//    data
		shiftCubeRight(vin.values.data(), vin.nlons, vin.nlats, vin.nlevs, shift);
	}

	printArray(lons);

	vin.writeOneShot("gpp_sl.nc");
	
//	vin.initMetaFromFile("data/gpp.intercept.2001-2010.nc");
//	vin.createNcInputStream(vector<string>(1,"data/gpp.intercept.2001-2010.nc"), glim, "none");
//	vin.printGrid();
////	fti.printValues();	
//	
//	
//	gVar vout;
//	vout.copyMeta(vin);
//	vout.createNcOutputStream();
//	fto.nlons = nlons;
//	fto.values.resize(fto.nlons*fto.nlats*fto.nlevs, 0);

//	NcFile_handle fto_handle;
//	fto_handle.open(filename+"_sl.nc","w", glimits);
//	fto_handle.writeCoords(fto);
//	NcVar* vVar = fto_handle.createVar(fto);
//	fto_handle.writeTimeValues(fto);
//	fto.printGrid();

//	for (int t=0; t<fti.ntimes; ++t){

//		fti_handle.readVar(fti,t);

//		for (int ilon=0; ilon<fti.nlons; ++ilon){
//			for (int ilat=0; ilat<fti.nlats; ++ilat){
//				float xlon = fto.lons[ilon]; if (xlon > 180) xlon -= 360;
//				float xlat = fto.lats[ilat];
//				for (int ilev=0; ilev<fti.nlevs; ++ilev){
//					fto(ilon,ilat,ilev) = fti.getCellValue(xlon, xlat, ilev);
//				}
//			}	
//		}	

//		fto_handle.writeVar(fto, vVar, t); // write data at time index ix
//		cout << t << "\n";	
//	}

//	fto_handle.close();
//	fti_handle.close();
//		
	cout << "> Successfully wrote fire NC file!!\n";
}






