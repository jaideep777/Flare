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

	// create a grid limits vector for convenience
	float glimits[] = {-180, 360, -90, 90};
	vector <float> glim(glimits, glimits+4);

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


	gVar vin;
	vin.createOneShot("data/gpp.intercept.2001-2010.nc");

	vector <float> &lons = vin.lons;
	printArray(lons);

	auto it = find_if(lons.begin(), lons.end(), [](float x){return x > 180;});	// find first lon > 180
	int shift = distance(lons.begin(), it);										// array to be shifted by that many elements
	for_each(it, lons.end(), [](float &x){x -= 360;});							// bring all lons > 180 t0 principle range
	
	if (it != lons.begin() && it != lons.end()){								// if shift is > 0 and < N, shift both
		shiftRight(lons.data(), lons.size(), shift);							//    lons and data
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






