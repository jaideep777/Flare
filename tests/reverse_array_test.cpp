#include <iostream>
#include <gsm.h>
#include <algorithm>
#include <cstring>
using namespace std;

// g++ -O3 -Wall -Wl,--no-as-needed -std=c++11 -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v3/include -L/home/jaideep/codes/libgsm_v3/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 reverse_array_test.cpp -pthread -l:libgsm.so.3 -lnetcdf_c++ 


int main(){

	vector <float> v(1000000);
	vector <float> w(1000000);
	for (int i=0; i<1000000; ++i) v[i]=i;
	
	// 
	clock_t start, end;
	for (int i=0; i<10; ++i) cout << v[i] << " ";
	cout << " ... ";
	for (int i=1000000-10; i<1000000; ++i) cout << v[i] << " ";
	cout << endl;
	start = clock();
	reverse(v.begin(), v.end());
	end = clock();
	cout << "Execution time for stl::reverse is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
	for (int i=0; i<10; ++i) cout << v[i] << " ";
	cout << " ... ";
	for (int i=1000000-10; i<1000000; ++i) cout << v[i] << " ";
	cout << endl;

	start = clock();
	reverseArray(v);
	end = clock();
	cout << "Execution time for brute reverse is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
	for (int i=0; i<10; ++i) cout << v[i] << " ";
	cout << " ... ";
	for (int i=1000000-10; i<1000000; ++i) cout << v[i] << " ";
	cout << endl;


	start = clock();
	memcpy(&w[0], &v[0], v.size()*sizeof(float));
	end = clock();
	cout << "Execution time for memcpy is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

	start = clock();
	copy(v.begin(), v.end(), w.begin());
	end = clock();
	cout << "Execution time for std::copy is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;



	start = clock();
	for (int i=0; i<1000; ++i) upper_bound(v.begin(), v.end(), 555555.5);
	end = clock();
	cout << "Execution time for stl::upper_bound is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

//	start = clock();
//	for (int i=0; i<1000; ++i) gsm_upper_bound(v, 555555.5, 0, v.size());
//	end = clock();
//	cout << "Execution time for gsm_upper_bound is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

	
	// test coarsegraining on a real large variable
	NcError err(NcError::silent_nonfatal);

	vector <float> ilim(4);
	ilim[0] = -180;
	ilim[1] = 180;
	ilim[2] = -90;
	ilim[3] = 90;
	
	gVar hires;
	hires.createOneShot("/media/jaideep/Totoro/Data/MODIS_50KM_06_16/converted_nc_files/LST_MOD11C3.A2006.nc", ilim);
	hires.printGrid();

	int nlons;
	vector <float> a = createCoord(ilim[0]+0.125, ilim[1]-0.125, 0.25, nlons);
	
//	start = clock();
//	for (int i=0; i<7200; ++i){
//		int my_upperbound = gsm_upper_bound(a, hires.lons[i]);
//	}
//	end = clock();
//	cout << "Execution time for gsm_upper_bound is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;


	start = clock();
	for (int i=0; i<7200; ++i){
		vector <float>::iterator stl_upperbound = upper_bound(a.begin(), a.end(), hires.lons[i]);
	}
	end = clock();
	cout << "Execution time for stl::upper_bound is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
	

//	vector <float> a = {0,1,2,3,4};
//	iota(a.begin(), a.end(), 5);
//	for_each(a.begin(), a.end(), [](float x){cout << x << " ";});
//	cout << endl;
//	
//	float missing_val= 7;
//	for_each(a.begin(), a.end(), [missing_val](float &x){ x = (x==missing_val)? missing_val:0;});
//	for_each(a.begin(), a.end(), [](float x){cout << x << " ";});
//	cout << endl;
	
	
	return 0;

}


