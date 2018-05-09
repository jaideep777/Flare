#include <iostream>
#include <gsm.h>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 buffered_vector_test.cpp -l:libgsm.so.2 -lnetcdf_c++ 

int main(){
	
	BufferedVector <float> v(5, 5.5);
	v.print("v");
	
	v.resize(2, 2.2);
	v.print("v");

	v.resize(8, 8.8);
	v.print("v");

	v[2] = 2;
	v.print("v");

	v.swap();
	v.print("v");
	

	BufferedVector <float> v2 = v;
	v2.resize(7);
	v.print("v");
	v2.print("v2");
	v[5] = 5;
	v.resize(6);
	
	BufferedVector <float> v3;
	v3 = v2 = v;
	v.print("v");
	v2.print("v2");
	v3.print("v3");

	v3.resize(0);
	v3.print("v3");


	return 0;
}


