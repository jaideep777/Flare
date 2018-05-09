#include <iostream>
#include <gsm.h>
#include <algorithm>
#include <cstring>
using namespace std;

// g++ -O3 -Wall -Wl,--no-as-needed -std=c++11 -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/Documents/libgsm_v2/include -L/home/jaideep/Documents/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 reverse_cube_test.cpp -pthread -l:libgsm.so.2 -lnetcdf_c++ 

//void reverseY(float v[], int nx, int ny){

//	for (int i=0; i<ny/2; ++i){
//		for (int j=0; j<nx; ++j){
//			swap(v[IX2(j,ny-i-1, nx)], v[IX2(j,i, nx)]);
////				temp = v[IX3(j,ny-i-1,k, nx,ny)];
////				v[IX3(j,ny-i-1,k, nx,ny)] = v[IX3(j,i,k, nx,ny)];
////				v[IX3(j,i,k, nx,ny)] = temp;
//		}
//	}

//}

int main(){

	vector <float> v(4*5*4);
	vector <float> w(4*5*4);
	for (int i=0; i<4*5*4; ++i) v[i]=i;
	
	printCube(&v[0], 4, 5, 4, 10);

	reverseCube(&v[0], 4, 5, 4);

	printCube(&v[0], 4, 5, 4, 10);

	return 0;

}




