#include <iostream>
#include <gsm.h>
#include <algorithm>
#include <cstring>
#include <thread>
#include <future>
using namespace std;

// g++ -O3 -Wall -Wl,--no-as-needed -std=c++11 -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/libgsm_v3/include -L/home/jaideep/codes/libgsm_v3/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 mutex_test.cpp -pthread -l:libgsm.so.3 -lnetcdf_c++ 

int increment_data(vector <float> &v){
//	usleep(5e6);
	for (int i=0; i<200; ++i)
	for (int i=0; i<v.size(); ++i) ++v[i];
	return 0;
}

void increment_membuf(vector <float> &v){
//	usleep(5e6);
	for (int i=0; i<200; ++i)
	for (int i=0; i<v.size(); ++i) ++v[i];
}
 

int main(){
	
	int n=10000000;
	vector <float> v(n);
	for (int i=0; i<n; ++i) v[i]=i;

	vector <float> v2 = v;
	
	// 
	clock_t start, end;
	start = clock();
	for (int i=0; i<200; ++i)
	for (int i=0; i<n; ++i) ++v[i];
	end = clock();
	cout << "Execution time for increment is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
	for (int i=0; i<200; ++i)
	for (int i=0; i<n; ++i) ++v2[i];
	end = clock();
	cout << "Execution time for increment is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

	start = clock();
	future <int> fut = async(launch::deferred, increment_data, ref(v));
	increment_data(ref(v2)); 
	int i = fut.get();
	end = clock();
	cout << "Execution time for increment is " << ((double) (end - start)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;

	for (int i=0; i<5; ++i) cout << v[i] << " ";
	cout << endl;
	for (int i=0; i<5; ++i) cout << v2[i] << " ";
	cout << endl;
	
	return 0;

}





