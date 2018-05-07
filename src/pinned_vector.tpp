#include <iostream>
#include "../include/gsm.h"
using namespace std;

template <typename T>
PinnedVector<T>::PinnedVector(){
	len = 0;
	data = NULL;
}

template <typename T>
PinnedVector<T>::PinnedVector(int n){
	data = new T[n];
	len = n;
}

template <typename T>
PinnedVector<T>::PinnedVector(int n, T initVal){
	data = new T[n];
	len = n;
	for (int i=0; i<n; ++i) data[i] = initVal;
}

template<typename T> // copy constructor
PinnedVector<T>::PinnedVector(const PinnedVector& pv){
	len = pv.len;
    data = new T[len];
    for(int i = 0; i < len; ++i) data[i] = pv.data[i];
}

template <typename T>
void PinnedVector<T>::print(string name){
	cout << name << ": " << len << " | ";
	for (int i=0; i<len; ++i){
		cout << data[i] << " ";
	}
	cout << endl;
}

template <typename T> // destructor
PinnedVector<T>::~PinnedVector(){
	if (len > 0) {
		delete [] data;
	}
}

template <typename T>
int PinnedVector<T>::size(){
	return len;
}

template <typename T>
void PinnedVector<T>::resize(int n, T fillVal){
	if (n <= 0){
		len = 0;
		delete [] data;
		data = NULL;
	}
	
	T * data_new = new T[n];
	if (n <= len){ // len >= n
		for (int i=0; i<n; ++i) data_new[i] = data[i];
	}
	else{	// len < n
		for (int i=0; i<len; ++i) data_new[i] = data[i];
		for (int i=len; i<n; ++i) data_new[i] = fillVal;
	}
	len = n;
	delete [] data;
	data = data_new;
}


template <typename T>
T& PinnedVector<T>::operator[](int i){
	return data[i];
}


template<typename T> //copy assignment operator
PinnedVector<T>& PinnedVector<T>::operator=(const PinnedVector& pv){
    delete [] data;
    data = new T[pv.len];
    len = pv.len;
    for(int i = 0; i < len; ++i) data[i] = pv.data[i];
    return *this;
}





