#include <iostream>
#include "../include/gsm.h"
using namespace std;

template <typename T>
BufferedVector<T>::BufferedVector(){
	len = 0;
	data = NULL;
	membuf = NULL;
}

template <typename T>
BufferedVector<T>::BufferedVector(int n){
	data = new T[n];
	membuf = new T[n];
	len = n;
}

template <typename T>
BufferedVector<T>::BufferedVector(int n, T initVal){
	data = new T[n];
	membuf = new T[n];
	len = n;
	for (int i=0; i<n; ++i) data[i] = initVal;
}

template<typename T> // copy constructor
BufferedVector<T>::BufferedVector(const BufferedVector& pv){
	len = pv.len;
    data = new T[len];
    membuf = new T[len];
    for(int i = 0; i < len; ++i) {
    	data[i] = pv.data[i];
    	membuf[i] = pv.membuf[i];
    }
}

template <typename T>
void BufferedVector<T>::print(string name){
	cout << name << " data-buf: " << len << " | ";
	for (int i=0; i<len; ++i){
		cout << data[i] << " ";
	}
	cout << endl;
	cout << name << "  mem-buf: " << len << " | ";
	for (int i=0; i<len; ++i){
		cout << membuf[i] << " ";
	}
	cout << endl;

}

template <typename T> // destructor
BufferedVector<T>::~BufferedVector(){
	if (len > 0) {
		delete [] data;
		delete [] membuf;
	}
}

template <typename T>
int BufferedVector<T>::size(){
	return len;
}

template <typename T>
void BufferedVector<T>::resize(int n, T fillVal){
	if (n <= 0){
		len = 0;
		delete [] data;
		delete [] membuf;
		data = NULL;
		membuf = NULL;
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
	
	delete [] membuf;
	membuf = new T[n];
}


template <typename T>
T& BufferedVector<T>::operator[](int i) const {
	return data[i];
}


template<typename T> //copy assignment operator
BufferedVector<T>& BufferedVector<T>::operator=(const BufferedVector& pv){
	delete [] data;
	delete [] membuf;
	data = new T[pv.len];
	membuf = new T[pv.len];
	len = pv.len;
	for(int i = 0; i < len; ++i) {
		data[i] = pv.data[i];
		membuf[i] = pv.membuf[i];
	}
	return *this;
}


template <typename T>
void BufferedVector<T>::swap(){
	T* temp;
	temp = data;
	data = membuf;
	membuf = temp;
}


