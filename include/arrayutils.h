#ifndef ARRAYUTILS_H
#define ARRAYUTILS_H

#include <vector>
#include <string>
#include <iostream>
#include "constants.h"
using namespace std;

/** \addtogroup utils 
	\{
*/

string int2str(int n);
float str2float(string s);
int str2int(string s);

// get 1D index of 3D point (ilon,ilat,ilev) from a cube (nlons,nlats,nlevs)
int IX3(int ix, int iy, int iz, int nx, int ny); 
// get 1D index of 2D point (ilon,ilat) from a cube (nlons,nlats)
int IX2(int ix, int iy, int nx);  

void printArray(float v[], int n, ostream &lfout=cout);	// print n elements of array on single line
void printArray(vector <float> &v, ostream &lfout = cout, string send="", int n=0);	// print float vector
void printArray2d(float v[], int rows, int columns);	// print 2d array with rows & columns
void printArray2d(vector <float> &v, int rows, int columns);	// print float vector 2d
void printCube(float v[], int nx, int ny, int nz=1, 
				float ignoreVal = std_missing_value); // print a data cube, ignore missing_values

//void setZero(vector <float> &v);	// set all of vector to zero
//void setValue(vector <float> &v, float value); // set all elements of vector to value

void reverseArray(vector <float> &orig);	// reverse given array
void reverseCube(float v[], int nx, int ny, 
				 int nz=1, int n4=1, int n5=1); // reverse a data cube along lat dimension

//vector <float> copyArray(vector <float> &v, int i2, int i1=0); // copy v[i1:i2] into new array (returned)

//int gsm_upper_bound (vector<float> &sorted_vec, float val, int first=0, int last=-1);

int ncIndexLo(vector <float> &v, float val); //!< lower bound, edge for outliers
int ncIndexHi(vector <float> &v, float val); //!< upper bound, edge for outliers
int lindexSW(vector <float> &v, float val);  //!< lower (S/W) bound, missing value for outliers
int indexC(vector <float> &v, float val);    //!< cell index by center, missing value for outliers
vector <float> max_vec(vector <float> &u, vector <float> &v);	//!< returns elementwise maximum 

// Array operations
float sum(vector <float> &v);	//!< Returns sum of vector
float avg(vector <float> &v);
//vector <float> dim_sum(vector <float> v, int idim, vector <float> dimsizes);
//vector <float> operator * (const vector <float> &x, const vector <float> &y);
//vector <float> operator / (const vector <float> &x, const float c);

/** \} */

#endif
