/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    libGSM - library for Gridded Spatial Modelling
    Copyright (C) 2016 Jaideep Joshi

	This file is part of libGSM.

    libGSM is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

	The author can be contacted at jaideep777@gmail.com 

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#include <algorithm>
#include <sstream>
#include <string>
#include "../include/gsm.h"


// convert from int to string
string int2str(int n){
	std::stringstream ss;
	ss << n;
	return ss.str();	
}

// convert from string to float
float str2float(string s){
	std::stringstream ss; ss.str(s);
	float f; ss >> f;
	return f;	
}

// convert from string to int
int str2int(string s){
	std::stringstream ss; ss.str(s);
	int f; ss >> f;
	return f;	
}


// get 1D index of 3D point (ix,iy,iz) from a cube (nx,ny,nz)
int IX3(int ix, int iy, int iz, int nx, int ny){
	return nx*ny*iz + nx*iy + ix;
}

// get 1D index of 2D point (ix,iy) from a cube (nx,ny)
int IX2(int ix, int iy, int nx){
	return nx*iy + ix;
}  


void printArray(float v[], int n, ostream &lfout){
	for (int i=0;i<n;++i) { lfout << v[i] << " ";}
	lfout << "\n"; 
}

void printArray(vector <float> &v, ostream &lfout, string send, int n){
	if (n==0) n = v.size();
	printArray(&v[0], n, lfout);
}

void printArray2d(float v[], int rows, int columns){
	for (int i=0; i<rows; ++i){
		for (int j=0;j<columns;++j) { cout << v[columns*i+j] << " ";}
		cout << "\n";
	}
	cout << "\n";
}

void printArray2d(vector <float> &v, int rows, int columns){
	printArray2d(&v[0], rows, columns);
}

void setValue(vector <float> &v, float value){
	for (int i=0; i<v.size(); ++i) v[i] = value;
}

void setZero(vector <float> &v){
	setValue(v, 0);
}

void reverseArray(vector <float> &orig){
	int b = orig.size();
	float swap;
	for(int a=0; a<--b; a++){ //increment a and decrement b until they meet eachother
		swap = orig[a];       //put what's in a into swap space
		orig[a] = orig[b];    //put what's in b into a
		orig[b] = swap;       //put what's in the swap (a) into b
	}
	//return orig;    //return the new (reversed) string (a pointer to it)
}

void printCube(vector <float> &v, int nx, int ny, int nz, float ignoreVal){
	if (nx*ny*nz != v.size()) {
		cout << "Error in printCube: dimensions mismatch!\n\n";
		return;
	}
	
	for (int k=0; k<nz; ++k){
		cout << "lev = " << k << ":\n";
		for (int i=0; i<ny; ++i){
			cout << "   ";
			for (int j=0; j<nx; ++j){
				if (v[IX3(j,i,k, nx,ny)] != ignoreVal) cout << v[IX3(j,i,k, nx,ny)] << " ";
				else cout << "--" << " ";
			}
			cout << "\n";
		}
		cout << "\n";	
	}	
}


void reverseCube(vector <float> &v, int nx, int ny, int nz, int n4, int n5){
	if (nx*ny*nz*n4*n5 != v.size()) {
		cout << "Error in reverseCube: dimensions mismatch!\n\n";
		return;
	}
	float temp;
	for (int k=0; k<nz; ++k){
		for (int i=0; i<ny/2; ++i){
			for (int j=0; j<nx; ++j){
				temp = v[IX3(j,ny-i-1,k, nx,ny)];
				v[IX3(j,ny-i-1,k, nx,ny)] = v[IX3(j,i,k, nx,ny)];
				v[IX3(j,i,k, nx,ny)] = temp;
			}
		}
	}	
}

bool ascComp(float a, float b){ return (a<b); }
bool dscComp(float a, float b){ return (a>b); }

// lower_bound = 1st element !< (>=) val
// upper_bound = last element !> (<=) val


// upper-bound (raw code from stl)
int gsm_upper_bound (vector<float> &sorted_vec, float val, int first, int last){
	if (last == -1 || last > sorted_vec.size()) last = sorted_vec.size();
	if (first < 0) first = 0;
	int count, step;
	int it;
	count = last - first;
	while (count > 0){
		int it = first; step=count/2; it += step;
		if (!(val< sorted_vec[it]))                 // or: if (!comp(val,*it)), for version (2)
			{ first=++it; count-=step+1;  }
		else count=step;
	}
	return first;
}

// returns lower bound in array for val. 
// if val is out of range, returns appropriate edge. does not return missing value
int ncIndexLo(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val >= v[v.size()-1]) return (v.size()-1);
		else if (val <= v[0]) return 0;
		else return (upper_bound(v.begin(), v.end(), val, ascComp) - v.begin() - 1);
	}
	else{
		if (val <= v[v.size()-1]) return (v.size()-1);
		else if (val >= v[0]) return 0;
		else return (lower_bound(v.begin(), v.end(), val, dscComp) - v.begin());
	}
}

// returns upper bound in array for val. 
// if val is out of range, returns appropriate edge. does not return missing value
int ncIndexHi(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val >= v[v.size()-1]) return (v.size()-1);
		else if (val <= v[0]) return 0;
		else return (lower_bound(v.begin(), v.end(), val, ascComp) - v.begin());
	}
	else{
		if (val <= v[v.size()-1]) return (v.size()-1);
		else if (val >= v[0]) return 0;
		else return (upper_bound(v.begin(), v.end(), val, dscComp) - v.begin() -1);
	}
}

// find index (m) of grid box such that P lies between m and m+1
// if val is out of range, returns missing value (-999)
//    1------2---P--3  P
//             ^G.C    ^outlier
int lindexSW(vector <float> &v, float val){
	bool asc = (v[1]>v[0])? true:false;
	if (asc){
		if (val > v[v.size()-1] || val < v[0]) return -999;	// if val exceeds edges, return -999
		else if (val == v[v.size()-1]) return v.size() -2;	// if val is on right edge, return 1 less
		else return (upper_bound(v.begin(), v.end(), val, ascComp) - v.begin() - 1);
	}
	else{
		// this case must not be used.. nonetheless, it may not cause trouble in interpolation
		CWARN << "lindexSW() invoked on descending vector!\n";
		if (val < v[v.size()-1] || val > v[0]) return -999;	// if val exceeds edges, return -999
		else if (val == v[v.size()-1]) return v.size() -2;	// if val is on right edge, return 1 less
		else return (upper_bound(v.begin(), v.end(), val, dscComp) - v.begin() -1);
	}
}

// find index (m) of grid box such that P lies closer to m than m+1 or m-1
// if val is out of range, returns missing value (-999)
// [---1---)[---2-P-)[---3-P-)     P
//               ^G.C    ^Sp.C   ^outlier
int indexC(vector <float> &v, float val){
	float dv = v[1]-v[0];
	bool asc = (dv>0)? true:false;
	if (asc){
		if (val >= v[v.size()-1]){
			if (val <= (v[v.size()-1] + (v[v.size()-1]-v[v.size()-2])/2)) return (v.size()-1);
			else return -999;
		}
		else if (val <= v[0]){
			if (val >= (v[0] - dv/2)) return 0;
			else return -999;
		} 
		else {
//			int m = distance(v.begin(), upper_bound(v.begin(), v.end(), val, ascComp)) - 1; // (last element <= val)
			int k = (val-v[0])/dv+1;
			int m = gsm_upper_bound(v, val, k-2, k+2) - 1; // (last element <= val) //TODO: deal with case where this is out of bounds
			return ((val - v[m]) < (v[m+1]-val))? m:(m+1);
		}
	}
	else{
		// this case must not be used.. 
		CERR << "indexC() invoked on descending vector!\n";
		return -999;
	}
}

// copy array[i1:i2] into new array
vector <float> copyArray(vector <float> &v, int i2, int i1){
	// if i1 > i2 swap
	if (i1 > i2){
		int temp = i1; i1 = i2; i2 = temp;
	}
	int n = i2-i1+1;
	vector <float> w;
	for (int i=0; i<n; ++i){
		w.push_back(v[i1+i]);
	}
	return w;
}


vector <float> max_vec(vector <float> &u, vector <float> &v){
	vector <float> temp(0);
	if (v.size() != u.size()){
		cout << "Error in max_vec: vectors not of same size!\n";
		return temp;
	}
	
	temp.resize(v.size());
	for (int i=0; i<temp.size(); ++i){
		temp[i] = (u[i] > v[i])? u[i]:v[i];
	}
	return temp;
}

// operations on vectors

float sum(vector <float> &v){
	float vsum = 0;
	for (int i=0; i<v.size(); ++i){
		vsum += v[i];
	}
	return vsum;
}

float avg(vector <float> &v){
	return sum(v)/v.size();
}

vector <float> dim_sum(vector <float> v, int idim, vector <float> dimsizes){

}


vector <float> operator * (const vector <float> &x, const vector <float> &y){
	vector <float> temp;
	if (x.size() != y.size()){
		CERR << "Vector multiplication: input vectors not compatible!!\n";
	}
	else{
		temp.resize(x.size());
		for (int i=0; i<temp.size(); ++i){
			temp[i] = x[i]*y[i];
		}
	}
	return temp;
}

vector <float> operator / (const vector <float> &x, const float c){
	vector <float> temp = x;
	for (int i=0; i<temp.size(); ++i){
		temp[i] /= c;
	}
	return temp;
}


//// index of cell containing given value
//// uses binary search, assumes array is in ascending order. therefore latSN required
//int cellIndex(float val, vector <float> &v){
//	int n = v.size();
//	int hi = n-1, lo=0;
//	float dxHi, dxLo;
//	int mid;
//	bool goUp=true;

////	if (v[n-1] < val) return n-1;
////	else if (v[0] > val) return 0;
//	bool asc = true;
//	if (v[1] < v[0]) asc = false;

//	while (1){
//		mid = (hi+lo)/2;

//		if (mid < n) dxHi = (v[mid+1] - v[mid])/2;
//		else dxHi = (v[mid] - v[mid-1])/2;
//		if (mid > 0) dxLo = (v[mid] - v[mid-1])/2;
//		else dxLo = (v[mid+1] - v[mid])/2;
//		
////		cout << lo << "\t" << mid << "\t" << hi << "\t(" 
////			 << v[mid]-dxLo << ", " << v[mid] << ", " << v[mid]+dxHi << ")\n"; // << val;
//		if ((v[mid] - dxLo) < val) goUp = asc;
//		if ((v[mid] + dxHi) > val) goUp = !asc;
//		if ((v[mid] - dxLo < val) && (v[mid] + dxHi) > val) return mid;
//		else if (lo == mid) return mid;

//		if (goUp) lo = mid; 
//		else hi = mid; 
//	}
//	return mid;	
//}

//// binary search: returns index of largest number smaller than value.
//inline int binsearch(float val, vector <float> &v){
//	int n = v.size();
//	int hi = n-1, lo=0, mid;
//	bool goUp = true;
//	bool asc = (v[1] > v[0])? true:false;

//	if (asc){
//		if (v[n-1] < val) return (n-1);
//		if (v[0] > val) return 0;
//	}
//	else{
//		if (v[n-1] > val) return (n-1);
//		if (v[0] < val) return 0;
//	}

//	// the actual search algo!
//	CDEBUG << "searching for: " << val << '\n';
//	while (1){
//		mid = (hi+lo)/2;
//		CDEBUG << lo << "\t" << mid << "\t" << hi << "\t(" 
//			 << v[mid-1] << ", " << v[mid] << ", " << v[mid+1] << ")\n"; // << val;
//		if (v[mid] < val) goUp = asc;
//		if (v[mid] > val) goUp = !asc;
//		if ((v[mid] < val) && (v[mid+1] > val)) break;
//		if (lo == mid) break;

//		if (goUp) lo = mid; 
//		else hi = mid; 
//	}
//	
//	return mid;
//}

//// index of largest number less than val, i.e. index of lower edge
//// uses binary search. v could be in ascending or descending order.
//int indexLow(float val, vector <float> &v){
//	bool asc = (v[1] > v[0])? true:false;
//	int mid = binsearch(val, v);
//	if (mid == 0 || mid == v.size()-1) return mid;
//	else return (asc)? mid:mid+1;
//}

//// index of smallest number greater than val, i.e. index of upper edge
//// uses binary search. v could be in ascending or descending order.
//// note: if v is ascending, indexHi is always = indexLow + 1, unless it exceeds bounds
//// there is no way to check if bounds have been exceeded, this must be checked seperately by the caller
//int indexHi(float val, vector <float> &v){
//	bool asc = (v[1] > v[0])? true:false;
//	int mid = binsearch(val, v);
//	if (mid == 0 || mid == v.size()-1) return mid;
//	else return (asc)? mid+1:mid;
//}


