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

#include <iostream>
#include <math.h>
#include "../include/gsm.h"


// create coordinate vector given x0, xf and number of points
vector <float> createCoord(float x0, float xf, int nx, float &dx){
	dx = (xf-x0)/(nx-1);		// use number of points to set dx
	vector <float> xvec(nx);	// create variable
	for (int i=0; i<nx; ++i) {xvec[i] = x0 + i*dx;}
	return xvec;
}

// create coordinate vector given x0, xf and resolution
vector <float> createCoord(float x0, float xf, float dx, int &nx){
	float xres = 1/dx;	
	nx = floor(fabs(xf-x0)*xres)+1;	// use resolution to set nx
	x0 += ((xf-x0) - (nx-1)*dx)/2;
	vector <float> xvec(nx);	// create variable
	for (int i=0; i<nx; ++i) {xvec[i] = x0 + i*dx;}
	return xvec;
}

// print a gridded variable along with 2d coordinates defined on x, y
void printVar(vector <float> &x, vector <float> &y, vector <float> &data){
	int nx = x.size(), ny = y.size();
	cout << "\nx->\t ";
	for (int i=0; i< nx; ++i) cout << x[i] << "\t ";
	cout << "\ny\t";
	for (int i=0; i< nx; ++i) cout << "-------\t";
	cout << '\n';
	for (int j=0; j< ny; ++j){
		cout << y[j] << "\t| ";
		for (int i=0; i< nx; ++i){
			cout << data[nx*j+i] << "\t ";
		}
		cout << '\n';
	}
}

// ****************************************************************
//	      f11     f21  
//	------(O)----- 2 ----- y[m]      SW  S  SE      lons --->// 
//		   |  P    |                  W  +  E     lats
//	       |       |                 NW  N  NE      |
//		   |       |                                v
//	------ 1 ----- 3 ----- y[m+1]      P ~ (x,y)         
//		  f12      f22                 O ~ [m,n] ~ uvs ~ 0
//		 x[m]    x[m+1]
//
// f = [f11, f12, f21, f22]
// ****************************************************************
// find indices of SW corner (0) of gridbox (points 1-4) containing point (x,y)
// --------------------------------------------------------------------------
vector <int> findGridBoxSW(float x,float y, vector <float> &lons, vector <float> &lats){
	vector <int> uvs(2);
	uvs[0] = lindexSW(lons, x);
	uvs[1] = lindexSW(lats, y);
	
	return uvs;
}

// ****************************************************************
//     +---v---+---v---+
//	   |       |       |
//     >  (O)  +   2   <   y[m]      SW  S  SE      lons --->
//     |     P |       |              W  +  E     lats
//	   ----+---+---+----             NW  N  NE      |
//	   |       |       |                            v
//	   >   1   +   3   <   y[m+1]      P ~ (x,y)         
//	   |       |       |               O ~ [m,n] ~ uvs ~ 0
//     +---^---+---^---+
//		 x[m]    x[m+1]
// f = [f11, f12, f21, f22]
// ****************************************************************
// find indices of center (0) of gridbox (point O) containing point (x,y)
// --------------------------------------------------------------------------
vector <int> findGridBoxC(float x,float y, vector <float> &lons, vector <float> &lats){
	vector <int> uvs(2,-999);
	uvs[0] = indexC(lons, x);
	uvs[1] = indexC(lats, y);
	
	return uvs;
}	


// During initialization, calculate the SW index for each point of model grid, 
// so that it can be resued at every step in bilinear
vector <int> bilIndices(vector <float> &lons, vector <float> &lats,
						vector <float> &mlons, vector <float> &mlats){

//	int gsm_nlats = mlats.size(), gsm_nlons = mlons.size();
	vector <int> indices(mlats.size()*mlons.size()*2);
	for (int ilat=0; ilat < mlats.size(); ++ilat){
		for (int ilon=0; ilon < mlons.size(); ++ilon){
			vector <int> uvs = findGridBoxSW(mlons[ilon], mlats[ilat], lons, lats);
			indices[IX2(ilon, ilat, mlons.size())*2 + 0] = uvs[0];	// there was a fatal mistake on this line
			indices[IX2(ilon, ilat, mlons.size())*2 + 1] = uvs[1];	// mlats.size() was being used!!!
		}
	}
	return indices;
}


inline float bilinear_mn(int m, int n, float x, float y, float iz,
				   vector <float> &lons, vector <float> &lats, 
				   vector <float> &data, float missingVal){
	// m, n are input

	// check if (x,y) lies in grid
	if (m < 0 || n < 0) return missingVal;

	// ------> uvs are OK.
	float x1 = lons[m], x2 = lons[m+1], y1 = lats[n], y2 = lats[n+1];
//	CDEBUG << "Lon bounds: (" << x1 << " < " << x << " < " << x2 << ")\n";
//	CDEBUG << "Lat bounds: (" << y1 << " < " << y << " < " << y2 << ")\n";
	// get f at nearest neighbors
	vector <float> f(4,missingVal);
	f[0] = data[IX3(m,   n,   iz,   lons.size(), lats.size())];
	f[1] = data[IX3(m,   n+1, iz,   lons.size(), lats.size())];
	f[2] = data[IX3(m+1, n,   iz,   lons.size(), lats.size())];
	f[3] = data[IX3(m+1, n+1, iz,   lons.size(), lats.size())];
	// check if any function value is missing.. if so, return missingval
	if (f[0] == missingVal || f[2] == missingVal || 
		f[1] == missingVal || f[3] == missingVal )   return missingVal;

	// ------> f is OK.
	// compute lterp value
	float xbil = (f[0] * ((x2-x )*(y2-y ))   // f11
				+ f[2] * ((x -x1)*(y2-y ))   // f21
				+ f[1] * ((x2-x )*(y -y1))   // f12
				+ f[3] * ((x -x1)*(y -y1)))  // f22
				/ ((x2-x1)*(y2-y1));		// x2=x1+dx, so no div-by-zero!
	//cout << "lterp value = " << xbilinear << '\n';
	return xbil;

}

// return bilinear interpolated value at point (x,y) from grid {lons,lats}
float bilinear(float x, float y, float iz, 
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &data, float missingVal){
	// get nearest neighbor indices
	vector <int> uvs = findGridBoxSW(x,y,lons, lats);
	int m = uvs[0], n = uvs[1];
	
	return bilinear_mn(m, n, x, y, iz, 
						lons, lats, 
						data, missingVal);
}

// return bilinear interpolated value at index (ilon,ilat) of grid {mlons, mlats} from grid {lons,lats}
// use bilIndices instead of finding them in situ
// {lons, lats} is input grid and {mlons, mlats} is output (typically model) grid
// data is of course, defined on {lons, lats}
float bilinear(int ilat, int ilon, int iz, vector <int> &indices,
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &mlons, vector <float> &mlats,
			   vector <float> &data, float missingVal){
	int nlons, nlats;
	
	// get nearest neighbor indices
	int m = indices[IX2(ilon, ilat, mlons.size())*2 + 0]; // this ID uses model nlons, nlats
	int n = indices[IX2(ilon, ilat, mlons.size())*2 + 1]; // because it is on model grid

	float x = mlons[ilon], y = mlats[ilat];

	return bilinear_mn(m, n, x, y, iz, 
						lons, lats, 
						data, missingVal);
}

inline float cellVal_mn(int m, int n, float x, float y, float iz, 
						vector <float> &lons, vector <float> &lats, 
						vector <float> &data, float missingVal){
	// m, n are input
	// check if (x,y) lies in grid
	if (m < 0 || n < 0) return missingVal;

	// ------> uvs are OK.
//	CDEBUG << "Cell found: (" << lons[m] << ", " << lats[n] << ")\n";
	return data[IX3(m,n,iz, lons.size(),lats.size())];
						
}

// return value at gridcell containing point (x,y) in grid {lons,lats}
float cellVal(float x, float y, float iz, 
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &data, float missingVal){
	// get index of cell containing (x,y)
	vector <int> uvs = findGridBoxC(x,y,lons, lats);
	int m = uvs[0], n = uvs[1];

	return cellVal_mn(m, n, x, y, iz, 
					  lons, lats, 
					  data, missingVal);
	
}

float cellVal(int ilat, int ilon, int iz, vector <int> &indices,
			  vector <float> &lons, vector <float> &lats, 
			  vector <float> &mlons, vector <float> &mlats,
			  vector <float> &data, float missingVal){
	// get nearest neighbor indices
	int m = indices[IX2(ilon, ilat, mlons.size())*2 + 0]; // this ID uses model nlons, nlats
	int n = indices[IX2(ilon, ilat, mlons.size())*2 + 1]; // because it is on model grid

	float x = mlons[ilon], y = mlats[ilat];

	return cellVal_mn(m, n, x, y, iz, 
					  lons, lats, 
					  data, missingVal);
}


// ------------------- OPERATIONS ON GVARS ---------------------------

// mask variable v using m as mask
// all values in v where (m <= val) are set to missing value
gVar mask(gVar &v, gVar &m, float val){
	if (v.nlons != m.nlons || v.nlats != m.nlats){
		CERR << "ERROR in masking: grids do not match.\n";
		gVar temp;
		return temp;
	}

	gVar temp(v); // must use deep copy because unmasked values from v must be retained.
	for (int ilev=0; ilev < v.nlevs; ++ilev){
		for (int ilat=0; ilat < v.nlats; ++ilat){
			for (int ilon=0; ilon < v.nlons; ++ilon){
				if (m(ilon,ilat,0) <= val || m(ilon,ilat,0) == m.missing_value) 
					temp(ilon,ilat,ilev) = v.missing_value;
			}
		}
	}
	return temp;
}



// copy values from v to out with cellRegridding (on grid of out)
// cellRegrid = regrid using cell values rather than interpolated ones
int cellRegridCube(gVar &v, gVar &out, vector <int> &indices){
	// get output lats lons from output var
	vector <float> xlons = out.lons, xlats = out.lats;
	// grids are OK. proceed with regridding...
	out.values.resize( xlons.size()*xlats.size()*v.nlevs );
	for (int ilev=0; ilev < v.nlevs; ++ilev){
		for (int ilat=0; ilat < xlats.size(); ++ilat){
			for (int ilon=0; ilon < xlons.size(); ++ilon){
				//float f = v.getValue(xlons[ilon], xlats[ilat], ilev);
				float f = cellVal(ilat, ilon, ilev, indices,
								   v.lons, v.lats, 
								   xlons, xlats,
								   v.values, v.missing_value);
				if (f != v.missing_value) out(ilon,ilat,ilev) = f;
				else out(ilon,ilat,ilev) = out.missing_value;
			}
		}
	}
	out.t = v.t; // copy the current time value because it is related to values
}

// copy values from v to out with interpolation (on grid of out)
int lterpCube(gVar &v, gVar &out, vector <int> &indices){
	// get output lats lons from output var
	vector <float> xlons = out.lons, xlats = out.lats;
	// grids are OK. proceed with regridding...
	out.values.resize( xlons.size()*xlats.size()*v.nlevs );
	for (int ilev=0; ilev < v.nlevs; ++ilev){
		for (int ilat=0; ilat < xlats.size(); ++ilat){
			for (int ilon=0; ilon < xlons.size(); ++ilon){
				//float f = v.getValue(xlons[ilon], xlats[ilat], ilev);
				float f = bilinear(ilat, ilon, ilev, indices,
								   v.lons, v.lats, 
								   xlons, xlats,
								   v.values, v.missing_value);
				if (f != v.missing_value) out(ilon,ilat,ilev) = f;
				else out(ilon,ilat,ilev) = out.missing_value;
			}
		}
	}
//	out.t = v.t; // copy the current time value because it is related to values
}


// interpolate variable v onto grid {xlons,xlats} (i.e. with metadata from v)
gVar lterp(gVar &v, vector <float> &xlons, vector <float> &xlats){
	gVar temp; temp.shallowCopy(v);
	temp.lats = xlats; temp.nlats = xlats.size();
	temp.lons = xlons; temp.nlons = xlons.size();
	temp.values.resize(temp.nlons*temp.nlats*temp.nlevs);
	for (int ilev=0; ilev < temp.nlevs; ++ilev){
		for (int ilat=0; ilat < temp.nlats; ++ilat){
			for (int ilon=0; ilon < temp.nlons; ++ilon){
				temp(ilon,ilat,ilev) = v.getValue(xlons[ilon], xlats[ilat], ilev);
			}
		}
	}
	//temp.values = tempvalues;
	return temp;
}


