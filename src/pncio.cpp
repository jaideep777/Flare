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

#include <math.h>
#include "../include/gsm.h"


vector <int> id2ix3(int id, int nx, int ny, int nz){
	// id = iz*ny*nx + iy*nx + ix; 
	vector <int> ids(3);
	int ix = id % nx;
	int iy = (id % (ny*nx) - ix)/nx;
	int iz = (id - iy*nx - ix)/nx/ny;
	ids[0] = ix;
	ids[1] = iy;
	ids[2] = iz;
	return ids;
}

// if file has only 1 variable, then ivar = ivar1
// ivar need to be specified only in case the file has multiple variables
int NcFile_handle::readVar_parallel(gVar &v, int itime, int iVar){
	if (mode != "r"){
		CERR << "ERROR in readVar: File not in read mode.\n";
		return 1;
	}

	CDEBUG << "readVar_it(" << v.varname << "): ";
	if (v.times.size() > 0) CDEBUGC << gt2string(v.ix2gt(itime)) << endl;
	else CDEBUGC << "2D map" << endl;
	
	// file is in read mode.. continue.
	if (iVar == -1) iVar = v.ivar1;

	// allocate space for values
	v.values.resize(v.nlevs*v.nlats*v.nlons); // no need to fill values 
	NcVar * vVar = dFile->get_var(iVar);
	
	// read data. cur is set at the SW corner of grid / grid-limits.
//	int nvals = v.nlevs*v.nlats*v.nlons;
//	int s1 = 0, e1 = nvals/2, s2 = nvals/2+1, e2 = nvals-1;
//	
//	vector <int> c1 = id2ix3(s1, v.nlons, v.nlats, v.nlevs);
//	vector <int> c2 = id2ix3(s2, v.nlons, v.nlats, v.nlevs);
	
	
	if (vVar->num_dims() == 4){ 
		vVar->set_cur(itime, 0, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner
		vVar->get(&v.values[0], 1, v.nlevs, v.nlats, v.nlons);
	}	
	else if (vVar->num_dims() == 3){
		cout << "Reading 3D file" << endl; 
		vVar->set_cur(itime, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner		
		vVar->get(&v.values[0], 1, v.nlats/2, v.nlons);
		cout << "ilat0-ilatf = " << ilat0 << ":" << ilatf << endl;
//		vVar->set_cur(itime, ilatf/2+1, ilon0); // set the starting time at itime and lat/lon/lev at SW corner		
//		vVar->get(&v.values[(ilatf/2+1)*v.nlons + ilon0], 1, v.nlats - v.nlats/2, v.nlons);
	}
	else if (vVar->num_dims() == 2){
		vVar->set_cur(ilat0, ilon0);
		vVar->get(&v.values[0], v.nlats, v.nlons);
		CWARN << "(" << v.varname << ") treating 2D Variable as lat-lon map..\n";
	}
	else{
		CERR << "Variables with only 2/3/4 dimensions are supported! Dims found: " 
			 << vVar->num_dims() << "\n";
		return 1;
	}
	
	// if lats are not in SN order, reverse the data along lats
	if (!latSN)	reverseCube(v.values, v.nlons, v.nlats, v.nlevs);

	// if either scale_factor or offset is present, convert data..
	if (v.scale_factor != 1 || v.add_offset != 0){
		for (int i=0; i< v.nlevs*v.nlats*v.nlons; ++i){
			if (v.values[i] != v.missing_value) 	// ignore missing values
				v.values[i] = v.values[i]*v.scale_factor + v.add_offset;	
		}
	}

	if (tVar){
		// set t to gday corresponding to times[itime]
		v.t = v.ix2gt(itime);
	}
	
	return 0;
}

