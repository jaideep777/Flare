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

#include <cmath>
#include <algorithm>
#include <sstream>
#include "../include/gvar.h"
#include "../include/ncio.h"
#include "../include/constants.h"
#include "../include/time_math.h"
#include "../include/arrayutils.h"
using namespace netCDF;
using namespace netCDF::exceptions;

/******************       class NcFile_handle      ***********************/

NcFile_handle::NcFile_handle(){
	levname = "lev"; levunits = "z";
	nlons = nlats = ntimes = 0; nlevs = 1;
	nvars = ncoords = 0;
	ilon0 = ilat0 = 0;
	wlonix = elonix = slatix = nlatix = 0;
	mplimited = false;
	dFile = NULL; // crucial - needed to check if dFile has been allocated
//	lonDim = latDim = levDim = tDim = NULL;
//	lonVar = latVar = levVar = tVar = NULL;
	firstVarID = -1;
}

// ------------------------- INIT ----------------------------
int NcFile_handle::open(string s, string m, float glimits[4]){
	fname = s;
	mode = m;
	CINFO << "Open file: " << s; gsm_log->flush();
	if (mode == "r"){ 
		dFile = new NcFile(s.c_str(), NcFile::read);
	}
	else if (mode == "w"){ 
		dFile = new NcFile(s.c_str(), NcFile::replace);
	}

	if (!dFile)	{CERR << "Failed to open File: " << s << endl; return 1;}
	
	if (mode == "r"){
		if (glimits[0] > 180) glimits[0] -= 360;
		if (glimits[1] > 180) glimits[1] -= 360; 

		if (glimits[1] < glimits[0]){
			CERR << "Failed to open file: " << s << " - Incorrect Lon bounds: [" << glimits[0] << ", " << glimits[1] << "]" << endl;
			return 1;
		}

		mplimited = (glimits[0] > -180 || glimits[1] < 180 || glimits[2] > -90 || glimits [3] < 90)? true:false; 
		wlon = glimits[0]; 
		elon = glimits[1]; 
		slat = glimits[2]; 
		nlat = glimits[3]; 
	}
		
	if (!dFile->isNull()) {CINFOC << "... Success!" << endl; return 0;}
	else return 1;
}

int NcFile_handle::close(){
	dFile->close();
	delete dFile;
	dFile = NULL;
	return 0;
}

NcFile_handle::~NcFile_handle(){
	//dFile->close();	including this statement gives seg fault!
//	delete dFile;
//	dFile = NULL;
}

//void NcFile_handle::setMapLimits(float xwlon, float xelon, float xslat, float xnlat){
//	if (xwlon > 180) xwlon -= 360;
//	if (xelon > 180) xelon -= 360; 
//	mplimited = (xwlon > -180 || xelon < 180 || xslat > -90 || xnlat < 90)? true:false;
//	wlon = xwlon;
//	elon = xelon;
//	slat = xslat;
//	nlat = xnlat;
//}


// ------------------------- READING ----------------------------

int NcFile_handle::readTime(gVar &v){
	if (!tVar.isNull()){
		clock_t start = clock(), end;
		CINFO << "  reading t... ";
		
		if (tVar.getDimCount() > 1) CERR << "Time var has more than 1 dimensions." << endl;
		
		ntimes = v.ntimes = tVar.getDim(0).getSize();	// otherwise constructor has init to 0

		v.times.resize(v.ntimes);
		tVar.getVar(&v.times[0]);
		
		string datestr; tVar.getAtt("units").getValues(datestr);


		string unit, junk, sdate, stime;
		stringstream ss;
		ss.clear(); ss.str(datestr);
		ss >> unit >> junk >> sdate >> stime;
		if (stime == "") stime = "0:0:0";
		
		v.tbase = ymd2gday(sdate) + hms2xhrs(stime); // note this time is in GMT
		if (unit == "hours") v.tscale = 1.0f;
		else if (unit == "days") v.tscale = 24.0f;
		else if (unit == "months") {
			CWARN << "WARNING: using months as time units! 365.2524 days/yr will be considered.\n";
			v.tscale = (365.2524/12.0)*24.0;
		}
		else {
			CERR << "ERROR setting base time in getCoords(): invalid time units!\n";
			return 1;
		}
		v.tstep = (v.times[v.times.size()-1] - v.times[0])/(v.times.size()-1)*v.tscale;    // average tstep in hours
		CINFOC << v.ntimes << " read: (" 
			   << gt2string(v.ix2gt(0)) << " --- " << gt2string(v.ix2gt(v.ntimes-1)) << ").";
		end = clock();
		CINFOC << " [" << double(end-start)/CLOCKS_PER_SEC*1000 << " ms]" << endl; 
	}
	// NOTE: For some reason, reading the time vector from .nc is slow. Takes ~1 ms. Rest are 0.01 ms

}


int NcFile_handle::getMeta(){
	clock_t start = clock(), end;

//	CINFO << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	// get number of variables
	CINFO << "> Reading metadata from file: " << fname << "\n";
	nvars = 0;
	nvars = dFile->getVarCount();
	CINFO << "    " << nvars << " variables found.\n";

	multimap<string,NcVar> vars_map = dFile->getVars();

	// Get the coordinate variables
	ncoords = 0;
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it ){
		
		string var = it->first;
		
		if (var == "lat" || var == "latitude" || var == "LAT" ){
			latVar = it->second;
			++ncoords;
		}

		else if (var == "lon" || var == "longitude" || var == "LON" ){
			lonVar = it->second;
			++ncoords;
		}

		else if (var == "lev" || var == "levels" || var == "LEV" ){
			levVar = it->second;
			++ncoords;
		}

		else if (var == "time" || var == "TIME" || var == "LEV" ){
			tVar = it->second;
			++ncoords;
		}

	}
	
	if (latVar.isNull()) CWARN << "Lat not found\n";
	if (lonVar.isNull()) CWARN << "Lon not found\n";
//	if (levVar.isNull()) CWARN << "Lev not found\n";
	if (tVar.isNull())   CWARN << "Time not found\n";

	CINFO << "    " << ncoords << " coordinates found" << endl;
		
	CINFO << "> Attempting to find 1st variable... ";
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it ){
		if (it->second.getDimCount() == ncoords) {
			firstVar = it->second; 
			break;
		}
	}
	firstVarID = firstVar.getId();
	CINFOC << firstVar.getName() << " (" << firstVarID << ")";
	end = clock();
	CINFOC << " [" << double(end-start)/CLOCKS_PER_SEC*1000 << " ms]" << endl; 

}

// THE MOST PAINFUL NCIO FUNCTION!!
// this function:
//	 reads most of the metadata - read # coords, coord values
//	 sets the lat lon limits
//	 sets the correct order of lats and lons according to model requirements
int NcFile_handle::readCoordData(gVar &v){
	 
	// read lats, lons metadata from file
	// read actual values only if rr is true

	// if lat order is found N-S in file, reverse it. 
	// in this model lat order will be S-N (i.e lats increase with index)
	// all data reading will depend crucially on the variable latSN
	clock_t start = clock(), end;
	CINFO << "> Reading coordinates.\n";
	if (!latVar.isNull()){
		CINFO << "  reading lats... ";
		nlats = v.nlats = latVar.getDim(0).getSize();	// otherwise constructor has init to 0
		v.lats.resize(v.nlats);
		latVar.getVar(&v.lats[0]);
		CINFOC << v.nlats << " read.\n";

		float dlat = v.lats[1] - v.lats[0];
		latSN = (dlat < 0)? false:true;
		
		if (mplimited){
			CINFO << "    trim lats... ";
			// set indices
			slatix = ncIndexLo(v.lats, slat);
			nlatix = ncIndexHi(v.lats, nlat);
			ilat0 = (nlatix > slatix)? slatix:nlatix;
			ilatf = (nlatix > slatix)? nlatix:slatix;
			// cut array
			//v.lats = copyArray(v.lats, nlatix, slatix);
			v.lats = vector<float> (v.lats.begin()+ilat0, v.lats.begin()+ilatf+1);
			v.nlats = fabs(nlatix - slatix) +1;
			CINFOC << v.nlats << " left: (" << v.lats[0] << " --- " << v.lats[v.nlats-1] << ").\n";
		}

		if (!latSN){
			CINFO << "    lats array is N-S. Reversing...";
			reverseArray(v.lats);
			CINFOC << " Done.\n";
		}
	}

	if (!lonVar.isNull()){
		CINFO << "  reading lons... ";
		nlons = v.nlons = lonVar.getDim(0).getSize();	// otherwise constructor has init to 0
		v.lons.resize(v.nlons);
		lonVar.getVar(&v.lons[0]);
		CINFOC << v.nlons << " read.\n";
		
		vector<float>& lons = v.lons;

		auto it = find_if(lons.begin(), lons.end(), [](float x){return x > 180;});	// find first lon > 180
		for_each(it, lons.end(), [](float &x){x -= 360;});							// bring all lons > 180 t0 principle range
		int shift = lons.size() - distance(lons.begin(), it);						// array to be shifted right by (n-it) elements

		CINFO << "    shift lons right by: " << shift << endl;

		if (it != lons.begin() && it != lons.end()){								// if shift is > 0 and < N, 
			shiftRight(lons.data(), lons.size(), shift);							//    shift lons 
		}

		int wlonix1 = ncIndexLo(lons, wlon);
		int elonix1 = ncIndexHi(lons, elon);

		wlonix = (wlonix1-shift+lons.size()) % lons.size();
		elonix = (elonix1-shift+lons.size()) % lons.size();

		splitRead = (wlonix > elonix);			// whether data will have to be read in 2 passes
		
		if (mplimited){
			CINFO << "    trim lons... ";

			v.lons = vector<float>(&v.lons[wlonix1], &v.lons[elonix1]+1);
			v.nlons = v.lons.size();
			CINFOC << v.nlons << " left: (" << v.lons[0] << " --- " << v.lons[v.nlons-1] << ").\n";
		}
		
		
	}

	if (!levVar.isNull()){
		CINFO << "  reading levs... ";
		nlevs = v.nlevs = levVar.getDim(0).getSize();	// otherwise constructor has init to 1
		v.levs.resize(v.nlevs);
		levVar.getVar(&v.levs[0]);
		CINFOC << v.nlevs << " read.\n";
	}

	// allocate space for values
	v.values.resize(v.nlevs*v.nlats*v.nlons); // no need to fill values 

	end = clock();
	CINFO << "  = [" << double(end-start)/CLOCKS_PER_SEC*1000 << " ms]" << endl; 
	return 0;
}


int NcFile_handle::readCoords(gVar &v){
	getMeta();
	v.ivar1 = firstVarID;
	readCoordData(v);
	readTime(v);
	CINFO << "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
	return 0;
}

//int NcFile_handle::getVarID(string varname){
//	NcVar * v;
//	for (int i=0; i< nvars; ++i){
//		v = dFile->getVar(i);
//		if (v->name() == varname) return i;
//	}
//	CWARN << "No variable named " << varname << " found.\n";
//}

// if file has only 1 variable, then ivar = ncoords since [0,ncoords-1] are coords
int NcFile_handle::readVarAtts(gVar &v, string vname){
	NcVar nVar;
	if (vname == "") nVar = firstVar;	// ivar not specified, default is -1, so set it to ivar1
	else  nVar = dFile->getVar(vname);

	string s; 
	float f;
	// name
	v.varname = nVar.getName();
	
	// units
	try{
		firstVar.getAtt("units").getValues(v.varunits);
	}
	catch(NcException &e){
		v.varunits = "NA";
	}
	
	// missing value
	try{
		firstVar.getAtt("missing_value").getValues(&v.missing_value);;
	}
	catch(NcException &e){
		try{
			firstVar.getAtt("_FillValue").getValues(&v.missing_value);
		}
		catch(NcException &e){
		}
	}

//	// scale factor
	try{
		firstVar.getAtt("scale_factor").getValues(&v.scale_factor);
	}
	catch(NcException &e){
		v.scale_factor = 1.f;
	}
	
	// offset
	try{
		firstVar.getAtt("add_offset").getValues(&v.add_offset);
	}
	catch(NcException &e){
		v.add_offset = 0.f;
	}

	return 0;
}


int NcFile_handle::read_data_block(NcVar* vVar, int ilon0, int ilat0, int ilev0, int _nlons, int _nlats, int _nlevs, int itime, vector<float>&values){

	vector <size_t> start, count;
	if (!tVar.isNull()){
		start.push_back(itime);
		count.push_back(1);
	}
	if (!levVar.isNull()){
		start.push_back(ilev0);
		count.push_back(_nlevs);
	}
	if (!latVar.isNull()){
		start.push_back(ilat0);
		count.push_back(_nlats);
	}
	if (!lonVar.isNull()){
		start.push_back(ilon0);
		count.push_back(_nlons);
	}

	vVar->getVar(start, count, values.data());

//	// read data. cur is set at the SW corner of grid / grid-limits.
//	if (vVar->num_dims() == 4){ 
//		vVar->set_cur(itime, ilev0, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner
//		vVar->get(&values[0], 1, _nlevs, _nlats, _nlons);
//	}	
//	else if (vVar->num_dims() == 3){ 
//		vVar->set_cur(itime, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner
//		vVar->get(&values[0], 1, _nlats, _nlons);
//	}
//	else if (vVar->num_dims() == 2){
//		vVar->set_cur(ilat0, ilon0);
//		vVar->get(&values[0], _nlats, _nlons);
//		CWARN << "(" << "" << ") treating 2D Variable as lat-lon map..\n";
//	}
//	else{
//		CERR << "Variables with only 2/3/4 dimensions are supported! Dims found: " 
//			 << vVar->num_dims() << "\n";
//		return 1;
//	}
	
	return 0;

}


// if file has only 1 variable, then ivar = ivar1
// ivar need to be specified only in case the file has multiple variables
int NcFile_handle::readVar(gVar &v, int itime, string vname){
	if (mode != "r"){
		CERR << "ERROR in readVar: File not in read mode.\n";
		return 1;
	}

	CDEBUG << "NcfileHandle::readVar(" << v.varname << ", t=" << itime <<  "): ";
	if (v.times.size() > 0) CDEBUGC << gt2string(v.ix2gt(itime));
	else CDEBUGC << "2D map" ;
	
	// file is in read mode.. continue.
	NcVar vVar;
	if (vname == "") vVar = firstVar;
	else vVar = dFile->getVar(vname);
	
	clock_t start = clock(), end;
	
	if (!splitRead){
		int nlons_trim = elonix-wlonix+1;
		v.values.resize(nlons_trim*v.nlats*v.nlevs);
		read_data_block(&vVar, wlonix, ilat0, 0, nlons_trim, v.nlats, v.nlevs, itime, v.values);
		if (v.nlons != nlons_trim) CERR << "Mismatch in calculation of nlons after trimming: " << v.nlons << " / " << nlons_trim << endl;
	}

	else {
		CDEBUGC << " in 2 passes: [0, " << elonix << "], [" << wlonix << ", " << nlons-1 << "] ";   
		int nlons1 = elonix-0+1;
		vector <float> seg1(nlons1*v.nlats*v.nlevs); 
		read_data_block(&vVar, 0, ilat0, 0, nlons1, v.nlats, v.nlevs, itime, seg1);

		int nlons2 = nlons-1-wlonix+1;
		vector <float> seg2(nlons2*v.nlats*v.nlevs); 
		read_data_block(&vVar, wlonix, ilat0, 0, nlons2, v.nlats, v.nlevs, itime, seg2);
		
		int nlons_trim = nlons1+nlons2;
		v.values.resize(nlons_trim*v.nlats*v.nlevs);
		
		for (int ilev=0; ilev < v.nlevs; ++ilev){
			for (int ilat=0; ilat < v.nlats; ++ilat){
				for (int ilon=0; ilon < nlons2; ++ilon){
					v.values[IX3(ilon,ilat,ilev, nlons_trim, v.nlats)] = seg2[IX3(ilon,ilat,ilev, nlons2, v.nlats)];
				}
			}
		}
		for (int ilev=0; ilev < v.nlevs; ++ilev){
			for (int ilat=0; ilat < v.nlats; ++ilat){
				for (int ilon=0; ilon < nlons1; ++ilon){
					v.values[IX3(ilon+nlons2,ilat,ilev, nlons_trim, v.nlats)] = seg1[IX3(ilon,ilat,ilev, nlons1, v.nlats)];
				}
			}
		}
		

		if (v.nlons != nlons_trim) CERR << "Mismatch in calculation of nlons after trimming: " << v.nlons << " / " << nlons_trim << endl;
	}

	end = clock();
	CDEBUGC << " [" << v.nlevs*v.nlons*v.nlats*sizeof(float)/1e6 << " Mb / " << double(end-start)/CLOCKS_PER_SEC*1000 << " ms @ " << v.nlevs*v.nlons*v.nlats*sizeof(float)/(double(end-start)/CLOCKS_PER_SEC)*1e-9 << " Gb/s]"<< endl;
	
	// if lats are not in SN order, reverse the data along lats
	if (!latSN)	reverseCube(&v.values[0], v.nlons, v.nlats, v.nlevs);	// TODO: This reversing can be moved to gVar, only as last step of getting data

	// if either scale_factor or offset is present, convert data..
	if (v.scale_factor != 1 || v.add_offset != 0){
		for (int i=0; i< v.nlevs*v.nlats*v.nlons; ++i){
			if (v.values[i] != v.missing_value) 	// ignore missing values
				v.values[i] = v.values[i]*v.scale_factor + v.add_offset;	
		}
	}

	if (!tVar.isNull()){
		// set t to gday corresponding to times[itime]
		v.t = v.ix2gt(itime);
	}
	
	return 0;
}

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	EXTREMELY PAINFUL FUNCTION! NEEDS TIME VARIABLES AS DOUBLE!
	This version of readVar reads variable values corresponding to global time 
	gt into gVar v.
	If exact time is not found in v's time values, then 
	the lower and upper bounds are read and the values calculated from the values
	corresponding to the bounds. (depending on mode)
	mode: 0 = hold, 1 = interpolate
	return values:
		0 = did not read from file because var was already up-to-date
		1 = read from tixlo
		2 = read from tixlo +1
		3 = interpolated
	   -1 = error		
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/	
int NcFile_handle::readVar_gt(gVar &v, double gt, int mode, string vname){
	// calculate index corresponding to gt


	int tixlo = v.gt2ix(gt);		// index for curret gt
	int tixlo_prev = v.gt2ix(v.t);	// index for gt last read into v

	CDEBUG << "readVar_gt(" << v.varname << "): " << gt2string(v.t) << " [ix="<<tixlo_prev << "] --> " << gt2string(gt) << " [ix=";
	CDEBUGC << tixlo << " = ";
	CDEBUGC << gt2string(v.ix2gt(tixlo)) << "] -> ";

	// in mode 0, if both these indices match, data is already in variable, so skip reading
	if (tixlo == tixlo_prev && mode == 0) { 
		CDEBUGC << "~ up to date. Not reading.\n";
		v.t = gt;	// set t to the actual gt and not to any index values
		return 0;	
	}
	// if times[index] = gt, then no need of interpolation
	// also, if mode = hold, values should be from tixlo anyway, so include that case here
	else if ( gt-v.ix2gt(tixlo) < t_tol || mode == 0){
		CDEBUGC << "~ Reading low, diff = " << (gt-v.ix2gt(tixlo))*24.0 << " hrs.\n";
		if (tixlo >= v.ntimes) CERR << "(" << v.varname << "): desired time out of bounds";
		readVar(v, tixlo, vname);
		v.t = gt;	// set t to the actual gt and not to any index values
		return 1;
	}
	else if ( v.ix2gt(tixlo+1)-gt < t_tol || mode == 0){	
		// we dont want to hold/interpolate for gt that is too close to the higher bound.. hence this case
		CDEBUGC << "~ Reading high, diff = " << (v.ix2gt(tixlo+1)-gt)*24.0 << " hrs.\n";
		readVar(v, tixlo+1, vname);
		v.t = gt;	// set t to the actual gt and not to any index values
		return 2;
	}
	else{	// read the values at next step and calculate.
	//TODO What happens in this mode if value is between 2 files?!
		CDEBUGC << "~ Interpolating, intermediate, diff = " << (gt-v.ix2gt(tixlo))*24.0 << " hrs.\n";
		gVar z; z.copyMeta(v);
		gVar w; w.copyMeta(v);	// make a new gVar to hold next step values, identical to v
		readVar(z, tixlo, vname);	// read this step in v
		readVar(w, tixlo+1, vname);	// read next step in w
		double t1 = z.ix2gt(tixlo), t2 = w.ix2gt(tixlo+1);
		if (fabs((t2-t1)*24.0 - v.tstep)*3600 > t_tol){	 // compare difference in seconds to t_tol (86 sec)
			CWARN << "WARNING in readVar_gt: calculated timestep does not match actual..\n";
			CWARN << "..difference is = " << fabs((t2-t1)*24.0 - z.tstep)*3600 << " sec, > 86 sec.\n";
		}	// just a check
		z = z + (w - z)*(gt-t1)/(t2-t1);	// linear interpolate
		v.copyValues(z);
		v.t = gt;	// set t to the actual gt and not to any index values
		return 3;
	}
	

	return -1;
}

// ------------------------- WRITING ----------------------------

int NcFile_handle::writeCoords(gVar &v, bool wr){
	if (mode != "w"){
		CERR << "ERROR in writeCoords: File not in write mode.\n";
		return 1;
	}
	// for each coord, add a dimension, variable, and add units and attributes
	// write coord values only if "wr" is true
	if (v.nlats > 0){
		latDim = dFile->addDim("lat", v.nlats);	// dim
		latVar = dFile->addVar("lat", ncFloat, latDim); // var
		latVar.putAtt("units", "degrees_north");	// attributes
//		latVar->putAtt("scale_factor", NcFloat, 1.0f);
//		latVar->putAtt("add_offset", NcFloat, 0.0f);
		if (wr) latVar.putVar(&v.lats[0]);	// write data
	}
	//CINFO << "Wrote lats to file." << endl;
	
	if (v.nlons > 0){
		lonDim = dFile->addDim("lon", v.nlons);
		lonVar = dFile->addVar("lon", ncFloat, lonDim);
		lonVar.putAtt("units", "degrees_east");
//		lonVar->add_att("scale_factor", 1.0f);
//		lonVar->add_att("add_offset", 0.0f);
		if (wr) lonVar.putVar(&v.lons[0]);
	}
	//CINFO << "Wrote lons to file." << endl;
	
	if (v.nlevs > 1){
		levDim = dFile->addDim(levname.c_str(), v.nlevs);
		levVar = dFile->addVar(levname.c_str(), ncFloat, levDim);
		levVar.putAtt("units", levunits.c_str());
//		levVar->add_att("scale_factor", 1.0f);
//		levVar->add_att("add_offset", 0.0f);
		if (wr) levVar.putVar(&v.levs[0]);
	}
	//CINFO << "Wrote levs to file." << endl;
	
	if (v.ntimes > 0){
		
		tDim = dFile->addDim("time");	// unlimited dimension
		tVar = dFile->addVar("time", ncFloat, tDim);
		string tunits = "";
		if (v.tscale == 1) tunits = "hours since " + gt2string(v.tbase);
		else if (v.tscale == 24) tunits = "days since " + gt2string(v.tbase);
		else if (fabs(v.tscale - (365.2524/12.0)*24.0) < 1){	// allow tolerance of 1 hr for hrs/month
			tunits = "months since " + gt2string(v.tbase);
			CWARN << "Using months as time units! 365.2524 days/yr are meant.\n";
		}
		else {
			CERR << "ERROR in writeCoords(): invalid time units! tscale = " << v.tscale << "\n";
		}
		// gday2ymd(int(v.tbase)) + " " + xhrs2hms(v.tbase - int(v.tbase));
		tVar.putAtt("units", tunits.c_str());
//		tVar->putAtt("scale_factor", 1.0f);
//		tVar->add_att("add_offset", 0.0f);
		// time vector not written now because ntimes may be determined later
	}
	CINFO << "> Written coordinates to file." << endl;
	return 0;
}	

// time coord is written seperately because several times number of times
// is determined later 
int NcFile_handle::writeTimeValues(gVar &v){
	tVar.putVar(&v.times[0]);
	CINFO << "> Written time vector to file." << endl;
}

NcVar NcFile_handle::createVar(gVar &v){
	if (mode != "w"){
		CERR << "ERROR in createVar: File not in write mode.\n";
		return NcVar();
	}
	if (v.varname == "" || v.varunits == ""){
		CERR << "ERROR in createVar: variable Name or Units not set.\n";
		return NcVar();
	}
	// file is in write mode.. continue.
	// include level dimension only if nlevs > 0.
	
	vector <NcDim> dims;
	if (!tDim.isNull())    dims.push_back(tDim);
	if (!levDim.isNull())  dims.push_back(levDim);
	if (!latDim.isNull())  dims.push_back(latDim);
	if (!lonDim.isNull())  dims.push_back(lonDim);
	
	
	NcVar vVar = dFile->addVar(v.varname, ncFloat, dims);
		
	// add variable attributes
	vVar.putAtt("units", v.varunits.c_str());
	vVar.putAtt("scale_factor", ncFloat, 1.0f);
	vVar.putAtt("add_offset", ncFloat, 0.0f);
	vVar.putAtt("_FillValue", ncFloat, v.missing_value);
	vVar.putAtt("missing_value", ncFloat, v.missing_value);
	
	// actual data not written here because ..
	// we need to write data into already created variable
	//vVar->put_rec(v.values, itime);
	CINFO << "> Successfully created variable in file." << endl;
	
	return vVar;
}

int NcFile_handle::writeVar(gVar &v, NcVar vVar, int itime){
	if (mode != "w"){
		CERR << "ERROR in writeVar: File not in write mode.\n";
		return 1;
	}
	clock_t startTime = clock(), end;
	CDEBUG << "NcFile_handle::writeVar (" << v.varname << ", t=" << itime << "): "; gsm_log->flush();

	vector <size_t> start, count;
	if (!tDim.isNull()){
		start.push_back(itime);
		count.push_back(1);
	}
	if (!levDim.isNull()){
		start.push_back(0);
		count.push_back(v.nlevs);
	}
	if (!latDim.isNull()){
		start.push_back(0);
		count.push_back(v.nlats);
	}
	if (!lonDim.isNull()){
		start.push_back(0);
		count.push_back(v.nlons);
	}

	vVar.putVar(start, count, v.values.data());

	end = clock();
	CDEBUGC << " [" << v.nlevs*v.nlons*v.nlats*sizeof(float)/1e6 << " Mb / " << double(end-startTime)/CLOCKS_PER_SEC*1000 << " ms @ " << v.nlevs*v.nlons*v.nlats*sizeof(float)/(double(end-startTime)/CLOCKS_PER_SEC)*1e-9 << " Gb/s]"<< endl;

	return 0;
}



