/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    libGSM - library for Gridded Spatial Modelling
    Copyright (C) 2016, 2018 Jaideep Joshi

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

#include "../include/gsm.h"
#include <cmath>
#include <algorithm>
#include <ctime>

gVar::gVar(){
	ntimes = 0; nlevs = 1; nlats = 0; nlons = 0;
	tbase = 0; tscale = 1; tstep = 1;
	varname = ""; varunits = "";
	missing_value = std_missing_value;
	ipvar = NULL;			// set all 4 pointers in class to Null
	ifile_handle = NULL;	// they will be allocated by filo IO functions
	ofile_handle = NULL;
	outNcVar = NULL;
	ivar1=-1;
	t=-1;
}

gVar::gVar(string name, string units, string tunits){
	ntimes = 0; nlevs = 1; nlats = 0; nlons = 0;
	tstep = t = 0;
	varname = name; varunits = units;
	missing_value = std_missing_value;
	ipvar = NULL;			// set all 4 pointers in class to Null
	ifile_handle = NULL;
	ofile_handle = NULL;
	outNcVar = NULL;
	ivar1=-1;
	t = -1;
	
	// read time unit string and set tbase and tscale
	string unit, junk, sdate, stime;
	stringstream ss;
	ss.clear(); ss.str(tunits);
	ss >> unit >> junk >> sdate >> stime;
	if (stime == "") stime = "0:0:0";
	
	tbase = ymd2gday(sdate) + hms2xhrs(stime); // note this time is in GMT
	if (unit == "hours") tscale = 1.0f;
	else if (unit == "days") tscale = 24.0f;
	else if (unit == "months") {
		CWARN << "Using months as time units! 365.2524 days/yr will be considered.\n";
		tscale = (365.2524/12.0)*24.0;
	}
	else {
		CERR << "ERROR setting base time in getCoords(): invalid time units!\n";
	}
}

// copy all data except values vector into self.
int gVar:: copyMeta(const gVar &v){
	ntimes = v.ntimes; nlevs = v.nlevs; nlats = v.nlats; nlons = v.nlons;
	times = v.times; levs = v.levs; lats = v.lats; lons = v.lons;
	tbase = v.tbase; tscale = v.tscale; tstep = v.tstep; //t = v.t;
	varname = v.varname; varunits = v.varunits;
	scale_factor = v.scale_factor; add_offset = v.add_offset;
	ncoords = v.ncoords; ivar1 = v.ivar1;
	missing_value = v.missing_value;
	gridlimits = v.gridlimits;
	
	values.resize(nlons*nlats*nlevs);
	
	t = v.t;
	return 0;
}

int gVar::initMetaFromFile(string filename){
	float glimits_globe[4] = {0, 360, -90, 90};
	ifile_handle = new NcFile_handle;
	int i = ifile_handle->open(filename, "r", glimits_globe);
	ifile_handle->readCoords(*this);
	ifile_handle->readVarAtts(*this);
	ifile_handle->close();
	delete ifile_handle; ifile_handle = NULL;
}


//// copy nciostream data 
//int gVar:: copyStreams(const gVar &v){
//	ifname = v.ifname; ofname = v.ofname;	
//	ifile_handle = v.ifile_handle; ofile_handle = v.ofile_handle;	
//	outNcVar = v.outNcVar;		
//	lwrite = v.lwrite; lwriteSP = v.lwriteSP;	
//	filenames = v.filenames;
//	lterp_indices = v.lterp_indices;
//	curr_file = v.curr_file;
//	ipvar = v.ipvar;
//}


// copy values vector and missing_value.
int gVar:: copyValues(const gVar &v){
	values = v.values;
	missing_value = v.missing_value;
	return 0;
}

int gVar::setCoords(vector <double> &t, vector <float> &le, vector <float> &la, vector <float> &lo){
	times = t; levs = le; lats = la; lons = lo;
	ntimes = t.size(); nlevs = le.size(); nlats = la.size(); nlons = lo.size();
	if (ntimes >2) tstep = (times[1] - times[0])*tscale;	// tstep in hours
	else tstep = 0;
	values.resize(nlons*nlats*nlevs);
}

int gVar::setTimeAtts(int xntimes, double xtbase, float xtscale){
	ntimes = xntimes;
	tbase = xtbase; 
	tscale = xtscale;	
}

int gVar::printGrid(ostream &lfout){
	lfout << "-------------------------------------------------------------------------\n";	
	lfout << "> Variable " << ivar1 << ": " << varname << " (" << varunits << ")\n";
	// grid info
	lfout << "> Grid:\n";
	lfout << "\t" << nlons << " lons: "; 
	if (nlons >0) lfout << lons[0] << " ... " << (lons[nlons-1]-lons[0])/(nlons-1) << " ... " << lons[nlons-1]; 
	lfout << "\n";
	lfout << "\t"  << nlats << " lats: ";
	if (nlats >0) lfout << lats[0] << " ... " << (lats[nlats-1]-lats[0])/(nlats-1) << " ... " << lats[nlats-1];
	lfout << "\n";
	lfout << "\t"  << nlevs << " levs.\n";
	lfout << "\t"  << ntimes << " times.\n";
	lfout << "\t" << "Current number of values: " << values.size() << '\n';
	lfout << "> Time:\n";
	lfout << "\tbase  gday (GMT) = " << gtstr6d(tbase) << ", i.e. " << gt2string(tbase) << '\n';
	if (ntimes > 0){
		double t0 = ix2gt(0);
		lfout << "\tfirst gday (GMT) = " << gtstr6d(t0) << ", i.e. " << gt2string(t0) << '\n';
		t0 = ix2gt(times.size()-1);
		lfout << "\tlast  gday (GMT) = " << gtstr6d(t0) << ", i.e. " << gt2string(t0) << '\n';
	}
	lfout << "\ttime step (tstep) = " << tstep << " hours.\n";
	lfout << "\thours per time unit (tscale) = " << tscale << "\n";
	lfout << "> Missing value = " << missing_value << "\n";
	lfout << "> Streams:\n";
	lfout << "\tInput: "; 
	if (ifile_handle != NULL) lfout << ifile_handle->dFile << " (ipvar = " << ipvar << ")\n";
	else  lfout << "--- (ipvar = " << ipvar << ")\n";
	lfout << "\tOutput: "; 
	if (ofile_handle != NULL) lfout << ofile_handle->dFile << "\n";
	else  lfout << "---\n";
	lfout << "-------------------------------------------------------------------------\n\n";
	lfout.flush();
	return 0;
}

int gVar::printGridIP(ostream &lfout){
	lfout << "Variable " << varname << " inputs from:\n";
	ipvar->printGrid(lfout);
}

int gVar::printValues(ostream &lfout){
	lfout << "> Values: (" << values.size() << ")\n";
	printCube(values, nlons, nlats, nlevs, missing_value);
}

// returns the index corresponding time just <= gt
int gVar::gt2ix(double gt){
//    CDEBUG << "t0 = " << times[0] << ", tbase = " << tbase << ", tscale = " << tscale << ", tstep = " << tstep << endl;
	return floor(((gt - tbase)*24.0 - times[0]*tscale)/tstep);	// essential to use floor. int() truncates towards 0! $%@#^$@#*@   
}

double gVar::ix2gt(int ix){	
	return tbase + times[ix]*tscale/24.0;
}

double gVar::ix2gt_IST(int ix){
	return tbase + 5.5/24.0 + times[ix]*tscale/24.0;	// middle term converts to IST
}


// *** Return interploated value of the variable at given coordinates  ***
float gVar::getValue(float xlon, float xlat, float ilev){
	return bilinear(xlon, xlat, ilev, lons, lats, values, missing_value);
}

// *** Return cell value of the variable at given coordinates  ***
float gVar::getCellValue(float xlon, float xlat, float ilev){
	return cellVal(xlon, xlat, ilev, lons, lats, values, missing_value);
}


// -------------------- FUNCTIONS ON SELF -------------------------

int gVar::fill(float f){
	for (int i=0; i<values.size(); ++i){
		if (values[i] != missing_value)
			values[i] = f;
	}
}

int gVar::sqrtVar(){
	for (int i=0; i<values.size(); ++i){
		if (values[i] != missing_value)
			values[i] = sqrt(values[i]);
	}
}

//int gVar::levMax(){
//	vector <float> vmax(nlons*nlats, -1e20);
//	
//	for (int ilev=2; ilev<veg.nlevs; ++ilev){
//		for (int ilat=0; ilat<veg.nlats; ++ilat){
//			for (int ilon=0; ilon<veg.nlons; ++ilon){
//				if (veg.values[ID(ilon, ilat, ilev)] > -1){
//					ft.values[ID(ilon, ilat, ilev)] = 1;
//				}
//			}
//		}
//	}
//}

// ------------------------- OPERATORS ----------------------------

// the element reference operators. They return a REFERENCE of the element desired
// return element by coord-indices  
float& gVar::operator () (int ilon, int ilat, int ilev){
	return values[nlons*nlats*ilev + nlons*ilat + ilon];
}

// return element by memory index
float& gVar::operator [] (int i){
	return values[i];
}


// all operators replace variable missing values with std missing value
// all operators operate only at places where both operands are non-missing
// 

gVar gVar::operator + (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		CERR << varname << " + " << v.varname << " : Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] + v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

gVar gVar::operator + (const float x){
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value)
			temp.values[i] = values[i] + x;
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};


gVar gVar::operator - (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		CERR << varname << " - " << v.varname << " : Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] - v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};


gVar gVar::operator - (const float x){
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value)
			temp.values[i] = values[i] - x;
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};


gVar gVar::operator * (const float x){
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value)
			temp.values[i] = values[i]*x;
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

gVar gVar::operator / (const float x){
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value)
			temp.values[i] = values[i]/x;
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

gVar gVar::operator * (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		CERR << varname << " * " << v.varname << " : Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] * v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};


gVar gVar::operator / (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		CERR << varname << " * " << v.varname << " : Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.copyMeta(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] / v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

 

// in all input stream functions, the fileIO members of ipvar are not touched 
// (in fact that cant be touched because they are private)
int gVar::createNcInputStream(vector <string> files, vector <float> glim){

	gridlimits = glim; //TODO should this be assigned to ipvar instead?
	if (glim[0] > lons[0] || glim[1] < lons[nlons-1] || glim[2] > lats[0] || glim[3] < lats[nlats-1]){
		CWARN << "Specified grid limits are narrower than the variable grid" << endl; 
	}
	
	filenames = files;	// check that at least 1 input file is specified
	if (filenames.size() < 1){
		CERR << "(" << varname << ") createInputStream: No input files specified\n"; 
		return 1;
	}
	curr_file = 0;
	
	ifile_handle = new NcFile_handle;	// create a file handle for input 
	ipvar = new gVar;					// allocate gVar for input
	
	loadInputFileMeta();	// read metadata from 1st input file
}

int gVar::loadInputFileMeta(){

	CDEBUG << "Attempting to load (" << varname << ") from file " << curr_file << ": " 	
		   << filenames[curr_file] << endl;

	// close any previously opened file
	if (ifile_handle->dFile != NULL) ifile_handle->close();

	// open new file
	int i = ifile_handle->open(filenames[curr_file], "r", &gridlimits[0]);
	if (i != 0) CERR << "NCFILE NOT VALID!!" << endl;

	// read metadata
	ifile_handle->readCoords(*ipvar);
	ifile_handle->readVarAtts(*ipvar);
	
	// calculate regridding indices
	lterp_indices = bilIndices(ipvar->lons, ipvar->lats, lons, lats);
	
}


int gVar::whichNextFile(double gt){
	double file_gtf = ipvar->ix2gt(ipvar->ntimes-1);	// first time
	double file_gt0 = ipvar->ix2gt(0);					// last time
	double file_dt = ipvar->tstep/24; 					// time  step in days
//	cout << varname << ": gt = " << gt2string(gt) << ", file limits = " << gt2string(file_gt0) << " --- " << gt2string(file_gtf+file_dt) << endl;
//	cout << "file_dt: " << file_dt << endl;
	if (gt >= file_gtf+file_dt) return curr_file+1;		// ---|---------||---------|======== <-- gt >= ix0
	else if (gt < file_gt0) return curr_file-1;			//   ixf   f1        f2   ix0            gt >= ixf+dt
	else return curr_file; 
}


int gVar::updateInputFile(double gt){
	
	int prev_file = curr_file;
	int next_file = whichNextFile(gt);

//	cout << prev_file << " " << curr_file << " " << next_file << endl;

	while (curr_file != next_file){

		if (next_file < 0 || next_file >= filenames.size()){
			CERR << "(" << varname << ") updateInputFile(): specified time out of range of given files (" << next_file << ").\n"; 
			return 1;
		}
		
		prev_file = curr_file;
		curr_file = next_file;
		loadInputFileMeta();
		next_file = whichNextFile(gt);
		
//		cout << prev_file << " " << curr_file << " " << next_file << endl;
		
		if (next_file == prev_file) {
			CERR << "(" << varname << ") updateInputFile(): Infinite loop (\n    " << filenames[curr_file] << "\n    " << filenames[next_file] 
				 << "\n    --------\n    " << filenames[min(curr_file, next_file)] << " - returned." << endl;
			curr_file = min(curr_file, next_file);
			break;
		}
	}

	return 0;
}



int gVar::closeNcInputStream(){
	if (ifile_handle->dFile != NULL) ifile_handle->close();
	delete ifile_handle; ifile_handle = NULL; 
	delete ipvar; ipvar = NULL;
}


int gVar::readVar_gt(double gt, int mode){
	if (ifile_handle == NULL){
		CERR << "gVar::readVar_gt(" << varname << "): NcInputStream not initialized" << endl;
		return 1;
	}
	 
	updateInputFile(gt);
	ifile_handle->readVar_gt(*ipvar, gt, mode, ipvar->ivar1);	// readCoords() would have set ivar1
	lterpCube(*ipvar, *this, lterp_indices);
	return 0;
}

int gVar::readVar_it(int tid){
	if (ifile_handle == NULL){
		CERR << "gVar::readVar_it(" << varname << "): NcInputStream not initialized" << endl;
		return 1;
	}

	ifile_handle->readVar(*ipvar, tid, ipvar->ivar1);
	lterpCube(*ipvar, *this, lterp_indices);
}

int gVar::createOneShot(string filename, vector<float> glim){
	ifname = filename;

	float glimits_globe[4] = {0, 360, -90, 90};
	if (glim.size() < 4) gridlimits = vector <float> (glimits_globe, glimits_globe+4);
	else gridlimits = glim;

	ifile_handle = new NcFile_handle;
	int i = ifile_handle->open(ifname, "r", &gridlimits[0]);
	ifile_handle->readCoords(*this);
	ifile_handle->readVarAtts(*this);
	ifile_handle->readVar(*this,0);
	ifile_handle->close();
	delete ifile_handle; ifile_handle = NULL;
}


int gVar::readOneShot(string filename, vector <float> glim){
	ipvar = new gVar();
	ipvar->createOneShot(filename, glim);
	lterp_indices = bilIndices(ipvar->lons, ipvar->lats, lons, lats);
	lterpCube(*ipvar, *this, lterp_indices);	// does not copy metadata
	delete ipvar; ipvar = NULL;
}




// output

int gVar::createNcOutputStream(string filename){
	ofname = filename;
	ofile_handle = new NcFile_handle;
	int i = ofile_handle->open(filename, "w", NULL);	// gridlimits are not required for writing
	ofile_handle->writeCoords(*this);
	if (ofile_handle->tVar) ofile_handle->writeTimeValues(*this);
	outNcVar = ofile_handle->createVar(*this);
	if (outNcVar->is_valid()) CDEBUG << "Succesfully created variable in file " << filename << endl;
}


int gVar::closeNcOutputStream(){
	ofile_handle->close();
	delete ofile_handle; ofile_handle = NULL;
}

int gVar::writeVar(int itime){
	return ofile_handle->writeVar(*this, outNcVar, itime);
}

int gVar::writeOneShot(string filename){
	createNcOutputStream(filename);
	writeVar(0);
	closeNcOutputStream();
}



// -----------------------------------------------------------------------
// Functions to read data and compute aggregates between specified times
// -----------------------------------------------------------------------
int gVar::readVar_reduce_mean(double gt1, double gt2){
	gVar temp; temp.copyMeta(*ipvar);
//	temp.values.resize(temp.nlons*temp.nlats*temp.nlevs);
	temp.fill(0);
	int count = 0;
	
	updateInputFile(gt1);	// this will give correct OR one previous file
	while (gt1 > ipvar->ix2gt(ipvar->times.size()-1)){	// increment curr_file as long as gt1 +outside file range
		++curr_file;
		loadInputFileMeta();
	}
	
	CDEBUG << "readVar_reduce_mean (" << varname << ") :" << gt2string(gt1) << " " << gt2string(gt2) << endl;
	while(1){
//		cout << (ipvar->times[0]) << " " << (ipvar->times[ipvar->times.size()-1]) << " "  << (gt1-ipvar->tbase)*24.0/ipvar->tscale <<  endl;
		int tstart = lower_bound(ipvar->times.begin(), ipvar->times.end(), (gt1-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin();	   // first elem >= gt1 
		int tend   = upper_bound(ipvar->times.begin(), ipvar->times.end(), (gt2-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin() -1;   // last elem <= gt2
//		cout << gt2string(ipvar->ix2gt(tstart)) << " " << gt2string(ipvar->ix2gt(tend)) << endl;

		if (tend < 0) break;

		for (int i=tstart; i<=tend; ++i){ 
			ifile_handle->readVar(*ipvar, i, ipvar->ivar1);	// readCoords() would have set ivar1
			temp = temp + *ipvar;
			++count;
		}

		if (tend >= ipvar->times.size()-1){ // if tend was the last time in file, then load next file and continue reading
			++curr_file;
			if (curr_file >= filenames.size()) break;
			loadInputFileMeta();
		}
		else break;
	}
	
	CDEBUG << "----------- Read " << count << " timesteps from " << varname << endl;
	temp = temp/count;	
	if (count > 0) lterpCube(temp, *this, lterp_indices);	// We want to preserve current values if no new values were read
}



gVar gVar::trend(double gt1, double gt2){

	
	gVar temp; temp.copyMeta(*ipvar);
	temp.fill(0);
	int count = 0;
	
	gVar Sxx = temp, 
		 Syy = temp, 
		 Sxy = temp;

	gVar sxy = temp, 
		 sx = temp, 
		 sy = temp, 
		 sxx = temp, 
		 syy = temp;

	gVar b1 = temp, 
		 s = temp,
		 t = temp; 
	
	updateInputFile(gt1);	// this will give correct OR one previous file
	while (gt1 > ipvar->ix2gt(ipvar->times.size()-1)){	// increment curr_file as long as gt1 +outside file range
		++curr_file;
		loadInputFileMeta();
	}

	clock_t start, end;
	double msecs;
	start = clock();
	
	CDEBUG << "readVar_reduce_mean (" << varname << ") :" << gt2string(gt1) << " " << gt2string(gt2) << endl;
	while(1){
//		cout << (ipvar->times[0]) << " " << (ipvar->times[ipvar->times.size()-1]) << " "  << (gt1-ipvar->tbase)*24.0/ipvar->tscale <<  endl;
		int tstart = lower_bound(ipvar->times.begin(), ipvar->times.end(), (gt1-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin();	   // first elem >= gt1 
		int tend   = upper_bound(ipvar->times.begin(), ipvar->times.end(), (gt2-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin() -1;   // last elem <= gt2
//		cout << gt2string(ipvar->ix2gt(tstart)) << " " << gt2string(ipvar->ix2gt(tend)) << endl;

		if (tend < 0) break;

		for (int i=tstart; i<=tend; ++i){ 
			ifile_handle->readVar(*ipvar, i, ipvar->ivar1);	// readCoords() would have set ivar1
//			double gt = ix2gt(i);
			
//			sxy = sxy + (*ipvar)*count;		// sum(xi*yi)
//			sx = sx + count;				// sum(xi)
//			sy = sy + (*ipvar);				// sum(yi)
//			sxx = sxx + count*count;		// sum(xi^2)
//			syy = syy + (*ipvar)*(*ipvar);	// sum(yi^2)
			
			for (int i=0; i<temp.values.size(); ++i){
				sxy[i] += (*ipvar)[i]*count;
				sx[i] += count;
				sy[i] += (*ipvar)[i];
				sxx[i] += count*count;
				syy[i] += (*ipvar)[i] * (*ipvar)[i];
			}
			
			++count;
		}

		if (tend >= ipvar->times.size()-1){ // if tend was the last time in file, then load next file and continue reading
			++curr_file;
			if (curr_file >= filenames.size()) break;
			else loadInputFileMeta();
		}
		else break;
	}
	
	CDEBUG << "----------- Read " << count << " timesteps from " << varname << endl;
	
	Sxy = sxy - sx*sy/count;
	Sxx = sxx - sx*sx/count;
	Syy = syy - sy*sy/count;
	
	b1 = Sxy/Sxx;
	s = (Syy - b1*Sxy)/(count-2);
	s.sqrtVar();
	
	gVar sqrt_Sxx = Sxx;
	sqrt_Sxx.sqrtVar();
	
	t = b1/(s/sqrt_Sxx);

	end = clock();
	msecs = ((double) (end - start)) * 1000 / CLOCKS_PER_SEC;
	
	cout << "Execution time is " << msecs/1000 << " sec" << endl;
//	if (count > 0) lterpCube(b1, *this, lterp_indices);	// We want to preserve current values if no new values were read
	return b1;
}




//int gVar::readVar_reduce_sd(double gt1, double gt2){
//	gVar s, ssq; 
//	s.copyMeta(*ipvar);
//	ssq.copyMeta(*ipvar);
//	s.values.resize(s.nlons*s.nlats*s.nlevs);
//	ssq.values.resize(ssq.nlons*ssq.nlats*ssq.nlevs);
//	s.fill(0); ssq.fill(0);
//	int count = 0;
//	for (double d = gt1; d < gt2; d += ipvar->tstep/24){
//		updateInputFile(d);
//		ifile_handle->readVar_gt(*ipvar, d, 0, ipvar->ivar1);	// read in mode 0 (Hold); readCoords() would have set ivar1
//		s = s + *ipvar;
//		ssq = ssq + (*ipvar)*(*ipvar);
//		++count;
//	}
//	if (count == 1) CERR << "Cannot compute SD with one data point" << endl;
//	gVar sd = ssq/(count-1) + s*s/count/(count-1);
//	sd.sqrtVar();
//	lterpCube(sd, *this, lterp_indices);
//}





