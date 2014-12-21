#include <math.h>
#include "../include/gsm.h"


/******************       class NcFile_handle      ***********************/

NcFile_handle::NcFile_handle(){
	levname = "lev"; levunits = "z";
	nlons = nlats = ntimes = 0; nlevs = 1;
	nvars = ncoords = 0;
	ilon0 = ilat0 = 0;
	wlonix = elonix = slatix = nlatix = 0;
	mplimited = false;
}

// ------------------------- INIT ----------------------------
int NcFile_handle::open(string s, string m, const float glimits[4]){
	fname = s;
	mode = m;
	if (mode == "r"){ 
		CDEBUG << "attaempting to open file: " << s << '\n';
		dFile = new NcFile(s.c_str(), NcFile::ReadOnly);
		if (!dFile)	{cout << "Error Opening File: " << s << "\n"; return 1;}
	}
	else if (mode == "w"){ 
		dFile = new NcFile(s.c_str(), NcFile::Replace);
		if (!dFile)	{cout << "Error Opening File: " << s << "\n"; return 1;}
	}
	
	mplimited = (glimits[0] > 0 || glimits[1] < 360 || glimits[2] > -90 || glimits [3] < 90)? true:false; 
	wlon = glimits[0]; //0.f;
	elon = glimits[1]; //360.f;
	slat = glimits[2]; //-90.f;
	nlat = glimits[3]; //90.f;
	
	if (dFile->is_valid()) return 0;
	else return 1;
}

int NcFile_handle::close(){
	dFile->close();
	return 0;
}

NcFile_handle::~NcFile_handle(){
	//dFile->close();	including this statement gives seg fault!
	delete dFile;
	dFile = NULL;
}

void NcFile_handle::setMapLimits(float xwlon, float xelon, float xslat, float xnlat){
	mplimited = true;
	wlon = xwlon;
	elon = xelon;
	slat = xslat;
	nlat = xnlat;
}


// ------------------------- READING ----------------------------

// THE MOST PAINFUL NCIO FUNCTION!!
// this function:
//	 reads most of the metadata - read # coords, coord values
//	 sets the lat lon limits
//	 sets the correct order of lats and lons according to model requirements
int NcFile_handle::readCoords(gVar &v, ostream &lfout, bool rr){

	// get number of variables
	lfout << "> Reading file: " << fname << '\n';
	nvars = 0;
	nvars = dFile->num_vars();
	lfout << '\t' << nvars << " variables found.\n";

	// Get the coordinate variables
	ncoords = 0;
	latVar = dFile->get_var("lat");
	if (!latVar) latVar = dFile->get_var("latitude");
	if (!latVar) latVar = dFile->get_var("LAT");
	if (!latVar) {CWARN << "lat not found!\n"; latVar = 0;}
	else ++ncoords;

	lonVar = dFile->get_var("lon");
	if (!lonVar) lonVar = dFile->get_var("longitude");
	if (!lonVar) lonVar = dFile->get_var("LON");
	if (!lonVar) {CWARN << "lon not found!\n"; lonVar = 0;}
	else ++ncoords;

	levVar = dFile->get_var("lev");
	if (!levVar) levVar = dFile->get_var("levels");
	if (!levVar) {/*CWARN << "lev not found!\n";*/ levVar = 0;}
	else ++ncoords;

	tVar = dFile->get_var("time");
	if (!tVar) tVar = dFile->get_var("TIME");
	if (!tVar) {CWARN << "time not found.\n"; tVar = 0;}
	else ++ncoords;
	lfout << "\t" << ncoords << " coordinates found.\n"; lfout.flush();
		
	int nbounds = 0;
	// several files have bounds for coordinates.. if found, ignore for now.	
	NcVar* latBnds = dFile->get_var("lat_bnds");
	if (!latBnds) latBnds = dFile->get_var("lat_bounds");
	else ++nbounds;

	NcVar* lonBnds = dFile->get_var("lon_bnds");
	if (!lonBnds) lonBnds = dFile->get_var("lon_bounds");
	else ++nbounds;

	NcVar* levBnds = dFile->get_var("lev_bnds");
	if (!levBnds) levBnds = dFile->get_var("lev_bounds");
	else ++nbounds;

	NcVar* tBnds = dFile->get_var("time_bnds");
	if (!tBnds) tBnds = dFile->get_var("time_bounds");
	else ++nbounds;

	lfout << "\t" << nbounds << " coordinate bounds found.\n";
	v.ncoords = ncoords; 
	
	v.ivar1 = ncoords+nbounds;	// assumes that 1st few variables are coords/bounds
								// will NOT always work. better set explicitly.

	// read lats, lons metadata from file
	// read actual values only is rr is true

	// if lat order is found N-S in file, reverse it. 
	// in this model lat order will be S-N (i.e lats increase with index)
	// all data reading will depend crucially on the variable latSN
	lfout << "> Reading coordinates.\n";
	if (latVar){
		CINFO << "reading lats... ";
		nlats = v.nlats = *(latVar->edges());	// otherwise constructor has init to 0
		v.lats.resize(v.nlats);
		if (rr) latVar->get(&v.lats[0], v.nlats);
		CINFOC << v.nlats << " read.\n";

		float dlat = v.lats[1] - v.lats[0];
		latSN = (dlat < 0)? false:true;
		
		if (mplimited){
			CINFO << "trimming lats array to match lat limits... ";
			// set indices
			slatix = ncIndexLo(v.lats, slat);
			nlatix = ncIndexHi(v.lats, nlat);
			ilat0 = (nlatix > slatix)? slatix:nlatix;
			// cut array
			v.lats = copyArray(v.lats, nlatix, slatix);
			v.nlats = fabs(nlatix - slatix) +1;
			CINFOC << v.nlats << " lats left.\n";
		}

		if (!latSN){
			CINFO << "lats array is N-S. Reversing...";
			reverseArray(v.lats);
			CINFOC << " Done.\n";
		}
	}

	if (lonVar){
		CINFO << "reading lons... ";
		nlons = v.nlons = *(lonVar->edges());	// otherwise constructor has init to 0
		v.lons.resize(v.nlons);
		if (rr) lonVar->get(&v.lons[0], v.nlons);
		CINFOC << v.nlons << " read.\n";

		if (mplimited){
			CINFO << "trimming lons array to match lon limits... ";
			// set indices
			wlonix = ncIndexLo(v.lons, wlon);
			elonix = ncIndexHi(v.lons, elon);
			ilon0 = (elonix > wlonix)? wlonix:elonix;
			// cut array
			v.lons = copyArray(v.lons, elonix, wlonix);
			v.nlons = elonix - wlonix +1;
			CINFOC << v.nlons << " lons left.\n";
		}
	}

	if (levVar){
		CINFO << "reading levs... ";
		nlevs = v.nlevs = *(levVar->edges());	// otherwise constructor has init to 1
		v.levs.resize(v.nlevs);
		if (rr) levVar->get(&v.levs[0], v.nlevs);
		CINFOC << v.nlevs << " read.\n";
	}
	
	if ( tVar){
		CINFO << "reading t... ";
		ntimes = v.ntimes = *(tVar->edges());	// otherwise constructor has init to 0
		v.times.resize(v.ntimes);
		if (rr) tVar->get(&v.times[0], v.ntimes);
		CINFOC << v.ntimes << " read.\n";
		// set itime0 from sim start time
		
	}

	// read and set time base settings
	// (this is from earlier function settbase() )
	NcAtt * a = tVar->get_att("units");
	string datestr = a->as_string(0);

	string unit, junk, sdate, stime;
	stringstream ss;
	ss.clear(); ss.str(datestr);
	ss >> unit >> junk >> sdate >> stime;
	if (stime == "") stime = "0:0:0";
	
	v.tbase = ymd2gday(sdate) + hms2xhrs(stime); // note this time is in GMT
	if (unit == "hours") v.tscale = 1.0f;
	else if (unit == "days") v.tscale = 24.0f;
	else if (unit == "months") {
		lfout << "WARNING: using months as time units! 365.2524 days/yr will be considered.\n";
		v.tscale = (365.2524/12.0)*24.0;
	}
	else {
		lfout << "ERROR setting base time in getCoords(): invalid time units!\n";
		return 1;
	}
	//lfout << "Completed read coords Function\n";
	v.tstep = (v.times[1] - v.times[0])*v.tscale;	// tstep in hours
	
	return 0;
}

int NcFile_handle::getVarID(string varname){
	NcVar * v;
	for (int i=0; i< nvars; ++i){
		v = dFile->get_var(i);
		if (v->name() == varname) return i;
	}
	CWARN << "No variable named " << varname << " found.\n";
}

// if file has only 1 variable, then ivar = ncoords since [0,ncoords-1] are coords
int NcFile_handle::readVarAtts(gVar &v, int ivar){
	if (ivar == -1) ivar = v.ivar1;	// ivar not specified, default is -1, so set it to ivar1
	NcVar * nVar = dFile->get_var(ivar);
	// name
	v.varname = nVar->name();
	// units
	NcAtt * A = nVar->get_att("units");
	if (A) v.varunits = A->as_string(0);
	delete A;	
	// missing value
	NcAtt * B = nVar->get_att("missing_value");
	if (B) v.missing_value = B->as_float(0);
	else {
		B = nVar->get_att("_FillValue");
		if (B) v.missing_value = B->as_float(0);
	}
	delete B;
	// scale factor
	v.scale_factor = 1.f;
	NcAtt *sc = nVar->get_att("scale_factor");
	if (sc) v.scale_factor = sc->as_float(0);
	delete sc;
	// offset
	v.add_offset = 0.0f;
	NcAtt *o = nVar->get_att("add_offset");
	if (o) v.add_offset = o->as_float(0);
	delete o;
	// apparently variables created with pointers need not be deleted..
	// deleting messes up the variable IDs
	return 0;
}

// if file has only 1 variable, then ivar = ivar1
// ivar need to be specified only in case the file has multiple variables
int NcFile_handle::readVar(gVar &v, int itime, int iVar){
	if (mode != "r"){
		cout << "ERROR in readVar: File not in read mode.\n";
		return 1;
	}

	// file is in read mode.. continue.
	if (iVar == -1) iVar = v.ivar1;
	// some crackpot nc files have coords last. so this check shouldnt be made. 
	// maybe default ivar should not be allowed!!
//	if (iVar < v.ncoords){
//		cout << "FATAL ERROR in readVar(): Attempted to read a coord variable.\n";
//		return 1;
//	}
	
	// allocate space for values
	v.values.resize(v.nlevs*v.nlats*v.nlons); // no need to fill values 
	NcVar * vVar = dFile->get_var(iVar);
	
	// read data. cur is set at the SW corner of grid / grid-limits.
	if (vVar->num_dims() == 4){ 
		vVar->set_cur(itime, 0, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner
		vVar->get(&v.values[0], 1, v.nlevs, v.nlats, v.nlons);
	}	
	else if (vVar->num_dims() == 3){ 
		vVar->set_cur(itime, ilat0, ilon0); // set the starting time at itime and lat/lon/lev at SW corner
		vVar->get(&v.values[0], 1, v.nlats, v.nlons);
	}
	else {
		cout << "\tVar Error: Cant handle 1D/2D variables yet..\n";
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
	// set t to gday corresponding to times[itime]
	v.t = v.ix2gt(itime);

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
int NcFile_handle::readVar_gt(gVar &v, double gt, int mode, int iVar){
	// calculate index corresponding to gt


	int tixlo = v.gt2ix(gt);		// index for curret gt
	int tixlo_prev = v.gt2ix(v.t);	// index for gt last read into v

	CDEBUG << "readVar_gt(): gt | ix | gt(ix) = " << gt2string(gt) << " | ";
	CDEBUGC << tixlo << " | ";
	CDEBUGC << gt2string(v.ix2gt(tixlo)) << " -> ";

	// in mode 0, if both these indices match, data is already in variable, so skip reading
	if (tixlo == tixlo_prev && mode == 0) { 
		v.t = gt;	// set t to the actual gt and not to any index values
		return 0;	
	}
	// if times[index] = gt, then no need of interpolation
	// also, if mode = hold, values should be from tixlo anyway, so include that case here
	else if ( gt-v.ix2gt(tixlo) < t_tol || mode == 0){
		CDEBUGC << "~low, diff = " << (gt-v.ix2gt(tixlo))*24.0 << " hrs. Read low.\n";
		readVar(v, tixlo, iVar);
		v.t = gt;	// set t to the actual gt and not to any index values
		return 1;
	}
	else if ( v.ix2gt(tixlo+1)-gt < t_tol || mode == 0){	
		// we dont want to hold/interpolate for gt that is too close to the higher bound.. hence this case
		CDEBUGC << "~high, diff = " << (v.ix2gt(tixlo+1)-gt)*24.0 << " hrs. Read Hi.\n";
		readVar(v, tixlo+1, iVar);
		v.t = gt;	// set t to the actual gt and not to any index values
		return 2;
	}
	else{	// read the values at next step and calculate.
		CDEBUGC << "intermediate, diff = " << (gt-v.ix2gt(tixlo))*24.0 << " hrs. Interpolate\n";
		gVar z; z.shallowCopy(v);
		gVar w; w.shallowCopy(v);	// make a new gVar to hold next step values, identical to v
		readVar(z, tixlo, iVar);	// read this step in v
		readVar(w, tixlo+1, iVar);	// read next step in w
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
		latDim = dFile->add_dim("lat", v.nlats);	// dim
		latVar = dFile->add_var("lat", ncFloat, latDim); // var
		latVar->add_att("units", "degrees_north");	// attributes
		latVar->add_att("scale_factor", 1.0f);
		latVar->add_att("add_offset", 0.0f);
		if (wr) latVar->put(&v.lats[0], v.nlats);	// write data
	}
	//CINFO << "Wrote lats to file." << endl;
	
	if (v.nlons > 0){
		lonDim = dFile->add_dim("lon", v.nlons);
		lonVar = dFile->add_var("lon", ncFloat, lonDim);
		lonVar->add_att("units", "degrees_east");
		lonVar->add_att("scale_factor", 1.0f);
		lonVar->add_att("add_offset", 0.0f);
		if (wr) lonVar->put(&v.lons[0], v.nlons);
	}
	//CINFO << "Wrote lons to file." << endl;
	
	if (v.nlevs > 1){
		levDim = dFile->add_dim(levname.c_str(), v.nlevs);
		levVar = dFile->add_var(levname.c_str(), ncFloat, levDim);
		levVar->add_att("units", levunits.c_str());
		levVar->add_att("scale_factor", 1.0f);
		levVar->add_att("add_offset", 0.0f);
		if (wr) levVar->put(&v.levs[0], v.nlevs);
	}
	//CINFO << "Wrote levs to file." << endl;
	
	if (v.ntimes > 0){
		tDim = dFile->add_dim("time");	// unlimited dimension
		tVar = dFile->add_var("time", ncFloat, tDim);
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
		tVar->add_att("units", tunits.c_str());
		tVar->add_att("scale_factor", 1.0f);
		tVar->add_att("add_offset", 0.0f);
		// time vector not written now because ntimes may be determined later
	}
	//CINFO << "Wrote time to file." << endl;
	return 0;
}	

// time coord is written seperately because several times number of times
// is determined later 
int NcFile_handle::writeTimeValues(gVar &v){
	tVar->put(&v.times[0], v.ntimes);
}

NcVar * NcFile_handle::createVar(gVar &v){
	if (mode != "w"){
		CERR << "ERROR in createVar: File not in write mode.\n";
		return NULL;
	}
	if (v.varname == "" || v.varunits == ""){
		CERR << "ERROR in createVar: variable Name or Units not set.\n";
		return NULL;
	}
	// file is in write mode.. continue.
	// include level dimension only if nlevs > 0.
	NcVar * vVar;
	if (v.nlevs <= 1) vVar = dFile->add_var(v.varname.c_str(), ncFloat, tDim,         latDim, lonDim);
	else			  vVar = dFile->add_var(v.varname.c_str(), ncFloat, tDim, levDim, latDim, lonDim);
	
	// add variable attributes
	vVar->add_att("units", v.varunits.c_str());
	vVar->add_att("scale_factor", 1.0f);
	vVar->add_att("add_offset", 0.0f);
	vVar->add_att("_FillValue", v.missing_value);
	vVar->add_att("mising_value", v.missing_value);
	
	// actual data not written here because ..
	// we need to write data into already created variable
	//vVar->put_rec(v.values, itime);
	
	return vVar;
}

int NcFile_handle::writeVar(gVar &v, NcVar* vVar, int itime){
	if (mode != "w"){
		cout << "ERROR in writeVar: File not in write mode.\n";
		return 1;
	}
	
	// actually write the data	
	vVar->put_rec(&v.values[0], itime);
	return 0;
}

