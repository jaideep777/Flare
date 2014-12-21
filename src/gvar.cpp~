#include "../include/gsm.h"


gVar::gVar(){
	ntimes = 0; nlevs = 1; nlats = 0; nlons = 0;
	tbase = 0; tscale = 1; tstep = 1;
	varname = ""; varunits = "";
	missing_value = std_missing_value;

	ifile_handle = new NcFile_handle();
	ofile_handle = new NcFile_handle();
}

gVar::gVar(string name, string units, string tunits){
	ntimes = 0; nlevs = 1; nlats = 0; nlons = 0;
	tstep = t = 0;
	varname = name; varunits = units;
	missing_value = std_missing_value;

	ifile_handle = new NcFile_handle;
	ofile_handle = new NcFile_handle;

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
	//cout << "Completed read coords Function\n";
}

// copy all data except values vector into self.
int gVar:: shallowCopy(const gVar &v){
	ntimes = v.ntimes; nlevs = v.nlevs; nlats = v.nlats; nlons = v.nlons;
	times = v.times; levs = v.levs; lats = v.lats; lons = v.lons;
	tbase = v.tbase; tscale = v.tscale; tstep = v.tstep; //t = v.t;
	varname = v.varname; varunits = v.varunits;
	scale_factor = v.scale_factor; add_offset = v.add_offset;
	ncoords = v.ncoords; ivar1 = v.ivar1;
	missing_value = v.missing_value;
	t = v.t;
	return 0;
}

// copy values vector and missing_value.
int gVar:: copyValues(const gVar &v){
	values = v.values;
	missing_value = v.missing_value;
	return 0;
}

int gVar::setCoords(vector <double> &t, vector <float> &le, vector <float> &la, vector <float> &lo){
	times = t; levs = le; lats = la; lons = lo;
	ntimes = t.size(); nlevs = le.size(); nlats = la.size(); nlons = lo.size();
	tstep = (times[1] - times[0])*tscale;	// tstep in hours
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
	lfout << "\t" << nlons << " lons: " << lons[0] << " ... " << lons[1]-lons[0] << " ... " << lons[nlons-1] << '\n';
	lfout << "\t"  << nlats << " lats: " << lats[0] << " ... " << lats[1]-lats[0] << " ... " << lats[nlats-1] << '\n';
	lfout << "\t"  << nlevs << " levs.\n";
	lfout << "\t"  << ntimes << " times.\n";
	lfout << "\t" << "Current number of values: " << values.size() << '\n';
	lfout << "> Time:\n";
	double tbaseIST = tbase + 5.5/24;
	lfout << "\tbase  gday (IST) = " << gtstr6d(tbaseIST) << ", i.e. " << gt2string(tbaseIST) << '\n';
	double t0 = ix2gt_IST(0);
	lfout << "\tfirst gday (IST) = " << gtstr6d(t0) << ", i.e. " << gt2string(t0) << '\n';
	t0 = ix2gt_IST(times.size()-1);
	lfout << "\tlast  gday (IST) = " << gtstr6d(t0) << ", i.e. " << gt2string(t0) << '\n';
	lfout << "\ttime step = " << tstep << " hours.\n";
	lfout << "\thours per time unit = " << tscale << "\n";
	lfout << "> Missing value = " << missing_value << "\n";
	lfout << "-------------------------------------------------------------------------\n\n";
	return 0;
}

int gVar::printValues(ostream &lfout){
	lfout << "> Values: (" << values.size() << ")\n";
	printCube(values, nlons, nlats, nlevs, missing_value);
}

// returns the index corresponding time just <= gt
int gVar::gt2ix(double gt){
	//CDEBUG << "tstep = " << tstep << "\n";
	return ((gt - tbase)*24.0 - times[0]*tscale)/tstep;
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

gVar gVar::operator + (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		cout << "ERROR summing variables: Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.shallowCopy(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] + v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

gVar gVar::operator - (const gVar &v){
	if (nlons != v.nlons || nlats != v.nlats || nlevs != v.nlevs){
		cout << "ERROR summing variables: Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.shallowCopy(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] - v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};

gVar gVar::operator * (const float x){
	gVar temp;
	temp.shallowCopy(*this);
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
	temp.shallowCopy(*this);
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
		cout << "ERROR summing variables: Grids not compatible.\n";
		gVar temp1;
		return temp1;
	}
	gVar temp;
	temp.shallowCopy(*this);
	temp.values.resize(nlevs*nlons*nlats, temp.missing_value);
	for (int i=0; i<nlevs*nlons*nlats; ++i){
		if (values[i] != missing_value && v.values[i] != v.missing_value)
			temp.values[i] = values[i] * v.values[i];
		else temp.values[i] = std_missing_value;
	}
	temp.missing_value = std_missing_value;
	return temp;
};


