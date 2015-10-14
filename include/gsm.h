#ifndef VECUTILS
#define VECUTILS

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <math.h>
#include <netcdfcpp.h>
using namespace std;


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  DEFS    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#define CDEBUG if (gsm_debug_on) cout << "<GSM debug> "
#define CDEBUGC if (gsm_debug_on) cout
#define CINFO if (gsm_info_on) cout << "<GSM info> "
#define CINFOC if (gsm_info_on) cout
#define CWARN if (gsm_warnings_on) cout << "GSM WARNING: "
#define CERR if (gsm_errors_on) cout << "GSM ERROR: "

const double t_tol = 1e-3;	// tolerance when comparing gt values
							// corresponds to ~86 sec
const float std_missing_value = 9.9e20;

extern bool gsm_info_on;
extern bool gsm_debug_on;
extern bool gsm_warnings_on;
extern bool gsm_errors_on;

const float glimits_globe[4] = {0, 360, -90, 90};
const float glimits_india[4] = {66.5, 100.5, 6.5, 38.5};
extern float glimits_custom[4];


/*~~~~~~~~~~~~~~~~~~~~~~~~~~   definition changing functions  ~~~~~~~~~~~~~~~~*/

//void setDebugFlag(bool f);
//void setInfoFlag(bool f);
//void setWarningFlag(bool f);
//void setErrorFlag(bool f);

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  TIME    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//	Definitions:

//	gday = global day = integer days since 1-3-0000 AD
//	xhrs = hours elapsed over and above integer days as fraction of 24 hrs. eg: 18 hrs = 0.75 xhrs
//	gt = global time = decimal days since 1-3-0000 AD (gday + xhrs)


// notes: 
// fractional day is on base 10 not base 24
// global time is in "days since 1-3-0000 AD"

// ------  date and time functions ------- 
string xhrs2hms(double dayf); // convert fractional day (dayf) to hh-mm-ss.s string
double hms2xhrs(string time); // convert hh-mm-ss(.s) string to fractional day 
int    ymd2gday(string date); // convert date string yyyy-mm-dd to global day
int    ymd2gday(int year, int month, int day); // convert y, m, d to global day
string gday2ymd(int g); // convert global time to readable date string yyyy-mm-dd
string gt2string(double gt); // convert gday (including day fraction) to full date-time
string gtstr6d(double gt);		// print gt upto 6 decimals

int gt2year(double gt);	// calculate year only (non-decimal) from gday
int gt2month(double gt);	// calculate current month from gday (NOTE: month ranges from 1-12 and NOT from 0-11)
int gt2day(double g);		// calculate day in month
int gt2daynum(double gt);	// calculate day of year
int gt2dayOfYear(double gt);	// calculate day of year

// ------  lat lon functions ------- 
float sex2decLL(string s);	// convert "ll mm ss D" to decimal lat/lon
string dec2sexLL(float lon); // convert decimal lat/lon to "ll mm ss D"



/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  ARRAYUTILS   ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

string int2str(int n);
float str2float(string s);
int str2int(string s);

// get 1D index of 3D point (ilon,ilat,ilev) from a cube (nlons,nlats,nlevs)
int IX3(int ix, int iy, int iz, int nx, int ny); 
// get 1D index of 2D point (ilon,ilat) from a cube (nlons,nlats)
int IX2(int ix, int iy, int nx);  

void printArray(float v[], int n);	// print n elements of array on single line
void printArray(vector <float> &v, int n=0);	// print float vector
void printArray2d(float v[], int rows, int columns);	// print 2d array with rows & columns
void printArray2d(vector <float> &v, int rows, int columns);	// print float vector 2d
void printCube(vector <float> &v, int nx, int ny, int nz=1, 
				float ignoreVal = std_missing_value); // print a data cube, ignore missing_values

void setZero(vector <float> &v);	// set all of vector to zero
void setValue(vector <float> &v, float value); // set all elements of vector to value

void reverseArray(vector <float> &orig);	// reverse given array
void reverseCube(vector <float> &v, int nx, int ny, 
				 int nz=1, int n4=1, int n5=1); // reverse a data cube along lat dimension

vector <float> copyArray(vector <float> &v, int i2, int i1=0); // copy v[i1:i2] into new array (returned)

int ncIndexLo(vector <float> &v, float val); // lower bound, edge for outliers
int ncIndexHi(vector <float> &v, float val); // upper bound, edge for outliers
int lindexSW(vector <float> &v, float val);  // lower (S/W) bound, missing value for outliers
int indexC(vector <float> &v, float val);    // cell index by center, missing value for outliers
vector <float> max_vec(vector <float> &u, vector <float> &v);	// return elementwise maximum 

// Array operations
float sum(vector <float> &v);
float avg(vector <float> &v);
vector <float> dim_sum(vector <float> v, int idim, vector <float> dimsizes);
vector <float> operator * (const vector <float> &x, const vector <float> &y);
vector <float> operator / (const vector <float> &x, const float c);


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GVAR    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

// all operations in gVar ignore fill values.
// lterp (in grid.h) may have to be rewritten because missing values grow around an existing one.
class NcFile_handle;

class gVar{
	public:
	int ntimes, nlevs, nlats, nlons;
	vector <float> levs, lats, lons;
	vector <double> times;
	vector <float> values;
	double tbase;
	float tscale, tstep;
	double t;	// the time for which values are currently held
	string varname, varunits;
	float scale_factor, add_offset;
	int ncoords, ivar1; // ivar1 is the index of 1st data variable
	float missing_value;
	
	string ifname, ofname;	// input and output filenames
	NcFile_handle *ifile_handle, *ofile_handle;	// ip and op file handles
	NcVar * outNcVar;		// output NcVar
	bool lwrite, lwriteSP;			// 'write to output' flag (nc, singlePointOutput)
	
	gVar();
	gVar(string name, string units, string tunits);
	
	int shallowCopy(const gVar &v);	// copy everything except values vector from v
	int copyValues(const gVar &v);	// copy values vector and missing_value from v

	int setCoords(vector <double> &t, vector <float> &le, vector <float> &la, vector <float> &lo);
	int setTimeAtts(int xntimes, double xtbase, float xtscale);

	int printGrid(ostream &lfout = std::cout);	// print succint grind info
	int printValues(ostream &lfout = std::cout); // print values
	
	int gt2ix(double gt);	// find index in time coord corresponding to global time (inc day fraction)
	double ix2gt(int ix);	// find GT fractional gday corresponding to given index 
	double ix2gt_IST(int ix); // find GT fractional gday for given index but add 5.5 hrs to convert to IST
	
	// get value by coords
	float getValue(float xlon, float xlat, float ilev = 0);		// get interpolated value at (xlon, xlat)
	float getCellValue(float xlon, float xlat, float ilev = 0);	// get cell-avg value for xlon, xlat) 
	
	// functions on gVars 
	int fill(float f);
	int sqrtVar();
	
	// operators
	float& operator () (int ilon, int ilat, int ilev);	// get reference to element by coord indices, like n-D array
	float& operator [] (int i);	// get reference to element by memory index, like n-D array
	gVar operator + (const gVar &v);
	gVar operator - (const gVar &v);
	gVar operator * (const gVar &v);
	gVar operator * (const float x);
	gVar operator / (const float x);
};

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  NCIO    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//	NcFile_handle class has resources to read/write Nc Files.
//	this class has no variables by itself and must input all relevant gVars
//	therefore destructor frees all pointers.
class NcFile_handle{
public:
	NcFile * dFile;
	string fname;
	string mode;
	static const int NC_ERR = 2;
	NcVar *tVar, *levVar, *latVar, *lonVar;
	NcDim *tDim, *levDim, *latDim, *lonDim;
	int ntimes, nlevs, nlats, nlons;
	int ncoords, nvars;
	string levname, levunits;
	bool latSN;	// if true, lats are SN, increase with index
	bool lonPos; // if true, lons range from 0-360 (desired), first index is lon = 0
	float wlon, elon, slat, nlat;	// map limits
	int wlonix, elonix, slatix, nlatix;	// map limits indices. 
	int ilon0, ilat0, itime0;
	// NOTE: above indices are defined on the native arrays and not reversed ones
	bool mplimited; // if true, map limits are set
	
	
	// constructor. Only initializes variables. Does not open file.
	NcFile_handle();	// open NC object (nc file)
	
	void setMapLimits(float xwlon, float xelon, float xslat, float xnlat);	// set lat-lon limits
	int  open(string s, string m, const float glimits[4] = glimits_globe); // open file s
	int close(); // close file
	~NcFile_handle();	// dFile must be deleted in destructor

	// reading functions
	int readCoords(gVar &v, ostream &lfout = cout, bool rr = true); // read coord metadata, read values if rr is true
	int getVarID(string varname);	// get variable id (ivar) from name
	int readVarAtts(gVar &v, int ivar = -1); // get variable name, units, missing_value etc.
	int readVar(gVar &v, int itime, int iVar = -1); // read variable with index iVar from file into "v"
	int readVar_gt(gVar &v, double gt, int mode, int iVar = -1); // read variable corresponding to GT gt
		// as specified by mode: 0 = hold values from previous step, 1 = interpolate
	
	// writing functions
	int writeCoords(gVar &v, bool wr = true); // write coord metadata, write values if wr is true
	NcVar * createVar(gVar &v); // create variable v into file matching gVar
	int writeVar(gVar &v, NcVar * vVar, int itime); // write values into vVar at time itime
	int writeTimeValues(gVar &v); // write time values
};
	

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GRID    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// create coordinate vector given x0, xf and number of points
vector <float> createCoord(float x0, float xf, int nx, float &dx);

// create coordinate vector given x0, xf and resolution
vector <float> createCoord(float x0, float xf, float dx, int &nx);

// print a gridded variable along with 2d coordinates defined on x, y
void printVar(vector <float> &x, vector <float> &y, vector <float> &data);


// find index (u,v) of gridbox containing point (x,y)
// ASSUMES THAT grid coordinates map to SW corners of gridboxes on ground.
//    uvs is the **SW corner** of gridbox that contains (x,y), i.e.
//    (x,y) must lie within [x0+mdx,x0+(m+1)dx] and ~ly y
// --------------------------------------------------------------------------
vector <int> findGridBoxSW(float x,float y, vector <float> &lons, vector <float> &lats);


// find index (u,v) of gridbox containing point (x,y)
// ASSUMES THAT coordinates map to centers of gridboxes on ground.
//    uvs is the **CENTER** of gridbox that contains (x,y), i.e.
//    (x,y) must lie within [x0 +/- mdx/2] and ~ly y
//    * Since dealing with India only, no need to deal with edges/poles
// --------------------------------------------------------------------------
vector <int> findGridBoxC(float x,float y, vector <float> &lons, vector <float> &lats);


// During initialization, calculate the SW index for each point of model grid, 
// so that it can be resued at every step in bilinear_mn
vector <int> bilIndices(vector <float> &lons, vector <float> &lats,
				 		vector <float> &mlons, vector <float> &mlats);


// return bilinear interpolated value at point (x,y) in grid {lons,lats}
float bilinear(float x, float y, float iz, 
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &data, float missingVal = std_missing_value);

// return bilinear interpolated value at point (x,y) in grid {lons,lats}
// use bilIndices instead of finding them in situ
float bilinear(int ilat, int ilon, int iz, vector <int> &indices,
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &mlons, vector <float> &mlats,
			   vector <float> &data, float missingVal = std_missing_value);

// return value at gridcell containing point (x,y) in grid {lons,lats}
float cellVal(float x, float y, float iz, 
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &data, float missingVal = std_missing_value);

// return value at gridcell containing point (x,y) in grid {lons,lats}
// use bilIndices instead of finding them in situ
float cellVal(int ilat, int ilon, int iz, vector <int> &indices,
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &mlons, vector <float> &mlats,
			   vector <float> &data, float missingVal = std_missing_value);

// mask variable v using m as mask
// for all points in v, if (m <= val) v = missing_value
gVar mask(gVar &v, gVar &m, float val = 0);

// *** Following regridding regridding functions copy metadata (with ofc, new lats/lons)
// interpolate variable v onto given grid (xlons, xlats)
gVar lterp(gVar &v, vector <float> &xlons, vector <float> &xlats); 
// ^ *** IMP! *** order of arguments was (v, xlats, xlons) in previous version!

// *** Following regridding functions DO NOT copy metadata.. 
// .. they assume (and check) compatibility
// .. input and output coords can be taken directly from input and output gVars
int lterpCube(gVar &v, gVar &out, vector <int> &indices);
int cellRegridCube(gVar &v, gVar &out, vector <int> &indices);



#endif

