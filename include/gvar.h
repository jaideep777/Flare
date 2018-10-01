#ifndef GVAR_H
#define GVAR_H

#include <iostream>
#include <vector>
#include <string>
#include <netcdfcpp.h>

using namespace std;


/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GVAR    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

// all operations in gVar ignore fill values.
// lterp (in grid.h) may have to be rewritten because missing values grow around an existing one.
class NcFile_handle;


/**
	\ingroup grid
	\brief Georeferenced Variable.
*/
class gVar{
	public:
	int ntimes, //!< Number of timesteps that this variable references (note that at any point, only one timestep is held in the variable).
		nlevs, 	//!< Number of levels in the data.
		nlats, 	//!< Number of latitudes (rows) in the data.
		nlons;	//!< Number of longitudes (columns) in the data.

	vector <float> levs, 	//!< Levels 
				   lats, 	//!< Latitudes associated with data rows
				   lons;	//!< Longitudes associated with data columns
	vector <double> times;	//!< Time vector of the variable (This is useful when reading/writing data to files).

	double tbase;			//!< Base time (values in the time vector are measured in units since this base time)

	float tscale, 	//!< Time unit in hours (hours/time-unit)
		  tstep; 	//!< Time-step in hours 

	double t;		//!< The time for which values are currently held
	
	string varname,		//!< Variable name 
		   varunits;	//!< Variable units
		   
	float scale_factor, add_offset;

	int ncoords,	//!< Number of coordinates in the associated inout file 
		ivar1; 		//!< index of 1st data variable in the associated inout file 
		
	float missing_value;	//!< Missing value (what value to treat as missing data)

	vector <float> gridlimits;	//!< Lat-Lon bounds 
	
	bool lwrite, lwriteSP;			// 'write to output' flag (nc, singlePointOutput)
	
	private:
	string ifname, ofname;	// input and output filenames
	NcFile_handle *ifile_handle, *ofile_handle;	// ip and op file handles
	NcVar * outNcVar;		// output NcVar
	vector <int> 		lterp_indices;			// interpolation indices (generated during init)
	vector <string>		filenames;					// list of filenames (genrated during init)
	int curr_file;
	gVar * ipvar;	// gVar to read data into
	string regriddingMethod;
	
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~
	public:
	vector <float> values;	//!< Data values. These are stored as a 1D array.
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~


	public:
	gVar();
	gVar(string name, string units, string tunits);

	int initMetaFromFile(string filename);
	
	int _copyMeta(const gVar & v);	// copy all metadata excecpt coords
	int copyMeta(const gVar &v);	// copy all metadata including coords and resize
	int copyMeta(const gVar &v, vector <float> &_lons, vector <float> &_lats, vector <float> &_levs); // copy meta from v but replace coords (useful in interpolations/coarsegraining functions)	
	
	int copyValues(const gVar &v);	// copy values vector and missing_value from v

	int setCoords(vector <double> &t, vector <float> &le, vector <float> &la, vector <float> &lo);
	int setTimeAtts(int xntimes, double xtbase, float xtscale);

	int printGrid(ostream &lfout = std::cout);	// print succint grind info
	int printGridIP(ostream &lfout = std::cout);	// print succint grind info
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
	float& operator [] (int i);	// get reference to element by memory index, like 1-D array
	gVar operator + (const gVar &v);
	gVar operator + (const float x);
	gVar operator - (const gVar &v);
	gVar operator - (const float x);
	gVar operator * (const gVar &v);
	gVar operator * (const float x);
	gVar operator / (const float x);
	gVar operator / (const gVar &v);

	void setRegriddingMethod(string m);
	int createNcInputStream(vector <string> files, vector <float> glim, string rm = "bilinear");
	int loadInputFileMeta();
	int whichNextFile(double gt);
	int updateInputFile(double gt);
	int closeNcInputStream();
	
	int readVar_gt(double gt, int mode); 
	int readVar_it(int tid);

	int readVar_reduce_mean(double gt1, double gt2);
//	int readVar_reduce_sd(double gt1, double gt2);
	gVar trend(double gt1, double gt2);
	gVar trend_gpu(double gt1, double gt2);
	
	// these 2 functions create a gVar in one shot by reading the first record in specified file
	// createOneShot uses file's coords, readOneShot uses variable's coords and interpolates data
	int createOneShot(string filename, vector<float> glim = vector <float> ());
	int readOneShot(string filename, vector<float> glim = vector <float> ());
	
	// writing functions
	int createNcOutputStream(string filename);
	int closeNcOutputStream();
	int writeVar(int itime); 

	int writeOneShot(string filename);
	
};

#endif


