#ifndef GVAR_H
#define GVAR_H

#include <iostream>
#include <vector>
#include <string>
#include <netcdf>

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
		   
	float scale_factor, add_offset;	// TODO: These should be moved to NcFile_handle, as these are file properties.

	int ncoords,	//!< Number of coordinates in the associated inout file 
		ivar1; 		//!< index of 1st data variable in the associated inout file 
		
	float missing_value;	//!< Missing value (what value to treat as missing data)

	vector <float> gridlimits;	//!< Lat-Lon bounds 
	
	bool lwrite, lwriteSP;			// 'write to output' flag (nc, singlePointOutput)
	
	private:
	string ifname, ofname;	// input and output filenames
	NcFile_handle *ifile_handle, *ofile_handle;	// ip and op file handles
	netCDF::NcVar outNcVar;		// output NcVar
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
	
	/**
		\brief Default Constructor	
	*/
	gVar();

	/**
		\brief Constructor with name and units	
	*/
	gVar(string name, string units, string tunits);

	/**
		\brief Set metadata from a NetCDF File.	
	*/
	int initMetaFromFile(string filename	//!< NetCDF file name
						);
	
	/**
		\brief Set Metadata (except coordinates) fom another Georeferenced variable.
		
		This function sets all metadata except the coordinates. This is useful in 
		functions like regridding where coordinates need to be modified. 	
	*/
	int _copyMeta(const gVar & v);	// copy all metadata excecpt coords

	/**
		\brief Set ALL metadata fom another Georeferenced variable.
		
		This function sets all metadata including the coordinates from the supplied variable.  	
	*/
	int copyMeta(const gVar &v);	// copy all metadata including coords and resize

	/**
		\brief Set metadata except cooridnates and set coordinates from specified values.
		
		This function sets all metadata except the coordinates from the supplied variable. 
		Coordinates are additionally set using the arguments supplied. 
		(Useful in interpolations/coarsegraining functions)	
	*/
	int copyMeta(const gVar &v, 			//!< Georeferenced variable from which to copy metadata
				 vector <float> &_lons, 	//!< New lons
				 vector <float> &_lats, 	//!< New lats
				 vector <float> &_levs		//!< New levels
				 ); 
	
	/**
		\brief Copy values from another variable.
	*/
	int copyValues(const gVar &v);	// copy values vector and missing_value from v


	/**
		\brief Set coordinates
	*/
	int setCoords(vector <double> &t,	//!< Times - Must be in "units since <base_date>"
				  vector <float> &le,	//!< Levels 
				  vector <float> &la,	//!< Lats - Must be in ascending order (-90 --> 90)
				  vector <float> &lo	//!< Lons  
				  );
				  
	/**
		\brief Set Time attributes (units, base date).
	*/
	int setTimeAtts(int xntimes, 		//!< Number of timesteps
					double xtbase, 		//!< Base time from which the values in the time vector are measured (must be days since 1 March 0000 AD). The ymd2gday() function may be used to calculate the base time in appropriate units.  
					float xtscale		//!< Hours per time-unit. For e.g., if time unit is days, xtscale = 24, if time unit is hours, xtscale = 1, etc. 
					);

	int printGrid(ostream &lfout = std::cout);	// print succint grind info
	int printGridIP(ostream &lfout = std::cout);	// print succint grind info
	int printValues(ostream &lfout = std::cout); // print values
	
	/** \name time-index conversions */
	//@{
	int gt2ix(double gt);	//!< find the index in variable's time vector that corresponds to global time gt (including day fraction). The highest index just <= gt is returned. 
	double ix2gt(int ix);	//!< Convert index ix in the variable's time vector to global time 
	double ix2gt_IST(int ix); //!< Convert index ix in the variable's time vector to global time +5.5 hours (Indian standard time) 
	//@}
	
	// get value by coords
	float getValue(float xlon, float xlat, float ilev = 0);		// get interpolated value at (xlon, xlat)
	float getCellValue(float xlon, float xlat, float ilev = 0);	// get cell-avg value for xlon, xlat) 
	
	// functions on gVars 
	int fill(float f);
	int sqrtVar();
	int logshift(float a);	//!< Performs log(a + self)
	
	/** \name operators */
	//@{
	float& operator () (int ilon, int ilat, int ilev);	//!< Get reference to the data element at specified coordinate indices.
	float& operator [] (int i);	//!< Get reference to the data element by directly accessing the 1D values array.
	gVar operator + (const gVar &v); //!< Add 2 gVars, returning missing_value when either operand is missing.
	gVar operator + (const float x); //!< Add a constant to gVar
	gVar operator - (const gVar &v); //!< Subtract 2 gVars, returning missing_value when either operand is missing.
	gVar operator - (const float x); //!< Subtract a constant from gVar
	gVar operator * (const gVar &v); //!< Multiply 2 gVars, returning missing_value when either operand is missing.
	gVar operator * (const float x); //!< Multiply gVar with a constant 
	gVar operator / (const float x); //!< Division, returning missing_value when either operand is missing.
	gVar operator / (const gVar &v); //!< Divide gVar with a constant
	//@} 
	

	private:
	int loadInputFileMeta();
	int whichNextFile(double gt);
	int updateInputFile(double gt, bool suppress_out_of_bounds_error_printing = false);
	
	public:
	/** \brief Set the regridding method to use when reading data from an input stream.
		\param m String specifying the regridding method. Currently, possible options 
		are "bilinear" and "none" (If "none" is specified, the lats-lons in the input 
		file must match exactly with those in the variable). In future releases, 
		coarseGrain and Nearest-Neighbour regridding will be supported.  		
	*/
	void setRegriddingMethod(string m);

	/** \name One-shot NetCDF reading and writing 
		These functions can be used to read or write a single time slice from/to a NetCDF file. 
		The functions automatically read/write all the necessary metadata. These functions are particularly 
		convenient quickly exchanging data from files and gVars.  
	*/
	//@{	 
	// these 2 functions create a gVar in one shot by reading the first record in specified file
	int createOneShot(string filename, vector<float> glim = vector <float> ()); 	//!< This function opens the specified file, reads the metadata and the first time slice, and closes the file. 
	int readOneShot(string filename, vector<float> glim = vector <float> ());		//!< This function opens the specified file, reads the first time slice and interpolates data into the variable's grid, closes the file. 
	int writeOneShot(string filename);
	//@}

	/** \name NetCDF input-output streams 
		These functions create input/output "streams" to repeatedly read time slices from one or more files.
		If the lats-lons in the file being read are different from those in the variable, the data is automatically interpolated 
		using the specified regridding method (the default regridding method is bilinear, but can be set using setRegriddingMethod()). 
	*/
	//@{	 
	/** \brief Create an input stream for reading data from one or more files. 
	
		Data is automatically interpolated using the specified regridding method. If data is spread over multiple files, all files must have the same coordinates.
	*/
	int createNcInputStream(vector <string> files, 		//!< A vector containing names of NetCDF files
							vector <float> glim, 		//!< Grid limits, to specify what subset of the data should be read. This should be in the order [west-lon, east-lon, south-lat, north-lat]       
							string rm = "bilinear"		//!< (optional) regridding method
						   ); 
	int createNcOutputStream(string filename);	//!< Create an output stream. Data will be written to file "filename".
	int closeNcInputStream();
	int closeNcOutputStream();

	int readVar_gt(double gt, int mode);	//!< Read time-slice for time gt from the input stream. If the stream comprises of multiple files, this function will automatically find the file which contains the time slice closest to gt and read data from that file.
	int readVar_it(int tid);				//!< Read time-slice from index tid from the NetCDF file currently opened in the input stream.
	int writeVar(int itime); 				//!< Write time-slice to index itime in the output stream
	int overwriteTime(vector <double> &times_new, string tunits);
	//@}

	
	/** \name High level operations using streams
		Before calling these functions, an input stream must be created using createNcInputStream()
	*/
	//@{	 
	int readVar_reduce_mean(double gt1, double gt2); //!< Calculate temporal mean of the variable during time interval gt1-gt2
//	int readVar_reduce_sd(double gt1, double gt2);
	gVar trend(double gt1, double gt2, gVar * P = NULL);		//!< Calculate trend (slope) of the variable during time interval gt1-gt2
	gVar yearlytrend(int year1, int year2, gVar * P = NULL);    //!< 

	gVar percentChange_yoy(int year1, int year2, gVar * P = NULL);    //!< 

	
	gVar trend_gpu(double gt1, double gt2);	//!< Calculate trend (slope) of the variable during time interval gt1-gt2 (Use GPU)
	//@}	
	

	
};

#endif


