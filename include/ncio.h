#ifndef NCIO_H
#define NCIO_H

#include <netcdfcpp.h>
#include "constants.h"
#include "gvar.h"
using namespace std;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  NCIO    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

//	NcFile_handle class has resources to read/write Nc Files.
//	this class has no variables by itself and must input all relevant gVars
//	therefore destructor frees all pointers.
/** \brief A handle for NetCDF files
*/
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
	int ilonf, ilatf;
	// NOTE: above indices are defined on the native arrays and not reversed ones
	bool mplimited; // if true, map limits are set
	
	int firstVarID;
	
	// constructor. Only initializes variables. Does not open file.
	NcFile_handle();	// open NC object (nc file)
	
	void setMapLimits(float xwlon, float xelon, float xslat, float xnlat);	// set lat-lon limits
	int  open(string s, string m, const float glimits[4]); // open file s
	int close(); // close file
	~NcFile_handle();	// dFile must be deleted in destructor

	// reading functions
	int getMeta();
	int readCoordData(gVar &v);
	int readTime(gVar &v);
	
	int readCoords(gVar &v); // read coord metadata, read values if rr is true
	int getVarID(string varname);	// get variable id (ivar) from name
	int readVarAtts(gVar &v, int ivar = -1); // get variable name, units, missing_value etc.
	int readVar(gVar &v, int itime, int iVar = -1); // read variable with index iVar from file into "v"
	int readVar_gt(gVar &v, double gt, int mode, int iVar = -1); // read variable corresponding to GT gt
		// as specified by mode: 0 = hold values from previous step, 1 = interpolate
	int readVar_parallel(gVar &v, int itime, int iVar = -1);
	
	// writing functions
	int writeCoords(gVar &v, bool wr = true); // write coord metadata, write values if wr is true
	NcVar * createVar(gVar &v); // create variable v into file matching gVar
	int writeVar(gVar &v, NcVar * vVar, int itime); // write values into vVar at time itime
	int writeTimeValues(gVar &v); // write time values
};

#endif

