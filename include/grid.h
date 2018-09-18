#ifndef GRID_H
#define GRID_H

#include <vector>
#include "constants.h"

using namespace std;

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~  GRID    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
// create coordinate vector given x0, xf and number of points
vector <float> createCoord(float x0, float xf, int nx, float &dx);

// create coordinate vector given x0, xf and resolution
vector <float> createCoord(double x0, double xf, double dx, int &nx);	// x0 and xf are centres of 1st and last gridbox
vector <float> createCoord_from_edges(double x0, double xf, double dx, int &nx); // x0 and xf are edges of 1st and last gridbox

// print a gridded variable along with 2d coordinates defined on x, y
void printVar(vector <float> &x, vector <float> &y, float * data);


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
			   float * data, float missingVal = std_missing_value);

// return bilinear interpolated value at point (x,y) in grid {lons,lats}
// use bilIndices instead of finding them in situ
float bilinear(int ilat, int ilon, int iz, vector <int> &indices,
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &mlons, vector <float> &mlats,
			   float * data, float missingVal = std_missing_value);

// return value at gridcell containing point (x,y) in grid {lons,lats}
float cellVal(float x, float y, float iz, 
			   vector <float> &lons, vector <float> &lats, 
			   float * data, float missingVal = std_missing_value);

// return value at gridcell containing point (x,y) in grid {lons,lats}
// use bilIndices instead of finding them in situ
float cellVal(int ilat, int ilon, int iz, vector <int> &indices,
			   vector <float> &lons, vector <float> &lats, 
			   vector <float> &mlons, vector <float> &mlats,
			   float * data, float missingVal = std_missing_value);

// mask variable v using m as mask
// for all points in v, if (m <= val) v = missing_value
gVar mask(gVar &v, gVar &m, float val = 0);

// *** Following regridding functions copy METADATA (except streans, and with new lats/lons)
// interpolate variable v onto given grid (xlons, xlats)
gVar lterp(gVar &v, vector <float> &xlons, vector <float> &xlats); 

// *** Following regridding functions copy ONLY VALUES.. 
// .. they assume (and check) compatibility
// .. input and output coords can be taken directly from input and output gVars
int lterpCube(gVar &v, gVar &out, vector <int> &indices);
int cellRegridCube(gVar &v, gVar &out, vector <int> &indices);

gVar coarseGrain_sum(gVar &hires, vector <float> &xlons, vector <float> &xlats);
gVar coarseGrain_mean(gVar &hires, vector <float> &xlons, vector <float> &xlats);

gVar binary(gVar v, float thresh=0);

#endif
