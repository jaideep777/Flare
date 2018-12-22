#include <iostream>
#include <vector>
#include <map>
#include <netcdf>
using namespace std;
using namespace netCDF;

int main(){
	
	NcFile f("data/gpp.intercept.2001-2010.nc", NcFile::read);
	if (!f.isNull() ) cout << "Success\n";
	
	NcVar lonVar = f.getVar("lon");
	if (lonVar.isNull()) cout << "Found lon\n";
	
	vector <NcDim> dims = lonVar.getDims();
	cout << "var has " << dims.size() << " dims\n";
	
	for (auto dim : dims) cout << dim.getSize() << " ";
	cout << "\n"; 	

	cout << "var has " << lonVar.getDim(0).getSize() << " values\n";

	vector <float> lons(dims[0].getSize());
	lonVar.getVar(lons.data());
	
	string s;
	lonVar.getAtt("units").getValues(s);
	cout << s << endl;
	
	multimap<string,NcVar> vars_map = f.getVars();
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it) cout << it->first << " "; cout << endl;
	
	map<string,NcGroup> coords_map = f.getCoordVars();
	for (auto it = coords_map.begin(); it != coords_map.end(); ++it) cout << it->first << " "; cout << endl;
	
//	for (auto it : coords_map) cout << it.first << " ";
	
	
	// Get the coordinate variables
	NcFile* dFile = &f;
	NcVar latVar, levVar, tVar;
	int ncoords = 0;
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it ){
		
		string var = it->first;
		
		if (it->first == "lat" || it->first == "latitude" || it->first == "LAT" ){
			latVar = it->second;
			++ncoords;
		}
//		else cout << "lat not found!\n";

		if (it->first == "lon" || it->first == "longitude" || it->first == "LON" ){
			lonVar = it->second;
			++ncoords;
		}
//		else cout << "lon not found!\n";


		if (it->first == "lev" || it->first == "levels" || it->first == "LEV" ){
			levVar = it->second;
			++ncoords;
		}
//		else cout << "lev not found!\n";


		if (it->first == "time" || it->first == "TIME" || it->first == "LEV" ){
			tVar = it->second;
			++ncoords;
		}
//		else cout << "time not found!\n";

	}
	
	if (latVar.isNull()) cout << "Lat not found\n";
	if (lonVar.isNull()) cout << "Lon not found\n";
	if (levVar.isNull()) cout << "Lev not found\n";
	if (tVar.isNull())   cout << "Time not found\n";
	
	
	NcVar firstVar;
	for (auto it = vars_map.begin(); it != vars_map.end(); ++it ){
		if (it->second.getDimCount() == ncoords) {firstVar = it->second; break;}
	}
	cout << firstVar.getName() << " (" << firstVar.getId() << ")";
	
	string a;
	NcVarAtt A = firstVar.getAtt("units");
	if (!A.isNull()) A.getValues(a);
	else a = "NA";
	cout << "var unit = " << a << endl;
	
	// missing value
	float b;
	A = firstVar.getAtt("missing_value");
	if (!A.isNull()) A.getValues(&b);
	else{
		A = firstVar.getAtt("_FillValue");
		if (!A.isNull()) A.getValues(&b);
	}
	cout << "missing Value = " << b << endl;

	
	return 0;
} 
