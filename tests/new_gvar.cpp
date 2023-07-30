#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <netcdf>
#include <chrono>
#include "../../tensorlib/include/tensor.h"

// print a std::vector via ofstream
// prints: size | v1 v2 v3 ...
template <class T>
std::ostream& operator << (std::ostream &os, const std::vector<T> &v) {
	//os << std::setprecision(12);
	os << v.size() << " | ";
	for (const auto &x : v) {
		os << x << ' ';
	}
	os << '\n';
	return os;
}

class NcFilePP : public netCDF::NcFile {
	public:
	std::multimap <std::string, netCDF::NcVar> vars_map;     // variable name --> NcVar map (for data variables)
	std::map      <std::string, netCDF::NcVar> coords_map;   // variable name --> NcVar map (for coordinate variables)
	std::map      <std::string, std::string> coordunits_map;
	std::map      <std::string, std::vector<float>> coordvalues_map;

	void readMeta(){
		// get all variable in the file in a name --> variable map
		vars_map = this->getVars();
		
		std::map<std::string, netCDF::NcGroup> coords_map_temp = this->getCoordVars();
		// fill coords map from obtained variables
		for (auto p : coords_map_temp){
			coords_map[p.first] = p.second.getVar(p.first);
		} 

		// remove coord vars from all variables, so only data variables are left
		for (auto p : coords_map){
			vars_map.erase(p.first);
		}

		// read values and units for coordinate vars
		for (auto p : coords_map){
			std::string s;
			try{
				p.second.getAtt("units").getValues(s);
				coordunits_map[p.first] = s;
			}
			catch(...){
				std::cout << "Warning: variable " << p.first << " has no unit.";
				coordunits_map[p.first] = "";
			}

			std::vector<float> coordvals(p.second.getDim(0).getSize());
			p.second.getVar(coordvals.data());
			coordvalues_map[p.first] = coordvals;
		}

	}

	void printMeta(){
		std::cout << "NC file >" << "\n";
		std::cout << "   coords:\n";
		for (auto p : coords_map) std::cout << "      " << p.first << " (" << coordunits_map[p.first] << ")\n";
		std::cout << "   vars (excluding coords):\n";
		for (auto p : vars_map) std::cout << "      " << p.first << "\n";
		std::cout << "   coord values:\n";
		for (auto p : coordvalues_map){
			std::cout << "      " << p.first << ": ";
			if (p.second.size() <= 6) std::cout << p.second;
			else{
				std::cout << p.second.size() << " | " << p.second[0] << " " << p.second[1] << " " << p.second[2] << " ... " 
				          << p.second[p.second.size()-3] << " " << p.second[p.second.size()-2] << " " << p.second[p.second.size()-1] << "\n";
			}
		}
		std::cout << "~~\n";
	}
};


template <class T>
class GeoCube : public Tensor<T> {
	public:
	std::string name;       // Variable name
	std::string unit;
	std::vector <std::string> dimnames;
	int unlim_idx = -1;
	int lon_idx, lat_idx;
	float scale_factor = 1.0, add_offset = 0.0; 
	
	std::string tunit = "";
	double tscale = 1;
	std::tm t_base = {};

	// data for standardizing dimension names (kept public to be editable)
	std::vector<std::string>   t_names_try = {"time"};
	std::vector<std::string> lev_names_try = {"lev", "level", "z"};
	std::vector<std::string> lat_names_try = {"lat", "latitude", "y"};
	std::vector<std::string> lon_names_try = {"lon", "longitude", "x"};

	private:
	netCDF::NcVar ncvar;
	std::vector <netCDF::NcDim>  ncdims;
	std::vector <size_t> dimsizes;
	std::vector <std::string> dimunits;

	std::map<std::string, std::string> renaming_map; // map used for renaming dimensions to standard names

	std::vector <size_t> start, count;

	public:
	void readMeta(NcFilePP &in_file, std::string varname = ""){
		// get the named variable, or the first variable in the file
		if (varname != "") ncvar = in_file.vars_map.find(varname)->second;
		else ncvar = in_file.vars_map.begin()->second;

		// get variable name and dimensions
		name = ncvar.getName();
		ncdims = ncvar.getDims();

		// get names and sizes of dimensions for this variable
		for (auto d : ncdims){
			dimnames.push_back(d.getName());
			dimsizes.push_back(d.getSize());
		}

		// ~~ standardize dimension names ~~
		// 1. create renaming map
		for (auto s :   t_names_try) renaming_map[s] = "time";
		for (auto s : lev_names_try) renaming_map[s] = "lev";
		for (auto s : lat_names_try) renaming_map[s] = "lat";
		for (auto s : lon_names_try) renaming_map[s] = "lon";

		for (auto& dname : dimnames){
			// 2. convert to lowecase
			std::transform(dname.begin(), dname.end(), dname.begin(),
						[](unsigned char c){ return std::tolower(c); });
			// 3. rename to standard name
			if (renaming_map.find(dname) != renaming_map.end()){
				dname = renaming_map[dname];
			}
		}

		// ~~ Get the dimension indices. i.e., which index in ncdim std::vector is lat, lon, etc
		lon_idx = std::find(dimnames.begin(), dimnames.end(), "lon") - dimnames.begin();
		lat_idx = std::find(dimnames.begin(), dimnames.end(), "lat") - dimnames.begin();
		if (lon_idx >= dimsizes.size() || lat_idx >= dimsizes.size()) throw std::runtime_error("Lat or Lon not found"); 
		
		// get unlimited dimension id
		for (int i=0; i<ncdims.size(); ++i) if (ncdims[i].isUnlimited()) unlim_idx = i;
	
		start.clear(); start.resize(ncdims.size(), 0);
		count = dimsizes;
		if (unlim_idx >=0) count[unlim_idx] = 1;

		// Get basic variable attributes
		// missing value
		try{ ncvar.getAtt("missing_value").getValues(&this->missing_value);}
		catch(netCDF::exceptions::NcException &e){
			try{ ncvar.getAtt("_FillValue").getValues(&this->missing_value);}
			catch(netCDF::exceptions::NcException &e){ std::cout << "Missing/Fill value not found. Setting to NaN\n";}
		}

		// scale factor
		try{ ncvar.getAtt("scale_factor").getValues(&scale_factor); }
		catch(netCDF::exceptions::NcException &e){ scale_factor = 1.f;}
		
		// offset
		try{ ncvar.getAtt("add_offset").getValues(&add_offset);}
		catch(netCDF::exceptions::NcException &e){ add_offset = 0.f;}

		if (std::find(dimnames.begin(), dimnames.end(), "time") != dimnames.end()) 
			parse_time_unit(in_file);
		else 
			std::cout << "Warning: Variable does not have a time dimension\n";

	}


	void print(bool b_values = false){
		std::cout << "Var (" << name << "):\n";
		std::cout << "   dim names: " << dimnames;
		std::cout << "   dim sizes: " << dimsizes;
		std::cout << "   lat/lon ids = " << lat_idx << " " << lon_idx << "\n";
		std::cout << "   unlimited dimension = ";
		if (unlim_idx < 0) std::cout << "NA\n";
		else std::cout << dimnames[unlim_idx] << " (" << unlim_idx << ")\n";
		std::cout << "   missing value = " << this->missing_value << "\n";
		std::cout << "   scale factor = " << scale_factor << "\n";
		std::cout << "   add offset = " << add_offset << "\n";
		std::cout << "   tbase = " << std::put_time(&t_base, "%Y-%m-%d %H:%M:%S %Z") << "\n";
		std::cout << "   tscale = " << tscale << " (" << tunit << ")" "\n";

		Tensor<T>::print(b_values);
	}

	
	void setMapLimits(float lon_w, float lon_e, float lat_s, float lat_n){

	}

	void setLatLonStarts(size_t lat_start, size_t lon_start){
		start[lat_idx] = lat_start;
		start[lon_idx] = lon_start;
	}

	void setLatLonCounts(size_t lat_count, size_t lon_count){
		count[lat_idx] = lat_count;
		count[lon_idx] = lon_count;
	}

	void readBlock(){
		std::cout << "Resizing tensor to: " << count;
		this->resize(count);
		ncvar.getVar(start, count, this->vec.data());
	}

	private:

	void parse_time_unit(NcFilePP &in_file){
		// parse time units
		std::string since;
		std::stringstream ss(in_file.coordunits_map["time"]);
		ss >> tunit >> since;

		if (since != "since") throw std::runtime_error("time unit is not in the required format (days/hours/etc since yyyy-mm-dd hh:mm:ss)");

		if      (tunit == "days")   tscale = 1;
		else if (tunit == "hours")  tscale = 1.0/24.0;
		else if (tunit == "months") tscale = 1.0*365.2425/12;

		std::stringstream ss1(in_file.coordunits_map["time"]);
		ss1 >> std::get_time(&t_base, std::string(tunit + " since %Y-%m-%d %H:%M:%S").c_str());
	
		// std::cout << std::put_time(&t_base, "%Y-%m-%d %H:%M:%S %Z") << '\n';
		// std::cout << "days since 1970-1-1: " << (std::chrono::duration_cast<std::chrono::hours>(std::chrono::system_clock::from_time_t(std::mktime(&t_base)).time_since_epoch())).count() << "\n";
		// std::cout << "julian days (since -4173/11/24: " << (std::chrono::duration_cast<std::chrono::hours>(std::chrono::system_clock::from_time_t(std::mktime(&t_base)).time_since_epoch()) + std::chrono::hours{58574100}).count() << "\n";
	}



};


int main(){

	NcFilePP in_file;
	in_file.open("data/gpp.2000-2015.nc", netCDF::NcFile::read);
	in_file.readMeta();
	in_file.printMeta();

	GeoCube<float> v;
	v.readMeta(in_file);
	v.setLatLonStarts(230, 161);
	v.setLatLonCounts(1, 2);

	v.readBlock();
	v.print(true);

	return 0;
} 

