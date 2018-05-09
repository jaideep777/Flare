#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
using namespace std;

#include <thread>

// g++ -Wl,--no-as-needed -std=c++11 -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 trend_parallel_io_compute.cpp -pthread -l:libgsm.so.2 -lnetcdf_c++ 

void readVar_threaded(NcFile_handle * h, gVar * g, int i){
	h->readVar(*g, i, g->ivar1);
}

gVar gVar::trend_par(double gt1, double gt2){

	
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
			thread t1(readVar_threaded, ifile_handle, ipvar, i);
			t1.join();
			
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




int main(){
	
	// ~~~~~~ Some NetCDF Essentials ~~~~~~~~
	// set NETCDF error behavior to non-fatal
	NcError err(NcError::silent_nonfatal);
	
	// specify log file for gsm
	ofstream gsml("gsm_log.txt");
	gsm_log = &gsml;

	// create a grid limits vector for convenience
	float glimits[] = {0, 360, -90, 90};
	vector <float> glim(glimits, glimits+4);
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~

	// create the coordinates for our georeferenced variable
//	int nlons, nlats, nlevs, ntimes;
//	vector <float> lons = createCoord(-180+.5/2,180-0.5/2,0.5,nlons);
//	vector <float> lats = createCoord(-90+.5/2,90-.5/2,0.5,nlats);
//	vector <float> levs = createCoord(1,1,1,nlevs);
//	vector <double> times(16*24); 
//	for (int i=0; i<times.size(); ++i) times[i]= ymd2gday("2001-1-1")+ i*365.2524/24 + 6 - ymd2gday("2000-1-1");

	string files[] = 
	{
		"/media/jaideep/WorkData/Fire_G/GPP_modis/gpp.2000-2015.nc",
	};

	vector <string> infiles(files, files+1); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.createNcInputStream(infiles, glim);
	hires.printGrid();

	gVar slope = hires.trend(ymd2gday("2000-1-1"), ymd2gday("2015-12-31"));
	
	slope.writeOneShot("npp.slope2.nc");
	
	return 0;

}


