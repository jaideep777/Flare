#include <iostream>
#include <gsm.h>
#include <netcdfcpp.h>
#include <vector>
#include <algorithm>
#include <cuda_runtime.h>
using namespace std;

// g++ -I/usr/local/netcdf-c-4.3.2/include -I/usr/local/netcdf-cxx-legacy/include -I/home/jaideep/codes/FIRE_CODES/libgsm_v2/include -L/home/jaideep/codes/FIRE_CODES/libgsm_v2/lib -L/usr/local/netcdf-cxx-legacy/lib -o 1 trend_test.cpp -l:libgsm.so.2 -lnetcdf_c++ 

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true){
   if (code != cudaSuccess){
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


__global__ void init_residuals_kernel(float * sxx, float * syy, float * sxy, float * sx, float * sy, int count, int nvals){
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < nvals){
		sxx[i] = sxy[i] = syy[i] = sx[i] = sy[i] = 0;
	}
}

__global__ void update_residuals_kernel(float * vals, float * sxx, float * syy, float * sxy, float * sx, float * sy, int count, int nvals){
	
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < nvals){
	
		sxy[i] += vals[i]*count;
		sx[i] += count;
		sy[i] += vals[i];
		sxx[i] += count*count;
		syy[i] += vals[i] * vals[i];
	
	}
	
}

__global__ void calc_metrics_kernel(float * b1, float * t, float * sxx, float * syy, float * sxy, float * sx, float * sy, int count, int nvals){
	
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	
	if (i < nvals){
		float Sxy = sxy[i] - sx[i]*sy[i]/count; 
		float Sxx = sxx[i] - sx[i]*sx[i]/count;
		b1[i] = Sxy/Sxx;	

		float Syy = syy[i] - sy[i]*sy[i]/count;
		float s = sqrt((Syy - b1[i]*Sxy)/(count-2));

		t[i] = b1[i]/(s/sqrt(Sxx));
	}
	
}


gVar gVar::trend_gpu(double gt1, double gt2){

	clock_t start, end;
	double msecs;
	start = clock();

	int nvals = nlons*nlats*nlevs;
	int blockSize = 128;
	int nblocks = ((nvals-1)/blockSize) + 1;

	float * sxx_dev, * sxy_dev, * syy_dev, *sx_dev, *sy_dev;
	float * b1_dev, *t_dev;
	float * var_dev;

	cudaMalloc((void**)& var_dev, nvals*sizeof(float));
	
	cudaMalloc((void**)& sxx_dev, nvals*sizeof(float));
	cudaMalloc((void**)& sxy_dev, nvals*sizeof(float));
	cudaMalloc((void**)& syy_dev, nvals*sizeof(float));
	cudaMalloc((void**)& sx_dev, nvals*sizeof(float));
	cudaMalloc((void**)& sy_dev, nvals*sizeof(float));
	
	cudaMalloc((void**)& b1_dev, nlons*nlats*nlevs*sizeof(float));
	cudaMalloc((void**)& t_dev,  nlons*nlats*nlevs*sizeof(float));
	
	gVar temp; temp.copyMeta(*ipvar);
	temp.fill(0);
	int count = 0;
	
	gVar b1 = temp, 
//		 s = temp,
		 t = temp; 
	
	updateInputFile(gt1);	// this will give correct OR one previous file
	while (gt1 > ipvar->ix2gt(ipvar->times.size()-1)){	// increment curr_file as long as gt1 +outside file range
		++curr_file;
		loadInputFileMeta();
	}	
	
	init_residuals_kernel <<< nblocks, blockSize >>> (sxx_dev, syy_dev, sxy_dev, sx_dev, sy_dev, 0, nvals);						

	CDEBUG << "readVar_reduce_mean (" << varname << ") :" << gt2string(gt1) << " " << gt2string(gt2) << endl;
	while(1){
//		cout << (ipvar->times[0]) << " " << (ipvar->times[ipvar->times.size()-1]) << " "  << (gt1-ipvar->tbase)*24.0/ipvar->tscale <<  endl;
		int tstart = lower_bound(ipvar->times.begin(), ipvar->times.end(), (gt1-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin();	   // first elem >= gt1 
		int tend   = upper_bound(ipvar->times.begin(), ipvar->times.end(), (gt2-ipvar->tbase)*24.0/ipvar->tscale) - ipvar->times.begin() -1;   // last elem <= gt2
//		cout << gt2string(ipvar->ix2gt(tstart)) << " " << gt2string(ipvar->ix2gt(tend)) << endl;

		if (tend < 0) break;

		for (int i=tstart; i<=tend; ++i){ 
			clock_t start1, end1;
			start1 = clock();
			ifile_handle->readVar(*ipvar, i, ipvar->ivar1);	// readCoords() would have set ivar1
			end1 = clock();
			cout << "time to read: " << ((double) (end1 - start1)) * 1000 / CLOCKS_PER_SEC << " ms" << endl;
			
			cudaMemcpy(var_dev, &(ipvar->values[0]), ipvar->values.size()*sizeof(float), cudaMemcpyHostToDevice);
			gpuErrchk(cudaGetLastError());

			update_residuals_kernel <<< nblocks, blockSize >>> (var_dev, sxx_dev, syy_dev, sxy_dev, sx_dev, sy_dev, count, nvals);						
			gpuErrchk(cudaGetLastError());
						
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
	
	calc_metrics_kernel <<< nblocks, blockSize >>> (b1_dev, t_dev, sxx_dev, syy_dev, sxy_dev, sx_dev, sy_dev, count, nvals);						
	
	cudaMemcpy(&(b1.values[0]), b1_dev, nvals*sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(&(t.values[0]),  t_dev,  nvals*sizeof(float), cudaMemcpyDeviceToHost);

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
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-01-01.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-01-09.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-01-17.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-01-25.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-02-02.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-02-10.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-02-18.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-02-26.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-03-06.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-03-14.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-03-22.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-03-30.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-04-07.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-04-15.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-04-23.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-05-01.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-05-09.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-05-17.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-05-25.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-06-02.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-06-10.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-06-18.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-06-26.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-07-04.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-07-12.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-07-20.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-07-28.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-08-05.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-08-13.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-08-21.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-08-29.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-09-06.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-09-14.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-09-22.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-09-30.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-10-08.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-10-16.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-10-24.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-11-01.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-11-09.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-11-17.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-11-25.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-12-03.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-12-11.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-12-19.nc",
		"/media/jaideep/HD-B1/LAI/LAI_new/MOD15A2H.A_LAI-2003-12-27.nc"
	};

	vector <string> infiles(files, files+46); 
	gVar hires;
	hires.initMetaFromFile(infiles[0]);
	hires.createNcInputStream(infiles, glim);
	hires.printGrid();

	gVar slope = hires.trend_gpu(ymd2gday("2003-1-1"), ymd2gday("2003-12-31"));
	
	slope.writeOneShot("npp.b1.nc");
	
	return 0;

}


