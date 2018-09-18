#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <curand_kernel.h>
#include <thrust/reduce.h>
#include <thrust/device_ptr.h>

using namespace std;

#include "../include/resource.h"
//#include "../headers/graphics.h"
#include "cuda_vector_math.cuh"
#include "cuda_device.h"



// extern curandGenerator_t generator_host;

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// KERNEL to set up RANDOM GENERATOR STATES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//__global__ void teRngStateSetup_kernel(int * rng_Seeds, curandState * rngStates, int nx, int ny){
//	//int tid = threadIdx.x;							// each block produces exactly the same random numbers
//	int tid_u = threadIdx.x + blockIdx.x*blockDim.x;	// each block produces different random numbers
//	if (tid_u >= nx*ny) return;
//	
//	curand_init (rng_Seeds[tid_u], 0, 0, &rngStates[tid_u]);
//}

//#define TE_PP_SEED time(NULL)

//void ResourceGrid::initRNG(){
//	srand(TE_PP_SEED);
//	for (int i=0; i<nx*ny; ++i) te_seeds_h[i] = rand(); 
//	cudaMemcpy( te_seeds_dev, te_seeds_h, sizeof(int)*nx*ny, cudaMemcpyHostToDevice);

//	int nt = 256, nb = (nx*ny-1)/nt+1;
//	teRngStateSetup_kernel <<< nb, nt>>> (te_seeds_dev, te_dev_XWstates, nx, ny);
//	getLastCudaError("RNG_kernel_launch");
//}



// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//          RESOURCE GRID
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void ResourceGrid::init(Initializer &I){
	nx = I.getScalar("nx"); 
	ny = I.getScalar("ny");
	D  = I.getScalar("D");
	dt = I.getScalar("dt");
	L  = I.getScalar("L");

	graphics = I.getScalar("graphicsQual")>0;

	res = new float[nx*ny];
	cudaMalloc((void**)&res_dev, sizeof(float)*nx*ny);
	cudaMalloc((void**)&res_new_dev, sizeof(float)*nx*ny);


	r = new float[nx*ny];
	K = new float[nx*ny];

	float r0 = I.getScalar("r");
	float K0 = I.getScalar("K");
	for (int i=0; i<nx*ny; ++i){
		r[i] = r0;
		K[i] = K0;
	}
	cudaMalloc((void**)&r_dev, sizeof(float)*nx*ny);
	cudaMalloc((void**)&K_dev, sizeof(float)*nx*ny);

	cudaMemcpy(r_dev, r, nx*ny*sizeof(float), cudaMemcpyHostToDevice);	
	cudaMemcpy(K_dev, K, nx*ny*sizeof(float), cudaMemcpyHostToDevice);	
	

	for (int i=0; i<nx*ny; ++i) res[i]=K0;
	res[ix2(128,128,256)] = K0;
	
	cudaMemcpy(res_dev, res, nx*ny*sizeof(float), cudaMemcpyHostToDevice);

//	if (graphics){
//		// create resource grid color-map
//		res_shape = ColorMap("res", false, 100, nx, 0, L);
//		float2 cmap_pos[res_shape.nVertices];
//		res_shape.createGridCentres(cmap_pos); 
//		res_shape.createShaders();
//		res_shape.createVBO(cmap_pos, res_shape.nVertices*sizeof(float2));	
//		res_shape.createColorBuffer();
//		res_shape.updateColors(res, nx*ny);
//	}

	cout << "total resource = " << sumResource() << endl;

}


void ResourceGrid::freeMemory(){
	delete [] res;
	delete [] r;
	delete [] K;
	cudaFree(res_dev);
	cudaFree(res_new_dev);
	cudaFree(r_dev);
	cudaFree(K_dev);
	
//	if (graphics){
//		res_shape.deleteShaders();
//		res_shape.deleteVBO();
//	}
}


//void ResourceGrid::graphics_updateArrays(){
//	cudaMemcpy(res, res_dev, nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
//	res_shape.updateColors(res, nx*ny, 0, 50);
//}



// =========================================================================================
//
//		Resource dynamics Kernels
//
// =========================================================================================


__global__ void diffusion_kernel(float * res, float * res_new, float D, int nx, int ny, float dt){

	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nx*ny) return;

	int ix = tid % nx;
	int iy = tid / nx;
	
	float grad_x = res[ix2(makePeriodicID(ix+1,nx),                iy,       nx)]
				 + res[ix2(makePeriodicID(ix-1,nx),                iy,       nx)]
				 - res[ix2(               ix,                      iy,       nx)]*2;
	float grad_y = res[ix2(               ix,       makePeriodicID(iy+1,ny), nx)]
				 + res[ix2(               ix,       makePeriodicID(iy-1,ny), nx)]
				 - res[ix2(               ix,                      iy,       nx)]*2;

	res_new[tid] = res[tid] + (D*grad_x+D*grad_y)*dt;	

}


void ResourceGrid::diffuse(){
	int nt = 256; int nb = (nx*ny-1)/nt + 1;

	diffusion_kernel <<<nb, nt>>> (res_dev, res_new_dev, D, nx, ny, dt);
	cudaMemcpy(res_dev, res_new_dev, nx*ny*sizeof(float), cudaMemcpyDeviceToDevice);
}


__global__ void resource_growth_kernel(float * res, float *r, float *Ke_all, float *K, float dt, int nx, int ny){
	
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nx*ny) return;

	float R = res[tid];
	R += ( r[tid]*R*(1-R/K[tid]) - Ke_all[tid]*R )*dt;  // resource is extracted after growth
	res[tid] = fmaxf(1e-6, R);	// resource should not go negative
}

void ResourceGrid::grow(float * ke_all_dev){
	int nt = 256; int nb = (nx*ny-1)/nt + 1;
	resource_growth_kernel <<<nb, nt>>> (res_dev, r_dev, ke_all_dev, K_dev, dt, nx, ny);
	getLastCudaError("resource growth kernel");
}


float ResourceGrid::sumResource(){
	thrust::device_ptr <float> arr_dev(res_dev);
	totalRes = thrust::reduce(arr_dev, arr_dev+nx*ny);
	return totalRes;
}



//void ResourceGrid::update(){
//	grow();
//	//diffuse();
//}



