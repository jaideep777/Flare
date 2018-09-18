#include <cmath>
#include <fstream>
#include <cstdlib>
#include <string>
#include <curand_kernel.h>
using namespace std;

#include "../include/turbulence.h"
//#include "../include/graphics.h"
#include "cuda_vector_math.cuh"
#include "cuda_device.h"
#include "../include/histogram.h"


extern curandGenerator_t generator_host;


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// KERNEL to set up RANDOM GENERATOR STATES
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
__global__ void teRngStateSetup_kernel(int * rng_Seeds, curandState * rngStates, int nx, int ny){
	//int tid = threadIdx.x;							// each block produces exactly the same random numbers
	int tid_u = threadIdx.x + blockIdx.x*blockDim.x;	// each block produces different random numbers
	if (tid_u >= nx*ny) return;
	
	curand_init (rng_Seeds[tid_u], 0, 0, &rngStates[tid_u]);
}

#define TE_PP_SEED 777

void TurbulenceEngine::initRNG(){
	srand(TE_PP_SEED);
	for (int i=0; i<nx*ny; ++i) te_seeds_h[i] = rand(); 
	cudaMemcpy( te_seeds_dev, te_seeds_h, sizeof(int)*nx*ny, cudaMemcpyHostToDevice);

	int nt = 256, nb = (nx*ny-1)/nt+1;
	teRngStateSetup_kernel <<< nb, nt>>> (te_seeds_dev, te_dev_XWstates, nx, ny);
	getLastCudaError("RNG_kernel_launch");
}

// =========================================================================================
//
//		Turbulence Kernels!!
//
// =========================================================================================


// =========================================================================================
//		Generate conjugate symetric noise matrix
// =========================================================================================
__global__ void te_generateNoise_kernel(float2* Zmat, int nx, int ny, curandState * rngStates){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nx*(ny/2+1)) return;
	
	int k = tid % nx;
	int m = tid / nx;	

//	if (m >= ny/2+1) return;
	
	float2 z1 = curand_normal2(&rngStates[tid]);
	Zmat[m*nx+k] = z1;

	// conjugate symetric element
	int kc = (nx-k)%nx, mc = (ny-m)%ny;
	Zmat[mc*nx+kc] =  make_float2(z1.x, -z1.y);

}

void TurbulenceEngine::generateNoise_gpu(){
	int nt = 256, nb = (nx*(ny/2+1)-1)/nt+1;
	te_generateNoise_kernel <<< nb, nt>>> (Z_dev, nx, ny, te_dev_XWstates);
	cudaMemcpy(Z, Z_dev, nx*ny*sizeof(float2), cudaMemcpyDeviceToHost);	

}


// =========================================================================================
//		evolve Psi in fourier domain for 1 time step
// =========================================================================================
__global__ void modifyPsi_kernel(float2 *Z_d, float2 *Psi_d, float *lambda_d, int nx, int ny,
								 float xi, float nu, float dt){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= nx*ny) return;
	
	int k = tid % nx;
	int m = tid / nx;	

	int kx = (k > nx/2)? k-nx : k;
	int ky = (m > ny/2)? m-ny : m;

	float normk = kx*kx+ky*ky;
	float sqrt_term;
	if (normk != 0){
		float sqrt_arg = xi*lambda_d[m*nx+k]*(1-exp(-2*nu*normk*dt))/2/nu/normk;
		if (sqrt_arg <0) sqrt_arg = 0;
		sqrt_term = sqrt(sqrt_arg);
	}
	else{ 
		sqrt_term = lambda_d[0]*sqrt(xi*dt);
	}	
	Psi_d[m*nx+k] = Psi_d[m*nx+k]*exp(-nu*normk*dt) + Z_d[m*nx+k]*sqrt_term;

}


void TurbulenceEngine::modifyPsi_gpu(){
	int nt = 256, nb = (nx*ny-1)/nt+1;
	modifyPsi_kernel <<<nb, nt >>> (Z_dev, Psi_dev, lambda_dev, nx, ny, xi, nu, dt);
//	cudaMemcpy(Psi, Psi_dev, nx*ny*sizeof(float2), cudaMemcpyDeviceToHost);
}


// =========================================================================================
//		calculate velocity field
// =========================================================================================
__global__ void calcVelField_kernel(float2* psi_d, float2* velfield, float L, int nx, int ny){
	int tid = blockIdx.x*blockDim.x + threadIdx.x;
	if (tid >= (nx-1)*(ny-1)) return;
	
	int ix = tid % (nx-1);
	int iy = tid / (nx-1);	
	
	velfield[iy*(nx-1)+ix].x =  (psi_d[(iy+1)*nx+ ix   ].x - psi_d[iy*nx+ix].x)/(2*L/nx);
	velfield[iy*(nx-1)+ix].y = -(psi_d[    iy*nx+(ix+1)].x - psi_d[iy*nx+ix].x)/(2*L/ny);
	
}

void TurbulenceEngine::calcVelocityField(){
	// calculate velocity field
	int nt = 256, nb = ((nx-1)*(ny-1)-1)/nt+1;
	calcVelField_kernel <<< nb, nt>>> (psi_dev, vel_field_dev, L, nx, ny);

//	cudaMemcpy(vel_field, vel_field_dev, (nx-1)*(ny-1)*sizeof(float2), cudaMemcpyDeviceToHost);
//	float u[(nx-1)*(ny-1)];
//	for (int i=0; i<(nx-1)*(ny-1); ++i) u[i] = sqrt(vel_field[i].x*vel_field[i].x + vel_field[i].y*vel_field[i].y);
//	float u_max = u[0], u_sum = u[0];
//	for (int i=1; i<(nx-1)*(ny-1); ++i) u_max = fmax(u_max, u[i]);
//	for (int i=1; i<(nx-1)*(ny-1); ++i) u_sum = u_sum + u[i];
//	float u_avg = u_sum/(nx-1)/(ny-1);
//	cout << u_max << " " << u_avg << endl;
}


// =========================================================================================
//
//		Turbulence Engine functions
//
// =========================================================================================


void TurbulenceEngine::init(Initializer &I){
	nx = I.getScalar("nx"); 
	ny = I.getScalar("ny");
	mu = I.getScalar("mu_t"); 
	nu = I.getScalar("nu_t"); 
	xi = I.getScalar("xi_t");
	lambda0 = I.getScalar("lambda0_t");
	dt = I.getScalar("dt_t");
	L = I.getScalar("L");
	xmin = 0; //-I.getScalar("arenaSize")/2;
	xmax = L; // I.getScalar("arenaSize")/2;
	ymin = 0; //-I.getScalar("arenaSize")/2;
	ymax = L; // I.getScalar("arenaSize")/2;

	nlevCol = 64;

	lambda = new float[nx*ny];
	Z = new cufftComplex[nx*ny];
	Psi = new cufftComplex[nx*ny];
	psi = new cufftComplex[nx*ny];
	vel_field = new float2[(nx-1)*(ny-1)];
	te_seeds_h = new int[nx*ny];

	cudaMalloc((void**)&lambda_dev, sizeof(float)*nx*ny);
	cudaMalloc((void**)&Z_dev, sizeof(cufftComplex)*nx*ny);
	cudaMalloc((void**)&Psi_dev, sizeof(cufftComplex)*nx*ny);
	cudaMalloc((void**)&psi_dev, sizeof(cufftComplex)*nx*ny);
	cudaMalloc((void**)&vel_field_dev, sizeof(float2)*(nx-1)*(ny-1));
	cudaMalloc((void**)&te_seeds_dev, nx*ny*sizeof(int));
	cudaMalloc((void**)&te_dev_XWstates, nx*ny*sizeof(curandState));
	getLastCudaError("alloc GPU arrays");

	// init RNG
	initRNG();

	// prepare to transform
	cout << "creating FFT plan config..." << endl;
	cufftPlan2d(&plan, nx, ny, CUFFT_C2C);
	getLastCudaError("create plan");

}


void TurbulenceEngine::freeMemory(){
	delete [] lambda;
	delete [] Z;
	delete [] Psi;
	delete [] psi;
	delete [] vel_field;
	delete [] te_seeds_h;

	cudaFree(lambda_dev);
	cudaFree(Z_dev);
	cudaFree(Psi_dev);
	cudaFree(psi_dev);
	cudaFree(vel_field_dev);
	cudaFree(te_seeds_dev);
	cudaFree(te_dev_XWstates);

}

void TurbulenceEngine::generateSpectrum(){
	cout << "Generating spectrum..." << endl;
//	ofstream fout("lambda.txt");
	for (int k=0; k<nx; ++k){
		for (int m=0; m<ny; ++m){
			int x = (k > nx/2)? k-nx : k;
			int y = (m > ny/2)? m-ny : m;
			lambda[m*nx+k] = lambda0*exp(-mu*sqrt(x*x+y*y));
//			fout << lambda[m*nx+k] << " ";
		}
//		fout << "\n";
	}
//	fout.close();
	
	cudaMemcpy(lambda_dev, lambda, nx*ny*sizeof(float), cudaMemcpyHostToDevice);
}


void TurbulenceEngine::calcEquilPsi(){
	// initial Psi
	for (int k=0; k<nx; ++k){
		for (int m=0; m<ny; ++m){
			int kx = (k > nx/2)? k-nx : k;
			int ky = (m > ny/2)? m-ny : m;

			float normk = kx*kx+ky*ky;
			if (normk != 0) {
				float num = sqrt(xi*lambda[m*nx+k]*(1-exp(-2*nu*normk*dt))/2/nu/normk);
				float den = 1-exp(-nu*normk*dt);
				Psi[m*nx+k] = Z[m*nx+k]*num/den;
			}
			else Psi[m*nx+k] = make_float2(0,0); //Z[m*nx+k];
//			cout << Psi[m*nx+k].x << " ";

			Psi[m*nx+k] = make_float2(1,0);
		}
//		cout << "\n";
	}	

	cudaMemcpy(Psi_dev, Psi, nx*ny*sizeof(float2), cudaMemcpyHostToDevice);
}


void TurbulenceEngine::transformPsi(){
	cufftExecC2C(plan, Psi_dev, psi_dev, CUFFT_INVERSE);
	getLastCudaError("fft");

	cudaDeviceSynchronize();
}


void TurbulenceEngine::normalize_psi(){
	cudaMemcpy(psi, psi_dev, nx*ny*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
	
	float tmin=psi[0].x, tmax=psi[0].x, tsum=0;
	for (int i=0; i<nx; ++i){
		for (int j=0; j<ny; ++j){
			tmin = min(tmin, psi[j*nx+i].x);
			tmax = max(tmax, psi[j*nx+i].x);
			tsum += psi[j*nx+i].x;
		}
	}
	float tmean = tsum/nx/ny;
	cout << tmin << " " << tmean << " " << tmax << endl;

	for (int i=0; i<nx; ++i){
		for (int j=0; j<ny; ++j){
			psi[j*nx+i] = (psi[j*nx+i] - tmean)/(tmax-tmin)*2;
		}
	}
	cudaMemcpy(psi_dev, psi, nx*ny*sizeof(cufftComplex), cudaMemcpyHostToDevice);
}


void TurbulenceEngine::update(){
	generateNoise_gpu();
	modifyPsi_gpu();
	transformPsi();
	calcVelocityField();
}

void TurbulenceEngine::updateColorMap(){
//	cudaMemcpy(psi, psi_dev, nx*ny*sizeof(cufftComplex), cudaMemcpyDeviceToHost);
//	float * tmp = new float[nx*ny];
//	for (int i=0; i<nx*ny; ++i) tmp[i] = psi[i].x;
//	glRenderer->setCmapColorBufferData(tmp, nx*ny, nlevCol);
//	//glutPostRedisplay();
//	delete [] tmp;
}



void TurbulenceEngine::printMap(string mapname, ofstream &fout){
	float * dat = new float[nx*ny];
	for (int iy=0; iy<ny; ++iy){
		for (int ix=0; ix<nx; ++ix){
			if (mapname == "psi") fout << psi[iy*nx+ix].x << " ";
			dat[iy*nx+ix] = psi[iy*nx+ix].x;
		}
		fout << '\n'; 
	}
	fout.flush();

	printSummary(dat, nx*ny);
	delete [] dat;
}


// Particle velocities due to turbulence

//	for (int t=0; t<1; ++t){
//		int ipx = (nx-1)*(ppos.x-xmin)/(xmax-xmin);
//		int ipy = (ny-1)*(ppos.y-ymin)/(ymax-ymin);
//		ppos += vel_field[ipy*(nx-1)+ipx]*dt*0.3;
//		makePeriodic(ppos.x, xmin, xmax);
//		makePeriodic(ppos.y, ymin, ymax);
//		pos_fout << ppos.x << " " << ppos.y << '\n';
//	}



