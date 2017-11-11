#include "gpu32.h"
#include "helper_cuda.h"

#define BLOCKSIZE 8

#define MAX_MAT_SIZE 2048
#define MAX_FREQS 1000
#define MAX_INTERFACES 10
#define MAX_WNUMS 100

#ifndef ORD_RICH
  #define ORD_RICH 1
#endif

#ifndef FTYPE	
  #define ftype float
#else
  #define ftype FTYPE
#endif

#include "residual/bisect_gpu.h"

//#define m_CopyToGPU(d, s,  bytes){ if(cudaSuccess != cudaMalloc((void**) &d, bytes)) boinc_temporary_exit(60); checkCudaErrors(cudaMemcpy((void*) d, (void*)s, bytes, cudaMemcpyHostToDevice)); }

#define m_CopyToGPU2(s, elements, type)\
	type* g_##s;\
	checkCudaErrors(cudaMalloc((void**) &g_##s, elements*sizeof(type)));\
	checkCudaErrors(cudaMemcpy((void*) g_##s, (void*)s, elements*sizeof(type), cudaMemcpyHostToDevice));

#define m_FreeHostAndGPU(s)\
	free(s);\
	checkCudaErrors(cudaFree(g_##s));

#define m_FreeGPU(s)\
	checkCudaErrors(cudaFree(g_##s));

__device__ void FillLocalArrays (
		const int tid,
		const ftype cb,
		const ftype rhob, 
		const int batch_sz, 
		const int cws_sz, 
		const ftype* cws, 
		ftype rhos[],
		ftype c1s[],
		ftype c2s[])
{
	// FIXME: magic constant!
	for (int i = 0; i < cws_sz +1; ++i)
		c1s[i] = 1500;

	for (unsigned i = 0; i < cws_sz - 1; i++)
	{
		c1s[i] = cws[i*batch_sz + tid];
		c2s[i] = cws[(i + 1)*batch_sz + tid];
		rhos[i] = 1;
	}

	c1s[cws_sz - 1] = cws[(cws_sz - 1)*batch_sz + tid];
	c2s[cws_sz - 1] = cws[(cws_sz - 1)*batch_sz + tid];
	rhos[cws_sz - 1] = 1;

	c1s[cws_sz] = cb;
	c2s[cws_sz] = cb;
	rhos[cws_sz] = rhob;
}

__device__ void FillDiagonals(
		const ftype omega,
		const ftype c[], 
		const int c_sz,
		const ftype rho[],
		const int interface_idcs[], 
		const int interface_idcs_sz,
		const ftype meshsizes[],
		ftype md[], 
		ftype ud[] /*sd*/, 
		int& mat_size)
{
	int N_points = c_sz;
	int layer_number = 0;

	ftype ld[MAX_MAT_SIZE];
	ftype dz = meshsizes[layer_number];
	for (int i = 0; i < N_points - 2; i++)
	{
		if ((layer_number < interface_idcs_sz) && (i == (interface_idcs[layer_number]-1)))
		{
			// special case of the point at the interface
			++layer_number;
			ftype dz_next = meshsizes[layer_number];
			ftype cp = c[i + 1];
			ftype dp = rho[i + 1];
			ftype cm = c[i];
			ftype dm = rho[i];
			ftype q = 1 / (dz_next * dm + dz * dp);

			ld[i] = 2 * q * dp / dz;
			// Magic!
			md[i] = -2 * q * (dz_next * dp + dz * dm) / (dz * dz_next) +
						omega * omega * q * (dz * dp * cp * cp + dz_next * dm * cm * cm) /
							(cp * cp * cm * cm);
			ud[i] = 2 * q * dm / dz_next;
			dz = dz_next;
		} 
		else 
		{
			// ordinary point
			ud[i] = (1 / (dz * dz)); 
			ld[i] = ud[i];
			md[i] = (-2 / (dz * dz) + omega * omega / (c[i + 1] * c[i + 1]));
		}
	}

	// TODO: merge me with ud-ld-md cycle and remove ld array!
	// Symmetrize the matrix
	for (int i = 0; i < N_points - 3; i++)
		ud[i] = sqrt(ud[i] * ld[i + 1]);
	mat_size = N_points - 2;
	// DIAGONALS!!!
}

__device__ Interval ComputeWavenumsLimits(
		const ftype omega, 
		const ftype c[], 
		const int c_sz)
{
	ftype cmin = c[0], cmax = c[0];
	for (int i = 0; i < c_sz; i++)
	{
		if (c[i] < cmin)
			cmin = c[i];
		if (c[i] > cmax)
			cmax = c[i];
	}
	ftype kappamax = omega / cmin;
	ftype kappamin = omega / cmax;
	return Interval {kappamin*kappamin, kappamax*kappamax};
}

__device__ void FillLayers(const int rr, 
		const int n_layers,
		const ftype* depths, 
		const ftype* rhos, 
		const ftype* c1s,
		const ftype* c2s, 
		const int* Ns_points, 
		ftype mesh[], 
		int interface_idcs[], int& interface_idcs_sz,
		ftype c[], int& c_sz, ftype rho[])
{
	c[0] = 0;
	rho[0] = 0;

	// TODO: Rewrite me, i am UGLY ((
	int n = 1; //total number of points
	ftype zp = 0;
	for (unsigned i = 0; i < n_layers; ++i)
	{
		int n_points_layer = Ns_points[i] * rr;
		ftype zc = depths[i];
		mesh[i] = (zc - zp) / n_points_layer; // dz

		c[n - 1] = c1s[i];
		rho[n - 1] = rhos[i];

		for (unsigned k = 1; k <= n_points_layer; ++k)
		{
			rho[n] = rhos[i];
			c[n] = (c1s[i] + (c2s[i] - c1s[i]) * k / n_points_layer);
			++n;
		}
		if (i < n_layers - 1)
			interface_idcs[i] = n - 1;
		zp = zc;
	}

	interface_idcs_sz = n_layers - 1;
	c_sz = n;
}

__device__ void ComputeWavenums(
		const ftype omega,
		const int n_layers,
		const int* Ns_points,
		const ftype* depths,
		const ftype rhos[],
		const ftype c1s[],
		const ftype c2s[],
		ftype wnums[],
		int& wnums_sz)
{
	// Strange things happen here...
	int  Ns_points_aligned [MAX_MAT_SIZE];
	for (int i = 0; i < n_layers; ++i)
		Ns_points_aligned[i] = 12 * (Ns_points[i] / 12);

	ftype coeff_extrap[4][4] = {
			{1,0,0,0},
			{-1, 2, 0, 0},
			{0.5, -4, 4.5, 0},
			{-1 / ftype(6), 4, -13.5, 32 / ftype(3)}};

	for (int rr = 1; rr <= ORD_RICH; ++rr)
	{
		ftype mesh [MAX_MAT_SIZE];
		int interface_idcs [MAX_INTERFACES]; 
		int interface_idcs_sz;
		ftype c [MAX_MAT_SIZE];
		int   c_sz;
		ftype rho [MAX_MAT_SIZE];
		FillLayers(rr, n_layers, depths, rhos, c1s, c2s, Ns_points_aligned, 
				mesh, interface_idcs, interface_idcs_sz, c, c_sz, rho);

		int mat_size;
		ftype md [MAX_MAT_SIZE];
		ftype sd [MAX_MAT_SIZE];
		FillDiagonals(omega, c, c_sz, rho, interface_idcs, interface_idcs_sz, mesh, 
				md, sd, mat_size);

		ftype wnums_rr [MAX_WNUMS];
		int wnums_rr_sz;
		Interval lim = ComputeWavenumsLimits(omega, c, c_sz);
		wnums_rr_sz = bisectGPU(md, sd, mat_size, lim.ll, lim.rl, wnums_rr);
		if (rr == 1) 
			wnums_sz = wnums_rr_sz;
		for (int i = 0; i < wnums_rr_sz; ++i)
			wnums[i] += (wnums_rr[i] * coeff_extrap[ORD_RICH-1][rr-1]);
	}
}

// This procedure computes MGV for a _single_ frequency
__device__ void ComputeModalGroupVelocities (
		const ftype freq,
		const int n_layers,
		const int* Ns_points,
		const ftype* depths,
		const ftype rhos[],
		const ftype c1s[],
		const ftype c2s[],
		ftype mgv[MAX_WNUMS],
		int& mgv_sz)
{
	ftype wnums1 [MAX_WNUMS] = {0}; int wnums1_sz;
	ftype wnums2 [MAX_WNUMS] = {0}; int wnums2_sz;
	// magic number for numerical differentiation procedure
	ftype deltaf = 0.05;
	ftype omega1 = 2 * LOCAL_M_PI * freq + deltaf;
	ftype omega2 = 2 * LOCAL_M_PI * freq;
	
	ComputeWavenums(omega1, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums1, wnums1_sz);
	ComputeWavenums(omega2, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums2, wnums2_sz);

	// Since with increase of omega the number of wave numbers
	// can only increase,  wnum2_sz <= wnum1_sz
	for (int i = 0; i < wnums2_sz; ++i)
		mgv[i] = (omega1 - omega2) / (sqrt(wnums1[i]) - sqrt(wnums2[i]));
	mgv_sz = wnums2_sz;
}

__global__ void EvalPoint_gpukernel(
		const int cws_sz, 
		const int dmaxsz,
		const ftype* cws, 
		const int* Ns_points,
		const ftype* depths,
		const ftype R, 
		const ftype tau, 
		const ftype rhob, 
		const ftype cb, 
		const ftype* freqs, 
		const int freqs_sz,
		const ftype* exp_delays,
		const int* exp_delays_sz,
		ftype* residual,
		int* n_res_global)
{
	int n_layers = cws_sz+1;
	
	const unsigned int tid = (blockIdx.x << BLOCKSIZE) + threadIdx.x;

	if (tid >= freqs_sz)
		return;
	ftype rhos[MAX_MAT_SIZE];
	ftype c1s[MAX_MAT_SIZE];
	ftype c2s[MAX_MAT_SIZE];
	FillLocalArrays(0, cb, rhob, 1, cws_sz, cws,  
			rhos, c1s, c2s);

	int n_residuals = 0;
	ftype residuals_local = 0;
	// Compute mgvs for all frequencies
	//assert (freqs_sz < MAX_FREQS);
	ftype calc_mgv[MAX_WNUMS];
	int calc_mgv_sz;
	ComputeModalGroupVelocities(freqs[tid], n_layers, Ns_points, depths, rhos, c1s, c2s, 
		calc_mgv, calc_mgv_sz);

	int min_size = calc_mgv_sz < exp_delays_sz[tid] ? 
		calc_mgv_sz : exp_delays_sz[tid];

	for (int j = 0; j < min_size; ++j) //iterate over modal velocities
	{
		ftype exp_delay = exp_delays[tid*dmaxsz + j];
		ftype calc_delay = R / calc_mgv[j];
		if (exp_delay > 0)
		{
			residuals_local += pow(exp_delay + tau - calc_delay, 2);
			++n_residuals;
		}
	}

	atomicAdd(residual, residuals_local);
	atomicAdd(n_res_global, n_residuals);

	//residual = sqrt(residuals_local/n_residuals);
}
void BisectResCalcGPU32::EvalPoint( //model
		const int dmaxsz,
		const int* Ns_points,
		const float* depths,
		const float* freqs, 
		const int freqs_sz,
		const float* exp_delays,
		const int* exp_delays_sz,
		//point
		const float R, 
		const float tau, 
		const float rhob, 
		const float cb, 
		const float* cws, 
		const int cws_sz, 
		//output
		float* residual,
		int* n_res_global)
{

	int n_layers = cws_sz+1;

	m_CopyToGPU2(cws, cws_sz, ftype);
	m_CopyToGPU2(freqs, freqs_sz, ftype);
	m_CopyToGPU2(exp_delays, dmaxsz*freqs_sz, ftype);
	m_CopyToGPU2(exp_delays_sz, freqs_sz, int);
	m_CopyToGPU2(depths, n_layers, ftype);
	m_CopyToGPU2(Ns_points, n_layers, int);
	m_CopyToGPU2(residual, 1, ftype);
	m_CopyToGPU2(n_res_global, 1, int);

	EvalPoint_gpukernel <<< freqs_sz/(1<<BLOCKSIZE), 1<<BLOCKSIZE >>> 
		(cws_sz, 
		 dmaxsz, 
		 g_cws, 
		 g_Ns_points, 
		 g_depths, 
		 R, 
		 tau, 
		 rhob, 
		 cb, 
		 g_freqs, 
		 freqs_sz,
		 g_exp_delays,
		 g_exp_delays_sz,
		 g_residual,
		 g_n_res_global);

	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	checkCudaErrors(err);
	//printf("\n Bla");


	checkCudaErrors(cudaMemcpy((void*) residual, (void*)g_residual, 
				sizeof(ftype), cudaMemcpyDeviceToHost));
	checkCudaErrors(cudaMemcpy((void*) n_res_global, (void*)g_n_res_global, 
				sizeof(int), cudaMemcpyDeviceToHost));
	//printf("\n Res_loc: %f %i", *residual, *n_res_global);

	m_FreeGPU(Ns_points);
	m_FreeGPU(depths);
	m_FreeGPU(exp_delays_sz);
	m_FreeGPU(exp_delays);
	m_FreeGPU(freqs);
	m_FreeGPU(cws);
	m_FreeGPU(residual);
	m_FreeGPU(n_res_global);
}

//BisectResCalcGPU32::BisectResCalcGPU32(std::string nm) { /*name_ = nm;*/ }

/*
BisectResCalcGPU32::BisectResCalcGPU32()
{ 
	
	// FIXME: Correct name assginment on constructor!!!!
	//size_t s = sizeof(ftype)*8;
	//name_ = "cpu"+std::to_string(s);
}

*/
