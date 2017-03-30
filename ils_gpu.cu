#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"
#include "assert.h"
#include "helper_cuda.h"

#define BLOCKSIZE 8

#define MAX_MAT_SIZE 2048
#define MAX_FREQS 1000
#define MAX_INTERFACES 10
#define MAX_WNUMS 100
#define ORD_RICH 1

#include "bisect_gpu.h"

//#define m_CopyToGPU(d, s,  bytes){ if(cudaSuccess != cudaMalloc((void**) &d, bytes)) boinc_temporary_exit(60); checkCudaErrors(cudaMemcpy((void*) d, (void*)s, bytes, cudaMemcpyHostToDevice)); }

#define m_CopyToGPU2(s, elements, type)\
	type* g_##s;\
	checkCudaErrors(cudaMalloc((void**) &g_##s, elements*sizeof(type)));\
	checkCudaErrors(cudaMemcpy((void*) g_##s, (void*)s, elements*sizeof(type), cudaMemcpyHostToDevice));

#define m_FreeHostAndGPU(s)\
	free(s);\
	checkCudaErrors(cudaFree(g_##s));\

__device__ void FillLocalArrays (
		const int tid,
		const float cb,
		const float rhob, 
		const int batch_sz, 
		const int cws_sz, 
		const float* cws, 
		float rhos[],
		float c1s[],
		float c2s[])
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
		const float omega,
		const float c[], 
		const int c_sz,
		const float rho[],
		const int interface_idcs[], 
		const int interface_idcs_sz,
		const float meshsizes[],
		float md[], 
		float ud[] /*sd*/, 
		int& mat_size)
{
	int N_points = c_sz;
	int layer_number = 0;

	float ld[MAX_MAT_SIZE];
	float dz = meshsizes[layer_number];
	for (int i = 0; i < N_points - 2; i++)
	{
		if ((layer_number < interface_idcs_sz) && (i == (interface_idcs[layer_number]-1)))
		{
			// special case of the point at the interface
			++layer_number;
			float dz_next = meshsizes[layer_number];
			float cp = c[i + 1];
			float dp = rho[i + 1];
			float cm = c[i];
			float dm = rho[i];
			float q = 1 / (dz_next * dm + dz * dp);

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
		const float omega, 
		const float c[], 
		const int c_sz)
{
	float cmin = c[0], cmax = c[0];
	for (int i = 0; i < c_sz; i++)
	{
		if (c[i] < cmin)
			cmin = c[i];
		if (c[i] > cmax)
			cmax = c[i];
	}
	float kappamax = omega / cmin;
	float kappamin = omega / cmax;
	return Interval {kappamin*kappamin, kappamax*kappamax};
}

__device__ void FillLayers(const int rr, 
		const int n_layers,
		const float* depths, 
		const float* rhos, 
		const float* c1s,
		const float* c2s, 
		const int* Ns_points, 
		float mesh[], 
		int interface_idcs[], int& interface_idcs_sz,
		float c[], int& c_sz, float rho[])
{
	c[0] = 0;
	rho[0] = 0;

	// TODO: Rewrite me, i am UGLY ((
	int n = 1; //total number of points
	float zp = 0;
	for (unsigned i = 0; i < n_layers; ++i)
	{
		int n_points_layer = Ns_points[i] / rr;
		float zc = depths[i];
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
		const float omega,
		const int n_layers,
		const int* Ns_points,
		const float* depths,
		const float rhos[],
		const float c1s[],
		const float c2s[],
		float wnums[],
		int& wnums_sz)
{
	// Strange things happen here...
	int  Ns_points_aligned [MAX_MAT_SIZE];
	for (int i = 0; i < n_layers; ++i)
		Ns_points_aligned[i] = 12 * (Ns_points[i] / 12);

	float coeff_extrap[4][4] = {
			{1,0,0,0},
			{-1, 2, 0, 0},
			{0.5, -4, 4.5, 0},
			{-1 / float(6), 4, -13.5, 32 / float(3)}};

	for (int rr = 1; rr <= ORD_RICH; ++rr)
	{
		float mesh [MAX_MAT_SIZE];
		int interface_idcs [MAX_INTERFACES]; 
		int interface_idcs_sz;
		float c [MAX_MAT_SIZE];
		int   c_sz;
		float rho [MAX_MAT_SIZE];
		FillLayers(rr, n_layers, depths, rhos, c1s, c2s, Ns_points_aligned, 
				mesh, interface_idcs, interface_idcs_sz, c, c_sz, rho);

		int mat_size;
		float md [MAX_MAT_SIZE];
		float sd [MAX_MAT_SIZE];
		FillDiagonals(omega, c, c_sz, rho, interface_idcs, interface_idcs_sz, mesh, 
				md, sd, mat_size);

		float wnums_rr [MAX_WNUMS];
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
		const float freq,
		const int n_layers,
		const int* Ns_points,
		const float* depths,
		const float rhos[],
		const float c1s[],
		const float c2s[],
		float mgv[MAX_WNUMS],
		int& mgv_sz)
{
	float wnums1 [MAX_WNUMS] = {0}; int wnums1_sz;
	float wnums2 [MAX_WNUMS] = {0}; int wnums2_sz;
	// magic number for numerical differentiation procedure
	float deltaf = 0.05;
	float omega1 = 2 * LOCAL_M_PI * freq + deltaf;
	float omega2 = 2 * LOCAL_M_PI * freq;
	
	ComputeWavenums(omega1, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums1, wnums1_sz);
	ComputeWavenums(omega2, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums2, wnums2_sz);

	// Since with increase of omega the number of wave numbers
	// can only increase,  wnum2_sz <= wnum1_sz
	for (int i = 0; i < wnums2_sz; ++i)
		mgv[i] = (omega1 - omega2) / (sqrt(wnums1[i]) - sqrt(wnums2[i]));
	mgv_sz = wnums2_sz;
}

__global__ void EvalPoints_gpukernel(
		const int batch_sz, 
		const int cws_sz, 
		const int dmaxsz,
		const float* cws, 
		const int* Ns_points,
		const float* depths,
		const float* R, 
		const float* tau, 
		const float* rhob, 
		const float* cb, 
		const float* freqs, 
		const int freqs_sz,
		//const float exp_delays[freqs_sz][dmaxsz], 
		const float* exp_delays,
		const int* exp_delays_sz,
		float* residuals)
{
	int n_layers = cws_sz+1;
	
	const unsigned int tid = (blockIdx.x << BLOCKSIZE) + threadIdx.x;

	float rhos[MAX_MAT_SIZE];
	float c1s[MAX_MAT_SIZE];
	float c2s[MAX_MAT_SIZE];
	FillLocalArrays(tid, cb[tid], rhob[tid], batch_sz, cws_sz, cws,  
			rhos, c1s, c2s);

	int n_residuals = 0;
	float residuals_local = 0;
	// Compute mgvs for all frequencies
	assert (freqs_sz < MAX_FREQS);
	for (int i = 0; i < freqs_sz; ++i)
	{
		float calc_mgv[MAX_WNUMS];
		int calc_mgv_sz;
		ComputeModalGroupVelocities(freqs[i], n_layers, Ns_points, depths, rhos, c1s, c2s, 
			calc_mgv, calc_mgv_sz);

		int min_size = calc_mgv_sz < exp_delays_sz[i] ? 
			calc_mgv_sz : exp_delays_sz[i];

		for (int j = 0; j < min_size; ++j) //iterate over modal velocities
		{
			float exp_delay = exp_delays[i*dmaxsz + j];
			float calc_delay = R[tid] / calc_mgv[j];
			if (exp_delay > 0)
			{
				residuals_local += pow(exp_delay + tau[tid] - calc_delay, 2);
				++n_residuals;
			}
		}

	}
	residuals[tid] = sqrt(residuals_local/n_residuals);
}

void EvalPointBatchGPU(
		std::vector <Point> &batch,
		const std::vector<double> &freqs_d,
		const std::vector<unsigned> &Ns_points_d,
		const std::vector<double> &depths_d,
		const std::vector<std::vector<double>> &modal_delays)
{
	// Transform AoS to SoA
	size_t sz = batch.size();
	size_t cws_sz = batch[0].cws.size();
	float *R =   (float*) malloc(sz*sizeof(float));
	float *tau = (float*) malloc(sz*sizeof(float));
	float *rhob= (float*) malloc(sz*sizeof(float));
	float *cb =  (float*) malloc(sz*sizeof(float));
	float *cws = (float*) malloc(sz*cws_sz*sizeof(float));
	for (size_t i = 0; i < sz; ++i)
	{
		R[i]    = batch[i].R;
		tau[i]  = batch[i].tau;
		rhob[i] = batch[i].rhob;
		cb[i]   = batch[i].cb;

		for (size_t j = 0; j < cws_sz; ++j)
			cws[j*sz + i] = batch[i].cws[j];
	}
	//TODO: stop converting this data every time
	assert (freqs_d.size() == modal_delays.size());
	
	std::cout << " copy const" << std::endl;
	// freqs array
	int freqs_sz = freqs_d.size();
	std::cout << " num freqs " << freqs_sz << std::endl;
	float *freqs = (float*) malloc(freqs_sz*sizeof(float));
	for (int i = 0; i < freqs_sz; ++i)
		freqs[i] = freqs_d[i];

	// exp_delays_sz
	int *exp_delays_sz = (int*) malloc(freqs_sz*sizeof(int));
	for (size_t i = 0; i < freqs_sz; ++i)
		exp_delays_sz[i] = modal_delays[i].size();

	// exp_delays 2d array
	int dmaxsz = 0;
	for (size_t i = 0; i < freqs_sz; ++i)
		dmaxsz = std::max(dmaxsz, exp_delays_sz[i]);
	float *exp_delays = (float*) malloc(dmaxsz*freqs_sz*sizeof(float));
	for (size_t i = 0; i < modal_delays.size(); ++i)
		for (size_t j = 0; j < modal_delays[i].size(); ++j)
			exp_delays[i*dmaxsz + j] = modal_delays[i][j];

	
	int n_layers = depths_d.size();
	float *depths = (float*) malloc(n_layers*sizeof(float));
	for (int i=0; i<n_layers; ++i)
		depths[i] = depths_d[i];

	int *Ns_points = (int*) malloc(n_layers*sizeof(int));
	for (int i=0; i<n_layers; ++i)
		Ns_points[i] = Ns_points_d[i];

	float *residuals = (float*) malloc(sz*sizeof(float));

	m_CopyToGPU2(R,   sz, float);
	m_CopyToGPU2(tau, sz, float);
	m_CopyToGPU2(rhob,sz, float);
	m_CopyToGPU2(cb,  sz, float);
	m_CopyToGPU2(cws, sz*cws_sz, float);
	m_CopyToGPU2(freqs, freqs_sz, float);
	m_CopyToGPU2(exp_delays, dmaxsz*freqs_sz, float);
	m_CopyToGPU2(exp_delays_sz, freqs_sz, int);
	m_CopyToGPU2(depths, n_layers, float);
	m_CopyToGPU2(Ns_points, n_layers, int);
	m_CopyToGPU2(residuals, sz, float);

	cudaEvent_t kernel_start, kernel_stop;
	cudaEventCreate(&kernel_start);
	cudaEventCreate(&kernel_stop);
	cudaEventRecord(kernel_start);

	EvalPoints_gpukernel <<< sz/(1<<BLOCKSIZE), 1<<BLOCKSIZE >>> 
		(sz, cws_sz, dmaxsz, g_cws, g_Ns_points, g_depths, g_R, g_tau, g_rhob, g_cb, g_freqs, freqs_sz,
			g_exp_delays, g_exp_delays_sz, g_residuals);
	#ifndef NDEBUG
	cudaThreadSynchronize();
	cudaError_t err = cudaGetLastError();
	checkCudaErrors(err);
	printf("\n Bla");
	#endif

	cudaThreadSynchronize();
	cudaEventRecord(kernel_stop);
	cudaEventSynchronize(kernel_stop);
	float runTime;
	cudaEventElapsedTime(&runTime, kernel_start, kernel_stop);


	float tm = runTime / 1000;
	printf("\n Time: %f", tm);
	
	checkCudaErrors(cudaMemcpy((void*) residuals, (void*)g_residuals, 
				sz*sizeof(float), cudaMemcpyDeviceToHost));
	for (size_t i = 0; i < sz; ++i)
		batch[i].residual = residuals[i];

	m_FreeHostAndGPU(Ns_points);
	m_FreeHostAndGPU(depths);
	m_FreeHostAndGPU(exp_delays_sz);
	m_FreeHostAndGPU(exp_delays);
	m_FreeHostAndGPU(freqs);
	m_FreeHostAndGPU(R);
	m_FreeHostAndGPU(tau);
	m_FreeHostAndGPU(rhob);
	m_FreeHostAndGPU(cb);
	m_FreeHostAndGPU(cws);
	m_FreeHostAndGPU(residuals);

}

Point sspemdd_sequential::generateRandomPoint()
{
	std::vector <unsigned> point_indexes;
	for (const auto &var :search_space)
	{
		std::uniform_int_distribution<int> uni(0, var.size()-1);
		point_indexes.push_back(uni(rng));
		//std::cout << " rand var " << randnum << std::endl;
	}
	Point point = fromPointIndexesToPoint(point_indexes);
	return std::move(point);
}

void sspemdd_sequential::ILSGPU(int ils_runs)
{
	std::cout << "Start ILS GPU" << std::endl;

	const size_t batch_size = 4*1280;
	Point global_record;

	std::cout << "Global record" << global_record.residual << std::endl;
	for (size_t i = 0; i < ils_runs; ++i)
	{
		std::vector <Point> batch;
		for (size_t j = 0; j < batch_size; ++j)
			batch.push_back(generateRandomPoint());
		std::cout << " start eval " << std::endl;
		EvalPointBatchGPU(batch, freqs, Ns_points, depths, modal_delays);
		Point best = *std::min_element(std::begin(batch), std::end(batch));
		std::cout << "Best of batch: " << best.residual << std::endl;
		if (best < global_record )
			global_record = best;
	}
	record_point = global_record;
	std::cout << "Global record" << global_record.residual << std::endl;
}
