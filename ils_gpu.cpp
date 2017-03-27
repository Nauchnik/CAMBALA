#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"
#include "assert.h"
#include "bisect_cpu.h"

#define MAX_MAT_SIZE 2048
#define MAX_FREQS 30
#define MAX_INTERFACES 10
#define MAX_WNUMS 100


void FillLocalArrays (
		const int tid,
		const float cb,
		const float rhob, 
		const int batch_sz, const int cws_sz, 
		//const float cws[cws_sz][batch_sz], 
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


void FillDiagonals(
		const float omega,
		const float c[], const int c_sz,
		const float rho[],
		const int interface_idcs[], const int interface_idcs_sz,
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


Interval ComputeWavenumsLimits(
		const float omega, 
		const float c[], const int c_sz)
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


void FillLayers(const int rr, 
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
	int n = 0; //total number of points
	float zp = 0;
	for (unsigned i = 0; i < n_layers; ++i)
	{
		int n_points_layer = Ns_points[i] / rr;
		float zc = depths[i];
		mesh[i] = (zc - zp) / n_points_layer; // dz

		++n;
		c[n - 1] = c1s[i];
		rho[n - 1] = rhos[i];

		for (unsigned k = 1; k <= n_points_layer; ++k)
		{
			++n;
			rho[n] = rhos[i];
			c[n] = (c1s[i] + (c2s[i] - c1s[i]) * k / n_points_layer);
		}
		if (i < n_layers - 1)
			interface_idcs[i] = n - 1;
		zp = zc;
	}
	interface_idcs_sz = n_layers - 1;
	c_sz = n;
}

void ComputeWavenums(
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

	int mat_size = n_layers;

	// Strange things happen here...
	int  Ns_points_aligned [MAX_MAT_SIZE];
	for (int i = 0; i < n_layers; ++i)
		Ns_points_aligned[i] = 12 * (Ns_points[i] / 12);

	float mesh [MAX_MAT_SIZE];
	int interface_idcs [MAX_INTERFACES]; 
	int interface_idcs_sz;
	float c [MAX_MAT_SIZE];
	int   c_sz;
	float rho [MAX_MAT_SIZE];
	FillLayers(1 /*rr*/, n_layers, depths, rhos, c1s, c2s, Ns_points_aligned, 
			mesh, interface_idcs, interface_idcs_sz, c, c_sz, rho);

	float md [MAX_MAT_SIZE];
	float sd [MAX_MAT_SIZE];
	FillDiagonals(omega, c, c_sz, rho, interface_idcs, interface_idcs_sz, mesh, 
			md, sd, mat_size);

	Interval lim = ComputeWavenumsLimits(omega, c, c_sz);
	wnums_sz = bisectCpu(md, sd, mat_size, lim.ll, lim.rl, wnums);
}



// This procedure computes MGV for a _single_ frequency
void ComputeModalGroupVelocities (
		const float freq,
		const int n_layers,
		const int* Ns_points,
		const float* depths,
		const float rhos[],
		const float c1s[],
		const float c2s[],
		float mgv[],
		int& mgv_sz)
{
	float wnums1 [MAX_WNUMS]; int wnums1_sz;
	float wnums2 [MAX_WNUMS]; int wnums2_sz;
	// magic number for numerical differentiation procedure
	float deltaf = 0.05;
	float omega1 = 2 * LOCAL_M_PI * freq + deltaf;
	float omega2 = 2 * LOCAL_M_PI * freq;
	
	ComputeWavenums(omega1, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums1, wnums1_sz);
	ComputeWavenums(omega2, n_layers, Ns_points, depths, rhos, c1s, c2s, wnums2, wnums2_sz);

	// Since with increase of omega the number of wave numbers
	// can only increase,  wnum2_sz <= wnum1_sz
	for (int i = 0; i < wnums2_sz; i++)
		mgv[i] = (omega1 - omega2) / (sqrt(wnums1[i]) - sqrt(wnums2[i]));
	mgv_sz = wnums2_sz;
}


float getResidual(
		const float R,
		const float* freqs, 
		const int freqs_sz, 
		const float* exp_delays,
		const int* exp_delays_sz,
		const int dmaxsz,
	       	const float calc_mgv[MAX_FREQS][MAX_WNUMS], 
		const int calc_mgv_sz[],
		const float tau)
{
	// Calculate avg distance between experimental and model delays
	// TODO: compare velocities instead of delays to speedup the
	// procedure
	float residual = 0;
	int n = 0;
	for (int i = 0; i < freqs_sz; ++i) //iterate over frequencies
	{
		for (int j = 0; j < calc_mgv_sz[i]; ++j) //iterate over modal velocities
		{
			if (exp_delays[i*dmaxsz + j] <= 0)
				continue;

			float calc_delay = 0;
			if (j < calc_mgv_sz[i]) // mode exists
				calc_delay = R / calc_mgv[i][j];
			else if ((i+1 < freqs_sz) && (j < calc_mgv_sz[i + 1])) // next freqs mode exists
				calc_delay = R / calc_mgv[i + 1][j];

			++n;
			// tau_comment: tau usage starts here 
			residual += pow(exp_delays[i*dmaxsz + j] + tau - calc_delay, 2);
		}
	}
	residual = sqrt(residual / n);
	return residual;
}

void EvalPointsCPU(
		const int batch_sz, 
		const int cws_sz, 
		const int dmaxsz,
		//const float cws[cws_sz][batch_sz], 
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
	float residual;
	for (int tid = 0; tid < batch_sz; ++tid)
	{
		float rhos[MAX_MAT_SIZE];
		float c1s[MAX_MAT_SIZE];
		float c2s[MAX_MAT_SIZE];
		FillLocalArrays(tid, cb[tid], rhob[tid], batch_sz, cws_sz, cws,  
				rhos, c1s, c2s);

		float calc_mgv[MAX_FREQS][MAX_WNUMS];
		int calc_mgv_sz[MAX_FREQS];
		// Compute mgvs for all frequencies
		for (int i = 0; i < freqs_sz; ++i)
			ComputeModalGroupVelocities(freqs[i], n_layers, Ns_points, depths, rhos, c1s, c2s, 
				calc_mgv[i], calc_mgv_sz[i]);

		residuals[tid] = getResidual(R[tid], freqs, freqs_sz, exp_delays, exp_delays_sz, dmaxsz, calc_mgv, calc_mgv_sz, tau[tid]);
	}
}

void EvalPointBatchCPU(
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
	
	// freqs array
	int freqs_sz = freqs_d.size();
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
		for (size_t j = 0; j < modal_delays.size(); ++j)
			exp_delays[i*dmaxsz + j] = modal_delays[i][j];

	
	int n_layers = depths_d.size();
	float *depths = (float*) malloc(n_layers*sizeof(float));
	for (int i=0; i<n_layers; ++i)
		depths[i] = depths_d[i];

	int *Ns_points = (int*) malloc(n_layers*sizeof(int));
	for (int i=0; i<n_layers; ++i)
		Ns_points[i] = Ns_points_d[i];

	float *residuals = (float*) malloc(sz*sizeof(float));
	EvalPointsCPU(sz, cws_sz, dmaxsz, cws, Ns_points, depths, R, tau, rhob, cb, freqs, freqs_sz,
			exp_delays, exp_delays_sz, residuals);

	for (size_t i = 0; i < sz; ++i)
		batch[i].residual = residuals[i];



	free(Ns_points);
	free(depths);
	free(exp_delays_sz);
	free(exp_delays);
	free(freqs);
	free(R);
	free(tau);
	free(rhob);
	free(cb);
	free(cws);
	free(residuals);

}

Point sspemdd_sequential::generateRandomPoint()
{
	std::vector <unsigned> point_indexes (search_space.size());
	for (size_t vi = 0; vi < search_space.size(); ++vi)
	{
		std::uniform_int_distribution<int> uni(0, search_space[vi].size());
		point_indexes[vi] = uni(rng);
	}
	Point point = fromPointIndexesToPoint(point_indexes);
	return std::move(point);
}


void sspemdd_sequential::ILSGPU(int ils_runs)
{
	std::cout << "Start ILS GPU" << std::endl;

	const size_t batch_size = 10;
	Point global_record;

	for (size_t i = 0; i < ils_runs; ++i)
	{
		std::vector <Point> batch;
		for (size_t j = 0; j < batch_size; ++i)
			batch.push_back(generateRandomPoint());
		EvalPointBatchCPU(batch, freqs, Ns_points, depths, modal_delays);
		Point best = *std::min_element(std::begin(batch), std::end(batch));
		std::cout << "Best of batch: " << best.residual << std::endl;
		if (global_record < best)
			global_record = best;
	}
	std::cout << "Global record" << global_record.residual << std::endl;
}
