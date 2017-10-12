#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"
#include "assert.h"
#include <cmath>


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

#include "bisect_cpu.h"


void FillLocalArrays (
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

void FillDiagonals(
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

Interval ComputeWavenumsLimits(
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

void FillLayers(const int rr, 
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

void ComputeWavenums(
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
void ComputeModalGroupVelocities (
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

void EvalPoint(
		const unsigned int tid,
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

	*residual += residuals_local;
	*n_res_global += n_residuals;

	//residual = sqrt(residuals_local/n_residuals);
}
void EvalPoints(
		const unsigned int tid,
		const int batch_sz, 
		const int cws_sz, 
		const int dmaxsz,
		const ftype* cws, 
		const int* Ns_points,
		const ftype* depths,
		const ftype* R, 
		const ftype* tau, 
		const ftype* rhob, 
		const ftype* cb, 
		const ftype* freqs, 
		const int freqs_sz,
		//const ftype exp_delays[freqs_sz][dmaxsz], 
		const ftype* exp_delays,
		const int* exp_delays_sz,
		ftype* residuals)
{
	int n_layers = cws_sz+1;
	

	ftype rhos[MAX_MAT_SIZE];
	ftype c1s[MAX_MAT_SIZE];
	ftype c2s[MAX_MAT_SIZE];
	FillLocalArrays(tid, cb[tid], rhob[tid], batch_sz, cws_sz, cws,  
			rhos, c1s, c2s);

	int n_residuals = 0;
	ftype residuals_local = 0;
	// Compute mgvs for all frequencies
	//assert (freqs_sz < MAX_FREQS);
	for (int i = 0; i < freqs_sz; ++i)
	{
		ftype calc_mgv[MAX_WNUMS];
		int calc_mgv_sz;
		ComputeModalGroupVelocities(freqs[i], n_layers, Ns_points, depths, rhos, c1s, c2s, 
			calc_mgv, calc_mgv_sz);

		int min_size = calc_mgv_sz < exp_delays_sz[i] ? 
			calc_mgv_sz : exp_delays_sz[i];

		for (int j = 0; j < min_size; ++j) //iterate over modal velocities
		{
			ftype exp_delay = exp_delays[i*dmaxsz + j];
			ftype calc_delay = R[tid] / calc_mgv[j];
			if (exp_delay > 0)
			{
				residuals_local += pow(exp_delay + tau[tid] - calc_delay, 2);
				++n_residuals;
			}
		}

	}
	residuals[tid] = sqrt(residuals_local/n_residuals);
}

void EvalPointCPU(
		search_space_point &point,
		const std::vector<double> &freqs_d,
		const std::vector<unsigned> &Ns_points_d,
		const std::vector<double> &depths_d,
		const std::vector<std::vector<double>> &modal_delays)
{
	// Transform AoS to SoA
	size_t cws_sz = point.cws.size();
	ftype *cws = (ftype*) malloc(cws_sz*sizeof(ftype));
	for (size_t i = 0; i < cws_sz; ++i)
		cws[i] = point.cws[i];
	//TODO: stop converting this data every time
	assert (freqs_d.size() == modal_delays.size());
	
	//std::cout << " copy const" << std::endl;
	// freqs array
	int freqs_sz = freqs_d.size();
	//std::cout << " num freqs " << freqs_sz << std::endl;
	ftype *freqs = (ftype*) malloc(freqs_sz*sizeof(ftype));
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
	ftype *exp_delays = (ftype*) malloc(dmaxsz*freqs_sz*sizeof(ftype));
	for (size_t i = 0; i < modal_delays.size(); ++i)
		for (size_t j = 0; j < modal_delays[i].size(); ++j)
			exp_delays[i*dmaxsz + j] = modal_delays[i][j];

	int n_layers = depths_d.size();
	ftype *depths = (ftype*) malloc(n_layers*sizeof(ftype));
	for (int i=0; i<n_layers; ++i)
		depths[i] = depths_d[i];

	int *Ns_points = (int*) malloc(n_layers*sizeof(int));
	for (int i=0; i<n_layers; ++i)
		Ns_points[i] = Ns_points_d[i];

	// output array
	ftype *residual = (ftype*) malloc(sizeof(ftype));
	residual[0] = 0;
	int *n_res_global = (int*) malloc(sizeof(int));
	n_res_global[0] = 0;


	for (int tid = 0; tid<freqs_sz; ++tid)
		EvalPoint 
			(
			 tid,
			 cws_sz, 
			 dmaxsz, 
			 cws, 
			 Ns_points, 
			 depths, 
			 point.R, 
			 point.tau, 
			 point.rhob, 
			 point.cb, 
			 freqs, 
			 freqs_sz,
			 exp_delays,
			 exp_delays_sz,
			 residual,
			 n_res_global);

	point.residual = std::sqrt(*residual / *n_res_global);

	free(Ns_points);
	free(depths);
	free(exp_delays_sz);
	free(exp_delays);
	free(freqs);
	free(cws);
	free(residual);
	free(n_res_global);
}
void EvalPointBatchCPU(
		std::vector <search_space_point> &batch,
		const std::vector<double> &freqs_d,
		const std::vector<unsigned> &Ns_points_d,
		const std::vector<double> &depths_d,
		const std::vector<std::vector<double>> &modal_delays)
{
	// Transform AoS to SoA
	size_t sz = batch.size();
	size_t cws_sz = batch[0].cws.size();
	ftype *R =   (ftype*) malloc(sz*sizeof(ftype));
	ftype *tau = (ftype*) malloc(sz*sizeof(ftype));
	ftype *rhob= (ftype*) malloc(sz*sizeof(ftype));
	ftype *cb =  (ftype*) malloc(sz*sizeof(ftype));
	ftype *cws = (ftype*) malloc(sz*cws_sz*sizeof(ftype));
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
	ftype *freqs = (ftype*) malloc(freqs_sz*sizeof(ftype));
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
	ftype *exp_delays = (ftype*) malloc(dmaxsz*freqs_sz*sizeof(ftype));
	for (size_t i = 0; i < modal_delays.size(); ++i)
		for (size_t j = 0; j < modal_delays[i].size(); ++j)
			exp_delays[i*dmaxsz + j] = modal_delays[i][j];

	
	int n_layers = depths_d.size();
	ftype *depths = (ftype*) malloc(n_layers*sizeof(ftype));
	for (int i=0; i<n_layers; ++i)
		depths[i] = depths_d[i];

	int *Ns_points = (int*) malloc(n_layers*sizeof(int));
	for (int i=0; i<n_layers; ++i)
		Ns_points[i] = Ns_points_d[i];

	ftype *residuals = (ftype*) malloc(sz*sizeof(ftype));

	for (int tid = 0; tid<freqs_sz; ++tid)
		EvalPoints
			(tid, sz, cws_sz, dmaxsz, cws, Ns_points, depths, R, tau, rhob, cb, freqs, freqs_sz,
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
