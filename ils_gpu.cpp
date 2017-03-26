#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"

Point sspemdd_sequential::generateRandomPoint()
{
	std::vector <size_t> point_indexes (search_space.size());
	for (size_t vi = 0; vi < search_space.size(); ++i)
		point_indexes[vi] = uni(rng);
	Point point = fromPointIndexesToPoint(point_indexes);
	return std::move(point);
}

void EvalPointBatchGPU(std::vector <Point> &batch)
{
	// Transform AoS to SoA
	size_t sz = batch.size();
	size_t sz_cws = batch[0].size();
	float *R = malloc(sz*sizeof(float));
	float *tau = malloc(sz*sizeof(float));
	float *rhob = malloc(sz*sizeof(float));
	float *cb = malloc(sz*sizeof(float));
	float (*cws)[sz] = malloc(sz*sz_cws*sizeof(float));
	for (size_t i = 0; i < sz; ++i)
	{
		R[i]    = batch[i].R;
		tau[i]  = batch[i].tau;
		rhob[i] = batch[i].rhob;
		cb[i]   = batch[i].cb;

		for (size_t j = 0; j < cws_sz; ++j)
			cws[j][i] = batch[i].cws[j];
	}

	float *residuals = malloc(sz*sizeof(float));
	EvalPointsGPU(sz, sz_cws, R, tau, rhob, cb, cws, residuals);

	for (size_t i = 0; i < sz; ++i)
		batch[i].residual = residuals[i];

	free(R);
	free(tau);
	free(rhob);
	free(cb);
	free(cws);
	free(residuals);

}

void sspemdd_sequential::ILSGPU()
{
	std::cout << "Start ILS GPU" << std::endl;

	const size_t batch_size = 10;
	unsigned global_record;

	for (size_t i = 0; i < ils_runs; ++i)
	{
		std::vector <Point> batch;
		for (size_t j = 0; j < batch_size; ++i)
			batch.push_back(generateRandomPoint());
		EvalPointBatchGPU(batch);
		Point best = std::min_element(batch.begin(), batch.end());
		std::cout << "Best of batch: " << best.residual << std::endl;
		if (global_record < best)
			global_record = best;
	}
	std::cout << "Global record" << global_record.residual << std::endl;
}


void FillLocalArrays (
		const int tid,
		const int batch_sz, const int cws_sz, 
		const float cws[cws_sz][batch_sz], 
		const float* cb,
		const float* rhob, 
		float rhos[],
		float c1s[],
		float c2s[])
{
	// FIXME: magic constant!
	for (int i = 0; i < cws_sz +1; ++i)
		c1s[i] = 1500;

	for (unsigned i = 0; i < cws_sz - 1; i++)
	{
		c1s[i] = cws[i][tid];
		c2s[i] = cws[i + 1][tid];
		rhos[i] = 1;
	}

	c1s[cws_sz - 1] = cws[cws_sz - 1][tid];
	c2s[cws_sz - 1] = cws[cws_sz - 1][tid];
	rhos[cws_sz - 1] = 1;

	c1s[cws_sz] = cb[tid];
	c2s[cws_sz] = cb[tid];
	rhos[cws_sz] = rhob[tid];
}


void FillDiagonals(
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

	float dz = meshsizes[layer_number];
	for (int i = 0; i < N_points - 2; i++)
	{
		if ((layer_number < interface_idcs_sz) && (i == (interface_idcs[layer_number]-1)))
		{
			/ special case of the point at the interface
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
						omeg * omeg * q * (dz * dp * cp * cp + dz_next * dm * cm * cm) /
							(cp * cp * cm * cm);
			ud[i] = 2 * q * dm / dz_next;
			dz = dz_next;
		} 
		else 
		{
			// ordinary point
			ud[i] = (1 / (dz * dz)); 
			ld[i] = ud[i];
			md[i] = (-2 / (dz * dz) + omeg * omeg / (c[i + 1] * c[i + 1]));
		}
	}

	// TODO: merge me with ud-ld-md cycle and remove ld array!
	// Symmetrize the matrix
	for (int i = 0; i < N_points - 3; i++)
		ud[i] = sqrt(ud[i] * ld(i + 1));
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
	float kappamax = omeg / cmin;
	float kappamin = omeg / cmax;
	return Interval {kappamin*kappamin, kappamax*kappamax};
}


void FillLayers(const int rr, 
		const float* depths, 
		const float* rhos, 
		const float* c1s,
		const float* c2s, 
		const float* Ns_points, 
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


int ComputeWnumbers(float omega,
		float wnums[])

	c1s
	Ns_points
	const float *rhos,
	const float *depths,
{
	int mat_size = depths_sz;

	// Strange things happen here...
	for (int i = 0; i < n_layers; ++i)
		Ns_points[ii] = 12 * (Ns_points[ii] / 12);

	float c [MAX_MAT_SIZE]; int c_sz;
	float rho [MAX_MAT_SIZE];
	float mesh [MAX_MAT_SIZE];
	int interface_idcs [MAX_INTERFACES]; int interface_idcs_sz;

	FillLayers(1 /*rr*/, depths, rhos, c1s, c2s, Ns_points, 
			mesh, interface_idcs, interface_idcs_sz, 
			c, c_sz, rho);
	int mat_size;
	float md [MAX_MAT_SIZE];
	float sd [MAX_MAT_SIZE];

	FillDiagonals(c, c_sz, rho, interface_idcs, interface_idcs_sz, mesh, 
			md, sd, mat_size);

	Interval lim = ComputeWavenumsLimits(omega, c, c_sz);
	return bisectCpu(md, sd, mat_size, lim.ll, lim.rl, wnums);
}



// This procedure computes MGV for a _single_ frequency
void ComputeModalGroupVelocities (float freq, )
{
	// magic number for numerical differentiation procedure
	float deltaf = 0.05;
	float omega1 = 2 * LOCAL_M_PI * freq + deltaf;
	float omega2 = 2 * LOCAL_M_PI * freq;
	int wnum1_sz = ComputeWnumbers(omega1, wnum1);
	int wnum2_sz = ComputeWnumbers(omega2, wnum2);

	// Since with increase of omega the number of wave numbers
	// can only increase,  wnum2_sz <= wnum1_sz
	for (int i = 0; i < wnum2_sz; j++)
		mgv[i] = (omega1 - omega2) / (sqrt(wnum1[i]) - sqrt(wnum2[i]));

}


float getResidual(
		const float* freqs, 
		const float* exp_delays,
	       	const float calc_mgv[][], const float calc_mgv_sz[],
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
			if (exp_delays[i][j] <= 0)
				continue;

			float calc_delay = 0;
			if (j < calc_mgv_sz[i]) // mode exists
				calc_delay = R / calc_mgv[i][j];
			else if ((i+1 < freqs_sz) && (j < calc_mgv[i + 1])) // next freqs mode exists
				calc_delay = R / calc_mgv[i + 1][j];

			++n;
			// tau_comment: tau usage starts here 
			residual += pow(exp_delays[i][j] + tau - calc_delay, 2);
		}
	}
	residual = sqrt(residual / n);
	return residual;
}

void EvalPointsGPU(const int size, const int size_cws, 
		const float* R, 
		const float* tau, 
		const float* rhob, 
		const float cws[size_cws][size], 
		const float* exp_delays,
		float* residuals)
{

	float residual;
	for (int tid = 0; tid < size; ++tid)
	{
		float rhos[MAX_MAT_SIZE];
		float c1s[MAX_MAT_SIZE];
		float c2s[MAX_MAT_SIZE];
		FillLocalArrays(tid, batch_sz, cws_sz, cws, cb, rhob, 
				rhos, c1s, c2s)

		float calc_mgv[MAX_FREQS][MAX_WNUMS];
		float calc_mgv_sz[MAX_FREQS];
		// Compute mgvs for all frequencies
		for (int i = 0; i < freqs_sz; ++i)
			ComputeModalGroupVelocities(freqs[i], rhos, c1s, c2s, 
				calc_mgv[i], calc_mgv_sz[i]);

		residuals[tid] = getResidual(freqs, exp_delays, calc_mgv, calc_mgv_sz, tau);
	}
}
