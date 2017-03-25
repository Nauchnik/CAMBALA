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
	EvalPointsCPU(sz, sz_cws, R, tau, rhob, cb, cws, residuals);

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


void PointToMatrixStrided (int batch_sz, int cws_sz, int mat_size,
		const float* R, 
		const float* tau, 
		const float* rhob, 
		const float cws[cws_sz][batch_sz], 
		float* md,
		float* sd)
{
	c1s.resize(n_layers_w + 1);
	for (auto &x : c1s)
		x = 1500;
	c2s.resize(n_layers_w + 1);

	for (unsigned jj = 0; jj < n_layers_w - 1; jj++)
	{
		c1s.at(jj) = point.cws.at(jj);
		c2s.at(jj) = point.cws.at(jj + 1);
		rhos.at(jj) = 1;
	}

	c1s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	c2s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	rhos.at(n_layers_w - 1) = 1;
	c1s.at(n_layers_w) = point.cb;
	c2s.at(n_layers_w) = point.cb;
	rhos.at(n_layers_w) = point.rhob;

}


void FillDiagonals()
{
	int N_points = (unsigned)c.size();
	unsigned layer_number = 0;
	double dz = meshsizes.at(layer_number);
	double dz_next = dz; //    ofstream myFile("thematrixdiags.txt");

	int next_interface_idx = (interface_idcs_sz > 0) 
		? interface_idcs[0] - 1
		: N_points;

	for (int i = 0; i < N_points - 2; i++)
	{
		// ordinary point
		ud[i] =  ld[i] = (1 / (dz * dz));
		md[i] = (-2 / (dz * dz) + omeg * omeg / (c[i + 1] * c[i + 1]));

		// special case of the point at the interface
		if (i == next_interface_idx)
		{	 
			// Magic!
			++layer_number;
			float cp = c[i + 1];
			float dp = rho[i + 1];
			float cm = c[i];
			float dm = rho[i];
			float q = 1 / (dz_next * dm + dz * dp);

			dz_next = meshsizes[layer_number];

			ld[i] = 2 * q * dp / dz;
			md[i] = -2 * q * (dz_next * dp + dz * dm) / (dz * dz_next) +
						omeg * omeg * q * (dz * dp * cp * cp + dz_next * dm * cm * cm) /
							(cp * cp * cm * cm);
			ud[i] = 2 * q * dm / dz_next;

			next_interface_idx = (interface_idcs_sz > layer_number) 
				? interface_idcs[layer_number] - 1
			       	: N_points;

			dz = dz_next;
		}
	}

	// Symmetrize the matrix
	for (int i = 0; i < N_points - 3; i++)
		ud[i] = sqrt(ud[i] * ld(i + 1));
	// DIAGONALS!!!

}


ComputeWavenumsLimits(float omega, c, float& ll, float& rl)
{
	float cmin = c[0], cmax = c[0];
	for (int ii = 0; ii < N_points; ii++)
	{
		if (c.at(ii) < cmin)
			cmin = c[ii];
		if (c.at(ii) > cmax)
			cmax = c[ii];
	}
	float kappamax = omeg / cmin;
	float kappamin = omeg / cmax;
	ll = kappamin * kappamin;
	rl = kappamin * kappamin;
}

void ComputeWnumbers(float omega,
		float* wnums,
		float* wnums_sz)

{
	int mat_size = depths_sz;


	FillDiagonals(md, sd);

	ComputeWavenumsLimits(omega, c, ll, ul);

	bisectCpu (md, sd, mat_size, ll, rl, eigenvalues);
}



// This procedure computes MGV for a _single_ frequency
void ComputeModalGroupVelocities (float freq, )
{
	// magic number for numerical differentiation procedure
	float deltaf = 0.05;
	float omega2 = 2 * LOCAL_M_PI * freq;
	float omega1 = omega1 + deltaf;
	ComputeWnumbers(omega1, wnum1, wnum1_sz);
	ComputeWnumbers(omega2, wnum2, wnum2_sz);

	// Since with increase of omega the number of wave numbers
	// can only increase,  wnum2 <= wnum1
	for (int i = 0; i < wnum2; j++)
		mgv[i] = (omega1 - omega2) / (wnum1[i] - wnum2[i]);

}

void ComputeMGVResidual ()
{
	float calc_mgv[MAX_FREQS][MAX_WNUMS];
	float calc_mgv_sz[MAX_FREQS];
	// Compute mgvs for all frequencies
	for (int i = 0; i < freqs_sz; ++i)
		ComputeModalGroupVelocities( freqs[i],
			R, tau, rhob, cws, cws_sz,
			calc_mgv[i], calc_mgv_sz[i]);

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
}

void EvalPointsStridedCPU(int size, int size_cws, 
		const float* R, 
		const float* tau, 
		const float* rhob, 
		const float cws[size_cws][size], 
		float* residuals)
{
	const int mat_size = ???????;

	float* md = malloc(size * mat_size * sizeof(float));
	float* sd = malloc(size * mat_size * sizeof(float));

	PointToMatrixStrided(size, size_cws, mat_size, R, tau, rhob, cws,
			/*outputs*/
			md, sd);

	float* eigenvalues = malloc(size * mat_size * sizeof(float));
	int*   eiv_size    = malloc(size * mat_size * sizeof(int));
	BisectCpuStrided(md, sd, mat_size,
		       /*outputs*/
			eigenvalues, eiv_sizes);

		

	ComputeMGVResidual(
		freqs_sz,
		exp_mode_numbers, exp_mode_numbers_sz,
		exp_delays, exp_delays_sz, 
		calc_mgv, calc_mgv_sz,
		/*outputs*/
		residual); 

	free (md);
	free (sd);
	free (eigenvalues);
	free (eiv_sizes);
}
