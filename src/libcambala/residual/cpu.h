#include "residual/interface.h"
#include "types.h"
#include <string>
#include "assert.h"
#include <cmath>
#define ELPP_STL_LOGGING
#include "easylogging++.h"

#define MAX_MAT_SIZE 2048
#define MAX_FREQS 1000
#define MAX_INTERFACES 10
#define MAX_WNUMS 100

#ifndef ORD_RICH
  #define ORD_RICH 1
#endif

#include "residual/bisect.h"

template <typename ftype> 
class BisectResCalcCPU : public ResidualCalculator
{
public:
	BisectResCalcCPU (std::string name); 
	BisectResCalcCPU () : name_("CLASSNAME") {}
	void LoadImmutableData (const Model& m);
	double CalculatePointResidual(Point& p);
	double CalculatePointResidual(const Model& m, Point& p);
	std::string getName();
private:
	const std::string name_;
	void FreeImmutableData();
	~BisectResCalcCPU  ();

	size_t freqs_sz_;
	ftype* freqs_ = NULL;

	int* exp_delays_sz_= NULL; //int used for compatibility
	ftype* exp_delays_= NULL; //size = freqs_sz_

	size_t depths_sz_; // n_layers
	ftype* depths_= NULL;

	int* Ns_points_= NULL; //int used for compatibility

	int n_layers_;
	int dmaxsz_;
};

template <typename ftype> 
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

template <typename ftype> 
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

template <typename ftype> 
Interval<ftype> ComputeWavenumsLimits(
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
	return Interval<ftype> {kappamin*kappamin, kappamax*kappamax};
}

template <typename ftype> 
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

template <typename ftype> 
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
		Interval<ftype> lim = ComputeWavenumsLimits(omega, c, c_sz);
		wnums_rr_sz = bisectGPU(md, sd, mat_size, lim.ll, lim.rl, wnums_rr);
		if (rr == 1) 
			wnums_sz = wnums_rr_sz;
		for (int i = 0; i < wnums_rr_sz; ++i)
			wnums[i] += (wnums_rr[i] * coeff_extrap[ORD_RICH-1][rr-1]);
	}
}

// This procedure computes MGV for a _single_ frequency
template <typename ftype> 
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

template <typename ftype> 
void EvalPoint( const unsigned int tid,
		//model
		const int dmaxsz,
		const int* Ns_points,
		const ftype* depths,
		const ftype* freqs, 
		const int freqs_sz,
		const ftype* exp_delays,
		const int* exp_delays_sz,
		//point
		const ftype R, 
		const ftype tau, 
		const ftype rhob, 
		const ftype cb, 
		const ftype* cws, 
		const int cws_sz, 
		//output
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


template <typename ftype> 
void BisectResCalcCPU<ftype>::FreeImmutableData()
{
	//TODO: add more safety checks!
	if (freqs_ != NULL)
	{
		free(freqs_);
		free(exp_delays_sz_);
		free(exp_delays_);
		free(depths_);
		free(Ns_points_);
	}
}

template <typename ftype> 
void BisectResCalcCPU<ftype>::LoadImmutableData (const Model& m)
{
	//LOG(DEBUG) << getName() << "->LoadImmutableData";
	FreeImmutableData();

	assert (m.freqs.size() == m.exp_delays.size());
	
	//std::cout << " copy const" << std::endl;
	// model.freqs
	freqs_sz_ = m.freqs.size();
	//std::cout << " num freqs " << freqs_sz << std::endl;
	freqs_ = (ftype*) malloc(freqs_sz_*sizeof(ftype));
	for (int i = 0; i < freqs_sz_; ++i)
		freqs_[i] = m.freqs[i];

	// model.exp_delays subvector sizes
	exp_delays_sz_ = (int*) malloc(freqs_sz_*sizeof(int));
	for (size_t i = 0; i < freqs_sz_; ++i)
		exp_delays_sz_[i] = m.exp_delays[i].size();


	// model.exp_delays 2d array
	dmaxsz_ = 0;
	for (size_t i = 0; i < freqs_sz_; ++i)
		dmaxsz_ = std::max(dmaxsz_, exp_delays_sz_[i]);
	exp_delays_ = (ftype*) malloc(dmaxsz_*freqs_sz_*sizeof(ftype));
	for (size_t i = 0; i < m.exp_delays.size(); ++i)
		for (size_t j = 0; j < m.exp_delays[i].size(); ++j)
			exp_delays_[i*dmaxsz_ + j] = m.exp_delays[i][j];

	// model.depths
	n_layers_ = m.depths.size();
	depths_ = (ftype*) malloc(n_layers_*sizeof(ftype));
	for (int i=0; i<n_layers_; ++i)
		depths_[i] = m.depths[i];

	// model.Ns_points
	Ns_points_ = (int*) malloc(n_layers_*sizeof(int));
	for (int i=0; i<n_layers_; ++i)
		Ns_points_[i] = m.Ns_points[i];
}

template <typename ftype> 
double BisectResCalcCPU<ftype>::CalculatePointResidual(const Model& m, Point& p)
{
	LoadImmutableData(m);
	return CalculatePointResidual(p);
}

template <typename ftype> 
double BisectResCalcCPU<ftype>::CalculatePointResidual(Point& p)
{
	if (freqs_ == NULL)
		exit(1);
	//LOG(DEBUG) << "Point: " << p.cws;
	// Transform AoS to SoA
	size_t cws_sz = p.cws.size();
	ftype *cws = (ftype*) malloc(cws_sz*sizeof(ftype));
	for (size_t i = 0; i < cws_sz; ++i)
		cws[i] = p.cws[i];

	// output array
	ftype *residual = (ftype*) malloc(sizeof(ftype));
	residual[0] = 0;
	int *n_res_global = (int*) malloc(sizeof(int));
	n_res_global[0] = 0;


	for (int tid = 0; tid<freqs_sz_; ++tid)
		EvalPoint <ftype> // <- needed for correct conversion of p.something
			(tid,
			 dmaxsz_, Ns_points_, depths_, freqs_, freqs_sz_, exp_delays_, exp_delays_sz_,
			 p.R, p.tau, p.rhob, p.cb, cws, cws_sz,
			 residual,
			 n_res_global);

	double residual_total = std::sqrt(*residual / *n_res_global);

	free(residual);
	free(cws);
	free(n_res_global);
	
	LOG(DEBUG) << "Residual: " << residual_total;
	return p.residual = residual_total;
}

template <typename ftype> 
BisectResCalcCPU<ftype>::BisectResCalcCPU(std::string name) : name_{name}
{
	LOG(DEBUG) << "Created residual calculator \""<< name_ << "\"." ;
}

template <typename ftype> 
std::string BisectResCalcCPU<ftype>::getName()
{
	return name_;
}

template <typename ftype> 
BisectResCalcCPU<ftype>::~BisectResCalcCPU ()
{
	//LOG(DEBUG) << "Destructor for "<< name_ << "\"." ;
	FreeImmutableData();
};
