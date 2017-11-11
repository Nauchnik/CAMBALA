#ifndef RESIDUAL_CPUH_H_
#define RESIDUAL_CPUH_H_
#include "residual/interface.h"
#include "types.h"
template <typename ftype> 
class BisectResCalcCPU : public ResidualCalculator
{

public:
	BisectResCalcCPU(std::string nm);
	BisectResCalcCPU();
	virtual ~BisectResCalcCPU  ();
private:
	void DoLoadModel (const Model& m);
	double DoComputeResidual(Point& p);
	const std::string name_;
	void FreeImmutableData();

	size_t freqs_sz_;
	ftype* freqs_ = NULL;

	int* exp_delays_sz_= NULL; //int used for compatibility
	ftype* exp_delays_= NULL; //size = freqs_sz_

	size_t depths_sz_; // n_layers
	ftype* depths_= NULL;

	int* Ns_points_= NULL; //int used for compatibility

	int n_layers_;
	int dmaxsz_;

protected:
	virtual void EvalPoint( //model
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
		int* n_res_global);
};
#endif
