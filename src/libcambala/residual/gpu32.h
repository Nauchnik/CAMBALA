#ifndef RESIDUAL_GPU32_H_
#define RESIDUAL_GPU32_H_
#include "residual/interface.h"
#include "types.h"
#include "residual/cpuh.h"

class BisectResCalcGPU32 : public BisectResCalcCPU <float>
{

public:
	BisectResCalcGPU32(std::string nm){};
	BisectResCalcGPU32(){};
	virtual ~BisectResCalcGPU32 (){};
private:
	void EvalPoint( //model
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
		int* n_res_global);
};
#endif
