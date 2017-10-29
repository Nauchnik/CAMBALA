#include "residual/interface.h"
#include "types.h"

class CLASSNAME : public ResidualCalculator
{
public:
	void LoadImmutableData (const Model& m);
	double CalculatePointResidual(Point& p);
	double CalculatePointResidual(const Model& m, Point& p);
private:
	void FreeImmutableData();
	inline ~CLASSNAME () {FreeImmutableData();};

	size_t freqs_sz_;
	ftype* freqs_;

	int* exp_delays_sz_; //int used for compatibility
	ftype* exp_delays_; //size = freqs_sz_

	size_t depths_sz_; // n_layers
	ftype* depths_;

	int* Ns_points_; //int used for compatibility

	int n_layers_;
	int dmaxsz_;
};
