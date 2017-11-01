#include "residual/interface.h"
#include "types.h"
#include <string>

class CLASSNAME : public ResidualCalculator
{
public:
	CLASSNAME(std::string name); 
	CLASSNAME() : name_("CLASSNAME") {}
	void LoadImmutableData (const Model& m);
	double CalculatePointResidual(Point& p);
	double CalculatePointResidual(const Model& m, Point& p);
	std::string getName();
private:
	const std::string name_;
	void FreeImmutableData();
	~CLASSNAME ();

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
