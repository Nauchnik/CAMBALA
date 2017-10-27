#ifndef RESIDUAL_SELECTOR_H_
#define RESIDUAL_SELECTOR_H_
#include <string>
#include "residual/interface.h"
#include "types.h"

class ResCalcSelector
{
public:
	ResCalcSelector();
	void addResidualCalculator(std::string name, ResCalc* rc);
	double computeResidual (const Model& m, const Point& p, 
			std::string calc_name = "default");
	inline void LoadModel(Model m){ m_ = m;};
private:
	Model m_;
};
#endif
