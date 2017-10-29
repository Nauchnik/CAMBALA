#ifndef RESIDUAL_SELECTOR_H_
#define RESIDUAL_SELECTOR_H_
#include <string>
#include <map>
#include "residual/interface.h"
#include "types.h"

class ResCalcSelector
{
public:
	ResCalcSelector();
	void AddResidualCalculator(std::string name, ResCalc* rc);
	double ComputeResidual (const Model& m, Point& p, 
			std::string calc_name = "default");
	double ComputeResidual (Point& p, 
			std::string calc_name = "default");
	inline void LoadModel(Model m){ m_ = m;};
private:
	Model m_;
	std::map <std::string, ResCalc*> calcs_;
};
#endif
