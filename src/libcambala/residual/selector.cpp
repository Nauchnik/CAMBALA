#include "residual/selector.h"



ResCalcSelector::ResCalcSelector()
{
	calcs_["default"] = nullptr;
	calcs_["precise"] = nullptr;
	calcs_["fast"] = nullptr;
}

void ResCalcSelector::AddResidualCalculator(std::string name, ResCalc* rc)
{
	//TODO: add more semantic checks
	calcs_[name] = rc;
	//add default calculators
	if( calcs_["default"]==nullptr)
	{
		calcs_["default"] = rc;
		calcs_["precise"] = rc;
		calcs_["fast"] = rc;
	}
}

double ResCalcSelector::ComputeResidual (const Model& m, Point& p, 
			std::string calc_name) 
{ 
	return calcs_[calc_name]->CalculatePointResidual(m, p);
}

double ResCalcSelector::ComputeResidual (Point& p, 
			std::string calc_name) 
{ 
	if (m_.depths.size() == 0)
		exit(1);

	return calcs_[calc_name]->CalculatePointResidual(m_, p);
}
