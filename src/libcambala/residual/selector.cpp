#include "residual/selector.h"

void ResCalcSelector::ResCalcSelector()
{
	calcs_["default"] = nullptr;
	calcs_["precise"] = nullptr;
	calcs_["fast"] = nullptr;
}

void ResCalcSelector::addResidualCalculator(std::string name, *ResCalc rc)
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

double ResCalcSelector::computeResidual (const Model& m, const Point& p, 
			std::string calc_name = "default") 
{ 
	return calcs_[calc_name]->CalculatePointResidual(m, p);
}
