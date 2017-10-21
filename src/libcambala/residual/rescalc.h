#ifndef RESIDUE_SELECTOR_H_
#define RESIDUE_SELECTOR_H_
#include <unordered_map>
#include <string>
#include "residual/calculator.h"


class ResCalcSelector
{
public:
	void addResidualCalculator(std::string name, *ResCalc rc);
	inline void setResidualCalculator(std::string name, *ResCalc rc) { calcs_[name] = rc; }
	double computeResidual (const Model& m, const Point& p, std::string calc_name);
private:
	std::unordered_map <std::string, *ResCalc> calcs_;
}

#endif
