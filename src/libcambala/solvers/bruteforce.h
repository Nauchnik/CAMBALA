#ifndef SOLVERS_BRUTEFORCE_H_
#define SOLVERS_BRUTEFORCE_H_
#include "solvers/interface.h"
#include <vector>
#include "solvers/discrete.h"

class BruteForce : public Solver
{
public:
	BruteForce();
	void Solve ();
	//inline void LoadData (std::vector<double> depths) { depths_ = depths; };
	Point getBestPoint();
	void SetResidualCalculatorSelector(ResCalcSelector* p_sel);
	void LoadSearchSpaceDims(SearchSpaceDims ssd);
private:
	//std::vector<double> depths_;
	Point recordPoint_;
	DiscreteSearchSpace ss_;
	ResCalcSelector* p_sel_ = NULL;
};
#endif
