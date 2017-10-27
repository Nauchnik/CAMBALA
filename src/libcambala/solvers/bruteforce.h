#ifndef SOLVERS_BRUTEFORCE_H_
#define SOLVERS_BRUTEFORCE_H_
#include "solvers/bruteforce.h"
#include "solvers/interface.h"
#include <vector>

class BruteForce : public Solver
{
public:
	void Solve ();
	void LoadData (std::vector<double> depths) { depths_ = depths; };
	Point getBestPoint() { return recordPoint_; };
	void SetResidualCalculatorSelector(ResCalcSelector* p_sel);
	void LoadSearchSpaceDims(SearchSpaceDims ssd);
private:
	std::vector<double> depths_;
	Point recordPoint_;
};
#endif
