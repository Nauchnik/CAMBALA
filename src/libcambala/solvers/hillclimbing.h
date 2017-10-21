#ifndef SOLVERS_HILLCLIMBING_H_
#define SOLVERS_HILLCLIMBING_H_
#include "solvers/interface.h"
#include <vector>

class HillClimbing : public Solver
{
public:
	void Solve ();
	void LoadData (std::vector<double> depths); { depths_ = depths; }
	Point getBestPoint() { return recordPoint_ };
private:
	std::vector<double> depths_;
	Point recordPoint_;
}
#endif
