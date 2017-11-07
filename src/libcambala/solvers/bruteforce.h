#ifndef SOLVERS_BRUTEFORCE_H_
#define SOLVERS_BRUTEFORCE_H_
#include "solvers/interface.h"
#include <vector>
#include "solvers/discrete.h"

class BruteForce : public Solver
{
public:
	BruteForce(std::string name);
	BruteForce();
	~BruteForce();
private:
	DiscreteSearchSpace ss_;

	void DoSolve();
	void DoLoadSearchSpaceDims(SearchSpaceDims ssd);
};
#endif
