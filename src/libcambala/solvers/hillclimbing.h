#ifndef SOLVERS_HILLCLIMBING_H_
#define SOLVERS_HILLCLIMBING_H_
#include "solvers/interface.h"
#include "solvers/discrete.h"
#include <vector>

class HillClimbing : public Solver
{
public:
	HillClimbing (std::string name); 
	HillClimbing (); 
	size_t ils_runs_ = 10;
	~HillClimbing();
private:
	DiscreteSearchSpace ss_;

	Point FindLocalMin(PointInds startInds);
	void CheckAxis(PointInds pointInds, size_t axisInd, Point& localRecord, bool forward);
	Point GetOrCompute(PointInds pointInds);

	void DoSolve();
	void DoLoadSearchSpaceDims(SearchSpaceDims ssd);
};
#endif
