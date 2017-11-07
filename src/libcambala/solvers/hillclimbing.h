#ifndef SOLVERS_HILLCLIMBING_H_
#define SOLVERS_HILLCLIMBING_H_
#include "solvers/interface.h"
#include "solvers/discrete.h"
#include <vector>

class HillClimbing : public Solver
{
public:
	void Solve ();
	void LoadData (std::vector<double> depths);
	Point getBestPoint();
	void SetResidualCalculatorSelector(ResCalcSelector* p_sel);
	void LoadSearchSpaceDims(SearchSpaceDims ssd);
	size_t ils_runs_ = 10;
private:
	std::vector<double> depths_;
	Point recordPoint_;
	DiscreteSearchSpace ss_;
	ResCalcSelector* p_sel_ = NULL;
	Point FindLocalMin(PointInds startInds);
	void UpdateRecord(Point p);
	void CheckAxis(PointInds pointInds, size_t axisInd, Point& localRecord, bool forward);
	Point GetOrCompute(PointInds pointInds);
};
#endif
