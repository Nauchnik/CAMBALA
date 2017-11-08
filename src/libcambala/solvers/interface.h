#ifndef SOLVERS_INTERFACE_H_
#define SOLVERS_INTERFACE_H_

// Using "Template Method" design pattern!

#include <vector>
#include "types.h"
#include "residual/interface.h"

class Solver
{
public:
	void Solve ();
	void Solve (Point p);
	Point getBestPoint();
	void SetResidualCalculator(ResCalc* rc);
	void LoadSearchSpaceDims(SearchSpaceDims ssd);
	std::string getName();
	int ResCalcCount_;

	virtual ~Solver();

protected:
	void UpdateRecord(Point P);
	ResCalc* rc_ = nullptr;
	Point recordPoint_;
	std::string name_;
	Point startingPoint_;

private:
	virtual void DoSolve() = 0;
	virtual void DoLoadSearchSpaceDims(SearchSpaceDims ssd) = 0;
};

#endif
