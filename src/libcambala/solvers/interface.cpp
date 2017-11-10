#include "solvers/interface.h"
#define ELPP_FEATURE_PERFORMANCE_TRACKING
#include "easylogging++.h"
#include "types.h"

void Solver::Solve ()
{
	startingPoint_ = Point(); // drop starting point to default
	TIMED_FUNC(timerObj);
	LOG(INFO) << "Solver " << getName() << " started";
	DoSolve();
	LOG(INFO) << "Solver " << getName() << " finished solving with result: " << getBestPoint().print();
	LOG(INFO) << "Solver " << getName() << " invoked ResCalc " << ResCalcCount_ << "times.";
}

void Solver::Solve (Point p)
{
	startingPoint_ = p;
	LOG(INFO) << "Using point " << startingPoint_.print() << " as starting point ";
	Solve();
}



void Solver::LoadSearchSpaceDims(SearchSpaceDims ssd)
{
	LOG(DEBUG) << "Solver " << getName() << " loading search space dimension data";
	LOG(DEBUG) << "Search space dimensions: " << ssd.print();
	DoLoadSearchSpaceDims(ssd);
}

void Solver::SetResidualCalculator(ResCalc* rc)
{
	LOG(INFO) << "Solver " << getName() << " attaching residual calculator "  ;
	rc_ = rc;
}

void Solver::UpdateRecord(Point p)
{
	recordPoint_ = p;
	LOG(DEBUG) << "New record found: " << recordPoint_.print();
}

Point Solver::getBestPoint() { return recordPoint_; }
Solver::~Solver() {}

std::string Solver::getName() {return name_;}
