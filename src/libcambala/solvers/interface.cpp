#include "solvers/interface.h"
#include "easylogging++.h"
#include "types.h"

void Solver::Solve ()
{
	LOG(INFO) << "Solver " << getName() << " started";
	DoSolve();
	LOG(INFO) << "Solver " << getName() << " finished solving with result: " << getBestPoint().print();
	LOG(INFO) << "Solver " << getName() << " invoked ResCalc " << ResCalcCount_ << "times.";
}
void Solver::LoadSearchSpaceDims(SearchSpaceDims ssd)
{
	LOG(DEBUG) << "Solver " << getName() << " loading search space dimension data";
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
