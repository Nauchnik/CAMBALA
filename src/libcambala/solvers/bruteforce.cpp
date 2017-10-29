#include "solvers/bruteforce.h"
#include "solvers/discrete.h"
#include <iostream>
#include <algorithm>
#include <vector>


void BruteForce::Solve()
{
	//cout << "findGlobalMinBruteForce()" << endl;


	std::vector<Point> Points_vec = ss_.getSearchSpacePointsVec();
	//cout << "Points_vec.size() " << Points_vec.size() << endl;

	for (auto &x : Points_vec)
		p_sel_->ComputeResidual(x);

	//recordPoint_ = *std::min_element(Points_vec.begin(), Points_vec.end());
}

void BruteForce::LoadSearchSpaceDims(SearchSpaceDims ssd)
{
	ss_ = DiscreteSearchSpace(ssd);
};

void BruteForce::SetResidualCalculatorSelector (ResCalcSelector* p_sel)
{
	p_sel_ = p_sel;
}

Point BruteForce::getBestPoint() { return recordPoint_; };
BruteForce::BruteForce(){};
