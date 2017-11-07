#include "solvers/bruteforce.h"
#include "solvers/discrete.h"
#include <iostream>
#include <algorithm>
#include <vector>
#define ELPP_STL_LOGGING
#include "easylogging++.h"


void BruteForce::DoSolve()
{
	LOG(DEBUG) << "Bruteforce.Solve()";
	//std::vector<Point> Points_vec = ss_.getSearchSpacePointsVec();
	//LOG(DEBUG) << "Bruteforce Points_vec.size " << Points_vec.size();

	/*
	for (auto &x : Points_vec)
		p_sel_->ComputeResidual(x);
		*/

	PointInds p(ss_.axes_.size(),0);
	while(ss_.IncreaseInd(p, 0))
	{
		Point x = ss_.Indexes2Point(p);
		//LOG(DEBUG) << "I: " << p;
		rc_->ComputeResidual(x);
		if (x.residual < recordPoint_.residual)
			recordPoint_ = x;
	}
}

void BruteForce::DoLoadSearchSpaceDims(SearchSpaceDims ssd)
{
	ss_ = DiscreteSearchSpace(ssd);
}


BruteForce::BruteForce()
{name_= "Bruteforce"; }
BruteForce::BruteForce(std::string name)
{name_ = name;}


BruteForce::~BruteForce() {}
