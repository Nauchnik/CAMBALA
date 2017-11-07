#include "solvers/hillclimbing.h"
#define ELPP_STL_LOGGING
#include "easylogging++.h"
#include <iostream>

void HillClimbing::Solve()
//void HillClimbing::Solve(Point startingPoint = Point())
{
	LOG(DEBUG) <<"findLocalMinHillClimbing";
	Point startingPoint = Point();
	PointInds startingPointInd = (startingPoint == Point())
		? ss_.getCenter()
		: ss_.Point2Indexes(startingPoint);

	// launch iterations of hill climbing
	for (size_t i = 0; i < ils_runs_; ++i)
	{
		LOG(INFO) << "iteration " << i << " of ILS";
		Point m = FindLocalMin(startingPointInd);
		if (m < recordPoint_)
			UpdateRecord(m);
		// permutate current global minimum point to obtain a new start point
		while (ss_.Checked(startingPointInd))
		{
			startingPointInd = ss_.PermutateInds(ss_.Point2Indexes(recordPoint_), 1.0/3.0);
		}
		LOG(DEBUG) << "new random cpi : " << startingPointInd;
	}
}

void HillClimbing::UpdateRecord(Point p)
{
	recordPoint_ = p;
	LOG(DEBUG) << "new record found: " << recordPoint_.print();
}

Point HillClimbing::FindLocalMin(PointInds startInds)
{
	// if changing of every variable will not lead to a record update, then a local min reached
	Point localRecord;
	for (;;)
	{
		PointInds roundInds = startInds;
		for (size_t i = 0; i < ss_.axes_.size(); ++i)
		{
			LOG(DEBUG) << i;
			CheckAxis (startInds, i, localRecord, /*forward*/1);
			CheckAxis (startInds, i, localRecord, /*backward*/0);
		}
		if (startInds == roundInds)
			break; // no changes = local minimum
	}
	return localRecord;
}

void HillClimbing::CheckAxis(PointInds pointInds, size_t axisInd, Point& localRecord, bool forward = true)
{
	size_t lastInd = ss_.axes_[axisInd].size()-1;
	if (!forward)
		--pointInds[axisInd]; //Exclude origin point when going backwards
	// Increment/decrement index while within axis size constraints
	for (size_t& i=pointInds[axisInd]; (i<=lastInd) && (i>=0); i+= (forward ? 1 : -1))
	{
		LOG(DEBUG) << "i: " << i << (forward ? " forward" : " backward");
		Point p = GetOrCompute(pointInds);
		if (p <= localRecord) // <= needed for algo to skip origin point
			localRecord = p;
		else
			break;
	}
}

Point HillClimbing::GetOrCompute(PointInds pointInds)
{
	Point p;
	if (ss_.Checked(pointInds))
		p = ss_.checked_[pointInds];
	else
	{
		p = ss_.Indexes2Point(pointInds);
		p_sel_->ComputeResidual(p);
		ss_.AddChecked(pointInds, p);
	}
	return p;
}

void HillClimbing::LoadSearchSpaceDims(SearchSpaceDims ssd)
{
	ss_ = DiscreteSearchSpace(ssd);
}

void HillClimbing::SetResidualCalculatorSelector (ResCalcSelector* p_sel)
{
	p_sel_ = p_sel;
}

Point HillClimbing::getBestPoint() { return recordPoint_; }
//HillClimbing::HillClimbing(){}
