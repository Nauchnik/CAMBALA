#include "solvers/discrete.h"
#include "rng.h"
#include "utils.h"
#include <algorithm>
#include <random>
#define ELPP_STL_LOGGING
#include "easylogging++.h"


PointInds DiscreteSearchSpace::getCenter()
{
	PointInds center;
	for (const auto& a: axes_)
		center.push_back((a.size()>1) ? (a.size()/2) : 0);
	return center;
}


Point DiscreteSearchSpace::Indexes2Point(PointInds indexes)
{
	Point point;
	point.cb   = axes_[0][indexes[0]];
	point.rhob = axes_[1][indexes[1]];
	point.R    = axes_[2][indexes[2]];
	point.tau  = axes_[3][indexes[3]];
	for (unsigned i = 4; i < axes_.size(); i++)
		point.cws.push_back(axes_[i][indexes[i]]);
	return point;
}


DiscreteSearchSpace::DiscreteSearchSpace(const SearchSpaceDims& ssd)
{
	axes_.push_back(getDimGrid(ssd.cb));
	axes_.push_back(getDimGrid(ssd.rhob));
	axes_.push_back(getDimGrid(ssd.R));
	axes_.push_back(getDimGrid(ssd.tau));

	// fill search_space_variables[4...] with cws
	for (auto d: ssd.cw)
		axes_.push_back(getDimGrid(d));
}

PointInds DiscreteSearchSpace::Point2Indexes(Point point)
{
	PointInds indexes(axes_.size());

	indexes[0] = std::find(axes_[0].begin(), axes_[0].end(), point.cb)   -axes_[0].begin();
	indexes[1] = std::find(axes_[1].begin(), axes_[1].end(), point.rhob) -axes_[1].begin();
	indexes[2] = std::find(axes_[2].begin(), axes_[2].end(), point.R)    -axes_[2].begin();
	indexes[3] = std::find(axes_[3].begin(), axes_[3].end(), point.tau)  -axes_[3].begin();
	for (size_t i = 4; i < axes_.size(); i++)
		indexes[i] = std::find(axes_[i].begin(), axes_[i].end(), point.cws[i - 4]) - axes_[i].begin();

	return indexes;
}

vector<Point> DiscreteSearchSpace::getSearchSpacePointsVec()
{
	//TODO: check for no duplicates
	vector <Point> all_points;
	PointInds p(axes_.size(),0);
	while(IncreaseInd(p, 0))
	{
		//LOG(DEBUG) << "P: " << p;
		all_points.push_back(Indexes2Point(p));
	}
	return all_points;
}

// former next_cartesian
// construct all combinations of search space parameters
// TODO: rewrite me as a template, if needed
bool DiscreteSearchSpace::IncreaseInd(PointInds& p, size_t k)
{
	bool last_element = (p[k] == (axes_[k].size()-1));
	if (!last_element)
	{
		p[k]++;
		return true;
	}
	else
	{
		bool last_axis = (k == (axes_.size()-1));
		if (last_axis)
			return false;
		p[k] = 0;
		return IncreaseInd(p, k+1);
	}
}

PointInds DiscreteSearchSpace::PermutateInds(PointInds p, float chance = 0.5)
{
	std::uniform_real_distribution<> rdis(0.0, 1.0);
	for (unsigned i = 0; i < axes_.size(); ++i)
	{
		if (rdis(Gen) < chance)
		{
			std::uniform_int_distribution<> idis (0, axes_[i].size()-1);
			p[i] = idis(Gen);
		}
	}
	return p;
}

void DiscreteSearchSpace::AddChecked(PointInds pd, Point p)
{
	// Add duplicates check
	checked_[pd] =  p;
}

bool DiscreteSearchSpace::Checked(PointInds p)
{
	return checked_.count(p)>0;
}






DiscreteSearchSpace::DiscreteSearchSpace() { }


/*
void CAMBALA_sequential::reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a)
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4, ...] - cws
	if (reduced_s_s_a.cb == false)
	{
		axes_[0].resize(1);
		axes_[0][0] = cb1;
	}
	if (reduced_s_s_a.rhob == false)
	{
		axes_[1].resize(1);
		axes_[1][0] = rhob1;
	}
	if (reduced_s_s_a.R == false)
	{
		axes_[2].resize(1);
		axes_[2][0] = R1;
	}
	if (reduced_s_s_a.tau == false)
	{
		axes_[3].resize(1);
		axes_[3][0] = tau1;
	}
	for (unsigned i=0; i < reduced_s_s_a.cws.size(); i++)
	{
		if (reduced_s_s_a.cws[i] == false)
		{
			axes_[4 + i].resize(1);
			axes_[4 + i][0] = cw1_arr[i];
		}
	}
}
*/
