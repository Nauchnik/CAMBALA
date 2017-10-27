#include <solvers/discrete.h>
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


void DiscreteSearchSpace::DiscreteSearchSpace(const SearchSpaceDims& ssd)
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
	// TODO: rewrite me for type correctness
	vector<vector<size_t>> search_space_indexes(axes_.size());
	for (unsigned i = 0; i < axes_.size(); i++)
		for (unsigned j = 0; j < axes_[i].size(); j++)
			search_space_indexes[i].push_back(j);

	vector<int> index_arr;
	vector<Point> points_vec;
	PointInds cur_point_indexes;
	while (CAMBALA_utils::next_cartesian(search_space_indexes, index_arr, cur_point_indexes))
		points_vec.push_back(Indexes2Point(cur_point_indexes));

	return points_vec;
}

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
