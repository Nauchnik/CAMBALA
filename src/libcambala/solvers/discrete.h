#ifndef SOLVERS_DISCRETE_H
#define SOLVERS_DISCRETE_H
using PointInds = vector <size_t>;

class DiscreteSearchSpace
{
public:
	//Cartesian axes
	vector<vector<double>> Axes_;// values of variables which form the search space
	void DiscreteSearchSpace (const SearchSpaceDims& ssd);
	PointInds Point2Indexes(Point point);
	vector<Point> getSearchSpacePointsVec();
}
#endif
