#ifndef SOLVERS_DISCRETE_H
#define SOLVERS_DISCRETE_H
#include <vector>
#include "types.h"
using namespace std;
using PointInds = vector <size_t>;

class DiscreteSearchSpace
{
public:
	//Cartesian axes
	vector<vector<double>> axes_;// values of variables which form the search space
	DiscreteSearchSpace (const SearchSpaceDims& ssd);
	DiscreteSearchSpace ();
	PointInds Point2Indexes(Point point);
	Point Indexes2Point(PointInds indexes);
	vector<Point> getSearchSpacePointsVec();
	bool IncreaseInd(PointInds& p, size_t k);
};

inline vector <double> getDimGrid(Dim d)
{
	vector<double> tmp_vec;
	for (size_t i = 0; (d.l + i*d.s) <= d.r; ++i)
		tmp_vec.push_back(d.l + i*d.s);
	
	return tmp_vec; // the compiler will optimize out the copy!
}
#endif
