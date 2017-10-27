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
	vector<vector<double>> Axes_;// values of variables which form the search space
	DiscreteSearchSpace (const SearchSpaceDims& ssd);
	PointInds Point2Indexes(Point point);
	vector<Point> getSearchSpacePointsVec();
};

inline vector <double> getDimGrid(Dim d)
{
	vector<double> tmp_vec;
	for (size_t i = 0; (d.l + i*d.s) <= d.r; ++i)
		tmp_vec.push_back(d.l + i*d.s);
	
	return tmp_vec; // the compiler will optimize out the copy!
}
#endif
