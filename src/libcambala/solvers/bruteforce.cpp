#include "solvers/bruteforce.h"
#include <iostream>
#include <algorithm>

void BruteForce::Solve()
{
	cout << "findGlobalMinBruteForce()" << endl;

	vector<Point> Points_vec = getSearchSpacePointsVec(depths_);
	cout << "Points_vec.size() " << Points_vec.size() << endl;

	for (auto &x : Points_vec)
		fillDataComputeResidual(x); // calculated residual is written to cur_point
	recordPoint_ = *std::min_element(Points_vec.begin(), Points_vec.end());
}
