
void CAMBALA_sequential::findGlobalMinBruteForce(vector<double> depths)
{
	cout << "findGlobalMinBruteForce()" << endl;

	vector<search_space_point> search_space_points_vec = getSearchSpacePointsVec(depths);
	cout << "search_space_points_vec.size() " << search_space_points_vec.size() << endl;

	for (auto &x : search_space_points_vec)
		fillDataComputeResidual(x); // calculated residual is written to cur_point
}
