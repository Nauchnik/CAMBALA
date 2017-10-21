
vector<vector<double>> CAMBALA_sequential::createDepthsArray(const double h)
{
	//layer_thickness_w
	//H_
	//h
	//d1_arr
	//
	vector <vector <double>> depths_vec;
	if (d1_arr.size() == 0)
	{
		size_t n_layers_w = cw1_arr.size();
		double layer_thickness_w = h / n_layers_w;
		vector <double> depths;
		for (unsigned jj = 1; jj <= n_layers_w; jj+)
			depths.push_back(layer_thickness_w*jj);
		depths.push_back(H_);
		depths_vec.push_back(depths);
	}
	else
	{
		vector<vector<double>> search_space_depths (d1_arr.size());
		for (unsigned i = 0; i < d2_arr.size(); i++)
		{
			double cur_val = d2_arr[i];
			for (;;)
			{
				search_space_depths[i].push_back(cur_val);
				cur_val -= d_step[i];
				if (cur_val < d1_arr[i])
					break;
			}
		}

		vector<int> index_arr;
		vector<double> tmp_depths;
		vector<vector<double>> ::iterator it;
		while (CAMBALA_utils::next_cartesian(search_space_depths, index_arr, tmp_depths))
		{
			vector<double> depths;
			double cur_treshold = tmp_depths[0] + 3;
			depths.push_back(tmp_depths[0]); // at least 1 water layer must exist
			for (unsigned i = 1; i < tmp_depths.size(); i++)
			{
				if (tmp_depths[i] >= cur_treshold)
				{
					depths.push_back(tmp_depths[i]);
					cur_treshold = tmp_depths[i] + 2;
				}
			}
			it = find(depths_vec.begin(), depths_vec.end(), depths);
			if (it == depths_vec.end())
				depths_vec.push_back(depths);
		}

		for (auto &x : depths_vec)
		{
			x.push_back(h);
			x.push_back(H_);
		}
	}

	ofstream ofile(depths_filename);
	for (auto &x : depths_vec) {
		for (auto &y : x)
			ofile << y << " ";
		ofile << endl;
	}
	ofile.close();
	cout << "depths_vec.size() " << depths_vec.size() << endl;

	return depths_vec;
}
