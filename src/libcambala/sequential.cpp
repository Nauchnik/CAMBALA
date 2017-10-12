#include "cambala_sequential.h"
#include "cambala_utils.h"
#include <iostream>
#include <complex>
#include <time.h>
#include <stdexcept>

CAMBALA_sequential::CAMBALA_sequential() :
	launch_type("bruteforce"),
	object_function_type("uniform"),
	output_filename("cambala_out"),
	depths_filename("cambala_depths_out"),
	H(0),
	nh(1),
	ncb(1),
	nrhob(1),
	nR(1),
	ntau(1),
	cb1(2000.0),
	cb2(2000.0),
	R1(3400.0),
	R2(3600.0),
	tau1(0.0),
	tau2(0.0),
	rhob1(2.0),
	rhob2(2.0),
	n_layers_w(1),
	iterated_local_search_runs(10),
	verbosity(1),
	N_total(1),
	isTimeDelayPrinting(false),
	ppm(0),
	rank(0)
{
	srand((unsigned)time(NULL));
	start_chrono_time = chrono::high_resolution_clock::now();
}
/*
A routine for computing delay residual.

Arguments:
1) Environment: five arrays of the same length: depth, c1s, c2s, rhos, Ns_points;
(each entry describes one layer as described in the comments to compute_wnumbers_extrap() )

2) Source-receive distance: R -- distance from the source to the receiver;

3) Experimental data: modal delays:
-- experimental_mode_numbers: number of modes for each frequency in the recorded signal
-- experimental_delays: experimental_delays[ii][jj] is the delay of jj+1-th mode for the frequency freqs[ii]

The routine computes the (uniform) residual (misfit) of experimental data and the "theoretical" delays for a given environment model.

It should be used as follows: for a set of environment models the residual should be computed. The minimal value of the residual indicates
the most "adequate" model.
*/

void CAMBALA_sequential::load_layers_data(
    string LayersFName,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points)
{
    ifstream Myfile(LayersFName);
    double d, c1,c2,rho;
    unsigned nsp;

    depths.clear();
    c1s.clear();
    c2s.clear();
    rhos.clear();
    Ns_points.clear();

	while (!(Myfile.eof()))
	{
		Myfile >> d;
		Myfile >> c1;
		Myfile >> c2;
		Myfile >> rho;
		Myfile >> nsp;

        depths.push_back(d);
        c1s.push_back(c1);
        c2s.push_back(c2);
        rhos.push_back(rho);
        Ns_points.push_back(nsp);

	}
	Myfile.close();

}

void CAMBALA_sequential::load_profile_deep_water(
    string ProfileFName,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points)
{
    ifstream Myfile(ProfileFName);
    double cp, cc, dc, dp;

    Myfile >> dp;
    Myfile >> cp;

    unsigned npc;

	while (!(Myfile.eof()))
	{
		Myfile >> dc;
		Myfile >> cc;

        depths.push_back(dc);
        c1s.push_back(cp);
        c2s.push_back(cc);
        rhos.push_back(1);

        npc = (unsigned)ppm*(dc - dp);
        Ns_points.push_back(npc);

        cp = cc;
        dp = dc;

	}
	Myfile.close();
}
//tau_comment: search and output,
//check for tau!
int CAMBALA_sequential::init(vector<double> depths)
{
	if ((!rank) && (verbosity > 0)) {
		cout << "init() started" << endl;
		cout << "depths: ";
		for (auto &x : depths)
			cout << x << " ";
		cout << endl;
	}

	n_layers_w = depths.size() - 1;

	if (!n_layers_w) {
		cerr << "n_layers_w == 0" << endl;
		return -1;
	}

	ncpl_arr = ncpl_init_arr;
	ncpl_arr.resize(n_layers_w);
	cw1_arr = cw1_init_arr;
	cw1_arr.resize(n_layers_w);
	cw2_arr = cw2_init_arr;
	cw2_arr.resize(n_layers_w);

	c1s.resize(n_layers_w + 1);
	for (auto &x : c1s)
		x = 1500;
	c2s.resize(n_layers_w + 1);
	for (auto &x : c2s)
		x = 1500;
	rhos.resize(n_layers_w + 1);
	for (auto &x : rhos)
		x = 1;
	Ns_points.resize(n_layers_w + 1);
	Ns_points[0] = (unsigned)round(ppm*depths[0]);
	for (unsigned i=1; i < depths.size(); i++ )
		Ns_points[i] = (unsigned)round(ppm*(depths[i] - depths[i-1]));

	N_total = nR*nrhob*ncb*ntau;
	for (auto &x : ncpl_arr)
		N_total *= (unsigned long long)x;

	if (!N_total) {
		cerr << "N_total == 0" << endl;
		cerr << "nR nrhob ncb ntau : " << endl;
		cerr << nR << " " << nrhob << " " << ncb << " " << ntau << endl;
		cerr << "ncpl_arr : " << endl;
		for (auto &x : ncpl_arr)
			cerr << x << " ";
		cerr << endl;
		cerr << "depths : " << endl;
		for (auto &x : depths)
			cerr << x << " ";
		cerr << endl;
		return -1;
	}

	if ( (!rank) && (verbosity > 0) )
		cout << "N_total " << N_total << endl;

	if (cw1_arr.size() != (depths.size() - 1)) {
		cerr << "cw1_arr.size() != (depths.size() - 1)";
		cerr << endl;
		exit(1);
	}
	loadValuesToSearchSpaceVariables();

	if ( (!rank) && (verbosity > 0) )
		cout << "init() finished" << endl;

	return 0;
}

int CAMBALA_sequential::createDepthsArray(const double h, vector<vector<double>> &depths_vec)
{
	if (launch_type == "bruteforce") {
		n_layers_w = cw1_init_arr.size();
		double layer_thickness_w = h / n_layers_w;
		vector<double> depths;
		for (unsigned jj = 1; jj <= n_layers_w; jj++)
			depths.push_back(layer_thickness_w*jj);
		depths.push_back(H);
		depths_vec.push_back(depths);
	}
	else if (launch_type == "ils")
	{
		if (d1_arr.size() == 0) {
			n_layers_w = cw1_arr.size();
			double layer_thickness_w = h / n_layers_w;
			vector<double> depths;
			for (unsigned jj = 1; jj <= n_layers_w; jj++)
				depths.push_back(layer_thickness_w*jj);
			depths.push_back(H);
			depths_vec.push_back(depths);
		}
		else
		{
			vector<vector<double>> search_space_depths;
			search_space_depths.resize(d1_arr.size());
			for (unsigned i = 0; i < d2_arr.size(); i++) {
				double cur_val = d2_arr[i];
				for (;;) {
					search_space_depths[i].push_back(cur_val);
					cur_val -= d_step[i];
					if (cur_val < d1_arr[i])
						break;
				}
			}

			vector<int> index_arr;
			vector<double> tmp_depths;
			vector<vector<double>> ::iterator it;
			double cur_treshold;
			while (CAMBALA_utils::next_cartesian(search_space_depths, index_arr, tmp_depths))
			{
				vector<double> depths;
				cur_treshold = tmp_depths[0] + 3;
				depths.push_back(tmp_depths[0]); // at least 1 water layer must exist
				for (unsigned i = 1; i < tmp_depths.size(); i++) {
					if (tmp_depths[i] >= cur_treshold) {
						depths.push_back(tmp_depths[i]);
						cur_treshold = tmp_depths[i] + 2;
					}
				}
				it = find(depths_vec.begin(), depths_vec.end(), depths);
				if (it == depths_vec.end())
					depths_vec.push_back(depths);
			}

			for (auto &x : depths_vec) {
				x.push_back(h);
				x.push_back(H);
			}
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

	return 0;
}

void CAMBALA_sequential::reportFinalResult()
{
	// fix final time
	chrono::high_resolution_clock::time_point t2;
	chrono::duration<double> time_span;
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - start_chrono_time);

	ofstream ofile(output_filename, ios_base::app);

	ofile << endl;
	ofile << "total solving time (chrono) " << time_span.count() << endl;
	ofile << "SEARCH ENDED!" << endl;
	ofile << "RESULTING VALUE:" << endl;
	ofile << "err = " << record_point.residual << ", parameters:" << endl;
	ofile << "c_b = " << record_point.cb << endl
			  << "tau = " << record_point.tau << endl
			  << "rho_b = " << record_point.rhob << endl
			  << "R = " << record_point.R << endl;
	ofile << "cws :" << endl;
	for (auto &x : record_point.cws)
		ofile << x << " ";
	ofile << endl;
	ofile << "depths " << endl;
	for (auto &x : record_point.depths)
		ofile << x << " ";
	ofile << endl;
	ofile << "total solving time " << time_span.count() << endl;

	ofile.close();
}

vector<search_space_point> CAMBALA_sequential::getSearchSpacePointsVec(vector<double> depths)
{
	vector<int> index_arr;
	vector<unsigned> cur_point_indexes;
	vector<vector<unsigned>> search_space_indexes;
	search_space_indexes.resize(search_space.size());
	for (unsigned i = 0; i < search_space.size(); i++)
		for (unsigned j = 0; j < search_space[i].size(); j++)
			search_space_indexes[i].push_back(j);

	vector<search_space_point> points_vec;
	while (CAMBALA_utils::next_cartesian(search_space_indexes, index_arr, cur_point_indexes))
		points_vec.push_back(fromPointIndexesToPoint(cur_point_indexes, depths));

	return points_vec;
}

void CAMBALA_sequential::reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a)
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4, ...] - cws
	if (reduced_s_s_a.cb == false) {
		search_space[0].resize(1);
		search_space[0][0] = cb1;
	}
	if (reduced_s_s_a.rhob == false) {
		search_space[1].resize(1);
		search_space[1][0] = rhob1;
	}
	if (reduced_s_s_a.R == false) {
		search_space[2].resize(1);
		search_space[2][0] = R1;
	}
	if (reduced_s_s_a.tau == false) {
		search_space[3].resize(1);
		search_space[3][0] = tau1;
	}
	for (unsigned i=0; i < reduced_s_s_a.cws.size(); i++) {
		if (reduced_s_s_a.cws[i] == false) {
			search_space[4 + i].resize(1);
			search_space[4 + i][0] = cw1_arr[i];
		}
	}
}


void CAMBALA_sequential::loadValuesToSearchSpaceVariables()
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws
	vector<double> tmp_vec;
	search_space.clear();

	// fill search_space_variables[0] with cb
	tmp_vec.resize(ncb);
	for (unsigned long long i = 0; i < ncb; i++)
		tmp_vec[i] = cb1 + (ncb == 1 ? 0 : i*(cb2 - cb1) / (ncb - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[1] with rhob
	tmp_vec.resize(nrhob);
	for (unsigned long long i = 0; i < nrhob; i++)
		tmp_vec[i] = rhob1 + (nrhob == 1 ? 0 : i*(rhob2 - rhob1) / (nrhob - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[2] with R
	tmp_vec.resize(nR);
	for (unsigned long long i = 0; i < nR; i++)
		tmp_vec[i] = R1 + (nR == 1 ? 0 : i*(R2 - R1) / (nR - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[3] with tau
	tmp_vec.resize(ntau);
	for (unsigned long long i = 0; i < ntau; i++)
		tmp_vec[i] = tau1 + (ntau == 1 ? 0 : i*(tau2 - tau1) / (ntau - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[4...] with cws
	for (unsigned long long i = 0; i < cw1_arr.size(); i++) {
		tmp_vec.resize(ncpl_arr[i]);
		for (unsigned long long j = 0; j < ncpl_arr[i]; j++)
			tmp_vec[j] = cw1_arr[i] + (ncpl_arr[i] == 1 ? 0 : j*(cw2_arr[i] - cw1_arr[i]) / (ncpl_arr[i] - 1));
		search_space.push_back(tmp_vec);
	}

	if (!rank)
		cout << "loadValuesToSearchSpaceVariables() finished" << endl;
}

search_space_point CAMBALA_sequential::fromPointIndexesToPoint( vector<unsigned> cur_point_indexes,
	                                                            vector<double> depths)
{
	search_space_point point;
	point.cb   = search_space[0][cur_point_indexes[0]];
	point.rhob = search_space[1][cur_point_indexes[1]];
	point.R    = search_space[2][cur_point_indexes[2]];
	point.tau  = search_space[3][cur_point_indexes[3]];
	for (unsigned i = 4; i < search_space.size(); i++)
		point.cws.push_back(search_space[i][cur_point_indexes[i]]);
	point.depths = depths;
	point.residual = START_HUGE_VALUE;
	return point;
}

vector<unsigned> CAMBALA_sequential::fromPointToPointIndexes( search_space_point point )
{
	vector<unsigned> cur_point_indexes;
	cur_point_indexes.resize(search_space.size());

	for (unsigned j = 0; j<search_space[0].size(); j++)
		if ( search_space[0][j] == point.cb ) {
			cur_point_indexes[0] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[1].size(); j++)
		if (search_space[1][j] == point.rhob) {
			cur_point_indexes[1] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[2].size(); j++)
		if (search_space[2][j] == point.R) {
			cur_point_indexes[2] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[3].size(); j++)
		if (search_space[3][j] == point.tau) {
			cur_point_indexes[3] = j;
			break;
		}
	for (unsigned i = 4; i < search_space.size(); i++) {
		for (unsigned j = 0; j < search_space[i].size(); j++)
			if (search_space[i][j] == point.cws[i - 4]) {
				cur_point_indexes[i] = j;
				break;
			}
	}

	return cur_point_indexes;
}

search_space_point CAMBALA_sequential::fromDoubleVecToPoint(vector<double> double_vec)
{
	search_space_point point;
	point.cb   = double_vec[0];
	point.rhob = double_vec[1];
	point.R    = double_vec[2];
	point.tau  = double_vec[3];
	for (unsigned i = 4; i < double_vec.size(); i++)
		point.cws.push_back(double_vec[i]);
	point.residual = START_HUGE_VALUE;
	return point;
}

// function for BOINC client application
search_space_point CAMBALA_sequential::fromStrToPoint(string str)
{
	search_space_point point;
	stringstream sstream;
	sstream << str;
	sstream >> point.residual >> point.cb >> point.rhob >> point.R >> point.tau;
	double val;
	while (sstream >> val)
		point.cws.push_back(val);
	return point;
}

void CAMBALA_sequential::fromPointToFile(const search_space_point &point, ofstream &ofile)
{
	ofile << point.residual << " " << point.cb << " " << point.rhob << " "
		<< point.R << " " << point.tau << " ";
	for (unsigned i = 0; i < point.cws.size(); i++)
		ofile << point.cws[i] << " ";
}

double CAMBALA_sequential::getRecordResidual()
{
	return record_point.residual;
}

void CAMBALA_sequential::getThreeValuesFromStr(string str, double &val1, double &val2, double &val3)
{
	val1 = val3 = -1;
	val2 = 1;
	string word1, word2, word3;
	for (auto &x : str)
		if (x == ':')
			x = ' ';
	stringstream sstream;
	sstream << str;
	sstream >> word1 >> word2 >> word3;
	istringstream(word1) >> val1;
	istringstream(word2) >> val2;
	istringstream(word3) >> val3;
	if (val3 == -1)
		val3 = val1;
}

int CAMBALA_sequential::readScenario(string scenarioFileName)
{
// read constant and variable values from a scenario file
	if ( (!rank) && (verbosity > 0) )
		cout << "scenarioFileName " << scenarioFileName << endl;
	ifstream scenarioFile(scenarioFileName.c_str());

	if (!scenarioFile.is_open()) {
		cerr << "scenarioFile with the name " << scenarioFileName << " wasn't openend" << endl;
		return -1;
	}

	string str, word, tmp_word;
	stringstream sstream;
	unsigned cw_index = 0, d_index = 0;
	double cur_val_step = 0, cur_val1 = 0, cur_val2 = 0;
	while (getline(scenarioFile, str)) {
		if ((str == "") || (str[0] == '%'))
			continue;
		sstream << str;
		sstream >> word;
		if (word.find("dtimes_file") != string::npos)
			sstream >> dtimesFileName;
		else if (word.find("spmag_file") != string::npos)
			sstream >> spmagFileName;
		else if (word == "H")
			sstream >> H;
		else if (word == "h") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			h1 = cur_val1;
			h2 = cur_val2;
			if (h1 == h2)
				nh = 1;
			else
				nh = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if ((word.size() >= 2) && (word[0] == 'c') && (word[1] == 'w')) {
			word = word.substr(2, word.size() - 2);
			istringstream(word) >> cw_index;
			if (cw1_init_arr.size() < cw_index + 1)
				cw1_init_arr.resize(cw_index + 1);
			if (cw2_init_arr.size() < cw_index + 1)
				cw2_init_arr.resize(cw_index + 1);
			if (ncpl_init_arr.size() < cw_index + 1)
				ncpl_init_arr.resize(cw_index + 1);
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			cw1_init_arr[cw_index] = cur_val1;
			cw2_init_arr[cw_index] = cur_val2;
			if (cur_val1 == cur_val2)
				ncpl_init_arr[cw_index] = 1;
			else
				ncpl_init_arr[cw_index] = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if ((word.size() == 2) && (word[0] == 'd') && (isdigit(word[1]))) {
			word = word.substr(1, word.size() - 1);
			istringstream(word) >> d_index;
			d_index--;
			if (d1_arr.size() < d_index + 1) {
				d1_arr.resize(d_index + 1);
				d2_arr.resize(d_index + 1);
				d_step.resize(d_index + 1);
			}
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			d1_arr[d_index] = cur_val1;
			d2_arr[d_index] = cur_val2;
			d_step[d_index] = cur_val_step;
		}
		else if (word == "R") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			R1 = cur_val1;
			R2 = cur_val2;
			if (R1 == R2)
				nR = 1;
			else
				nR = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "rhob") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			rhob1 = cur_val1;
			rhob2 = cur_val2;
			if (rhob1 == rhob2)
				nrhob = 1;
			else
				nrhob = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "cb") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			cb1 = cur_val1;
			cb2 = cur_val2;
			if (cb1 == cb2)
				ncb = 1;
			else
				ncb = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "tau") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			tau1 = cur_val1;
			tau2 = cur_val2;
			if (tau1 == tau2)
				ntau = 1;
			else
				ntau = abs(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "ils_iterations")
			sstream >> iterated_local_search_runs;
		else if (word == "function_type")
			sstream >> object_function_type;
		else if (word == "ppm")
			sstream >> ppm;
		sstream.str(""); sstream.clear();
	}
	
	if (ppm == 0) { // if ppm wasn't set directly
		if ((object_function_type == "uniform") || (object_function_type == "weighted"))
			ppm = 2;
		else if ((object_function_type == "uniform2") || (object_function_type == "weighted2"))
			ppm = 1;
		else {
			cerr << "unknown object_function_type " << object_function_type << endl;
			exit(1);
		}
	}

	if (!cw1_init_arr.size()) {
		cerr << "!cw1_init_arr.size()" << endl;
		return -1;
	}
	if (!H) {
		cerr << "!H" << endl;
		return -1;
	}
	cw1_arr = cw1_init_arr;
	cw2_arr = cw2_init_arr;

	input_params_sstream << "Parameters :" << endl;
	input_params_sstream << "launch_type " << launch_type << endl;
	input_params_sstream << "object_function_type " << object_function_type << endl;
	input_params_sstream << "ppm " << ppm << endl;
	input_params_sstream << "iterated_local_search_runs " << iterated_local_search_runs << endl;
	input_params_sstream << "cw1_init_arr :" << endl;
	for (auto &x : cw1_init_arr)
		input_params_sstream << x << " ";
	input_params_sstream << endl;
	input_params_sstream << "cw2_init_arr :" << endl;
	for (auto &x : cw2_init_arr)
		input_params_sstream << x << " ";
	input_params_sstream << endl;
	input_params_sstream << "ncpl_init_arr :" << endl;
	for (auto &x : ncpl_init_arr)
		input_params_sstream << x << " ";
	input_params_sstream << endl;
	input_params_sstream << "nR " << nR << endl;
	input_params_sstream << "R1 " << R1 << endl;
	input_params_sstream << "R2 " << R2 << endl;
	input_params_sstream << "ntau " << ntau << endl;
	input_params_sstream << "tau1 " << tau1 << endl;
	input_params_sstream << "tau2 " << tau2 << endl;
	input_params_sstream << "nrhob " << nrhob << endl;
	input_params_sstream << "rhob1 " << rhob1 << endl;
	input_params_sstream << "rhob2 " << rhob2 << endl;
	input_params_sstream << "ncb " << ncb << endl;
	input_params_sstream << "cb1 " << cb1 << endl;
	input_params_sstream << "cb2 " << cb2 << endl;
	input_params_sstream << "dtimes_file " << dtimesFileName << endl;
	input_params_sstream << "spmag_file " << spmagFileName << endl;

	if (!rank) {
		ofstream ofile(output_filename, ios_base::out);
		ofile << input_params_sstream.str();
		ofile.close();
	}

	if (((rank == 0) || (rank == 1)) && (verbosity > 0))
		cout << "readScenario() finished" << endl;

	return 0;
}

int CAMBALA_sequential::readInputDataFromFiles()
{
	ifstream dtimesFile(dtimesFileName.c_str());
	if (!dtimesFile.is_open()) {
		cerr << "dtimesFile " << dtimesFileName << " wasn't opened" << endl;
		return -1;
	}
	stringstream myLineStream;
	string myLine;
	double buff;
	vector<double> buffvect;
	// reading the "experimental" delay time data from a file
	while (getline(dtimesFile, myLine)) {
		myLine.erase(remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete windows endline symbol for correct reading
		myLineStream << myLine;
		myLineStream >> buff;
		freqs.push_back(buff);

		buffvect.clear();
		while (!myLineStream.eof()) {
			myLineStream >> buff;
			buffvect.push_back(buff);
		}

		mode_numbers.push_back((unsigned)buffvect.size());
		modal_delays.push_back(buffvect);
		myLineStream.str(""); myLineStream.clear();
	}
	dtimesFile.close();

	if (object_function_type.find("weighted") != string::npos) {
		weight_coeffs.clear();
		ifstream spmagFile(spmagFileName.c_str());
		if (spmagFile.is_open()) {
			buffvect.clear();
			while (getline(spmagFile, myLine)) {
				myLine.erase(remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete windows endline symbol for correct reading
				myLineStream << myLine;
				myLineStream >> buff;

				buffvect.clear();
				while (!myLineStream.eof()) {
					myLineStream >> buff;
					buffvect.push_back(buff);
				}

				weight_coeffs.push_back(buffvect);
				myLineStream.str(""); myLineStream.clear();
			}
			spmagFile.close();

			if ((!rank) && (verbosity > 0)) {
				cout << "weight_coeffs.size() " << weight_coeffs.size() << endl;
				cout << "weight_coeffs first 10 lines : " << endl;
				for (unsigned i = 0; i < 10; i++) {
					for (auto &x : weight_coeffs[i])
						cout << x << " ";
					cout << endl;
				}
			}
		}
		else {
			cerr << "spmagfile " << spmagFileName << " wasn't opened" << endl;
			return -1;
		}
	}

	if ( (!rank) && (verbosity > 0) )
		cout << "readInputDataFromFiles() finished " << endl;
	return 0;
}

search_space_point CAMBALA_sequential::getNonRandomStartPoint( vector<double> depths )
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws

	search_space_point point;

	point.cb = (cb1 == cb2) ? cb1 : (cb2 - cb1) / 2;
	point.rhob = (rhob1 == rhob2) ? rhob1 : (rhob2 - rhob1) / 2;
	point.R = (R1 == R2) ? R1 : (R2 - R1) / 2;
	point.tau = (tau1 == tau2) ? tau1 : (tau2 - tau1) / 2;

	point.cws.resize(cw1_arr.size());
	for ( unsigned i = 0; i < point.cws.size(); i++ ) {
		if (cw1_arr[i] == cw2_arr[i])
			point.cws[i] = cw1_arr[i];
		else if ((START_CW_VALUE >= cw1_arr[i]) && (START_CW_VALUE <= cw2_arr[i])) // if cws should be modified
			point.cws[i] = START_CW_VALUE;
		else
			point.cws[i] = (cw2_arr[i] - cw1_arr[i]) / 2;
	}

	point.depths = depths;

	return point;
}

void CAMBALA_sequential::printDelayTime(double R, vector<unsigned> mode_numbers, vector<vector<double>> modal_group_velocities)
{
	string ofile_name = "delayTimeOutput_" + object_function_type + ".txt";
	ofstream ofile(ofile_name);
	for (unsigned ii = 0; ii < freqs.size(); ii++) {
		unsigned mnumb = mode_numbers.at(ii);
		ofile << freqs.at(ii) << "\t";
		for (unsigned jj = 0; jj < mnumb; jj++)
			ofile << R / modal_group_velocities[ii][jj] << "\t";
		ofile << endl;
	}
	ofile.close();
}

double CAMBALA_sequential::directPointCalc( search_space_point point )
{
	isTimeDelayPrinting = true;
	return fillDataComputeResidual(point);
}




void CAMBALA_sequential::PrintPoint(const search_space_point& point)
{
	cout << endl;
	cout << endl << "New residual minimum:" << endl;
	cout << "err = " << point.residual << ", parameters:" << endl;
	cout << "c_b = " << point.cb
		<< ", rho_b = " << point.rhob
		<< ", tau = " << point.tau
		<< ", R = " << point.R << endl;
	cout << "cws_min :" << endl;
	for (auto x : point.cws)
		cout << x << " ";
	cout << endl;
	cout << "depths " << endl;
	for (auto x : point.depths)
		cout << x << " ";
	cout << endl;
	cout << "Ns_points " << endl;
	for (auto x : Ns_points)
		cout << x << " ";
	cout << endl;
	cout << endl;
}


void CAMBALA_sequential::solve()
{
	search_space_point record_point;
	record_point.cb       = START_HUGE_VALUE;
	record_point.rhob     = START_HUGE_VALUE;
	record_point.R        = START_HUGE_VALUE;
	record_point.tau      = START_HUGE_VALUE;
	record_point.residual = START_HUGE_VALUE;

	for (unsigned i = 0; i < nh; i++)
	{
		double cur_h = h1 + (nh == 1 ? 0 : i*(h2 - h1) / (nh - 1));
		vector<vector<double>> depths_vec;
		createDepthsArray(cur_h, depths_vec);
		if (launch_type == "ils")
		{
			for (unsigned j = 0; j < depths_vec.size(); j++)
			{
				init(depths_vec[j]);
				//CAMBALA_seq.findGlobalMinBruteForce(depths_vec[i]);
				findLocalMinHillClimbing(depths_vec[j]);
				cout << "Processed " << j + 1 << " out of " << depths_vec.size() << " depths" << endl;
				// Moved from ComputeResidual
				if (point.residual < record_point.residual)
				{
					record_point = point;
					if (verbosity > 0)
						PrintPoint(record_point);
				}
			}
		}
		else
		{
			init(depths_vec[0]);
			findGlobalMinBruteForce(depths_vec[0]);
		}
		cout << "Processed " << i + 1 << " out of " << nh << " h (max depths)" << endl;
