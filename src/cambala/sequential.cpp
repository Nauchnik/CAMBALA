#include "sequential.h"
#include <iostream>
#include <complex>
#include <time.h>
#include <stdexcept>
#include "utils.h"
#include "compute.h"

using namespace CAMBALA_compute;
using namespace CAMBALA_utils;

CAMBALA_sequential::CAMBALA_sequential() :
	launch_type("bruteforce"),
	object_function_type("uniform"),
	output_filename("cambala_out"),
	depths_filename("cambala_depths_out"),
	dtimesFileName(""),
	spmagFileName(""),
	H(0),
	nh(1),
	n_layers_w(1),
	iterated_local_search_runs(10),
	init_iterated_local_search_runs(10),
	verbosity(1),
	N_total(1),
	isTimeDelayPrinting(false),
	ppm(0),
	is_mpi(false),
	rank(0)
{
	record_point.cb       = START_HUGE_VALUE;
	record_point.rhob     = START_HUGE_VALUE;
	record_point.R        = START_HUGE_VALUE;
	record_point.tau      = START_HUGE_VALUE;
	record_point.residual = START_HUGE_VALUE;
	srand((unsigned)time(NULL));
	start_chrono_time = chrono::high_resolution_clock::now();
#ifdef _MPI
	is_mpi = true;
#endif
}

// tau_comment: search and output,
// check for tau!
int CAMBALA_sequential::init(vector<double> depths)
{
	if ((!rank) && (verbosity > 0)) {
		cout << "init() started" << endl;
		cout << "depths: \n" << doubleVecToStr(depths) << endl;
	}

	n_layers_w = depths.size() - 1;

	if (!n_layers_w) {
		cerr << "n_layers_w == 0" << endl;
		return -1;
	}

	cw_vec_vec = cw_init_vec_vec;
	cw_vec_vec.resize(n_layers_w);

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
	
	N_total = R_vec.size()*rhob_vec.size()*cb_vec.size()*tau_vec.size();
	for (auto &cw_vec : cw_vec_vec)
		N_total *= (unsigned long long)cw_vec.size();

	if (!N_total) {
		cerr << "N_total == 0" << endl;
		cerr << "depths : " << doubleVecToStr(depths) << endl;
		return -1;
	}

	if ( (!rank) && (verbosity > 0) )
		cout << "N_total " << N_total << endl;

	if (cw_vec_vec.size() != (depths.size() - 1)) {
		cerr << "cw_vec_vec.size() != (depths.size() - 1)";
		cerr << endl;
		exit(1);
	}
	loadValuesToSearchSpaceVariables();

	if ( (!rank) && (verbosity > 0) )
		cout << "init() finished" << endl;

	return 0;
}

void CAMBALA_sequential::writeOutputData(stringstream &sstream)
{
	ofstream ofile(output_filename, ios_base::app);
	ofile << sstream.rdbuf();
	sstream.clear(); sstream.str("");
	ofile.close(); ofile.clear();
}

vector<vector<double>> CAMBALA_sequential::createDepthsArray()
{
	vector<vector<double>> depths_vec;
	if (d1_vec.size() == 0) { // static depths mode
		n_layers_w = cw_init_vec_vec.size();
		double layer_thickness_w = h / n_layers_w;
		vector<double> depths;
		for (unsigned jj = 1; jj <= n_layers_w; jj++)
			depths.push_back(layer_thickness_w*jj);
		depths.push_back(H);
		depths_vec.push_back(depths);
	}
	else { // dynamic depths mode
		vector<vector<double>> search_space_depths;
		search_space_depths.resize(d1_vec.size());
		for (unsigned i = 0; i < d2_vec.size(); i++) {
			double cur_val = d2_vec[i];
			for (;;) {
				search_space_depths[i].push_back(cur_val);
				cur_val -= d_step_vec[i];
				if (cur_val < d1_vec[i])
					break;
			}
		}

		vector<int> index_arr;
		vector<double> tmp_depths;
		vector<vector<double>> ::iterator it;
		double cur_treshold;
		while (next_cartesian(search_space_depths, index_arr, tmp_depths)) {
			vector<double> depths;
			cur_treshold = tmp_depths[0] + 3;
			depths.push_back(tmp_depths[0]); // at least 1 water layer must exist
			for (unsigned i = 1; i < tmp_depths.size(); i++) {
				if (tmp_depths[i] >= cur_treshold) {
					depths.push_back(tmp_depths[i]);
					cur_treshold = tmp_depths[i] + 2;
				}
			}
			if (depths.size() < d1_vec.size()) // skip short depths
				continue; 
			it = find(depths_vec.begin(), depths_vec.end(), depths);
			if (it == depths_vec.end())
				depths_vec.push_back(depths);
		}

		for (auto &x : depths_vec) {
			x.push_back(h);
			x.push_back(H);
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

void CAMBALA_sequential::reportFinalResult()
{
	// fix final time
	chrono::high_resolution_clock::time_point t2;
	chrono::duration<double> time_span;
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - start_chrono_time);

	stringstream sstream;
	
	sstream << endl;
	sstream << "SEARCH ENDED!" << endl;
	sstream << "SOLVING TIME " << time_span.count() << endl;
	sstream << "RESULTING VALUE:" << endl;
	sstream << strPointData(record_point);
	
	cout << sstream.str();

	writeOutputData(sstream);
}

void CAMBALA_sequential::findGlobalMinBruteForce(vector<double> depths)
{
	cout << "Brute force mode" << endl;

	vector<search_space_point> search_space_points_vec = getSearchSpacePointsVec(depths);
	cout << search_space_points_vec.size() << " points in the search space" << endl;

	int processed_points = 0;
	for (auto &x : search_space_points_vec) {
		fillDataComputeResidual(x); // calculated residual is written to cur_point
		processed_points++;
		if ((verbosity>0) && (processed_points % 10 == 0) ) {
			cout << processed_points << " out of " << search_space_points_vec.size() << " points have been processed \n";
		}
	}
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
	while (next_cartesian(search_space_indexes, index_arr, cur_point_indexes))
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
		search_space[0][0] = cb_vec[0];
	}
	if (reduced_s_s_a.rhob == false) {
		search_space[1].resize(1);
		search_space[1][0] = rhob_vec[0];
	}
	if (reduced_s_s_a.R == false) {
		search_space[2].resize(1);
		search_space[2][0] = R_vec[0];
	}
	if (reduced_s_s_a.tau == false) {
		search_space[3].resize(1);
		search_space[3][0] = tau_vec[0];
	}
	for (unsigned i=0; i < reduced_s_s_a.cws.size(); i++) {
		if (reduced_s_s_a.cws[i] == false) {
			search_space[4 + i].resize(1);
			search_space[4 + i][0] = cw_vec_vec[i][0];
		}
	}
}

double CAMBALA_sequential::fillDataComputeResidual( search_space_point &point )
{ // finally specify sound speed in water
  // the parameters are transformed into the arrays c1s, c2s, rhos
	if (verbosity > 1)
		cout << "fillDataComputeResidual()" << endl;
	if (point.cws.size() != point.depths.size() - 1) {
		cerr << "point.cws.size() != point.depths.size() - 1" << endl;
		cerr << point.cws.size() << " " << point.depths.size() - 1 << endl;
		exit(1);
	}
	for (unsigned jj = 0; jj < n_layers_w - 1; jj++) {
		c1s.at(jj) = point.cws.at(jj);
		c2s.at(jj) = point.cws.at(jj + 1);
		rhos.at(jj) = 1;
	}
	c1s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	c2s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	rhos.at(n_layers_w - 1) = 1;
	c1s.at(n_layers_w) = point.cb;
	c2s.at(n_layers_w) = point.cb;
	rhos.at(n_layers_w) = point.rhob;
	vector<double> depths = point.depths;
	if (depths.size() == 0) {
		cerr << "depths.size() == 0" << endl;
		exit(-1);
	}

	if (verbosity > 1) {
		/*for (unsigned jj = 0; jj <= n_layers_w; jj++)
			cout << "Layer #" << jj + 1 << ": c=" << c1s.at(jj) << "..." << c2s.at(jj) << "; rho=" << rhos.at(jj) << "; np=" << Ns_points.at(jj) << endl;
		cout << residual << endl << endl;*/
		cout << "depths : ";
		for (auto &x : depths)
			cout << x << " ";
		cout << endl;
		cout << "Ns_points : ";
		for (auto &x : Ns_points)
			cout << x << " ";
		cout << endl;
	}

	if (object_function_type == "uniform") {
		point.residual = compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, mode_numbers);
	}
	else if (object_function_type == "uniform2") {
		point.residual = compute_modal_delays_residual_uniform2(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, mode_numbers);
	}
	else if (object_function_type == "Wan_uniform") {
		point.residual = compute_modal_delays_residual_LWan(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, mode_numbers);
	}
	else if (object_function_type == "weighted") {
		point.residual = compute_modal_delays_residual_weighted(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, weight_coeffs, mode_numbers);
	}
	else if (object_function_type == "weighted2") {
		point.residual = compute_modal_delays_residual_weighted2(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, weight_coeffs, mode_numbers);
	}
	else {
		cerr << "unknown object_function_type " << object_function_type << endl;
		exit(1);
	}

	if ( verbosity > 1 )
		cout << "point.residual " << point.residual << endl;

	if (point < record_point) {
		record_point = point;
		if ((!is_mpi) && (verbosity > 0)) {
			cout << endl;
			cout << endl << "New residual minimum:" << endl;
			cout << strPointData(record_point);
		}
	}

	return point.residual;
}

void CAMBALA_sequential::updateRecordPoint(const search_space_point point)
{
	if ( (isCorrectCalculatedPoint(point)) && (point < record_point) )
		record_point = point;
}

void CAMBALA_sequential::loadValuesToSearchSpaceVariables()
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws
	search_space.clear();

	// fill search_space_variables[0] with cb
	search_space.push_back(cb_vec);
	// fill search_space_variables[1] with rhob
	search_space.push_back(rhob_vec);
	// fill search_space_variables[2] with R
	search_space.push_back(R_vec);
	// fill search_space_variables[3] with tau
	search_space.push_back(tau_vec);
	// fill search_space_variables[4...] with cws
	for (auto cw_vec : cw_vec_vec)
		search_space.push_back(cw_vec);
}

search_space_point CAMBALA_sequential::findLocalMinHillClimbing(vector<double> depths)
{
	//if (verbosity > 1)
	//	cout << "findLocalMinHillClimbing" << endl;
	unsigned u_val = 1;
	for (unsigned i = 0; i < depths.size() - 2; i++)
		u_val *= (unsigned)depths[i];
	srand(u_val);

	// choose random point in the search space
	/*for (unsigned i = 0; i < search_space.size(); i++) // i stands for variable_index
	local_record_point_indexes[i] = rand() % search_space[i].size(); // get random index
	cur_point_indexes = local_record_point_indexes;*/

	search_space_point local_record_point = getNonRandomStartPoint(depths);
	vector<unsigned> local_record_point_indexes = fromPointToPointIndexes( local_record_point );
	fillDataComputeResidual(local_record_point); // calculated residual is written to cur_point

	bool isCheckRequired = false;
	for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
		if (search_space[i].size() > 1) {
			isCheckRequired = true;
			break;
		}
	}
	if ( (!isCheckRequired) && (verbosity > 0) ) {
		cout << "1 element in a search space, fast exit" << endl;
		return local_record_point;
	}

	checked_points.reserve(N_total);
	checked_points.push_back(local_record_point);
	unsigned skipped_points = 0;
	bool isContinueDimension;
	vector<unsigned> cur_point_indexes;
	search_space_point cur_point;
	// launch iterations of hill climbing
	for (unsigned run_index = 0; run_index < iterated_local_search_runs; run_index++) {
		bool isLocalMin;
		//if (verbosity > 1)
		//	cout << "iteration " << run_index << " of ILS" << endl;
		do { // do while local min not reached
			isLocalMin = true; // if changing of every variable will not lead to a record updata, then a local min reached
			for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
				if (search_space[i].size() == 1) {
					//cout << "one value of a variable, skip it" << endl;
					continue;
				}
				//cout << "variable_index " << variable_index << endl;
				cur_point_indexes = local_record_point_indexes;
				vector<unsigned> point_indexes_before_increase = cur_point_indexes;
				unsigned index_from = cur_point_indexes[i]; // don't check index twice
				//if (verbosity > 0)
				//	cout << "index_from " << index_from << endl;
				bool isDecreaseTurn = false;
				bool isTriggerIncToDec = false;
				do { // change value of a variabble while it leads to updating of a record
					double old_record_residual = local_record_point.residual;
					if (isDecreaseTurn) {
						if (isTriggerIncToDec) { // move to a point before increase
							cur_point_indexes = point_indexes_before_increase;
							isTriggerIncToDec = false;
						}
						if (cur_point_indexes[i] == 0)
							cur_point_indexes[i] = search_space[i].size() - 1;
						else
							cur_point_indexes[i]--;
					}
					else {
						cur_point_indexes[i]++;
						if (cur_point_indexes[i] == search_space[i].size())
							cur_point_indexes[i] = 0;
					}
					if (cur_point_indexes[i] == index_from) {
						//if (verbosity > 0)
						//	cout << "cur_point_indexes[i] == index_from. Break iteration." << endl;
						break;
					}
					/*if (verbosity > 0) {
						cout << "checking index " << cur_point_indexes[i] <<
							", max index " << search_space[i].size() - 1 << endl;
						cout << "cur_point_indexes" << endl;
						for (unsigned j = 0; j < cur_point_indexes.size(); j++)
							cout << cur_point_indexes[j] << " ";
						cout << endl;
					}*/
					cur_point = fromPointIndexesToPoint(cur_point_indexes, depths);
					if (find(checked_points.begin(), checked_points.end(), cur_point) != checked_points.end()) {
						skipped_points++;
						continue;
					}
					double d_val = fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
					checked_points.push_back(cur_point);
					isContinueDimension = false;
					if (d_val < old_record_residual) { // new record was found
						local_record_point.residual = d_val;
						local_record_point_indexes = cur_point_indexes;
						isLocalMin = false;
						isContinueDimension = true;
					}
					if ((!isContinueDimension) && (!isDecreaseTurn)) {
						isDecreaseTurn = true; // try to decrease current value
						//if (verbosity > 0)
						//	cout << "isDecreaseTurn " << isDecreaseTurn << endl;
						isContinueDimension = true;
						isTriggerIncToDec = true;
						continue;
					}
				} while (isContinueDimension);
			}
		} while (!isLocalMin);

		//cout << endl << "*** local minimum in hill climbing" << endl;
		//cout << "local record of residual " << record_point.residual << endl;
		//cout << "-----" << endl;
		//if (verbosity > 0)
		//	cout << "new random cur_point_indexes : " << endl;

		for(;;) {
			// prmutate current global minimum point to obtain a new start point
			for (unsigned i = 0; i < search_space.size(); i++) {
				if (search_space[i].size() == 1) {
					cur_point_indexes[i] = 0;
					if (verbosity > 0)
						cout << cur_point_indexes[i] << " ";
					continue;
				}
				unsigned rand_numb = rand();
				if (rand_numb % 3 == 0)
					cur_point_indexes[i] = local_record_point_indexes[i];
				else
					cur_point_indexes[i] = (rand_numb % search_space[i].size());
			}
			cur_point = fromPointIndexesToPoint(cur_point_indexes, depths);
			// if a new point doesn't exist in the check-list, then use this point as a new start
			if (find(checked_points.begin(), checked_points.end(), cur_point) == checked_points.end())
				break;
		}
		/*if (verbosity > 0) {
			cout << "new random point" << endl;
			for (unsigned j = 0; j < cur_point_indexes.size(); j++)
				cout << cur_point_indexes[j] << " ";
			cout << endl;
		}*/

		fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
		checked_points.push_back(cur_point);
		// new start point
		local_record_point = cur_point;
		local_record_point_indexes = cur_point_indexes;

		cout << "checked_points size " << checked_points.size() << endl;
		cout << "skipped_points " << skipped_points << endl;
		cout << "---" << endl;
	}

	return local_record_point;
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

double CAMBALA_sequential::getRecordResidual()
{
	return record_point.residual;
}

int CAMBALA_sequential::readScenario(string scenarioFileName)
{
// read constant and variable values from a scenario file
	if ( (!rank) && (verbosity > 0) )
		cout << "scenario file " << scenarioFileName << endl;
	ifstream scenarioFile(scenarioFileName.c_str());

	if (!scenarioFile.is_open()) {
		cerr << "scenarioFile with the name " << scenarioFileName << " wasn't openend" << endl;
		return -1;
	}

	if ((!rank) && (verbosity > 0))
		cout << "scenario file was opened" << endl;

	string str, word, tmp_word;
	stringstream sstream;
	unsigned cw_index = 0, d_index = 0;
	while (getline(scenarioFile, str)) {
		if ((str == "") || (str[0] == '%'))
			continue;
		sstream << str;
		sstream >> word;
		if (word.find("dtimes_file") != string::npos)
			sstream >> dtimesFileName;
		else if (word.find("spmag_file") != string::npos)
			sstream >> spmagFileName;
		else if (word == "function_type")
			sstream >> object_function_type;
		else if (word.find("launch_type") != string::npos)
			sstream >> launch_type;
		else if (word == "H")
			sstream >> H;
		else if (word == "h") {
			sstream >> h;
			nh = 1;
		}
		else if ((word.size() >= 3) && (word[0] == 'c') && (word[1] == 'w') && (isdigit(word[2]))) {
			word = word.substr(2, word.size() - 2);
			istringstream(word) >> cw_index;
			if (cw_vec_vec.size() < cw_index + 1)
				cw_init_vec_vec.resize(cw_index + 1);
			sstream >> word;
			fillArrayStep(word, cw_init_vec_vec[cw_index]);
		}
		else if ((word.size() == 2) && (word[0] == 'd') && (isdigit(word[1]))) {
			word = word.substr(1, word.size() - 1); // read depths
			istringstream(word) >> d_index;
			d_index--;

			if (d1_vec.size() < d_index + 1) {
				d1_vec.resize(d_index + 1);
				d2_vec.resize(d_index + 1);
				d_step_vec.resize(d_index + 1);	
			}
				sstream >> word;
			// here we need step, not count, so setParameterWithCount() is not required^M
			getThreeValuesFromStr(word, d1_vec[d_index], d_step_vec[d_index], d2_vec[d_index]);
		}
		else if (word == "R") {
			sstream >> word;
			fillArrayStep(word, R_vec);
		}
		else if (word == "rhob") {
			sstream >> word;
			fillArrayStep(word, rhob_vec);
		}
		else if (word == "cb") {
			sstream >> word;
			fillArrayStep(word, cb_vec);
		}
		else if (word == "tau") {
			sstream >> word;
			fillArrayStep(word, tau_vec);
		}
		else if (word == "ppm")
			sstream >> ppm;
		else if (word == "ils_iterations")
			sstream >> init_iterated_local_search_runs;
		else if (word == "verbosity")
			sstream >> verbosity;
		sstream.str(""); sstream.clear();
	}
	
	if (ppm == 0) { // if ppm wasn't set directly
		if ((object_function_type == "uniform") || (object_function_type == "Wan_uniform") || (object_function_type == "weighted"))
			ppm = 2;
		else if ((object_function_type == "uniform2") || (object_function_type == "weighted2"))
			ppm = 1;
		else {
			cerr << "unknown object_function_type " << object_function_type << endl;
			exit(1);
		}
	}
	
	if (dtimesFileName == "") {
		cerr << "dtimesFileName == "" \n";
		exit(-1);
	}
	if (!cw_init_vec_vec.size()) {
		cerr << "cw_init_vec_vec is empty \n";
		exit(-1);
	}
	if (H <= 0) {
		cerr << "incorrect H " << H << endl;
		exit(-1);
	}

	input_params_sstream << "Input parameters :" << endl;
	input_params_sstream << "launch_type " << launch_type << endl;
	input_params_sstream << "object_function_type " << object_function_type << endl;
	input_params_sstream << "ppm " << ppm << endl;
	input_params_sstream << "init_iterated_local_search_runs " << init_iterated_local_search_runs << endl;
	input_params_sstream << "cw_init_vec_vec :" << endl;
	for (auto &x : cw_init_vec_vec)
		input_params_sstream << doubleVecToStr(x) << endl;
	input_params_sstream << "d1_vec :" << endl;
	input_params_sstream << doubleVecToStr(d1_vec) << endl;
	input_params_sstream << "d2_vec :" << endl;
	input_params_sstream << doubleVecToStr(d2_vec) << endl;
	input_params_sstream << "d_step_vec :" << endl;
	input_params_sstream << doubleVecToStr(d_step_vec) << endl;
	input_params_sstream << "R_vec :" << endl;
	input_params_sstream << doubleVecToStr(R_vec) << endl;
	input_params_sstream << "cb_vec :" << endl;
	input_params_sstream << doubleVecToStr(cb_vec) << endl;
	input_params_sstream << "rhob_vec :" << endl;
	input_params_sstream << doubleVecToStr(rhob_vec) << endl;
	input_params_sstream << "tau_vec :" << endl;
	input_params_sstream << doubleVecToStr(tau_vec) << endl;
	input_params_sstream << "dtimes_file " << dtimesFileName << endl;
	input_params_sstream << "spmag_file " << spmagFileName << endl;
	input_params_sstream << "launch_type " << launch_type << endl;

	if (!rank) {
		cout << input_params_sstream.str() << endl;
		writeOutputData(input_params_sstream);
	}

	if (((rank == 0) || (rank == 1)) && (verbosity > 0))
		cout << "readScenario() finished" << endl;

	return 0;
}

int CAMBALA_sequential::fillArrayStep(const string word, vector<double> &vec)
{
	double left_bound_val = 0, step = 0, right_bound_val = 0;
	getThreeValuesFromStr(word, left_bound_val, step, right_bound_val);
	vec.push_back(left_bound_val);
	if (left_bound_val == right_bound_val)
		return 0;
	double val = left_bound_val;
	for (;;) {
		val += step;
		if (val > right_bound_val)
			break;
		vec.push_back(val);
	}
	if (find(vec.begin(),vec.end(),right_bound_val) == vec.end())
		vec.push_back(right_bound_val);
	return 0;
}

int CAMBALA_sequential::readInputDataFromFiles()
{
	ifstream dtimesFile(dtimesFileName.c_str());
	if (!dtimesFile.is_open()) {
		cerr << "dtimesFile " << dtimesFileName << " wasn't opened \n";
		exit(-1);
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
		ifstream spmagFile(spmagFileName);
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
			exit(-1);
		}
	}

	if ( (!rank) && (verbosity > 0) )
		cout << "readInputDataFromFiles() finished " << endl;
	return 0;
}

search_space_point CAMBALA_sequential::getNonRandomStartPoint( vector<double> depths )
{
	search_space_point point;
	point.cb = getMidVecValue(cb_vec);
	point.rhob = getMidVecValue(rhob_vec);
	point.R = getMidVecValue(R_vec);
	point.tau = getMidVecValue(tau_vec);
	point.cws.resize(cw_vec_vec.size());
	for (unsigned i=0; i<cw_vec_vec.size(); i++)
		point.cws[i] = getMidVecValue(cw_vec_vec[i]);
	point.depths = depths;

	return point;
}

double CAMBALA_sequential::directPointCalc( search_space_point point )
{
	if (verbosity > 1)
		cout << "start directPointCalc()" << endl;
	isTimeDelayPrinting = true;
	return fillDataComputeResidual(point);
}

void CAMBALA_sequential::solve()
{
	// clear out file
	ofstream ofile(output_filename, ios_base::out);
	ofile.close();

	vector<vector<double>> depths_vec = createDepthsArray();
	if (launch_type == "ils") {
		for (unsigned j = 0; j < depths_vec.size(); j++) {
			init(depths_vec[j]);
			//CAMBALA_seq.findGlobalMinBruteForce(depths_vec[i]);
			iterated_local_search_runs = init_iterated_local_search_runs;
			findLocalMinHillClimbing(depths_vec[j]);
			cout << "Processed " << j + 1 << " out of " << depths_vec.size() << " depths" << endl;
		}
	}
	else {
		init(depths_vec[0]);
		findGlobalMinBruteForce(depths_vec[0]);
	}
}
