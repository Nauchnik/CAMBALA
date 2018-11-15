#ifndef SEQUENTIAL_H
#define SEQUENTIAL_H

#include <vector>
#include <complex>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <math.h>
#include "linalg.h"
#include "point.h"

using namespace std;
using namespace CAMBALA_point;

const unsigned MAX_DEPTHS_VECTORS = 1000000;

class CAMBALA_sequential
{
public:
	CAMBALA_sequential();
	double H;
	unsigned long long nh;
	double h;
	vector<double> cb_vec;
	vector<double> R_vec;
	vector<double> tau_vec;
	vector<double> rhob_vec;
	vector<vector<double>> cw_vec_vec;
	vector<vector<double>> cw_init_vec_vec;
	vector<double> d1_vec;
	vector<double> d2_vec;
	vector<double> d_step_vec;
	string object_function_type;
	unsigned ppm;
	unsigned long long n_layers_w;
	unsigned long long init_iterated_local_search_runs;
	unsigned long long iterated_local_search_runs;
	vector<unsigned> mode_numbers;
	vector<vector<double>> modal_delays;
	vector<double> freqs;
	vector<double> c1s;
	vector<double> c2s;
	vector<double> rhos;
	vector<unsigned> Ns_points;
	vector<vector<double>> weight_coeffs;
	int verbosity;
	string dtimesFileName;
	string spmagFileName;
	int rank;
	bool isTimeDelayPrinting;
	string launch_type; // bruteforce | ils
	stringstream input_params_sstream;
	string output_filename;
	string depths_filename;
	bool is_mpi;
	
	// Oleg's functions
	vector<vector<double>> search_space; // values of variables which form a search space
	int readScenario(string scenarioFileName);
	int fillArrayStep(const string word, vector<double> &vec);
	int readInputDataFromFiles();
	void writeOutputData(stringstream &sstream);
	int init(vector<double> depths);
	void solve();
	vector<vector<double>> createDepthsArray();
	void loadValuesToSearchSpaceVariables();
	double getRecordResidual();
	double fillDataComputeResidual( search_space_point &point );
	void updateRecordPoint(const search_space_point point);
	search_space_point getNonRandomStartPoint(vector<double> depths);
	vector<search_space_point> getSearchSpacePointsVec(vector<double> depths);
	search_space_point findLocalMinHillClimbing(vector<double> depths);
	search_space_point fromPointIndexesToPoint(vector<unsigned> cur_point_indexes, vector<double> depths);
	vector<unsigned> fromPointToPointIndexes(search_space_point point);
	void findGlobalMinBruteForce(vector<double> depths);
	void reportFinalResult();
	void reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a);
	double directPointCalc( search_space_point point );

protected:
	unsigned long long N_total;
	search_space_point record_point;
	chrono::high_resolution_clock::time_point start_chrono_time;
	// hill climbing
	vector<search_space_point> checked_points;
};

#endif
