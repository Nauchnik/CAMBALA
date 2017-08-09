/*****************************************************************************************
// SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data -- Copyright(c) 2015-2016
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#ifdef _MPI
#include <mpi.h>
#endif

#include "sspemdd_sequential.h"
#include "sspemdd_parallel.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

int main(int argc, char **argv)
{
	unsigned ncpl = 0; // search mesh within each water layer

	std::string scenarioFileName = "";
	sspemdd_sequential sspemdd_seq;
	unsigned iterated_local_search_runs = 10;
	int verbosity = 0;

#ifdef _DEBUG
	argc = 2;
	//argv[1] = "test_hydro_r_uniform260_3layers.txt";
	argv[1] = "test_hydro_r_uniform260_4layers_ssp.txt";
	//argv[1] = "39_hydro_r_uniform260.txt";
	//argv[1] = "true_scenario_2.txt";
	argv[2] = "1"; // iterated_local_search_runs
	verbosity = 2;
#endif

	if (argc >= 2)
		scenarioFileName = argv[1];
	else {
		cout << "Usage : scenarioFileName [iterated_local_search_runs] [verbosity] [-test]" << endl;
		exit(1);
	}
	if (argc >= 3)
		iterated_local_search_runs = atoi(argv[2]);
	if (argc >= 4)
		verbosity = atoi(argv[3]);
	bool isTestLaunch = false;
	if (argc >= 5) {
		string test_str= argv[4];
		if (test_str == "-test")
			isTestLaunch = true;
	}
	
#ifndef _MPI
	// sequential mode
	cout << "iterated_local_search_runs " << iterated_local_search_runs << endl;
	cout << "verbosity " << verbosity << endl;

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

	sspemdd_seq.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_seq.verbosity = verbosity;
	// read scenario, modal_delays, mode_numbers and freqs, then determine the search space
	sspemdd_seq.readScenario(scenarioFileName);
	sspemdd_seq.readInputDataFromFiles();
	vector<vector<double>> depths_vec;
	sspemdd_seq.createDepthsArray(depths_vec);

	if (isTestLaunch) {
		search_space_point point;
		// true values
		point.cb = 1700;
		point.R = 7000;
		point.rhob = 1.7;
		point.tau = 0;
		point.cws = { 1500, 1498, 1493, 1472, 1462 };
		point.depths = { 10, 20, 30, 40, 50 };
		sspemdd_seq.init(point.depths);
		sspemdd_seq.object_function_type = "uniform2";
		sspemdd_seq.directPointCalc( point );
		sspemdd_seq.reportFinalResult();
		cout << "true value test" << endl;
		return 0;
	}
	
	/*ofstream depths_file("depths.txt");
	depths_file << depths_vec.size() << " depths for the scenario " << scenarioFileName << endl;
	for (auto &x : depths_vec) {
		for (auto &y : x)
			depths_file << y << " ";
		depths_file << endl;
	}*/

	/*cout << "Unit test" << endl;
	vector<double> test_depths = { 9, 16, 36, 50, 300 };
	sspemdd_seq.init(test_depths);
	search_space_point test_point;
	test_point.cb = 1700;
	test_point.rhob = 1.7;
	test_point.tau = 0;
	test_point.R = 7000;
	test_point.cws = { 1459, 1487, 1488, 1467 };
	test_point.depths = test_depths;
	double test_residual = sspemdd_seq.fillDataComputeResidual(test_point);
	cout << "*** test_residual" << test_residual << endl;
	cout << "*** canonical residual " << 0.0113626 << endl;
	cout << endl;*/
	
	for (unsigned i = 0; i < depths_vec.size(); i++) {
		sspemdd_seq.init(depths_vec[i]);
		//sspemdd_seq.findGlobalMinBruteForce(depths_vec[i]);
		sspemdd_seq.findLocalMinHillClimbing(depths_vec[i]);
		cout << "Processed " << i + 1 << " out of " << depths_vec.size() << " depths" << endl;
	}

	sspemdd_seq.reportFinalResult();

	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	cout << "main() total time " << time_span.count() << " s" << endl;

#else
	int rank = 0;
	int corecount = 1;

	// parallel mode
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (!rank) {
		cout << "iterated_local_search_runs " << iterated_local_search_runs << endl;
		cout << "verbosity " << verbosity << endl;
	}

	sspemdd_parallel sspemdd_par;
	sspemdd_par.rank = rank;
	sspemdd_par.corecount = corecount;

	sspemdd_par.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_par.verbosity = verbosity;

	sspemdd_par.readScenario(scenarioFileName);
	sspemdd_par.readInputDataFromFiles();
	
	double cur_time = MPI_Wtime();
	sspemdd_par.MPI_main();
	cout << "MPI_main() total time " << MPI_Wtime() - cur_time << " s" << endl;
#endif
	
	return 0;
}