/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#ifdef _MPI
#include <mpi.h>
#endif

#include "cambala_sequential.h"
#include "cambala_parallel.h"

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
	CAMBALA_sequential CAMBALA_seq;
	int verbosity = 0;

#ifdef _DEBUG
	argc = 3;
	argv[1] = "./scenarios/803_bottom_r_full_uniform.txt";
	//argv[1] = "./boinc_CAMBALA_app/in";
	//argv[1] = "39_hydro_r_uniform260.txt";
	//argv[1] = "true_scenario_2.txt";
	argv[2] = "1";
	//argv[3] = "-test";
	verbosity = 2;
#endif

	if (argc >= 2)
		scenarioFileName = argv[1];
	else {
		cout << "Usage : scenarioFileName [verbosity] [-test]" << endl;
		exit(1);
	}
	if (argc >= 3)
		verbosity = atoi(argv[2]);
	bool isTestLaunch = false;
	if (argc >= 4) {
		string test_str= argv[3];
		if (test_str == "-test")
			isTestLaunch = true;
	}
	
#ifndef _MPI
	// sequential mode%
	cout << "verbosity " << verbosity << endl;

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

	CAMBALA_seq.verbosity = verbosity;
	// read scenario, modal_delays, mode_numbers and freqs, then determine the search space
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();

	if (isTestLaunch) {
		search_space_point point;
		// true values
		point.cb = 1700;
		point.R = 7000;
		point.rhob = 1.7;
		point.tau = 0;
		point.cws = { 1500, 1498, 1493, 1472, 1462 };
		point.depths = { 10, 20, 30, 40, 50, 300 };
		CAMBALA_seq.init(point.depths);
		CAMBALA_seq.object_function_type = "uniform2";
		CAMBALA_seq.directPointCalc( point );
		CAMBALA_seq.reportFinalResult();
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
	
	CAMBALA_seq.solve();
	CAMBALA_seq.reportFinalResult();

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

	if (!rank)
		cout << "verbosity " << verbosity << endl;

	CAMBALA_parallel CAMBALA_par;
	CAMBALA_par.rank = rank;
	CAMBALA_par.corecount = corecount;
	CAMBALA_par.verbosity = verbosity;

	CAMBALA_par.readScenario(scenarioFileName);
	CAMBALA_par.readInputDataFromFiles();
	
	double cur_time = MPI_Wtime();
	CAMBALA_par.MPI_main();
	cout << "MPI_main() total time " << MPI_Wtime() - cur_time << " s" << endl;
#endif
	
	return 0;
}