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
		cout << "Usage : scenarioFileName [iterated_local_search_runs] [verbosity]" << endl;
		exit(1);
	}
	if (argc >= 3)
		iterated_local_search_runs = atoi(argv[2]);
	if (argc >= 4)
		verbosity = atoi(argv[3]);

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

	/*ofstream depths_file("depths.txt");
	depths_file << depths_vec.size() << " depths for the scenario " << scenarioFileName << endl;
	for (auto &x : depths_vec) {
		for (auto &y : x)
			depths_file << y << " ";
		depths_file << endl;
	}*/

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