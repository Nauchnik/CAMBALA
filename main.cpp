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

	std::string dtimesFileName = "50_ac_modes_R7km_dtimes.txt";
	std::string spmagFileName = "50_ac_modes_R7km_spmag.txt";
	int launchType = 11; // 0 - deafult; 1 - fixed cw1=1490; 2 - cw1>cw2>...>cwn; 3 - both 1 and 2
	sspemdd_sequential sspemdd_seq;
	unsigned iterated_local_search_runs = 10;
	int verbosity = 0;

#ifdef _DEBUG
	argc = 5;
	argv[1] = "50_ac_modes_R7km_dtimes.txt";
	//argv[1] = "dtimes_synth_thcline_hf.txt";
	argv[2] = "11"; // launchType
	argv[3] = "21"; // ncpl
	argv[4] = "10"; // iterated_local_search_runs
#endif

	std::cout << "spmagFileName " << spmagFileName << std::endl;
	
	if (argc >= 2) {
		dtimesFileName = argv[1];
		std::cout << "dtimesFileName " << dtimesFileName << std::endl;
	}
	else {
		std::cout << "Usage : thcline_file [launchType] [ncpl] [iterated_local_search_runs] [verbosity]" << std::endl;
		return 1;
	}
	
	if (argc >= 3) {
		launchType = atoi(argv[2]);
		std::cout << "launchType " << launchType << std::endl;
	}

	if (argc >= 4) {
		ncpl = atoi(argv[3]);
		std::cout << "ncpl " << ncpl << std::endl;
	}

	if (argc >= 5) {
		iterated_local_search_runs = atoi(argv[4]);
		std::cout << "iterated_local_search_runs " << iterated_local_search_runs << std::endl;
	}

	if (argc >= 6) {
		verbosity = atoi(argv[5]);
		std::cout << "verbosity " << verbosity << std::endl;
	}

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;
	
#ifndef _MPI
	// sequential mode

	sspemdd_seq.ncpl = ncpl;
	sspemdd_seq.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_seq.verbosity = verbosity;
	// read modal_delays, mode_numbers and freqs, then determine the search space
	sspemdd_seq.readInputDataFromFiles(dtimesFileName, spmagFileName, launchType);
	sspemdd_seq.init();
	
	//sspemdd_seq.findGlobalMinBruteForce();
	sspemdd_seq.findLocalMinHillClimbing();
	sspemdd_seq.report_final_result();

	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "main() total time " << time_span.count() << " s" << std::endl;

#else
	int rank = 0;
	int corecount = 1;
	
	// parallel mode
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sspemdd_parallel sspemdd_par;
	sspemdd_par.rank = rank;
	sspemdd_par.corecount = corecount;

	// read modal_delays, mode_numbers and freqs, then determine the search space
	sspemdd_par.readDtimesFromFile(dtimesFileName, launchType);
	sspemdd_par.readSpmagFromFile(spmagFileName);
	
	sspemdd_par.ncpl = ncpl;
	sspemdd_par.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_par.verbosity = verbosity;
	
	sspemdd_par.init();
	double cur_time = MPI_Wtime();
	sspemdd_par.MPI_main();
	std::cout << "MPI_main() total time " << MPI_Wtime() - cur_time << " s" << std::endl;
#endif
	
	return 0;
}