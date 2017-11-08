/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
// Vadim Bulavintsev (Delft University of Technology)
*****************************************************************************************/

#ifdef _MPI
#include <mpi.h>
#include "cambala_mpi.h"
#endif

#include "residual/cpu.h"
#include "cambala.h"


#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#define ELPP_FEATURE_PERFORMANCE_TRACKING
#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP


using std::cout;
using std::endl;

int main(int argc, char *argv[])
{
	TIMED_FUNC(timerObj);
	//unsigned ncpl = 0; // search mesh within each water layer

	std::string scenarioFileName = "";
	int verbosity = 0;

#ifdef _DEBUG
	argc = 3;
	argv[1] = "./scenarios_5/504_1_bottom_r_full_weighted2.txt";
	//argv[1] = "./boinc_CAMBALA_app/in";
	//argv[1] = "39_hydro_r_uniform260.txt";
	//argv[1] = "true_scenario_2.txt";
	argv[2] = "1";
	//argv[3] = "-test";
	verbosity = 2;
#endif

	if (argc >= 2)
	{
		scenarioFileName = argv[1];
	}
	else 
	{
		cout << "Usage : scenarioFileName [verbosity] [-test]" << endl;
		exit(1);
	}
	if (argc >= 3)
		verbosity = atoi(argv[2]);
	bool isTestLaunch = false;
	if (argc >= 4)
	{
		string test_str= argv[3];
		if (test_str == "-test")
			isTestLaunch = true;
	}
	
#ifndef _MPI
	// sequential mode%
	//cout << "verbosity " << verbosity << endl;


	// read scenario, modal_delays, mode_numbers, freqs and parameters of the search space
	Scenario scenario(scenarioFileName);
	CAMBALA cambala;
	

	ResCalc* cpu64 = new BisectResCalcCPU <double> ("cpu64");
	cambala.AddResidualCalculator("cpu64", cpu64);

	ResCalc* cpu32 = new BisectResCalcCPU <float> ("cpu32");
	cambala.AddResidualCalculator("cpu32", cpu32);

	cambala.calcs_["fast"] = cpu32;
	cambala.calcs_["precise"] = cpu64;
	//TIMED_SCOPE(timerBlkObj, "CambalaSolve");
	cambala.Solve(scenario);
	cambala.reportFinalResult();
	delete cpu64;
	delete cpu32;


	/*ofstream depths_file("depths.txt");
	depths_file << depths_vec.size() << " depths for the scenario " << scenarioFileName << endl;
	for (auto &x : depths_vec) {
		for (auto &y : x)
			depths_file << y << " ";
		depths_file << endl;
	}*/
	



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
