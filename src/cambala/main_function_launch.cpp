/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS),
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#include "sequential.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

double cpu_time();

int main(int argc, char **argv)
{
	double start_time = cpu_time();

	string scenarioFileName;
	search_space_point point;
	// input data
	point.cb = 1730;
	point.R = 6080;
	point.rhob = 2.4;
	point.tau = 0.063;
	point.cws = { 1452, 1450, 1448 };
	point.depths = { 5, 20, 37, 500 };

	if (argc == 1) {
        cerr << "Usage: program scenario_file \n";
        exit(1);
	}

    scenarioFileName = argv[1];

	CAMBALA_sequential CAMBALA_seq;
	CAMBALA_seq.verbosity = 2;
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();
	CAMBALA_seq.init(point.depths);
	double direct_point_calc_time = cpu_time();
	CAMBALA_seq.directPointCalc( point );
	direct_point_calc_time = cpu_time() - direct_point_calc_time;
	cout << "direct_point_calc_time " << direct_point_calc_time << endl;
	CAMBALA_seq.reportCurrentResult(true);

	cout << "test done" << endl;
	cout << "final time " << cpu_time() - start_time << endl;

	return 0;
}

double cpu_time() {
	return (double)clock() / CLOCKS_PER_SEC;
}
