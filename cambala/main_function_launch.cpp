/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#include "cambala_sequential.h"

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
	point.cb = 1715;
	point.R = 6993;
	point.rhob = 1.56;
	point.tau = 0;
	point.cws = { 1500, 1498, 1493, 1472, 1462 };
	point.depths = { 10, 20, 30, 40, 50, 300 };

	if (argc > 1)
		scenarioFileName = argv[1];
	else {
		//scenarioFileName = "./504_bottom_r_full_weighted.txt";
		scenarioFileName = "./504_1_bottom_r_full_weighted2.txt";
	}
	
	CAMBALA_sequential CAMBALA_seq;
	CAMBALA_seq.verbosity = 2;
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();
	CAMBALA_seq.init(point.depths);
	double direct_point_calc_time = cpu_time();
	CAMBALA_seq.directPointCalc( point );
	direct_point_calc_time = cpu_time() - direct_point_calc_time;
	cout << "direct_point_calc_time " << direct_point_calc_time << endl;
	CAMBALA_seq.reportFinalResult();

	cout << "test done" << endl;
	cout << "final time " << cpu_time() - start_time << endl;
	
	return 0;
}

double cpu_time() {
	return (double)clock() / CLOCKS_PER_SEC;
}