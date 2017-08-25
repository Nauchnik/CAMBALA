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

int main(int argc, char **argv)
{
	string scenarioFileName; 
	search_space_point point;
	// input data
	point.cb = 1715;
	point.R = 6993;
	point.rhob = 1.56;
	point.tau = 0;
	point.cws = { 1500, 1498, 1493, 1472, 1462 };
	point.depths = { 10, 20, 30, 40, 50, 300 };

	//scenarioFileName = "./504_bottom_r_full_weighted.txt";
	scenarioFileName = "./504_1_bottom_r_full_weighted2.txt";
	//
	CAMBALA_sequential CAMBALA_seq;
	CAMBALA_seq.verbosity = 2;
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();

	CAMBALA_seq.init(point.depths);
	CAMBALA_seq.directPointCalc( point );
	CAMBALA_seq.reportFinalResult();

	cout << "test done" << endl;
	
	return 0;
}