/*****************************************************************************************
// CAMBALA: Coupled Acoustic Modes -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS),
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#include "../cambala/sequential.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

int main(int argc, char **argv)
{
	if (argc < 3) {
		cerr << "Usage: prog scenario smac-params, e.g. -d1 '2' -d2 '5'";
		exit(-1);
	}
	
	string scenarioFileName = argv[1];
	search_space_point point;
	point.residual = 1e+308;
	for (unsigned i=2; i<argc-1; i++) {
		string word = argv[i];
		if ((word.size() >= 4) && (word[0] == '-') && (word[1] == 'c') && (word[2] == 'w')) {
			word = word.substr(3, word.size() - 3);
			unsigned cw_index;
			istringstream(word) >> cw_index;
			if (point.cws.size() < cw_index + 1)
				point.cws.resize(cw_index + 1);
			string str_val = argv[i + 1];
			stringstream sstream;
			//sstream << str_val.substr(1, str_val.size() - 2);
			sstream << str_val;
			sstream >> point.cws[cw_index];
		}
		else if ((word.size() >= 3) && (word[0] == '-') && (word[1] == 'd')) {
			word = word.substr(2, word.size() - 2);
			unsigned d_index;
			istringstream(word) >> d_index;
			d_index--;
			if (point.depths.size() < d_index + 1)
				point.depths.resize(d_index + 1);
			string str_val = argv[i + 1];
			stringstream sstream;
			//sstream << str_val.substr(1, str_val.size() - 2);
			sstream << str_val;
			sstream >> point.depths[d_index];
		}
	}

	cout << "cws :\n";
	for (unsigned i = 0; i < point.cws.size(); i++)
		cout << point.cws[i] << " ";
	cout << endl;
	cout << "depths :\n";
	for (unsigned i = 0; i < point.depths.size(); i++)
		cout << point.depths[i] << " ";
	cout << endl;
	
	CAMBALA_sequential CAMBALA_seq;
	CAMBALA_seq.verbosity = 0;
	CAMBALA_seq.readScenario(scenarioFileName);
	CAMBALA_seq.readInputDataFromFiles();
	point.cb   = CAMBALA_seq.cb1;
	point.R    = CAMBALA_seq.R1;
	point.rhob = CAMBALA_seq.rhob1;
	point.tau  = CAMBALA_seq.tau1;
	vector<vector<double>> depths_vec = CAMBALA_seq.createDepthsArray();
	cout << "depths_vec size " << depths_vec.size() << endl;

	point.depths.push_back(CAMBALA_seq.h);
	point.depths.push_back(CAMBALA_seq.H);
	cout << "full depths :\n";
	for (unsigned i = 0; i < point.depths.size(); i++)
		cout << point.depths[i] << " ";
	cout << endl;
	
	double val = 0;
	vector<vector<double>>::iterator it = find(depths_vec.begin(), depths_vec.end(), point.depths);
	if (it == depths_vec.end())
		val = 1e100;
	else {
		CAMBALA_seq.init(point.depths);
		val = CAMBALA_seq.fillDataComputeResidual(point);
	}
	cout << "Result for SMAC: SUCCESS, 0, 0, " << val << ", 0";
	
	return 0;
}