#ifndef SCENARIO_H_
#define SCENARIO_H_
#include <vector>
#include <string>
#include "types.h"

using namespace std;

class Scenario
{
public:
	double H_;
	double h_;
	string dtimesFileName_;
	string spmagFileName_;
	SearchSpaceDims ssd_;
	vector <Dim> depthsDim_;
	size_t iterated_local_search_runs_;
	string object_function_type_;
	size_t ppm_;
	vector <double> freqs_; // loaded from 1st column of dtimesFileName_
	vector <vector<double>>modal_delays_; // loaded from dtimesFilename_
	vector <vector<double>> spmag_; // loaded from smpagFilename_
	int readFile(string scenarioFileName);
	void readInputDataFromFiles();
	void print();
};

#endif
