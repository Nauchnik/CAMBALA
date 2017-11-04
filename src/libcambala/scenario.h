#ifndef SCENARIO_H_
#define SCENARIO_H_
#include <vector>
#include <string>
#include "types.h"

using std::string;
using std::vector;

class Scenario
{
public:
	double H_;
	Dim hDim_;
	string dtimesFileName_;
	string spmagFileName_;
	SearchSpaceDims ssd_;
	vector <Dim> depthsDim_;
	size_t iterated_local_search_runs_ = 10;
	string object_function_type_ = "uniform";
	size_t ppm_ = 0;
	vector <double> freqs_; // loaded from 1st column of dtimesFileName_
	vector <vector<double>>modal_delays_; // loaded from dtimesFilename_
	vector <vector<double>> spmag_; // loaded from smpagFilename_
	int readFile(string scenarioFileName);
	void readInputDataFromFiles();
	void print();
	Scenario(string scenarioFileName);
};

#endif
