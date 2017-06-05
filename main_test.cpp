/*****************************************************************************************
// SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#include "sspemdd_sequential.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

int main()
{
	string scenarioFileName = "test_scenario.txt";
	vector<double> freqs = { 20, 30, 40, 50, 60, 70 };
	vector<vector<double>> depths_vec;

	sspemdd_sequential sspemdd;
	sspemdd.verbosity = 0;
	sspemdd.readScenario(scenarioFileName);
	sspemdd.createDepthsArray(depths_vec);
	sspemdd.init(depths_vec[0]);

	vector<double> depths = depths_vec[0];
	vector<double> c1s = sspemdd.c1s;
	vector<double> c2s = sspemdd.c2s;
	vector<double> rhos = sspemdd.rhos;
	vector<unsigned> Ns_points = sspemdd.Ns_points;
	unsigned long long n_layers_w = sspemdd.n_layers_w;
	vector<double> cws = sspemdd.cw1_arr;
	double cb = sspemdd.cb1;
	double rhob = sspemdd.rhob1;

	for (unsigned jj = 0; jj < n_layers_w - 1; jj++) {
		c1s.at(jj) = cws.at(jj);
		c2s.at(jj) = cws.at(jj + 1);
		rhos.at(jj) = 1;
	}
	c1s.at(n_layers_w - 1) = cws.at(n_layers_w - 1);
	c2s.at(n_layers_w - 1) = cws.at(n_layers_w - 1);
	rhos.at(n_layers_w - 1) = 1;
	c1s.at(n_layers_w) = cb;
	c2s.at(n_layers_w) = cb;
	rhos.at(n_layers_w) = rhob;
	double deltaf = 0.05; // ???
	unsigned int ordRich = 1; // ???
	vector<double> result_vec;
	ofstream ofile("output.txt");
	for (unsigned i = 0; i < freqs.size(); i++) {
		double omeg1 = 2 * LOCAL_M_PI*(freqs[i] + deltaf / 2);
		result_vec = sspemdd.compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, 1, ordRich);
		for (unsigned j = 0; j < result_vec.size(); j++)
			ofile << result_vec[j] << " ";
		ofile << endl;
	}
	ofile.close();
	
	cout << "test passed" << endl;
	
	return 0;
}