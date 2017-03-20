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

int main()
{
	std::string scenarioFileName = "test_scenario.txt";
	std::vector<double> freqs = { 20, 30, 40, 50, 60, 70 };

	sspemdd_sequential sspemdd;
	sspemdd.verbosity = 0;
	sspemdd.readScenario(scenarioFileName);
	sspemdd.init();

	std::vector<double> depths = sspemdd.depths;
	std::vector<double> c1s = sspemdd.c1s;
	std::vector<double> c2s = sspemdd.c2s;
	std::vector<double> rhos = sspemdd.rhos;
	std::vector<unsigned> Ns_points = sspemdd.Ns_points;
	unsigned long long n_layers_w = sspemdd.n_layers_w;
	std::vector<double> cws = sspemdd.cw1_arr;
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
	std::vector<double> result_vec;
	std::ofstream ofile("output.txt");
	for (unsigned i = 0; i < freqs.size(); i++) {
		double omeg1 = 2 * LOCAL_M_PI*(freqs[i] + deltaf / 2);
		result_vec = sspemdd.compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, 1, ordRich);
		for (unsigned j = 0; j < result_vec.size(); j++)
			ofile << result_vec[j] << " ";
		ofile << std::endl;
	}
	ofile.close();
	
	std::cout << "test passed" << std::endl;
	
	return 0;
}