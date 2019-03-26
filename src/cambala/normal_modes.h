#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>

using namespace std;

//const double LOCAL_M_PI = 3.14159265358979323846;
//const complex<double> Iu(0.0, 1.0);

class NormalModes
{
public:
	NormalModes();

	// media parameters
	vector<double> M_c1s;
	vector<double> M_c2s;
	vector<double> M_rhos;
	vector<unsigned> M_Ns_points;
	vector<double> M_depths;

	// source/receiver parameters
	vector<double> zr;
	vector<double> R;
	vector<double> freq;

	//
	double iModesSubset;
	unsigned ordRich;
	double omeg; // sound frequency

	// input functions
	

	// output functions

	// computation functions
	vector<double> compute_wnumbers_extrap_lin_dz();
	vector<double> compute_wnumbers(vector<double> &c, vector<double> &rho,
		vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	
private:

	// parameters of wave numbers and mode functions
	vector<double> khs;
	vector<vector<double>> mfunctions_zr;
};

#endif
