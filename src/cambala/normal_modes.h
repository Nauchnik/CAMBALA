#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>

const double M_PI = 3.14159265358979323846;

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
	vector<double> freqs;

	//
	double iModesSubset;
	unsigned ordRich;
	double omeg; // sound frequency

	// input functions
	

	// output functions
	

	// main computation functions
	vector<double> compute_wnumbers_extrap_lin_dz(double &omeg);
	vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho,
								    vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	double compute_modal_delays_residual_uniform(double R, double tau, vector<vector<double>> &experimental_delays, 
		                            vector<unsigned> &experimental_mode_numbers);
	
private:
	// parameters of wave numbers and mode functions
	vector<double> khs;
	vector<vector<double>> mfunctions_zr;

	// additional computation functions
	double RK4(double omeg2, double kh2, double deltah, double c1, double c2, unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	double Layer_an_exp(double omeg2, double kh2, double deltah, double c, unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	void load_layers_data(string LayersFName);
	void load_profile_deep_water(string ProfileFName, const unsigned ppm);
	int compute_modal_grop_velocities(double deltaf, vector<vector<double>> &modal_group_velocities,
		vector<unsigned> &mode_numbers);
};

#endif
