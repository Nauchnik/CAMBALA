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
	vector<double> R_vec;
	vector<double> freqs;


	// parameters of wave numbers and mode functions
	vector<double> khs;
	vector<vector<double>> mfunctions_zr;
	vector<vector<double>> modal_group_velocities;
	vector<unsigned> mode_numbers;

	//
	double iModesSubset;
	unsigned ordRich;
	vector<vector<double>> weight_coeffs;
	
	// input functions
	

	// output functions
	
	
	// main computation functions

	vector<double> compute_wnumbers_extrap_lin_dz(double &omeg);
	vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho,
								    vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	vector<double> compute_wnumbers_extrap2(double &omeg);
	double compute_modal_delays_residual_uniform(const double R, const double tau, vector<vector<double>> &experimental_delays,
		                            vector<unsigned> &experimental_mode_numbers);
	int compute_modal_grop_velocities(const double deltaf);
	int compute_modal_grop_velocities2(const double deltaf);
	double compute_modal_delays_residual_LWan1(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan_weighted(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_uniform2(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);

private:
	// additional computation functions
	double RK4(const double omeg2, const double kh2, const double deltah, const double c1, const double c2, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	double Layer_an_exp(const double omeg2, const double kh2, const double deltah, const double c, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	void load_layers_data(const string LayersFName);
	void load_profile_deep_water(const string ProfileFName, const unsigned ppm);
};

#endif
