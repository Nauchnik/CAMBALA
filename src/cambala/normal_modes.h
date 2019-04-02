#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>

using namespace std;

const double M_PI = 3.14159265358979323846;
const complex<double> M_Iu(0.0, 1.0);

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
	
	// input/output functions
	void load_layers_data(const string LayersFName);
	void load_profile_deep_water(const string ProfileFName, const unsigned ppm);
	void printDelayTime(const double R);
	
	// main computation functions
	double compute_modal_delays_residual_uniform(const double R, const double tau, vector<vector<double>> &experimental_delays,
		                            vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan1(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan_weighted(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_uniform2(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_weighted(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_weighted2(const double R, const double tau, vector<vector<double>> &experimental_delays, 
		vector<unsigned> &experimental_mode_numbers
	);

	vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho,
		vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	vector<double> compute_wnumbers_extrap2(double &omeg);
	vector<double> compute_wnumbers_extrap_lin_dz(double &omeg);
	void compute_wmode1(double &omeg, vector<unsigned> &Ns_points_m, const double kh, vector<double> &phi, vector<double> &dphi);
	double compute_wmode_vg(double &omeg, vector<unsigned> &Ns_points_m, const double kh, vector<double> &phi);
	void compute_wmode(double &omeg, const double kh, vector<double> &phi, vector<double> &dphi);
	vector<complex<double>> compute_cpl_pressure(const double f, vector<double> &Rr);
	void compute_mfunctions_zr(double &omeg, vector<vector<double>> &mfunctions_zr);
	void compute_all_mfunctions(double &omeg);
	int compute_modal_grop_velocities(const double deltaf);
	int compute_modal_grop_velocities2(const double deltaf);
	int compute_wnumbers_bb(const double deltaf, const unsigned flOnlyTrapped);

private:
	// additional computation functions
	double RK4(const double omeg2, const double kh2, const double deltah, const double c1, const double c2, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	double Layer_an_exp(const double omeg2, const double kh2, const double deltah, const double c, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
};

#endif
