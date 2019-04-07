#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>

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
	vector<double> M_depths;
	vector<unsigned> M_Ns_points;

	// source/receiver parameters
	vector<double> zr;
	vector<double> Rs;
	vector<double> freqs;

	// parameters of wave numbers and mode functions
	double iModesSubset;
	unsigned ppm;
	unsigned ordRich;
	vector<unsigned> mode_numbers;
	vector<vector<double>> modal_group_velocities;
	vector<double> khs;
	vector<vector<double>> mfunctions_zr;
	
	// input/output functions
	void read_data(const string scenarioFileName);
	void write_result(const string resultFileName);
	
	// main computation functions
	vector<double> compute_wnumbers(const double omeg, vector<double> &c, vector<double> &rho,
		vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	vector<double> compute_wnumbers_extrap2(const double omeg);
	vector<double> compute_wnumbers_extrap_lin_dz(const double omeg);
	void compute_wmode1(const double omeg, const double kh, vector<double> &phi, vector<double> &dphi);
	double compute_wmode_vg(const double omeg, const double kh, vector<double> &phi);
	void compute_wmode(const double omeg, const double kh, vector<double> &phi, vector<double> &dphi);
	vector<complex<double>> compute_cpl_pressure(const double f, vector<double> &Rr);
	void compute_mfunctions_zr(const double omeg);
	void compute_all_mfunctions(const double omeg);
	int compute_modal_grop_velocities(const double deltaf);
	int compute_modal_grop_velocities2(const double deltaf);
	int compute_wnumbers_bb(const double deltaf, const unsigned flOnlyTrapped);
	
	// inversion functions
	double compute_modal_delays_residual_uniform(const double R, const double tau, vector<vector<double>> &experimental_delays,
		                            vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_uniform2(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan1(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_weighted(const double R, const double tau, 
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers,
		vector<vector<double>> weight_coeffs);
	double compute_modal_delays_residual_weighted2(const double R, const double tau, vector<vector<double>> &experimental_delays, 
		vector<unsigned> &experimental_mode_numbers, vector<vector<double>> weight_coeffs);
	double compute_modal_delays_residual_LWan_weighted(const double R, const double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers,
		vector<vector<double>> weight_coeffs);
	
private:
	// additional computation functions
	double RK4(const double omeg2, const double kh2, const double deltah, const double c1, const double c2, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	double Layer_an_exp(const double omeg2, const double kh2, const double deltah, const double c, const unsigned Np,
		vector<double> &phi0, vector<double> &dphi0);
	vector<double> fillArrayStep(const string word);
	void getThreeValuesFromStr(string str, double &val1, double &val2, double &val3);
	vector<double> parseArrayBrackets(string word);
	vector<double> parseVector(stringstream &sstream);
};

#endif
