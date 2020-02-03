#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <iostream>
#include <fstream>
#include <complex>
#include <algorithm>
#include <vector>
#include <string>

using namespace std;

#ifndef M_PI
const double M_PI = 3.14159265358979323846;
#endif

#ifndef M_Iu
const complex<double> M_Iu(0.0, 1.0);
#endif

const int NCV_COEF = 6; // coef for ncv - parameter that controls the convergence speed of arpack

class NormalModes
{
public:
	NormalModes();

	// media parameters
	vector<double> M_c1s;
	vector<double> M_c2s;
	vector<double> M_rhos;
	vector<double> M_betas;
	vector<double> M_depths;
	vector<unsigned> M_Ns_points;

	// source/receiver parameters
	vector<double> zr;
	vector<double> Rs;
	double f; // frequency

	// parameters of wave numbers and mode functions
	double iModesSubset;
	unsigned ppm;
	unsigned ordRich;
	vector<unsigned> mode_numbers;
	vector<vector<double>> modal_group_velocities;
	vector<double> khs;
	vector<vector<double>> mfunctions_zr;
	vector<double> mattenuation;
	
	// input/output
	string wnumbers_out_file_name = "kj_wedge_att.txt";
	string mfunctions_out_file_name = "phizr_wedge.txt";
	string modal_group_velocities_out_file_name = "vgr.txt";
	void read_scenario(const string scenarioFileName);
	void write_result(const string resultFileName);
	void write_wnumbers();
	void write_mfunctions_zr();
	void write_modal_group_velocities();
	void print_wnumbers();
	void print_mfunctions_zr();
	void print_modal_group_velocities();
	
	// main computation functions
	vector<double> compute_wnumbers(const double omeg, vector<double> &c, vector<double> &rho,
		vector<unsigned> &interface_idcs, vector<double> &meshsizes);
	vector<double> compute_wnumbers_extrap2(const double omeg);
	vector<double> compute_wnumbers_extrap_lin_dz(const double omeg);
	void compute_wmode1(const double omeg, const double kh, vector<double> &phi, vector<double> &dphi);
	double compute_wmode_vg(const double omeg, const double kh, vector<double> &phi);
	void compute_wmode(const double omeg, const double kh, vector<double> &phi, vector<double> &dphi);
	vector<complex<double>> compute_cpl_pressure(const double f, vector<double> &Rr);
	void compute_mfunctions_zr(double omeg = -1);
	void compute_all_mfunctions(const double omeg);
	void compute_khs(double omeg = -1);
	void compute_mattenuation(double omeg = -1);
	void compute_modal_group_velocities_fixed_freq();
	int compute_modal_grop_velocities(const double deltaf, vector<double> freqs);
	int compute_modal_grop_velocities2(const double deltaf, vector<double> freqs);
	int compute_wnumbers_bb(const double deltaf, const unsigned flOnlyTrapped, vector<double> freqs);
	void compute_mfunctions_airy(double omeg = -1, double eps = 1e-10);
	
	// inversion functions
	double compute_modal_delays_residual_uniform(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_uniform2(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan1(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_LWan(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);
	double compute_modal_delays_residual_weighted(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers,
		vector<vector<double>> weight_coeffs);
	double compute_modal_delays_residual_weighted2(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers, 
		vector<vector<double>> weight_coeffs);
	double compute_modal_delays_residual_LWan_weighted(const double R, const double tau, vector<double> freqs,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers,
		vector<vector<double>> weight_coeffs);

	// meta data and functions
	vector<vector<double>> all_depths;
	void compute_for_all_depths();
	
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
	vector<vector<double>> parseTwoDimVector(stringstream &sstream);
	double integrate(const vector<double>&, const double&, unsigned, unsigned);
	//
	string eigen_type;                // "alglib" or "arpack"
	int verbosity;                    // 0 - silent, 1- short, 2 - full
	int arpack_required_eigen_values; // number of eigen values required by arpack
};

#endif
