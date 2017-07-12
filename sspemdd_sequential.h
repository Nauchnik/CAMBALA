#ifndef SSPEMDD_SEQUENTIAL_H
#define SSPEMDD_SEQUENTIAL_H

#include <vector>
#include <complex>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <math.h>
#include "linalg.h"

const double LOCAL_M_PI = 3.14159265358979323846;
const std::complex<double> Iu(0.0,1.0);
const double START_HUGE_VALUE = 1e100;

using namespace std;

struct search_space_point
{
	double R;
	double tau;
	double rhob;
	double cb;
	vector<double> cws;
	double residual;
	vector<double> depths;
	
	bool operator==(const search_space_point& a) const
	{
		return (R == a.R && tau == a.tau && rhob == a.rhob && cb == a.cb && cws == a.cws && depths == a.depths);
	}
};

struct reduced_search_space_attribute
{
	bool R;
	bool tau;
	bool rhob;
	bool cb;
	vector<bool> cws;
};

class sspemdd_sequential
{
public:
	sspemdd_sequential();
	double h;
	double H;
	unsigned long long ncb;
	unsigned long long nrhob;
	unsigned long long nR;
	unsigned long long ntau;  //tau_comment: added tau to the class declaration
	double cb1;
	double cb2;
	vector<double> cw1_arr;
	vector<double> cw2_arr;
	vector<double> cw1_init_arr;
	vector<double> cw2_init_arr;
	vector<double> d1_arr;
	vector<double> d2_arr;
	vector<double> d_step;
	vector<unsigned long long> ncpl_init_arr;
	vector<unsigned long long> ncpl_arr;
	std::string object_function_type;
	double R1;
	double R2;
	double tau1;    //tau_comment: added tau to the class declaration
	double tau2;    //tau_comment: added tau to the class declaration
	double rhob1;
	double rhob2;
	unsigned long long n_layers_w;
	unsigned long long iterated_local_search_runs;
	vector<unsigned> mode_numbers;
	vector<vector<double>> modal_delays;
	vector<double> freqs;
	vector<double> c1s;
	vector<double> c2s;
	vector<double> rhos;
	vector<unsigned> Ns_points;
	vector<vector<double>> weight_coeffs;
	int verbosity;
	std::string dtimesFileName;
	std::string spmagFileName;
	int rank;

	// functions by Oleg
	vector<vector<double>> search_space; // values of variables which form a search space
	int readScenario(std::string scenarioFileName);
	int readInputDataFromFiles();
	int init(vector<double> depths);
	int createDepthsArray(vector<vector<double>> &depths_vec);
	double getRecordResidual();
	double fillDataComputeResidual( search_space_point &point );
	vector<search_space_point> getSearchSpacePointsVec(vector<double> depths);
	void findGlobalMinBruteForce(vector<double> depths);
	void loadValuesToSearchSpaceVariables();
	void findLocalMinHillClimbing(vector<double> depths);
	void reportFinalResult();
	void getThreeValuesFromStr(std::string str, double &val1, double &val2, double &val3);
	void reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a);

	// functions by Pavel
	//tau_comment: added tau to the arguments of compute_modal_delays_residual_uniform()

    void load_layers_data(std::string LayersFName,
	vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,
	vector<unsigned> &Ns_points);


	void load_profile_deep_water(std::string ProfileFName,unsigned ppm,
	vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,
	vector<unsigned> &Ns_points);

	double compute_modal_delays_residual_uniform(vector<double> &freqs, vector<double> &depths,
		vector<double> &c1s, vector<double> &c2s, vector<double> &rhos,
		vector<unsigned> &Ns_points, double R, double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);

    double compute_modal_delays_residual_uniform2(vector<double> &freqs, vector<double> &depths,
		vector<double> &c1s, vector<double> &c2s, vector<double> &rhos,
		vector<unsigned> &Ns_points, double R, double tau,
		vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);

	double compute_modal_delays_residual_weighted(vector<double> &freqs,
		vector<double> &depths, vector<double> &c1s, vector<double> &c2s,
		vector<double> &rhos, vector<unsigned> &Ns_points, double R,
		double tau, vector<vector<double>> &experimental_delays,
		vector<vector<double>> &weight_coeffs, vector<unsigned> &experimental_mode_numbers
	);

	vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho,
		vector<unsigned> &interface_idcs, vector<double> &meshsizes,
		double iModesSubset);

	vector<double> compute_wnumbers_extrap(double &omeg, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos,
		vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

    vector<double> compute_wnumbers_extrap2(double &omeg, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos,
		vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

	vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos,
		vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

    void compute_wmode(double &omeg, vector<double> &depths, vector<double> &c1s,
        vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, double kh,	vector<double> &phi,
	vector<double> &dphi);

    void compute_wmode1(double &omeg, vector<double> &depths, vector<double> &c1s,
        vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, double kh,	vector<double> &phi,
	vector<double> &dphi);

	vector<std::complex<double>> compute_cpl_pressure(double f, vector<double> &depths, vector<double> &c1s,
	vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, vector<double> &Rr,
	vector<double> &zr, double iModesSubset, unsigned ordRich);

	void compute_all_mfunctions(double &omeg, vector<double> &depths, vector<double> &c1s,
	vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, vector<double> &khs);

	void compute_mfunctions_zr(double &omeg, vector<double> &depths, vector<double> &c1s, vector<double> &c2s,
	vector<double> &rhos, vector<unsigned> &Ns_points, vector<double> &khs, vector<double> &zr, vector<vector<double>> &mfunctions_zr);

	double compute_wmode_vg(double &omeg, vector<double> &depths, vector<double> &c1s, vector<double> &c2s,
    vector<double> &rhos, vector<unsigned> &Ns_points, double kh, vector<double> &phi);

	int compute_modal_grop_velocities(vector<double> &freqs, double deltaf, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points,
		double iModesSubset, unsigned &ordRich, vector<vector<double>> &modal_group_velocities,
		vector<unsigned> &mode_numbers);

    int compute_modal_grop_velocities2(vector<double> &freqs, double deltaf, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points,
		double iModesSubset, unsigned &ordRich, vector<vector<double>> &modal_group_velocities,
		vector<unsigned> &mode_numbers);

	int compute_wnumbers_bb(vector<double> &freqs, double deltaf, vector<double> &depths, vector<double> &c1s,
		vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points,
		unsigned flOnlyTrapped, unsigned &ordRich, vector<vector<double>> &modal_group_velocities,
		vector<unsigned> &mode_numbers);

protected:
	unsigned long long N_total;
	search_space_point record_point;
	std::chrono::high_resolution_clock::time_point start_chrono_time;
	// hill climbing
	search_space_point fromPointIndexesToPoint( vector<unsigned> cur_point_indexes, vector<double> depths);
	search_space_point fromDoubleVecToPoint(vector<double> double_vec);
	vector<search_space_point> checked_points;
};

#endif
