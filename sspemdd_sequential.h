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

struct search_space_point
{
	double R;
	double tau;
	double rhob;
	double cb;
	std::vector<double> cws;
	double residual;
	std::vector<double> depths;
	
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
	std::vector<bool> cws;
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
	std::vector<double> cw1_arr;
	std::vector<double> cw2_arr;
	std::vector<double> d1_arr;
	std::vector<double> d2_arr;
	std::vector<double> d_step;
	std::vector<unsigned long long> ncpl_arr;
	std::string object_function_type;
	double R1;
	double R2;
	double tau1;    //tau_comment: added tau to the class declaration
	double tau2;    //tau_comment: added tau to the class declaration
	double rhob1;
	double rhob2;
	unsigned long long n_layers_w;
	unsigned long long iterated_local_search_runs;
	std::vector<unsigned> mode_numbers;
	std::vector<std::vector<double>> modal_delays;
	std::vector<double> freqs;
	std::vector<double> c1s;
	std::vector<double> c2s;
	std::vector<double> rhos;
	std::vector<unsigned> Ns_points;
	std::vector<std::vector<double>> weight_coeffs;
	int verbosity;
	std::string dtimesFileName;
	std::string spmagFileName;
	int rank;

	// functions by Oleg
	std::vector<std::vector<double>> search_space; // values of variables which form a search space
	int readScenario(std::string scenarioFileName);
	int readInputDataFromFiles();
	int init(std::vector<double> depths);
	int createDepthsArray(std::vector<std::vector<double>> &depths_vec);
	double getRecordResidual();
	double fillDataComputeResidual(search_space_point &point );
	std::vector<search_space_point> getSearchSpacePointsVec();
	void findGlobalMinBruteForce(std::vector<double> depths);
	void loadValuesToSearchSpaceVariables();
	void findLocalMinHillClimbing(std::vector<double> depths);
	void reportFinalResult();
	void getThreeValuesFromStr(std::string str, double &val1, double &val2, double &val3);
	void reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a);
	void updateGlobalMin();

	// functions by Pavel
	//tau_comment: added tau to the arguments of compute_modal_delays_residual_uniform()

    void load_layers_data(std::string LayersFName,
	std::vector<double> &depths,std::vector<double> &c1s,std::vector<double> &c2s,std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points);


	void load_profile_deep_water(std::string ProfileFName,unsigned ppm,
	std::vector<double> &depths,std::vector<double> &c1s,std::vector<double> &c2s,std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points);

	double compute_modal_delays_residual_uniform(std::vector<double> &freqs, std::vector<double> &depths,
		std::vector<double> &c1s, std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, double R, double tau,
		std::vector<std::vector<double>> &experimental_delays, std::vector<unsigned> &experimental_mode_numbers);

    double compute_modal_delays_residual_uniform2(std::vector<double> &freqs, std::vector<double> &depths,
		std::vector<double> &c1s, std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, double R, double tau,
		std::vector<std::vector<double>> &experimental_delays, std::vector<unsigned> &experimental_mode_numbers);


	double compute_modal_delays_residual_weighted(std::vector<double> &freqs,
		std::vector<double> &depths, std::vector<double> &c1s, std::vector<double> &c2s,
		std::vector<double> &rhos, std::vector<unsigned> &Ns_points, double R,
		double tau, std::vector<std::vector<double>> &experimental_delays,
		std::vector<std::vector<double>> &weight_coeffs, std::vector<unsigned> &experimental_mode_numbers
	);

	std::vector<double> compute_wnumbers(double &omeg, std::vector<double> &c, std::vector<double> &rho,
		std::vector<unsigned> &interface_idcs, std::vector<double> &meshsizes,
		double iModesSubset);

	std::vector<double> compute_wnumbers_extrap(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

    std::vector<double> compute_wnumbers_extrap2(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

	std::vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, double iModesSubset, unsigned &ordRich);

    void compute_wmode(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
        std::vector<double> &c2s,std::vector<double> &rhos,std::vector<unsigned> &Ns_points, double kh,	std::vector<double> &phi,
	std::vector<double> &dphi);

    void compute_wmode1(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
        std::vector<double> &c2s,std::vector<double> &rhos,std::vector<unsigned> &Ns_points, double kh,	std::vector<double> &phi,
	std::vector<double> &dphi);


	std::vector<std::complex<double>> compute_cpl_pressure(double f, std::vector<double> &depths, std::vector<double> &c1s,
	std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points, std::vector<double> &Rr,
	std::vector<double> &zr, double iModesSubset, unsigned ordRich);

	void compute_all_mfunctions(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
	std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points, std::vector<double> &khs);

	void compute_mfunctions_zr(double &omeg, std::vector<double> &depths, std::vector<double> &c1s, std::vector<double> &c2s,
	std::vector<double> &rhos, std::vector<unsigned> &Ns_points, std::vector<double> &khs, std::vector<double> &zr, std::vector<std::vector<double>> &mfunctions_zr);

	double compute_wmode_vg(double &omeg, std::vector<double> &depths, std::vector<double> &c1s, std::vector<double> &c2s,
    std::vector<double> &rhos, std::vector<unsigned> &Ns_points, double kh, std::vector<double> &phi);

	int compute_modal_grop_velocities(std::vector<double> &freqs, double deltaf, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		double iModesSubset, unsigned &ordRich, std::vector<std::vector<double>> &modal_group_velocities,
		std::vector<unsigned> &mode_numbers);

    int compute_modal_grop_velocities2(std::vector<double> &freqs, double deltaf, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		double iModesSubset, unsigned &ordRich, std::vector<std::vector<double>> &modal_group_velocities,
		std::vector<unsigned> &mode_numbers);

	int compute_wnumbers_bb(std::vector<double> &freqs, double deltaf, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		unsigned flOnlyTrapped, unsigned &ordRich, std::vector<std::vector<double>> &modal_group_velocities,
		std::vector<unsigned> &mode_numbers);

protected:
	unsigned long long N_total;
	search_space_point record_point;
	search_space_point global_record_point;
	std::chrono::high_resolution_clock::time_point start_chrono_time;
	// hill climbing
	search_space_point fromPointIndexesToPoint( std::vector<unsigned> cur_point_indexes );
	search_space_point fromDoubleVecToPoint(std::vector<double> double_vec);
	std::vector<search_space_point> checked_points;
};

#endif
