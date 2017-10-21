#ifndef CAMBALA_SEQUENTIAL_H
#define CAMBALA_SEQUENTIAL_H

#include <vector>
#include <complex>
#include <fstream>
#include <algorithm>
#include <chrono>
#include <math.h>
#include "linalg.h"

using namespace std;


class CAMBALA_sequential
{
public:
	CAMBALA_sequential();
	vector<double> cw1_init_arr;
	vector<double> cw2_init_arr;
	vector<double> d1_arr;
	vector<double> d2_arr;
	vector<double> d_step;
	vector<unsigned long long> ncpl_init_arr;
	unsigned ppm;
	double h1_;
	double h2_;
	double R1_;
	double R2_;
	double tau1_;    //tau_comment: added tau to the class declaration
	double tau2_;    //tau_comment: added tau to the class declaration
	double rhob1_;
	double rhob2_;
	unsigned long long n_layers_w;
	unsigned long long iterated_local_search_runs;
	int verbosity;
	int rank;
	bool isTimeDelayPrinting;
	string launch_type; // bruteforce | ils
	string output_filename;
	string depths_filename;
	
	// Oleg's functions
	//vector<vector<double>> search_space; // values of variables which form a search space
	int readScenario(string scenarioFileName);
	int readInputDataFromFiles();
	int init(vector<double> depths);
	void solve();
	int createDepthsArray(const double h, vector<vector<double>> &depths_vec);
	void loadValuesToSearchSpaceVariables();
	double getRecordResidual();
	Point getNonRandomStartPoint( vector<double> depths );
	double fillDataComputeResidual( Point &point );
	vector<Point> getSearchSpacePointsVec(vector<double> depths);
	Point findLocalMinHillClimbing(vector<double> depths);
	void findGlobalMinBruteForce(vector<double> depths);
	void reportFinalResult();
	void getThreeValuesFromStr(string str, double &val1, double &val2, double &val3);
	void reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a);
	double directPointCalc( Point point );
	void printDelayTime(double R, vector<unsigned> mode_numbers, vector<vector<double>> modal_group_velocities);

	//
	Point fromPointIndexesToPoint(vector<unsigned> cur_point_indexes, vector<double> depths);
	vector<unsigned> fromPointToPointIndexes(Point point);
	Point fromDoubleVecToPoint(vector<double> double_vec);
	Point fromStrToPoint(string str); // BOINC client application
	void fromPointToFile(const Point &point, ofstream &ofile); // BOINC client application

	// functions by Pavel
	// tau_comment: added tau to the arguments of compute_modal_delays_residual_uniform()

    void load_layers_data(string LayersFName,
	vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,
	vector<unsigned> &Ns_points);

	void load_profile_deep_water(string ProfileFName, vector<double> &depths, vector<double> &c1s, 
		vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points);

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

    double compute_modal_delays_residual_weighted2(vector<double> &freqs,
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

	vector<complex<double>> compute_cpl_pressure(double f, vector<double> &depths, vector<double> &c1s,
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
	Point record_point;
	chrono::high_resolution_clock::time_point start_chrono_time;
	// hill climbing
	vector<Point> checked_points;
};

#endif
