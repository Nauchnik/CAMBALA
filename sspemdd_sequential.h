#ifndef SSPEMDD_SEQUENTIAL_H
#define SSPEMDD_SEQUENTIAL_H

#include "linalg.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <vector>
#include <tuple>
#include <random>
#include <limits>

const double LOCAL_M_PI = 3.14159265358979323846;
const double START_HUGE_VALUE = 1e100;

struct Point
{
	double R;
	double tau;
	double rhob;
	double cb;
	std::vector<double> cws;
	double residual = std::numeric_limits<double>::infinity();

	bool operator==(const Point &a) const
	{
		return (R == a.R && tau == a.tau && rhob == a.rhob && cb == a.cb && cws == a.cws);
	}
	bool operator<(const Point& r)
	{
		return (residual < r.residual);
		/*
		return ((residual < r.residual) || 
			// strict weak ordering
			(std::tie(R, tau, rhob, cb, cws) < (std::tie(r.R, r.tau, r.rhob, r.cb, r.cws))));
			*/
	}
};

class sspemdd_sequential
{
  public:
	sspemdd_sequential();
	std::string object_function_type;
	double h;
	double H;
	unsigned long long ncb;
	unsigned long long nrhob;
	unsigned long long nR;
	unsigned long long ntau; // tau_comment: added tau to the class declaration
	double cb1;
	double cb2;
	double R1;
	double R2;
	double tau1; // tau_comment: added tau to the class declaration
	double tau2; // tau_comment: added tau to the class declaration
	double rhob1;
	double rhob2;
	unsigned long long n_layers_w;
	unsigned long long iterated_local_search_runs;
	std::vector<unsigned> mode_numbers;
	std::vector<std::vector<double> > modal_delays;
	std::vector<double> freqs;
	std::vector<double> depths;
	std::vector<double> c1s;
	std::vector<double> c2s;
	std::vector<double> rhos;
	std::vector<unsigned> Ns_points;
	std::vector<std::vector<double> > weight_coeffs;
	int verbosity;
	std::string dtimesFileName;
	std::string spmagFileName;
	int rank;
	std::vector<double> cw1_arr;
	std::vector<double> cw2_arr;
	std::vector<unsigned long long> ncpl_arr;

	std::mt19937 rng;

	// functions by Oleg
	std::vector<std::vector<double>> search_space; // values of variables which form a search space
	void readScenario(const std::string &scenarioFileName);
	void readInputDataFromFiles();
	void init();
	double getRecordResidual();
	double fill_data_compute_residual(Point &point);
	void findGlobalMinBruteForce();
	void loadValuesToSearchSpaceVariables();
	void findLocalMinHillClimbing();
	void report_final_result();
	void getThreeValuesFromStr(std::string str, double &val1, double &val2, double &val3);

	// methods by Vader
	Point generateRandomPoint();
	void ILSGPU(int ils_runs=1);

	// functions by Pavel
	// tau_comment: added tau to the arguments of compute_modal_delays_residual_uniform()
	double compute_modal_delays_residual_uniform(
		std::vector<double> &freqs, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		double R, double tau, std::vector<std::vector<double>> &experimental_delays,
		std::vector<unsigned> &experimental_mode_numbers);

	double compute_modal_delays_residual_weighted(
		std::vector<double> &freqs, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		double R, double tau, std::vector<std::vector<double>> &experimental_delays,
		std::vector<std::vector<double>> &weight_coeffs,
		std::vector<unsigned> &experimental_mode_numbers);

	std::vector<double> compute_wnumbers(double &omeg, std::vector<double> &c,
										 std::vector<double> &rho,
										 std::vector<unsigned> &interface_idcs,
										 std::vector<double> &meshsizes, unsigned flOnlyTrapped);

	std::vector<double> compute_wnumbers_extrap(double &omeg, std::vector<double> &depths,
												std::vector<double> &c1s, std::vector<double> &c2s,
												std::vector<double> &rhos,
												std::vector<unsigned> &Ns_points,
												unsigned flOnlyTrapped, unsigned &ordRich);

	std::vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, std::vector<double> &depths,
													   std::vector<double> &c1s,
													   std::vector<double> &c2s,
													   std::vector<double> &rhos,
													   std::vector<unsigned> &Ns_points,
													   unsigned flOnlyTrapped, unsigned &ordRich);

	int compute_modal_grop_velocities(std::vector<double> &freqs, double deltaf,
									  std::vector<double> &depths, std::vector<double> &c1s,
									  std::vector<double> &c2s, std::vector<double> &rhos,
									  std::vector<unsigned> &Ns_points, unsigned flOnlyTrapped,
									  unsigned &ordRich,
									  std::vector<std::vector<double>> &modal_group_velocities,
									  std::vector<unsigned> &mode_numbers);

	int compute_wnumbers_bb(std::vector<double> &freqs, double deltaf, std::vector<double> &depths,
							std::vector<double> &c1s, std::vector<double> &c2s,
							std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
							unsigned flOnlyTrapped, unsigned &ordRich,
							std::vector<std::vector<double>> &modal_group_velocities,
							std::vector<unsigned> &mode_numbers);

  protected:
	unsigned long long N_total;
	Point record_point;
	std::chrono::high_resolution_clock::time_point start_chrono_time;
	// hill climbing
	Point fromPointIndexesToPoint(std::vector<unsigned> cur_point_indexes);
	Point fromDoubleVecToPoint(std::vector<double> double_vec);
	std::vector<Point> checked_points;
};

#endif
