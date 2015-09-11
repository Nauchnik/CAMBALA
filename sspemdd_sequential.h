#ifndef SSPEMDD_SEQUENTIAL_H
#define SSPEMDD_SEQUENTIAL_H

#include <vector>
#include <algorithm>
#include "linalg.h"

const double LOCAL_M_PI = 3.14159265358979323846;

class sspemdd_sequential
{
public:
	std::vector<double> compute_wnumbers(double &omeg, std::vector<double> &c, std::vector<double> &rho,
		std::vector<unsigned> &interface_idcs, std::vector<double> &meshsizes,
		unsigned flOnlyTrapped);

	std::vector<double> compute_wnumbers_extrap(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, unsigned flOnlyTrapped, unsigned &ordRich);

	std::vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos,
		std::vector<unsigned> &Ns_points, unsigned flOnlyTrapped, unsigned &ordRich);

	int compute_modal_grop_velocities(std::vector<double> &freqs, double deltaf, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		unsigned flOnlyTrapped, unsigned &ordRich, std::vector<std::vector<double>> &modal_group_velocities,
		std::vector<unsigned> &mode_numbers);

	double compute_modal_delays_residual_uniform(std::vector<double> &freqs, std::vector<double> &depths, std::vector<double> &c1s,
		std::vector<double> &c2s, std::vector<double> &rhos, std::vector<unsigned> &Ns_points,
		double R, std::vector<std::vector<double>> &experimental_delays,
		std::vector<unsigned> &experimental_mode_numbers);
};

#endif