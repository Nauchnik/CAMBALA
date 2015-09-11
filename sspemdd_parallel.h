#ifndef SSPEMDD_PARALLEL_H
#define SSPEMDD_PARALLEL_H

#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include "sspemdd_sequential.h"

class sspemdd_parallel
{
public:
	sspemdd_parallel();
	~sspemdd_parallel();
	int rank;
	int corecount;
	double mpi_start_time;
	unsigned task_array_len;
	unsigned result_array_len;
	unsigned ncb;
	unsigned nrhob;
	unsigned nR;
	unsigned ncpl;
	double residual;
	double cb1;
	double cb2;
	double cw1;
	double cw2;
	double R1;
	double R2;
	double rhob1;
	double rhob2;
	double cb_cur;
	double R_cur;
	double rhob_cur;
	unsigned n_layers_w;
	double res_min;
	double cb_min;
	double rhob_min;
	double R_min;
	std::vector<std::vector<double>> cws_all_cartesians;
	std::vector<double> cws_cur;
	std::vector<double> cws_min;
	std::vector<double> cws_fixed;
	std::vector<double> c1s;
	std::vector<double> c2s;
	std::vector<double> rhos;
	std::vector<unsigned> Ns_points;
	std::vector<double> depths;
	std::vector<double> freqs;
	std::vector<std::vector<double>> modal_group_velocities;
	std::vector<unsigned> mode_numbers;
	std::vector<std::vector<double>> modal_delays;

	double *task_array;
	double *result_array;
	void control_process();
	void computing_process();
};

#endif