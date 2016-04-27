/*****************************************************************************************
// SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data -- Copyright(c) 2015-2016
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS), 
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#ifdef _MPI
#include <mpi.h>
#endif

#include "sspemdd_sequential.h"
#include "sspemdd_parallel.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

// layer data
struct layer
{
	std::vector<double> zend; // vector of finite depths
	double cbeg;			  // velocity of a sound at the beginning of a layer
	double cend;			  // velocity of a sound at the end of a layer
	double dbeg;              // density at the beginning of a layer
	double dend;              // density at the end of a layer
};

struct sspemdd_params
{
	double cb1;
	double cb2;
	unsigned ncb;
	double rhob1;
	double rhob2;
	unsigned nrhob;
	double R1;
	double R2;
	unsigned nR;
	double tau1;
	double tau2;
	unsigned ntau;
	double cw1;
	double cw2;
	unsigned ncpl;
	unsigned iterated_local_search_runs;
	unsigned n_layers_w;
	int launchType;
	bool isHomogeneousWaterLayer;
	std::vector<double> depths;
	std::vector<double> c1s;
	std::vector<double> c2s;
	std::vector<double> rhos;
	std::vector<unsigned>Ns_points;
	std::vector<unsigned>mode_numbers;
	std::vector<std::vector<double>> modal_delays;
	std::vector<double> freqs;
	int verbosity;
};

int main(int argc, char **argv)
{
	//SEARCH SPACE
	double cb1 = 2000;
	double cb2 = 2000;
	unsigned ncb = 1;

	double rhob1 = 2;
	double rhob2 = 2;
	unsigned nrhob = 1;

	double R1 = 3400;
	double R2 = 3600;
	unsigned nR = 41;

	double tau1 = 0;
	double tau2 = 0;
	unsigned ntau = 1;

	double cw1 = 1450;
	double cw2 = 1500;
	unsigned ncpl = 0; // search mesh within each water layer

	//std::vector<std::vector<double>> modal_group_velocities;
	std::vector<unsigned> mode_numbers;
	std::vector<std::vector<double>> modal_delays;
	std::vector<double> freqs;
	double buff;
	std::vector<double> buffvect;
	std::string myLine, myFileName = "8000_extracted.txt";
	int launchType = 0; // 0 - deafult; 1 - fixed cw1=1490; 2 - cw1>cw2>...>cwn; 3 - both 1 and 2
	sspemdd_sequential sspemdd_seq;
	unsigned iterated_local_search_runs = 10;
	int verbosity = 0;

#ifdef _DEBUG
	argc = 5;
	argv[1] = "8000_extracted.txt";
	//argv[1] = "dtimes_synth_thcline_hf.txt";
	argv[2] = "0"; // launchType
	argv[3] = "1"; // ncpl
	argv[4] = "10"; // iterated_local_search_runs
#endif
	
	if (argc >= 2) {
		myFileName = argv[1];
		std::cout << "myFileName " << myFileName << std::endl;
	}
	else {
		std::cout << "Usage : thcline_file [launchType] [ncpl] [iterated_local_search_runs] [verbosity]" << std::endl;
		return 1;
	}
	
	if (argc >= 3)
		launchType = atoi(argv[2]);

	if (argc >= 4) {
		ncpl = atoi(argv[3]);
		std::cout << "ncpl " << ncpl << std::endl;
	}

	if (argc >= 5) {
		iterated_local_search_runs = atoi(argv[4]);
		std::cout << "iterated_local_search_runs " << iterated_local_search_runs << std::endl;
	}

	if (argc >= 6) {
		verbosity = atoi(argv[5]);
		std::cout << "verbosity " << verbosity << std::endl;
	}
	
	std::ifstream myFileSynth(myFileName.c_str());
	std::stringstream myLineStream;
	// reading the "experimental" delay time data from a file
	while (std::getline(myFileSynth, myLine)) {
		myLine.erase(std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete windows endline symbol for correct reading
		myLineStream << myLine;
		myLineStream >> buff;
		freqs.push_back(buff);

		buffvect.clear();
		while (!myLineStream.eof()) {
			myLineStream >> buff;
			buffvect.push_back(buff);
			mode_numbers.push_back((unsigned)buffvect.size());
		}

		modal_delays.push_back(buffvect);
		myLineStream.str(""); myLineStream.clear();
	}
	myFileSynth.close();

	/*  //SYMTHETIC TEST #1
		//OLD version: bottom halfspace + Range


		vector<double> depths{ 90, 600 };
		vector<double> c1s{ 1500, 2000 };
		vector<double> c2s{ 1500, 2000 };
		vector<double> rhos{ 1, 2 };
		vector<unsigned> Ns_points{ 180, 1020 };

		// END TEST #1
		*/

		// SYNTHETIC TEST #2
		// new version: ssp in water + bottom halfspace + Range

		// environment model
		// waveguide depth: H
		// water column depth: h
		// n_layers_w in the water,

	bool isHomogeneousWaterLayer = false;

	unsigned n_layers_w;
	if (myFileName == "8000_extracted.txt") {
		isHomogeneousWaterLayer = true;
		n_layers_w = 1;
	}
	else
		n_layers_w = 5;

	unsigned ppm = 2;
	double h = 90;
	double H = 600;
	double layer_thickness_w = h / n_layers_w;
	//int layer_np = round(h/n_layers_w);

	std::vector<double> depths;
	for (unsigned jj = 1; jj <= n_layers_w; jj++) 
		depths.push_back(layer_thickness_w*jj);
	depths.push_back(H);

	std::vector<double> c1s(n_layers_w + 1, 1500);
	std::vector<double> c2s(n_layers_w + 1, 1500);
	std::vector<double> rhos(n_layers_w + 1, 1);
	std::vector<unsigned> Ns_points(n_layers_w + 1, (unsigned)round(ppm*layer_thickness_w));
	Ns_points.at(n_layers_w) = (unsigned)round(ppm*(H - h));

	// END TEST #

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

	//TEST BLOCK! PLEASE DONT REMOVE, COMMENT IF NECESSARY
	/*std::cout << "TEST BLOCK" << std::endl;
	std::vector<double> c1s_t    { 1490, 1490, 1480, 1465, 1460, 2000};
	std::vector<double> c2s_t    { 1490, 1480, 1465, 1460, 1460, 2000};
	std::vector<double> rhos_t   { 1, 1, 1, 1, 1, 2};
	std::vector<double> depths_t   { 18, 36, 54, 72, 90, 600};
	std::vector<unsigned> Ns_points_t  { 36, 36, 36, 36, 36, 1020};
	double deltaf_t = 0.5;
	std::vector<std::vector<double>> modal_group_velocities_t;
	std::vector<unsigned> mode_numbers_t;
	double R_t = 3500;
	double residual_t;
	unsigned rord = 3;

	sspemdd_seq.compute_modal_grop_velocities( freqs, deltaf_t, depths_t, c1s_t, c2s_t, rhos_t, Ns_points_t, 1, rord, modal_group_velocities_t, mode_numbers_t );

	for (unsigned ii=0; ii<freqs.size();  ii++) {
		for (unsigned jj=0; jj<mode_numbers_t.at(ii);  jj++) {
			std::cout << modal_group_velocities_t[ii][jj] << "; ";
		}
		std::cout << std::endl;
	}

	residual_t = sspemdd_seq.compute_modal_delays_residual_uniform(freqs, depths_t, c1s_t, c2s_t, rhos_t, Ns_points, R_t, modal_delays, mode_numbers);

	std::cout << " TEST comparison." << std::endl << "RESIDUAL: " << residual_t << std::endl << std::endl;;
	*/
	//END OF TEST BLOCK!

	// set other parameters of the search space
	if ((launchType >= 4) && (launchType <= 6)) {
		nR = 1;
		R1 = 3500;
		R2 = 3500;
	}
	
	if ((launchType == 5) || (launchType == 7)) {
		cb1 = cb2 = 3000;
		ncb = 1;
		rhob1 = rhob2 = 3;
		nrhob = 1;
	}
	else if ((launchType == 6) || (launchType == 8)) {
		cb1 = cb2 = 4000;
		ncb = 1;
		rhob1 = rhob2 = 4;
		nrhob = 1;
	}

	if (isHomogeneousWaterLayer) {
		cb1 = 1600;
		cb2 = 1900;
		ncb = 61;
		rhob1 = 1.4;
		rhob2 = 2.0;
		nrhob = 9;
		tau1 = -0.015;
		tau2 = 0.015;
		ntau = 301;
		R1 = R2 = 8000;
		nR = 1;
		cw1 = cw2 = 1500;
		ncpl = 1;
	}
	
	std::cout << "Input parameters :" << std::endl;
	std::cout << "ncpl " << ncpl << std::endl;
	std::cout << "n_layers_w " << n_layers_w << std::endl;
	std::cout << "nR " << nR << std::endl;
	std::cout << "ntau " << ntau << std::endl;
	std::cout << "nrhob " << nrhob << std::endl;
	std::cout << "ncb " << ncb << std::endl;
	
#ifndef _MPI
	// sequential mode

	sspemdd_seq.cb1 = cb1;
	sspemdd_seq.cb2 = cb2;
	sspemdd_seq.ncb = ncb;
	sspemdd_seq.rhob1 = rhob1;
	sspemdd_seq.rhob2 = rhob2;
	sspemdd_seq.nrhob = nrhob;
	sspemdd_seq.R1 = R1;
	sspemdd_seq.R2 = R2;
	sspemdd_seq.nR = nR;
	sspemdd_seq.tau1 = tau1;
	sspemdd_seq.tau2 = tau2;
	sspemdd_seq.ntau = ntau;
	sspemdd_seq.cw1 = cw1;
	sspemdd_seq.cw2 = cw2;
	sspemdd_seq.ncpl = ncpl;
	sspemdd_seq.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_seq.n_layers_w = n_layers_w;
	sspemdd_seq.launchType = launchType;
	sspemdd_seq.isHomogeneousWaterLayer = isHomogeneousWaterLayer;
	sspemdd_seq.depths = depths;
	sspemdd_seq.c1s = c1s;
	sspemdd_seq.c2s = c2s;
	sspemdd_seq.rhos = rhos;
	sspemdd_seq.Ns_points = Ns_points;
	sspemdd_seq.mode_numbers = mode_numbers;
	sspemdd_seq.modal_delays = modal_delays;
	sspemdd_seq.freqs = freqs;
	sspemdd_seq.verbosity = verbosity;

	sspemdd_seq.init();
	
	//sspemdd_seq.findGlobalMinBruteForce();
	sspemdd_seq.findLocalMinHillClimbing();
	sspemdd_seq.report_final_result();

	// fix final time
	/*t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	
	std::cout << std::endl;
	std::cout << "total solving time " << time_span.count() << std::endl;*/
#else
	int rank = 0;
	int corecount = 1;
	
	// parallel mode
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	sspemdd_parallel sspemdd_par;
	sspemdd_par.rank = rank;
	sspemdd_par.corecount = corecount;

	sspemdd_par.cb1 = cb1;
	sspemdd_par.cb2 = cb2;
	sspemdd_par.ncb = ncb;
	sspemdd_par.rhob1 = rhob1;
	sspemdd_par.rhob2 = rhob2;
	sspemdd_par.nrhob = nrhob;
	sspemdd_par.R1 = R1;
	sspemdd_par.R2 = R2;
	sspemdd_par.nR = nR;
	sspemdd_par.tau1 = tau1;
	sspemdd_par.tau2 = tau2;
	sspemdd_par.ntau = ntau;
	sspemdd_par.cw1 = cw1;
	sspemdd_par.cw2 = cw2;
	sspemdd_par.ncpl = ncpl;
	sspemdd_par.iterated_local_search_runs = iterated_local_search_runs;
	sspemdd_par.n_layers_w = n_layers_w;
	sspemdd_par.launchType = launchType;
	sspemdd_par.isHomogeneousWaterLayer = isHomogeneousWaterLayer;
	sspemdd_par.depths = depths;
	sspemdd_par.c1s = c1s;
	sspemdd_par.c2s = c2s;
	sspemdd_par.rhos = rhos;
	sspemdd_par.Ns_points = Ns_points;
	sspemdd_par.mode_numbers = mode_numbers;
	sspemdd_par.modal_delays = modal_delays;
	sspemdd_par.freqs = freqs;
	sspemdd_par.verbosity = verbosity;
	
	sspemdd_par.init();
	sspemdd_par.MPI_main();
#endif

	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	std::cout << "main() total time " << time_span.count() << " s" << std::endl;
	
	return 0;
}