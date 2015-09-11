/*************************************************************************
// SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data -- Copyright(c) 2015
// Pavel Petrov(V.I.Il'ichev Pacific Oceanological Institute), 
// Oleg Zaikin(V.M.Matrosov Institute for System Dynamics and Control Theory)
*************************************************************************/

#ifdef _MPI
#include <mpi.h>
#endif

#include "sspemdd_sequential.h"
#include "sspemdd_parallel.h"
#include "sspemdd_utils.h"

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
	double cbeg;         // velocity of a sound at the beginning of a layer
	double cend;         // velocity of a sound at the end of a layer
	double dbeg;         // density at the beginning of a layer
	double dend;         // density at the end of a layer
};

int main(int argc, char **argv)
{
	unsigned rord = 3;
	std::vector<std::vector<double>> modal_group_velocities;
	std::vector<unsigned> mode_numbers;
	std::vector<std::vector<double>> modal_delays;
	std::vector<double> freqs;
	double buff;
	std::vector<double> buffvect;
	std::string myLine, myFileName = "dtimes_synth_thcline_hhf.txt";
	int launchType = 0; // 0 - deafult, 1: fixed cw1=1490, 2: cw1>cw2>...>cwn, 3 : both 1 and 2
	sspemdd_sequential sspemdd_seq;
	
	if (argc >= 2) {
		myFileName = argv[1];
		std::cout << "myFileName " << myFileName << std::endl;
	}

	if (argc >= 3)
		launchType = atoi(argv[2]);
	
	std::ifstream myFileSynth(myFileName.c_str());

	// reading the "experimental" delay time data from a file
	while (std::getline(myFileSynth, myLine)){
		myLine.erase(std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete window symbol for correct reading
		std::stringstream myLineStream(myLine);
		myLineStream >> buff;
		freqs.push_back(buff);
		buffvect.clear();

		while (!myLineStream.eof()) {
			myLineStream >> buff;
			buffvect.push_back(buff);
			mode_numbers.push_back((unsigned)buffvect.size());
		}

		modal_delays.push_back(buffvect);
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

	unsigned ppm = 2;
	double h = 90;
	double H = 600;
	unsigned n_layers_w = 5;
	double layer_thickness_w = h / n_layers_w;
	//int layer_np = round(h/n_layers_w);

	std::vector<double> depths;
	for (unsigned jj = 1; jj <= n_layers_w; jj++){ depths.push_back(layer_thickness_w*jj); }
	depths.push_back(H);

	std::vector<double> c1s(n_layers_w + 1, 1500);
	std::vector<double> c2s(n_layers_w + 1, 1500);
	std::vector<double> rhos(n_layers_w + 1, 1);
	std::vector<unsigned> Ns_points(n_layers_w + 1, (unsigned)round(ppm*layer_thickness_w));
	Ns_points.at(n_layers_w) = (unsigned)round(ppm*(H - h));

	// END TEST #2

	double residual = 1e50; // start huge value
	double R = 3500;

	//SEARCH SPACE
	double cb1 = 2000;
	double cb2 = 2000;
	double cb_min = 1e50, cb_cur;
	unsigned ncb = 0;

	double rhob1 = 2;
	double rhob2 = 2;
	double rhob_min = 1e50, rhob_cur;
	unsigned nrhob = 0;

	double R1 = 3400;
	double R2 = 3600;
	double R_min = 1e50, R_cur = 0;
	unsigned nR = 0;
	
	double cw1 = 1450;
	double cw2 = 1500;
	unsigned ncpl; // search mesh within each water layer
	std::vector<double> cws_cur(n_layers_w, 1500);
	std::vector<double> cws_min(n_layers_w, 1500);
	std::vector<double> cws_fixed{ 1490, 1490, 1480, 1465, 1460 }; // use when ncpl = 1
	
	if (cws_fixed.size() != n_layers_w) {
		std::cerr << "cws_fixed.size() != n_layers_w" << std::endl;
		std::cerr << cws_fixed.size() << " != " << n_layers_w << std::endl;
		exit(1);
	}
	
	double res_min = 1e50;

	// fix start time
	std::chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;

#ifdef _MPI
	// parallel mode
	int rank;
	int corecount;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	sspemdd_parallel sspemdd_par;
	sspemdd_par.rank = rank;
	sspemdd_par.corecount = corecount;
	sspemdd_par.res_min = res_min;
	sspemdd_par.cb_min = cb_min;
	sspemdd_par.rhob_min = rhob_min;
	sspemdd_par.R_min = R_min;
	sspemdd_par.ncb = 1;
	sspemdd_par.nrhob = 1;
	sspemdd_par.nR = 41;
	sspemdd_par.ncpl = 10;
	sspemdd_par.residual = residual;
	sspemdd_par.cb1 = cb1;
	sspemdd_par.cb2 = cb2;
	sspemdd_par.cw1 = cw1;
	sspemdd_par.cw2 = cw2;
	sspemdd_par.R1 = R1;
	sspemdd_par.R2 = R2;
	sspemdd_par.rhob1 = rhob1;
	sspemdd_par.rhob2 = rhob2;
	sspemdd_par.n_layers_w = n_layers_w;
	sspemdd_par.cws_cur = cws_cur;
	sspemdd_par.cws_min = cws_min;
	sspemdd_par.cws_fixed = cws_fixed; // use when ncpl = 1
	sspemdd_par.c1s = c1s;
	sspemdd_par.c2s = c2s;
	sspemdd_par.rhos = rhos;
	sspemdd_par.Ns_points = Ns_points;
	sspemdd_par.depths = depths;
	sspemdd_par.freqs = freqs;
	sspemdd_par.modal_group_velocities = modal_group_velocities;
	sspemdd_par.mode_numbers = mode_numbers;
	sspemdd_par.modal_delays = modal_delays;
#else
	// sequential mode
	// small values for fast checking of correctness
	ncb = 1;
	nrhob = 1;
	nR = 2;
	ncpl = 2;

 /*   //TEST BLOCK! PLEASE DONT REMOVE, COMMENT IF NECESSARY
    vector<double> c1s_t    { 1490, 1490, 1480, 1465, 1460, 2000};
    vector<double> c2s_t    { 1490, 1480, 1465, 1460, 1460, 2000};

	vector<double> rhos_t   { 1, 1, 1, 1, 1, 2};
	vector<double> depths_t   { 18, 36, 54, 72, 90, 600};
	vector<unsigned> Ns_points_t  { 36, 36, 36, 36, 36, 1020};
	double deltaf_t = 0.5;
	vector<vector<double>> modal_group_velocities_t;
	vector<unsigned> mode_numbers_t;
	double R_t = 3500;
	double residual_t;

    compute_modal_grop_velocities( freqs, deltaf_t, depths_t, c1s_t, c2s_t, rhos_t, Ns_points_t, 1, rord, modal_group_velocities_t, mode_numbers_t );

    for (unsigned ii=0; ii<freqs.size();  ii++) {

        for (unsigned jj=0; jj<mode_numbers_t.at(ii);  jj++) {
			cout << modal_group_velocities_t[ii][jj] << "; ";
        }
        cout << endl;
    }

    residual_t = wnumbers_obj.compute_modal_delays_residual_uniform(freqs, depths_t, c1s_t, c2s_t, rhos_t, Ns_points, R_t, modal_delays, mode_numbers);

    cout << " TEST comparison." << endl << "RESIDUAL: " << residual_t << endl << endl;;

    ncb = 1;
    nrhob = 1;
    nR = 1;
    ncpl = 6;
    cb1 = 2000;
    cb2 = 2000;

    rhob1 = 2;
    rhob2 = 2;

    R1 = 3500;
    R2 = 3500;

    cw1 = 1460;
    cw2 = 1490;

    //END OF TEST BLOCK! */
#endif
	
	if (argc >= 4) {
		ncpl = atoi(argv[3]);
		std::cout << "new ncpl " << ncpl << std::endl;
	}
	
	std::cout << "Input parameters :" << std::endl;
	std::cout << "ncpl " << ncpl << std::endl;
	std::cout << "n_layers_w " << n_layers_w << std::endl;
	std::cout << "nR " << nR << std::endl;
	std::cout << "nrhob " << nrhob << std::endl;
	std::cout << "ncb " << ncb << std::endl;
	// Calculate total number of points in a search space
	unsigned N_total = (unsigned)round(pow(ncpl, n_layers_w))*nR*nrhob*ncb;
	std::cout << "Formula for calculating total number of points in a search space: (ncpl^n_layers_w)*nR*nrhob*ncb" << std::endl;
	std::cout << "N_total " << N_total << std::endl;
	std::cout << "launchType " << launchType << std::endl;
	
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
	
	// make cws_all_cartesians - all cartesians of all speeds in water
	std::vector<std::vector<double>> cws_vii; // all variants for every depth
	std::vector<int> index_arr;
	std::vector<double> cws_vi;
	std::vector<std::vector<double>> cws_all_cartesians;
	bool isAdding;
	if (ncpl == 1)
		cws_all_cartesians.push_back(cws_fixed);
	else {
		for (unsigned ncpl_cur = 0; ncpl_cur < ncpl; ncpl_cur++)
			cws_vi.push_back(cw1 + ncpl_cur*(cw2 - cw1) / (ncpl - 1));
		for (unsigned i = 0; i < n_layers_w; i++)
			cws_vii.push_back(cws_vi);
		while (SSPEMDD_utils::next_cartesian(cws_vii, index_arr, cws_vi)) {
			switch (launchType){
			case 1 :
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 2:
				isAdding = true;
				for (unsigned cws_vi_index = 0; cws_vi_index < cws_vi.size()-1; cws_vi_index++)
					if (cws_vi[cws_vi_index] <= cws_vi[cws_vi_index+1])
						isAdding = false;
				if (isAdding)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 3:
				if ( cws_vi[0] != 1490 )
					isAdding = false;
				else {
					isAdding = true;
					for (unsigned cws_vi_index = 0; cws_vi_index < cws_vi.size() - 1; cws_vi_index++)
						if (cws_vi[cws_vi_index] <= cws_vi[cws_vi_index + 1])
							isAdding = false;
				}
				if (isAdding)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 7: // like case 1
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 8: // like case 1
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			default:
				cws_all_cartesians.push_back(cws_vi);
				break;
			}
		}
	}
	
	std::cout << "cws_all_cartesians.size() " << cws_all_cartesians.size() << std::endl;
	
#ifdef _MPI
	sspemdd_par.cws_all_cartesians = cws_all_cartesians;
#endif

#ifndef _MPI
	// sequential mode

	// inverting for bottom halfspace parameters + sound speed in water
	for (unsigned cur_ncb = 0; cur_ncb < ncb; cur_ncb++)
		for (unsigned cur_nrhob = 0; cur_nrhob < nrhob; cur_nrhob++)
			for (unsigned cur_nR = 0; cur_nR < nR; cur_nR++) {
				// specify bottom parameters;
				if (ncb > 1) {cb_cur = cb1 + cur_ncb*(cb2 - cb1) / (ncb-1);}
				else { cb_cur = cb1; }
				if (nrhob > 1) {rhob_cur = rhob1 + cur_nrhob  *(rhob2 - rhob1) / (nrhob-1);}
                else { rhob_cur = rhob1; }
				// specify range
				if (nR > 1) {R_cur = R1 + cur_nR*(R2 - R1) / (nR-1);}
				else {R_cur = R1;}
				for (auto &cws_cur : cws_all_cartesians) { // and finally specify sound speed in water
					// the parameters are transformed into the arrays c1s, c2s, rhos
					for (unsigned jj = 0; jj < n_layers_w - 1; jj++){
						c1s.at(jj) = cws_cur.at(jj);
						c2s.at(jj) = cws_cur.at(jj + 1);
						rhos.at(jj) = 1;
					}
					c1s.at(n_layers_w - 1) = cws_cur.at(n_layers_w - 1);
					c2s.at(n_layers_w - 1) = cws_cur.at(n_layers_w - 1);
					rhos.at(n_layers_w - 1) = 1;
					c1s.at(n_layers_w) = cb_cur;
					c2s.at(n_layers_w) = cb_cur;
					rhos.at(n_layers_w) = rhob_cur;

					for (unsigned jj = 0; jj <= n_layers_w; jj++){
						std::cout << "Layer #" << jj + 1 << ": c=" << c1s.at(jj) << "..." << c2s.at(jj) << "; rho=" << rhos.at(jj) << "; np=" << Ns_points.at(jj) << std::endl;
					}
					
					residual = sspemdd_seq.compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points, R_cur, modal_delays, mode_numbers);
					std::cout << residual << std::endl << std::endl;

					if (residual < res_min) {
						res_min  = residual;
						cb_min   = cb_cur;
						rhob_min = rhob_cur;
						R_min    = R_cur;
						cws_min  = cws_cur;
						std::cout << std::endl;
						std::cout << std::endl << "New residual minimum:" << std::endl;
						std::cout << "err=" << res_min << ", parameters:" << std::endl;
						std::cout << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << std::endl;
						std::cout << "cws_min :" << std::endl;
						for (auto &x : cws_min)
							std::cout << x << " ";
					}
				}
			}

	// fix final time
	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	std::cout << "SEARCH ENDED!" << std::endl;
	std::cout << "RESULTING VALUE:" << std::endl;
	std::cout << "err=" << res_min << ", parameters:" << std::endl;
	std::cout << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << std::endl;
	std::cout << "time " << time_span.count() << std::endl;

#else
	// MPI mode
	if (rank == 0)
		sspemdd_par.control_process();
	else if (rank > 0)
		sspemdd_par.computing_process();
#endif
	return 0;
}
