// +-----------------------------------------------------------------------------------+
// | Client application for the volunteer computing project Acoustics@home             |
// +-----------------------------------------------------------------------------------+
// | Pacific Oceanological Institute, Institute for System Dynamics and Control Theory |
// +-----------------------------------------------------------------------------------+
// | Authors: Pavel Petrov, Oleg Zaikin                                                |
// +-----------------------------------------------------------------------------------+

#ifdef _MPI
#include <mpi.h>
#endif

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "linalg.h"

const double LOCAL_M_PI = 3.14159265358979323846;

using namespace std;

// layer data
struct layer
{
	vector<double> zend; // vector of finite depths
	double cbeg;         // velocity of a sound at the beginning of a layer
	double cend;         // velocity of a sound at the end of a layer
	double dbeg;         // density at the beginning of a layer
	double dend;         // density at the end of a layer
};

template< typename T >
bool next_cartesian(std::vector<T> &vii, std::vector<int> &index_arr, T &cur_vi)
{
	if (index_arr.size() == 0) { // init
		index_arr.resize(vii.size());
		//for( auto &x : index_arr )
		//	x = 0;
		for (std::vector<int> ::iterator it = index_arr.begin(); it != index_arr.end(); ++it)
			*it = 0;
	}
	if (index_arr[0] == -1)
		return false;
	// get current value
	cur_vi.resize(vii.size());
	for (unsigned i = 0; i < index_arr.size(); ++i)
		cur_vi[i] = vii[i][index_arr[i]];
	// check if last iteration
	bool IsLastValue = true;
	for (unsigned i = 0; i < index_arr.size(); ++i) {
		if (index_arr[i] != vii[i].size() - 1) {
			IsLastValue = false;
			break;
		}
	}
	if (IsLastValue)
		index_arr[0] = -1; // condition of stopping
	else {
		// find last changable row to increase its value
		unsigned last_changable = (unsigned)index_arr.size() - 1;
		while (last_changable != -1){
			if (index_arr[last_changable] < vii[last_changable].size() - 1)
				break;
			--last_changable;
		}
		index_arr[last_changable]++;
		for (unsigned i = last_changable + 1; i < index_arr.size(); ++i)
			index_arr[i] = 0;
	}

	return true;
}

vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho, vector<unsigned> &interface_idcs, vector<double> &meshsizes,unsigned flOnlyTrapped);
vector<double> compute_wnumbers_extrap(double &omeg, vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, unsigned flOnlyTrapped,unsigned &ordRich);
vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, unsigned flOnlyTrapped,unsigned &ordRich);
int compute_modal_grop_velocities( vector<double> &freqs, double deltaf, vector<double> &depths, vector<double> &c1s, vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, unsigned flOnlyTrapped, unsigned &ordRich, vector<vector<double>> &modal_group_velocities, vector<unsigned> &mode_numbers );
double compute_modal_delays_residual_uniform( vector<double> &freqs, vector<double> &depths, vector<double> &c1s, vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, double R, vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);

int main(int argc, char **argv)
{
	/*
		double freq = 50;
		double c_w = 1500;
		double c_b = 2000;
		double rho_w = 1;
		double rho_b = 2;
		double dz = 1;
		unsigned nz = 501;
		unsigned ib = 90;    //at POINT = 89 we have water, at POINT = 90 we have bottom
		//   ii = 89,                  at ii = 90

		cout.precision(15);

		double omeg = 2*LOCAL_M_PI*freq;

		vector<double> input_c;
		vector<double> input_rho;
		vector<double> input_mesh { dz,dz };
		vector<unsigned> input_interf_idcs { ib };
		vector<double> out_wnum;

		for ( unsigned ii = 0; ii<nz; ii++ ){
		if (ii<ib){                   //
		input_c.push_back(c_w);
		input_rho.push_back(rho_w);
		}
		else {
		input_c.push_back(c_b);
		input_rho.push_back(rho_b);
		}

		}

		cout << "freq = " << freq << endl;
		cout << "omega = " << omeg << endl;

		out_wnum = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,1);

		for (unsigned ii=0; ii<out_wnum.size();  ii++) {

		cout << ii << "->" << out_wnum.at(ii) << endl;


		}

		cout << "NEW: extrapolation" << endl;

		vector<double> depths {90,600};
		vector<double> c1s  {1500,2000};
		vector<double> c2s  {1500,2000};
		vector<double> rhos  {1,2};
		vector<unsigned> Ns_points  {180,1020};
		unsigned rord = 3;
		vector<double> freqs {20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
		vector<vector<double>> modal_group_velocities;
		vector<unsigned> mode_numbers;
		double deltaf = 0.5;

		out_wnum = compute_wnumbers_extrap_lin_dz(omeg,depths,c1s,c2s,rhos,Ns_points,1,rord);

		cout << "Extrapolated ev:" << endl;
		for (unsigned ii=0; ii<out_wnum.size();  ii++) {

		cout << ii << "->" << out_wnum.at(ii) << endl;


		}

		compute_modal_grop_velocities( freqs, deltaf, depths,c1s, c2s, rhos, Ns_points, 1, rord, modal_group_velocities, mode_numbers );

		ofstream myFile("mgv.txt");
		//    for (int ii=0; ii<N_points-2; ii++){
		//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
		//    }
		//

		for (unsigned ii=0; ii<freqs.size();  ii++) {
		cout << "f=" << freqs.at(ii) << endl;

		for (unsigned jj=0; jj<mode_numbers.at(ii);  jj++) {

		cout << modal_group_velocities[ii][jj] << endl;
		myFile << std::fixed << std::setprecision(16) << modal_group_velocities[ii][jj] << "\t";

		}
		myFile << endl;
		}

		myFile.close();

		*/

	// reading the "experimental" delay time data from a file:

	unsigned rord = 3;
	vector<vector<double>> modal_group_velocities;
	vector<unsigned> mode_numbers;
	vector<vector<double>> modal_delays;
	vector<double> freqs;
	freqs.clear();
	double buff;
	vector<double> buffvect;
	mode_numbers.clear();
	string myLine, myFileName = "dtimes_synth_thcline_hhf.txt";
	int launchType = 0; // 0 - deafult, 1: fixed cw1=1490, 2: cw1>cw2>...>cwn, 3 : both 1 and 2
	
	if (argc >= 2) {
		myFileName = argv[1];
		cout << "myFileName " << myFileName << endl;
	}

	if (argc >= 3)
		launchType = atoi(argv[2]);
	
	ifstream myFileSynth(myFileName.c_str());

	while (getline(myFileSynth, myLine)){
		myLine.erase(std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete window symbol for correct reading
		stringstream myLineStream(myLine);
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

	vector<double> depths;
	for (unsigned jj = 1; jj <= n_layers_w; jj++){ depths.push_back(layer_thickness_w*jj); }
	depths.push_back(H);

	vector<double> c1s(n_layers_w + 1, 1500);
	vector<double> c2s(n_layers_w + 1, 1500);
	vector<double> rhos(n_layers_w + 1, 1);
	vector<unsigned> Ns_points(n_layers_w + 1, (unsigned)round(ppm*layer_thickness_w));
	Ns_points.at(n_layers_w) = (unsigned)round(ppm*(H - h));

	// END TEST #2

	double residual = 1e50; // start huge value
	double R = 3500;

	//SEARCH SPACE
	double cb1 = 2000;
	double cb2 = 2000;
	double cb_min = 1e50, cb_cur;
	unsigned ncb;

	double rhob1 = 2;
	double rhob2 = 2;
	double rhob_min = 1e50, rhob_cur;
	unsigned nrhob;

	double R1 = 3400;
	double R2 = 3600;
	double R_min = 1e50, R_cur;
	unsigned nR;

	double cw1 = 1450;
	double cw2 = 1500;
	unsigned ncpl; // search mesh within each water layer
	vector<double> cws_cur(n_layers_w, 1500);
	vector<double> cws_min(n_layers_w, 1500);
	vector<double> cws_fixed{ 1490, 1490, 1480, 1465, 1460 }; // use when ncpl = 1
	
	if (cws_fixed.size() != n_layers_w) {
		cerr << "cws_fixed.size() != n_layers_w" << endl;
		cerr << cws_fixed.size() << " != " << n_layers_w << endl;
		exit(1);
	}
	
	double res_min = 1e50;

	// fix start time
	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	chrono::high_resolution_clock::time_point t2;
	chrono::duration<double> time_span;

#ifdef _MPI
	// parallel mode
	int corecount, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &corecount);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	double mpi_start_time = MPI_Wtime();
	// large values
	ncb = 1;
	nrhob = 1;
	nR = 41;
	ncpl = 10;
#else
	// sequential mode
	// small values for fast checking of correctness
	ncb = 1;
	nrhob = 1;
	nR = 41;
	ncpl = 10;

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


    residual_t = compute_modal_delays_residual_uniform(freqs, depths_t, c1s_t, c2s_t, rhos_t, Ns_points, R_t, modal_delays, mode_numbers);

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
		cout << "new ncpl " << ncpl << endl;
	}
	
	unsigned N_total = (unsigned)round(pow(ncpl, n_layers_w))*nR*nrhob*ncb;
	cout << "N_total " << N_total << endl;
	cout << "launchType " << launchType << endl;

	if (launchType == 4)
		nR = 1;
	
	// make cws_all_cartesians - all cartesians of all speeds in water
	vector<vector<double>> cws_vii; // all variants for every depth
	vector<int> index_arr;
	vector<double> cws_vi;
	vector<vector<double>> cws_all_cartesians;
	bool isAdding;
	if (ncpl == 1)
		cws_all_cartesians.push_back(cws_fixed);
	else {
		for (unsigned ncpl_cur = 0; ncpl_cur < ncpl; ncpl_cur++)
			cws_vi.push_back(cw1 + ncpl_cur*(cw2 - cw1) / (ncpl - 1));
		for (unsigned i = 0; i < n_layers_w; i++)
			cws_vii.push_back(cws_vi);
		while (next_cartesian(cws_vii, index_arr, cws_vi)) {
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
			default:
				cws_all_cartesians.push_back(cws_vi);
				break;
			}
			
		}
	}

	cout << "cws_all_cartesians.size() " << cws_all_cartesians.size() << endl;

#ifndef _MPI
	// sequential mode

	// inverting for bottom halfspace parameters + sound speed in water!
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
						cout << "Layer #" << jj + 1 << ": c=" << c1s.at(jj) << "..." << c2s.at(jj) << "; rho=" << rhos.at(jj) << "; np=" << Ns_points.at(jj) << endl;
					}

					residual = compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points, R_cur, modal_delays, mode_numbers);
					cout << residual << endl << endl;

					if (residual < res_min) {
						res_min  = residual;
						cb_min   = cb_cur;
						rhob_min = rhob_cur;
						R_min    = R_cur;
						cws_min  = cws_cur;
						cout << endl;
						cout << endl << "New residual minimum:" << endl;
						cout << "err=" << res_min << ", parameters:" << endl;
						cout << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << endl;
						cout << "cws_min :" << endl;
						for (auto &x : cws_min)
							cout << x << " ";
					}
				}
			}

	// fix final time
	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	cout << "SEARCH ENDED!" << endl;
	cout << "RESULTING VALUE:" << endl;
	cout << "err=" << res_min << ", parameters:" << endl;
	cout << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min <<  endl;
	cout << "time " << time_span.count() << endl;

#else
	// MPI mode
	const int result_array_len = 9;
	double result_array[result_array_len];
	const int task_array_len = n_layers_w;
	double task_array[task_array_len];
	
	if ( rank == 0 ) {
		// sequential mode
		cout << "Start residual is: " << residual << endl;

		// brute force minimum search
		cout << "BRUTE FORCE MINIMUM SEARCH" << endl;
		cout << "Search space:" << endl;
		cout << cb1 << "< c_b <" << cb2 << ", step " << (cb2 - cb1) / ncb << endl;
		cout << R1 << "< Range <" << R2 << ", step " << (R2 - R1) / nR << endl;
		cout << cw1 << "< cws <" << cw2 << " , step " << (cw2 - cw1) / ncpl << endl;
		cout << rhob1 << "< rho_b <" << rhob2 << ", step " << (rhob2 - rhob1) / nrhob << endl;

		stringstream sstream_out;
		sstream_out << "MPI control process" << std::endl;
		std::chrono::high_resolution_clock::time_point t1, t2, finding_new_bkv_start_time, now_time;
		std::chrono::duration<double> time_span;

		sstream_out << "cws_all_cartesians.size() " << cws_all_cartesians.size() << std::endl;
		if (cws_all_cartesians.size() < corecount) {
			std::cerr << "cws_all_cartesians.size() < corecount" << std::endl;
			std::cerr << cws_all_cartesians.size() << " < " << corecount << std::endl;
			exit(1);
		}

		unsigned send_task_count = 0;
		unsigned processed_task_count = 0;
		
		// sending first part of tasks
		for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
			for (int j=0; j < task_array_len; j++)
				task_array[j] = cws_all_cartesians[send_task_count][j];
			MPI_Send(task_array, task_array_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
			send_task_count++;
		}
		sstream_out << "send_task_count " << send_task_count << std::endl;
		std::ofstream ofile( "mpi_out" );
		ofile << sstream_out.rdbuf();
		sstream_out.clear(); sstream_out.str("");
		ofile.close(); ofile.clear();

		// get results and send new tasks on idle computing processes
		while ( processed_task_count < cws_all_cartesians.size() ) {
			MPI_Recv( result_array, result_array_len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
			processed_task_count++;

			residual   = result_array[0];
			cb_cur     = result_array[1];
			rhob_cur   = result_array[2];
			R_cur      = result_array[3];
			cws_cur[0] = result_array[4];
			cws_cur[1] = result_array[5];
			cws_cur[2] = result_array[6];
			cws_cur[3] = result_array[7];
			cws_cur[4] = result_array[8];

			/*sstream_out << "recv residual " << residual << std::endl;
			sstream_out << "recv cb_cur "   << cb_cur   << std::endl;
			sstream_out << "recv rhob_cur " << rhob_cur << std::endl;
			sstream_out << "recv R_cur "    << R_cur    << std::endl;*/

			if (residual < res_min) {
				res_min = residual;
				cb_min = cb_cur;
				rhob_min = rhob_cur;
				R_min = R_cur;
				cws_min = cws_cur;
				sstream_out << endl << "New residual minimum:" << endl;
				sstream_out << "err=" << res_min << ", parameters:" << endl;
				sstream_out << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << endl;\
				sstream_out << "cws_min :" << endl;
				for ( auto &x : cws_min )
					sstream_out << x << " ";
				sstream_out << endl;
				sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
			}
			// if free tasks for sending
			if ( send_task_count < cws_all_cartesians.size() ) {
				for (int j=0; j < task_array_len; j++)
					task_array[j] = cws_all_cartesians[send_task_count][j];
				MPI_Send( task_array, task_array_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
				send_task_count++;
				sstream_out << "send_task_count " << send_task_count << std::endl;
			}
			
			ofile.open("mpi_out", std::ios_base::app);
			ofile << sstream_out.rdbuf();
			sstream_out.clear(); sstream_out.str("");
			ofile.close(); ofile.clear();
		}
		
		// send stop-messages
		for (int j=0; j < task_array_len; j++)
			task_array[j] = -1;
		for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++)
			MPI_Send(task_array, task_array_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		
		sstream_out << endl << "SEARCH ENDED!" << endl;
		sstream_out << "RESULTING VALUE:" << endl;
		sstream_out << "err=" << res_min << ", parameters:" << endl;
		sstream_out << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << endl;
		sstream_out << "cws_min :" << endl;
		for ( auto &x : cws_min )
			sstream_out << x << " ";
		sstream_out << endl;
		sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << endl;

		ofile.open("mpi_out", std::ios_base::app);
		ofile << sstream_out.rdbuf();
		sstream_out.clear(); sstream_out.str("");
		ofile.close(); ofile.clear();

		MPI_Finalize();
	}
	else if (rank > 0) {
		double local_res_min = res_min;
		double local_cb_min;
		double local_rhob_min;
		double local_R_min;
		vector<double> local_cws_min(n_layers_w, 1500);

		vector<vector<double>> cws_vii; // all variants for every depth
		vector<int> index_arr;
		vector<double> cws_vi;
		vector<vector<double>> cws_all_cartesians;
		for ( int ncpl_cur = 0; ncpl_cur <= ncpl; ncpl_cur++ )
			cws_vi.push_back(cw1 + ncpl_cur*(cw2 - cw1) / ncpl);
		if ( rank == 1 ) {
			cout << "cws_vi" << endl;
			for ( auto &x : cws_vi )
				cout << x << " ";
		}
		for (unsigned i=0; i < n_layers_w; i++)
			cws_vii.push_back(cws_vi);
		while (next_cartesian(cws_vii, index_arr, cws_vi))
			cws_all_cartesians.push_back(cws_vi);

		cout << "cws_all_cartesians.size() " << cws_all_cartesians.size() << endl;

		for(;;) {
			MPI_Recv(task_array, task_array_len, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			// if stop-message then finalize
			if ( ( task_array[0] == -1 ) && ( task_array[2] == -1 ) && ( task_array[3] == -1 ) ) {
				MPI_Finalize();
				break;
			}
			
			for(int j=0; j < task_array_len;j++)
				cws_cur[j] = task_array[j];

			if ( rank == 1 )
				for(int j=0; j < task_array_len;j++)
					cout << "recv cws_cur " << cws_cur[j] << endl;
			
			// inverting for bottom halfspace parameters + sound speed in water!
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

						residual = compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points, R_cur, modal_delays, mode_numbers);

						if (residual < local_res_min) {
							local_res_min  = residual;
							local_cb_min   = cb_cur;
							local_rhob_min = rhob_cur;
							local_R_min    = R_cur;
							local_cws_min  = cws_cur;
						}
					}
			
			// send current local minimum to the control process
			result_array[0] = local_res_min;
			result_array[1] = local_cb_min ;
			result_array[2] = local_rhob_min;
			result_array[3] = local_R_min;
			result_array[4] = local_cws_min[0];
			result_array[5] = local_cws_min[1];
			result_array[6] = local_cws_min[2];
			result_array[7] = local_cws_min[3];
			result_array[8] = local_cws_min[4];

			MPI_Send(result_array, result_array_len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
		}
	}
#endif
	return 0;
}

/*
A routine for computing delay residual.
Arguments:
1) Environment: five arrays of the same length: depth, c1s, c2s, rhos, Ns_points;
(each entry describes one layer as described in the comments to compute_wnumbers_extrap() )

2) Source-receive distance: R -- distance from the source to the receiver;

3) Experimental data: modal delays:
    -- experimental_mode_numbers: number of modes for each frequency in the recorded signal
    -- experimental_delays: experimental_delays[ii][jj] is the delay of jj+1-th mode for the frequency freqs[ii]

The routine computes the (uniform) residual (misfit) of experimental data and the "theoretical" delays for a given environment model.

It should be used as follows: for a set of environment models the residual should be computed. The minimal value of the residual indicates
the most "adequate" model.
*/
double compute_modal_delays_residual_uniform( vector<double> &freqs,
                                     vector<double> &depths,
                                     vector<double> &c1s,
                                     vector<double> &c2s,
                                     vector<double> &rhos,
                                     vector<unsigned> &Ns_points,
                                     double R,
                                     vector<vector<double>> &experimental_delays,
                                     vector<unsigned> &experimental_mode_numbers
                                     )
{
    unsigned rord = 3;
    unsigned flTrappedOnly = 1;
    double deltaf = 0.5;
    double residual = 0;
    unsigned mnumb;
    double mdelay;

    vector<vector<double>> modal_group_velocities;
    vector<unsigned> mode_numbers;

    compute_modal_grop_velocities( freqs, deltaf, depths,c1s, c2s, rhos, Ns_points, flTrappedOnly, rord, modal_group_velocities, mode_numbers );

    for (unsigned ii=0; ii<freqs.size();  ii++) {
		mnumb = min(mode_numbers.at(ii), experimental_mode_numbers.at(ii) );
        for (unsigned jj=0; jj<mnumb;  jj++) {
			if (experimental_delays[ii][jj]>0) {
				mdelay =  R/modal_group_velocities[ii][jj];
                residual = residual + pow(experimental_delays[ii][jj]-mdelay,2);
            }
        }
    }

	residual = sqrt(residual);

	return residual;
}


int compute_modal_grop_velocities(      vector<double> &freqs,
                                        double deltaf,
                                        vector<double> &depths,
                                        vector<double> &c1s,
                                        vector<double> &c2s,
                                        vector<double> &rhos,
                                        vector<unsigned> &Ns_points,
                                        unsigned flOnlyTrapped,
                                        unsigned &ordRich,
                                        vector<vector<double>> &modal_group_velocities,
                                        vector<unsigned> &mode_numbers
                                        )
{
    mode_numbers.clear();
    modal_group_velocities.clear();

    vector<double> out_wnum1;
    vector<double> out_wnum2;
    vector<double> mgv_ii;
    unsigned nwnum;
	unsigned nfr = (unsigned)freqs.size();
    double omeg1, omeg2;

    for (unsigned ii=0; ii<nfr; ii++) {
        out_wnum1.clear();
        out_wnum2.clear();
        mgv_ii.clear();
		omeg1 = 2*LOCAL_M_PI*(freqs.at(ii) + deltaf / 2);
        out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1,depths,c1s,c2s,rhos,Ns_points,1,ordRich);
		nwnum = (unsigned)out_wnum1.size();

        /*
        cout << "f=" << freqs.at(ii) << "Hz" << endl;

        for (unsigned jj=0; jj < nwnum; jj++ )
        {
            cout << "k_" << jj+1 << "=" << out_wnum1.at(jj) << endl;
        }
        */

		omeg2 = 2*LOCAL_M_PI*(freqs.at(ii) - deltaf / 2);
        out_wnum2 = compute_wnumbers_extrap_lin_dz(omeg2,depths,c1s,c2s,rhos,Ns_points,1,ordRich);
        nwnum = min(nwnum, (unsigned)out_wnum2.size());

        for (unsigned jj=0; jj < nwnum; jj++ )
        {
            mgv_ii.push_back(  (omeg1 - omeg2)/(out_wnum1.at(jj) - out_wnum2.at(jj) )  );
        }

        modal_group_velocities.push_back( mgv_ii );
        mode_numbers.push_back( nwnum );
    }

	return 0;
}


/*
General considerations:
1) It is better to choose mesh in such a way that mesh size is the same for all z. This gives better accuracy!
(the reason: if the diagonal elements in the upper and lower diagonals vary with row index j, then accuracy is crippled, the
 closer are u_j to each other -- the better, ideally they should be equal)
2) compute_wnumbers_extrap_lin_dz() makes the number of points within each layer a multiple of 12;
it is better to set all the numbers to multiple of 12 in advance
3) Richardson extrapolation of the order 3 gives reasonable accuracy
*/
vector<double> compute_wnumbers_extrap(       double &omeg, // sound frequency
											  vector<double> &depths,
											  vector<double> &c1s,
											  vector<double> &c2s,
											  vector<double> &rhos,
											  vector<unsigned> &Ns_points,
											  unsigned flOnlyTrapped,
											  unsigned &ordRich)
/*  subroutine for computing wavenumbers for a given waveguide structure
    the computation is performed by the FD method for certain meshsize,
    Richardson extrapolation is used to improve the

Layer structure:

depths_{i-1}----c=c1s_i--
*
*
... <-Ns_points_i       rho = rhos_i
*
*
depths_{i}------c=c2s_i--

Other parameters:
omeg -- cyclic frequency, omeg = 2*Pi*f;
flOnlyTrapped -- flag to determine the mode subset: set to 0 to compute all propagating modes, i.e. such that k^2>=0, otherwise only trapped modes are computed
ordRich -- order of the Richardson extrapolation;

the top of the first layer is z=0
*/
{
    vector<double> coeff_extrap;
	switch (ordRich) {
	case 1 :
		coeff_extrap.assign({ 1 });
		break;
	case 2 :
		coeff_extrap.assign({ -1, 2 });
		break;
	case 4 :
		coeff_extrap.assign({ -1 / double(6), 4, -13.5, 32 / double(3) });
	default :
		ordRich = 3;
		coeff_extrap.assign({ 0.5, -4, 4.5 });
		//coeff_extrap.assign({0.1, -0.6, 1.5});
		break;
	}
//
//    cout << "Richardson coeffs" << endl;
//    for (unsigned ii=0; ii<coeff_extrap.size() ; ii++ ){
//        cout << coeff_extrap.at(ii) << endl;
//    }

    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh;
    vector<unsigned> input_interf_idcs;
    vector<double> out_wnum2;
    vector<double> wnum_extrapR;
    double zc = 0;
    double zp = 0;
    double dz = 0;
    unsigned m_wnum = 1000;

	unsigned n_layers = (unsigned)depths.size();
    unsigned n_points_total = 0;
    unsigned n_points_layer = 0;

// outer loop for Richardson coefficient rr
    for (unsigned rr = 1; rr <= ordRich; rr++){
        input_c.clear();
        input_rho.clear();
        input_interf_idcs.clear();
        input_mesh.clear();
        out_wnum2.clear();

        input_c.push_back(0);
        input_rho.push_back(0);
        n_points_total = 1;
        zp = 0;

        for (unsigned ll = 0; ll<n_layers; ll++){
            zc = depths.at(ll);
            n_points_layer = Ns_points.at(ll)*rr;
            dz = (zc - zp)/(  n_points_layer  );
            input_mesh.push_back(  dz  );
            input_c.at(n_points_total-1) = c1s.at(ll) ;
            input_rho.at(n_points_total-1) = rhos.at(ll) ;

            n_points_total = n_points_total + n_points_layer;

            for (unsigned kk=1; kk<= n_points_layer; kk++) {
                input_rho.push_back(rhos.at(ll));
                input_c.push_back( c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk/n_points_layer );
            }

            if (ll < n_layers - 1) {
                input_interf_idcs.push_back(n_points_total-1);
            }
            zp = zc;
        }

        //cout << "rr=" << rr << endl;

        out_wnum2 = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,flOnlyTrapped);
		m_wnum = min(m_wnum, (unsigned)out_wnum2.size());

        if (rr == 1) { wnum_extrapR.assign(m_wnum,0);}

        for (unsigned mm=0; mm<m_wnum; mm++ ) {
            wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr-1);
        }

        /*
            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
            }
        */
    }

	for (unsigned mm=0; mm<m_wnum; mm++ ) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}

	return wnum_extrapR;
}



vector<double> compute_wnumbers_extrap_lin_dz( double &omeg, // sound frequency
											   vector<double> &depths,
											   vector<double> &c1s,
											   vector<double> &c2s,
											   vector<double> &rhos,
											   vector<unsigned> &Ns_points,
											   unsigned flOnlyTrapped,
											   unsigned &ordRich )
/*  subroutine for computing wavenumbers for a given waveguide structure
    the computation is performed by the FD method for certain meshsize,
    Richardson extrapolation is used to improve the

Layer structure:

depths_{i-1}----c=c1s_i--
*
*
... <-Ns_points_i       rho = rhos_i
*
*
depths_{i}------c=c2s_i--

Other parameters:
omeg -- cyclic frequency, omeg = 2*Pi*f;
flOnlyTrapped -- flag to determine the mode subset: set to 0 to compute all propagating modes, i.e. such that k^2>=0, otherwise only trapped modes are computed
ordRich -- order of the Richardson extrapolation;

the top of the first layer is z=0
	*/
{
    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh;
    vector<unsigned> input_interf_idcs;
    vector<double> out_wnum2;
    vector<double> wnum_extrapR;
    double zc = 0;
    double zp = 0;
    double dz = 0;
    unsigned m_wnum = 1000;

	unsigned n_layers = (unsigned)depths.size();
    unsigned n_points_total = 0;
    unsigned n_points_layer = 0;

    vector<double> coeff_extrap;
	switch (ordRich) {
	case 1 :
		coeff_extrap.assign({1});
		break;
	case 2 :
		coeff_extrap.assign({1.333333333333333, -0.333333333333333});
		break;
	case 4 :
		// coeff_extrap.assign({1.595325630252102, -0.788449059052564, 0.216346153846154, -0.023222725045691});
		coeff_extrap.assign({1.6, -0.8, 0.228571428571429, -0.028571428571429});
		break;
	default:
		ordRich = 3;
		coeff_extrap.assign({1.5, -0.6, 0.1});
		//coeff_extrap.assign({0.1, -0.6, 1.5});
		break;
	}

    // number of points in each layer is multiple of 12
    // this allows us to use nz_ii = nz/ii, ii = 1,2,3,4
    for (unsigned ii=0; ii < n_layers; ii++ ){
        Ns_points.at(ii) = 12*(Ns_points.at(ii)/12);
    }

//    cout << "Richardson coeffs" << endl;
//    for (int ii=0; ii<coeff_extrap.size() ; ii++ ){
//        cout << coeff_extrap.at(ii) << endl;
//    }

// outer loop for Richardson coefficient rr
    for (unsigned rr = 1; rr <= ordRich; rr++){
        input_c.clear();
        input_rho.clear();
        input_interf_idcs.clear();
        input_mesh.clear();
        out_wnum2.clear();

        input_c.push_back(0);
        input_rho.push_back(0);
        n_points_total = 1;
        zp = 0;

        for (unsigned ll = 0; ll<n_layers; ll++){
            zc = depths.at(ll);
            n_points_layer = Ns_points.at(ll)/rr;
            dz = (zc - zp)/(  n_points_layer  );

//            cout << "np=" << n_points_layer << "  " << "dz=" << dz <<endl;

            input_mesh.push_back(  dz  );
            input_c.at(n_points_total-1) = c1s.at(ll) ;
            input_rho.at(n_points_total-1) = rhos.at(ll) ;

            n_points_total = n_points_total + n_points_layer;

            for (unsigned kk=1; kk<= n_points_layer; kk++) {
                input_rho.push_back(rhos.at(ll));
                input_c.push_back( c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk/n_points_layer );
            }

            if (ll < n_layers - 1) {
                input_interf_idcs.push_back(n_points_total-1);
            }
            zp = zc;
        }

//        cout << "rr=" << rr << endl;

        out_wnum2 = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,flOnlyTrapped);
		m_wnum = min(m_wnum, (unsigned)out_wnum2.size());

        if (rr == 1) { wnum_extrapR.assign(m_wnum,0);}

        for (unsigned mm=0; mm<m_wnum; mm++ ) {
            wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr-1);
        }

//            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
//
//                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
//            }
//            cout << endl;
    }

	for (unsigned mm=0; mm<m_wnum; mm++ ) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}
	return wnum_extrapR;
}


vector<double> compute_wnumbers( double &omeg, // sound frequency
								 vector<double> &c,
								 vector<double> &rho,
								 vector<unsigned> &interface_idcs,
								 vector<double> &meshsizes,
								 unsigned flOnlyTrapped                 // set flOnlyTrapped = 0 to compute all propagating modes, i.e. such that k^2>=0
											  )
{
	// prepare the three diagonals of the matrix to be passed to the EIG function
    // for the c = c_j, j=0... N_points
    // interfaces are at z = z_m,  interface_idcs = {m}, if empty then we have NO interfaces
    // mesh size in the j-th layer is meshsizes.at(j); this vector has at least one element,
    // for the k-th interface interface_idcs.at(k-1) we have meshsizes meshsizes.at(k-1) and meshsizes.at(k) below and above the interface respectively
    // for c(interface_idcs.at(k-1)) the value of c is the one BELOW the k-th interface
    //(i.e. for the water-bottom interface at the boundary we take the value from the bottom)

	vector<double> md;
	vector<double> ud;
	vector<double> ld;
    vector<double> wnumbers2;

	int N_points = (unsigned)c.size();
    unsigned layer_number = 0;
    double dz = meshsizes.at(layer_number);
    double dz_next = dz;//    ofstream myFile("thematrixdiags.txt");
//    for (int ii=0; ii<N_points-2; ii++){
//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
//    }
//    myFile.close();
    double q = 0;
    double cp, cm, dp, dm, cmin, cmax, kappamax, kappamin;
    int next_interface_idx;

    if ( interface_idcs.size() > 0 )
    {
        next_interface_idx = interface_idcs.at(0)-1;
    }
    else
    {
        next_interface_idx = N_points;
    }

    cmin = c.at(0);
    cmax = c.at(0);

    for( int ii=0; ii < N_points; ii++ ){
        if (c.at(ii) < cmin) { cmin = c.at(ii); }
        if (c.at(ii) > cmax) { cmax = c.at(ii); }
    }

    kappamax = omeg/cmin;
    kappamin = omeg/cmax;

    if (flOnlyTrapped == 0 )
        kappamin = 0;

    for(int ii=0; ii < N_points-2; ii++ ){
        // ordinary point
        ud.push_back( 1/(dz*dz) );
        ld.push_back( 1/(dz*dz) );
        md.push_back(-2/(dz*dz)  + omeg*omeg/(c.at(ii+1)*c.at(ii+1)) );

        // special case of the point at the interface

        if (ii == next_interface_idx) {         //ii -- z(ii+1), z(0) = 0
            layer_number = layer_number + 1;    // вообще ii=89 -- вода, в ii=90 -дно,
                                                // здесь ii = 89 -- интерфейс, уже дно
            cp = c.at(ii+1);
            dp = rho.at(ii+1);
            cm = c.at(ii);
            dm = rho.at(ii);
            q = 1/( dz_next*dm + dz*dp );

            dz_next = meshsizes.at(layer_number);

            ld.at(ii) = 2*q*dp/dz;
            md.at(ii) = -2*q*( dz_next*dp + dz*dm )/(dz*dz_next) + omeg*omeg*q*( dz*dp*cp*cp + dz_next*dm*cm*cm )/( cp*cp*cm*cm ) ;
            ud.at(ii) = 2*q*dm/dz_next;

            if ( interface_idcs.size() > layer_number )
            {
                next_interface_idx = interface_idcs.at(layer_number) - 1;
            }
            else
            {
                next_interface_idx = N_points;
            }

            dz = dz_next;
        }
    }

    // HERE WE CALL THE EIG ROUTINE!!!
    //input: diagonals ld, md, ud + interval [0 k_max]
    //output: wnumbers2 = wave numbers squared

    alglib::real_2d_array eigenvectors; // V - собств вектор
	alglib::real_1d_array eigenvalues; // Lm -собств знач
	eigenvalues.setlength(N_points-2);
    alglib::ae_int_t eigen_count = 0;

    alglib::real_1d_array main_diag, second_diag;
    main_diag.setlength(N_points-2);
    second_diag.setlength(N_points-3);

    for ( int ii=0; ii < N_points-3; ii++ ) {
        second_diag[ii] = sqrt(ud.at(ii)*ld.at(ii+1));
        main_diag[ii] = md.at(ii);
	}
    main_diag[N_points-3] = md.at(N_points-3);

//    ofstream myFile("thematrixdiags.txt");
//    for (int ii=0; ii<N_points-2; ii++){
//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
//    }
//    myFile.close();

//    ofstream myFile1("thematrixdiags_sym.txt");
//    for (int ii=0; ii<N_points-3; ii++){
//        myFile1 << std::fixed << std::setprecision(16) << main_diag[ii] << "  " << second_diag[ii] << endl;
//    }
//    myFile1.close();

    alglib::smatrixtdevdr(main_diag,  second_diag,  N_points-2,  0, kappamin*kappamin,  kappamax*kappamax, eigen_count,eigenvectors);

    for (int ii=0; ii<eigen_count; ii++) {
        wnumbers2.push_back( main_diag[eigen_count-ii-1] );
    }

    return wnumbers2;
}
