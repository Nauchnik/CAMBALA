#include "sspemdd_parallel.h"
#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

sspemdd_parallel::sspemdd_parallel() :
	ncb(0),
	nrhob(0),
	nR(0),
	task_array_len ( 0 ),
	result_array_len ( 0 ),
	corecount ( 0 ),
	rank ( 0 ),
	residual ( 0.0 ),
	cb1 ( 0.0 ),
	cb2 (0.0),
	cw1 (0.0),
	cw2 (0.0),
	R1 (0.0),
	R2 (0.0),
	rhob1 (0.0),
	rhob2 (0.0),
	cb_cur (0.0),
	R_cur (0.0),
	rhob_cur(0.0),
	n_layers_w(0),
	res_min(1e100),
	cb_min(0.0),
	rhob_min(0.0),
	R_min(0.0)
{}

sspemdd_parallel::~sspemdd_parallel()
{}

void sspemdd_parallel::control_process()
{
#ifdef _MPI
	MPI_Status status;
	std::stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	sstream_out << "Start residual is: " << residual << std::endl;

	// brute force minimum search
	sstream_out << "BRUTE FORCE MINIMUM SEARCH" << std::endl;
	sstream_out << "Search space:" << std::endl;
	sstream_out << cb1 << " < c_b < " << cb2 << std::endl;
	sstream_out << R1 << " < Range < " << R2 << std::endl;
	sstream_out << cw1 << " < cws < " << cw2 << std::endl;
	sstream_out << rhob1 << "< rho_b < " << rhob2 << std::endl;
	
	sstream_out << "MPI control process" << std::endl;
	sstream_out << "cws_all_cartesians.size() " << cws_all_cartesians.size() << std::endl;
	sstream_out << "cws_all_cartesians[0]" << std::endl;
	for (unsigned j = 0; j < cws_all_cartesians[0].size(); j++)
		sstream_out << cws_all_cartesians[0][j] << " ";
	sstream_out << std::endl;
	if (cws_all_cartesians.size() < (unsigned)corecount) {
		std::cerr << "cws_all_cartesians.size() < corecount" << std::endl;
		std::cerr << cws_all_cartesians.size() << " < " << corecount << std::endl;
		exit(1);
	}

	unsigned send_task_count = 0;
	unsigned processed_task_count = 0;

	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		for (int j = 0; j < (int)task_array_len; j++)
			task_array[j] = cws_all_cartesians[send_task_count][j];
		MPI_Send(task_array, task_array_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		/*sstream_out << "sending cws as task_array " << std::endl;
		for (int j = 0; j < (int)task_array_len; j++)
			sstream_out << task_array[j] << " ";
		sstream_out << std::endl;*/
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << std::endl;
	std::ofstream ofile("mpi_out");
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	double stop_message = -1;
	cws_cur.resize(n_layers_w);

	// get results and send new tasks on idle computing processes
	while (processed_task_count < cws_all_cartesians.size()) {
		MPI_Recv(result_array, result_array_len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		sstream_out << "processed_task_count " << processed_task_count << std::endl;
		
		residual   = result_array[0];
		cb_cur     = result_array[1];
		rhob_cur   = result_array[2];
		R_cur      = result_array[3];
		cws_cur[0] = result_array[4];
		cws_cur[1] = result_array[5];
		cws_cur[2] = result_array[6];
		cws_cur[3] = result_array[7];
		cws_cur[4] = result_array[8];
		
		/*sstream_out << "recv new result " << std::endl;
		sstream_out << "recv residual "   << residual   << std::endl;
		sstream_out << "recv cb_cur "     << cb_cur     << std::endl;
		sstream_out << "recv rhob_cur "   << rhob_cur   << std::endl;
		sstream_out << "recv R_cur "      << R_cur      << std::endl;
		sstream_out << "recv cws_cur[0] " << cws_cur[0] << std::endl;
		sstream_out << "recv cws_cur[1] " << cws_cur[1] << std::endl;
		sstream_out << "recv cws_cur[2] " << cws_cur[2] << std::endl;
		sstream_out << "recv cws_cur[3] " << cws_cur[3] << std::endl;
		sstream_out << "recv cws_cur[4] " << cws_cur[4] << std::endl;*/
		
		if (residual < res_min) {
			res_min = residual;
			cb_min = cb_cur;
			rhob_min = rhob_cur;
			R_min = R_cur;
			cws_min = cws_cur;
			sstream_out << std::endl << "New residual minimum:" << std::endl;
			sstream_out << "err=" << res_min << ", parameters:" << std::endl;
			sstream_out << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << std::endl; \
				sstream_out << "cws_min :" << std::endl;
			for (auto &x : cws_min)
				sstream_out << x << " ";
			sstream_out << std::endl;
			sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << std::endl;
		}
		// if free tasks for sending
		if (send_task_count < cws_all_cartesians.size()) {
			for (unsigned j = 0; j < task_array_len; j++)
				task_array[j] = cws_all_cartesians[send_task_count][j];
			MPI_Send(task_array, task_array_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			send_task_count++;
			sstream_out << "send_task_count " << send_task_count << std::endl;
		}
		else {
			// send stop-messages
			for (unsigned j = 0; j < task_array_len; j++)
				task_array[j] = -1;
				MPI_Send(task_array, task_array_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}
		
		ofile.open("mpi_out", std::ios_base::app);
		ofile << sstream_out.rdbuf();
		sstream_out.clear(); sstream_out.str("");
		ofile.close(); ofile.clear();
	}
	
	sstream_out << std::endl << "SEARCH ENDED!" << std::endl;
	sstream_out << "RESULTING VALUE:" << std::endl;
	sstream_out << "err=" << res_min << ", parameters:" << std::endl;
	sstream_out << "c_b=" << cb_min << ", rho_b=" << rhob_min << ", R=" << R_min << std::endl;
	sstream_out << "cws_min :" << std::endl;
	for (auto &x : cws_min)
		sstream_out << x << " ";
	sstream_out << std::endl;
	sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << std::endl;

	ofile.open("mpi_out", std::ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();

	MPI_Finalize();
#endif
}

void sspemdd_parallel::computing_process()
{
#ifdef _MPI
	MPI_Status status;
	double local_res_min = res_min;
	double local_cb_min;
	double local_rhob_min;
	double local_R_min;
	std::vector<double> local_cws_min(n_layers_w, 1500);
	sspemdd_sequential sspemdd_seq;
	cws_cur.resize(n_layers_w);
	
	if ( rank == 1 )
		std::cout << "computing process, rank " << rank << std::endl;
	
	for (;;) {
		MPI_Recv(task_array, task_array_len, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		
		// if stop-message then finalize
		if (task_array[0] == -1) {
			std::cout << "rank " << rank << " received stop-message" << std::endl;
			MPI_Finalize();
			break;
		}
		
		for (unsigned j = 0; j < task_array_len; j++)
			cws_cur[j] = task_array[j];
		/*if (rank == 1) {
			for (unsigned j = 0; j < task_array_len; j++)
				std::cout << "recv cws_cur " << cws_cur[j] << std::endl;
		}*/

		// inverting for bottom halfspace parameters + sound speed in water!
		for (unsigned cur_ncb = 0; cur_ncb < ncb; cur_ncb++)
			for (unsigned cur_nrhob = 0; cur_nrhob < nrhob; cur_nrhob++)
				for (unsigned cur_nR = 0; cur_nR < nR; cur_nR++) {
					// specify bottom parameters;
					if (ncb > 1) { cb_cur = cb1 + cur_ncb*(cb2 - cb1) / (ncb - 1); }
					else { cb_cur = cb1; }
					if (nrhob > 1) { rhob_cur = rhob1 + cur_nrhob  *(rhob2 - rhob1) / (nrhob - 1); }
					else { rhob_cur = rhob1; }
					// specify range
					if (nR > 1) { R_cur = R1 + cur_nR*(R2 - R1) / (nR - 1); }
					else { R_cur = R1; }

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

					residual = sspemdd_seq.compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points, R_cur, modal_delays, mode_numbers);

					if (residual < local_res_min) {
						local_res_min = residual;
						local_cb_min = cb_cur;
						local_rhob_min = rhob_cur;
						local_R_min = R_cur;
						local_cws_min = cws_cur;
						/*if (rank == 1) {
							std::cout << "on rank 1 residual < local_res_min : " 
								<< residual << " < " << local_res_min << std::endl;
						}*/
					}
				}

		// send current local minimum to the control process
		result_array[0] = local_res_min;
		result_array[1] = local_cb_min;
		result_array[2] = local_rhob_min;
		result_array[3] = local_R_min;
		result_array[4] = local_cws_min[0];
		result_array[5] = local_cws_min[1];
		result_array[6] = local_cws_min[2];
		result_array[7] = local_cws_min[3];
		result_array[8] = local_cws_min[4];
		/*std::cout << "sending local_res_min " << local_res_min << std::endl;
		std::cout << "sending local_cb_min " << local_cb_min << std::endl;
		std::cout << "sending local_rhob_min " << local_rhob_min << std::endl;
		std::cout << "sending local_R_min " << local_R_min << std::endl;
		std::cout << " sending local_cws_min[0]" << local_cws_min[0] << std::endl;
		std::cout << " sending local_cws_min[1]" << local_cws_min[1] << std::endl;
		std::cout << " sending local_cws_min[2]" << local_cws_min[2] << std::endl;
		std::cout << " sending local_cws_min[3]" << local_cws_min[3] << std::endl;
		std::cout << " sending local_cws_min[4]" << local_cws_min[4] << std::endl;*/
		
		MPI_Send(result_array, result_array_len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
#endif
}

void sspemdd_parallel::allocateArrays()
{
	task_array_len = n_layers_w;
	result_array_len = 9;
	task_array = new double[task_array_len];
	result_array = new double[result_array_len];
}

void sspemdd_parallel::deallocateArrays()
{
	std::cout << "deallocateArrays() start" << std::endl;
	delete[] task_array;
	delete[] result_array;
	std::cout << "deallocateArrays() end" << std::endl;
}