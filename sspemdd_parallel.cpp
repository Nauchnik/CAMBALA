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
	task_len ( 0 ),
	result_len ( 0 ),
	corecount ( 0 ),
	rank ( 0 )
{}

sspemdd_parallel::~sspemdd_parallel()
{}

void sspemdd_parallel::MPI_main()
{
	allocateArrays();

	if (rank == 0)
		control_process();
	else if (rank > 0)
		computing_process();
	
	deallocateArrays();
}

void sspemdd_parallel::control_process()
{
#ifdef _MPI
	MPI_Status status;
	std::stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	sstream_out << "MPI control process" << std::endl;
	sstream_out << "Start residual is: " << residual << std::endl;
	sstream_out << "Search space:" << std::endl;
	sstream_out << cb1 << " < c_b < " << cb2 << std::endl;
	sstream_out << R1 << " < Range < " << R2 << std::endl;
	sstream_out << cw1 << " < cws < " << cw2 << std::endl;
	sstream_out << rhob1 << "< rho_b < " << rhob2 << std::endl;
	sstream_out << tau1 << "< tau < " << tau2 << std::endl;
	sstream_out << task_array_len << "< task_array_len < " << task_array_len << std::endl;
	
	std::vector<int> index_arr;
	std::vector<std::vector<unsigned>> search_space_indexes;
	std::vector<unsigned> cur_point_indexes;
	search_space_indexes.resize(search_space.size());
	for (unsigned variable_index = 0; variable_index < search_space.size(); variable_index++)
		for (unsigned j = 0; j < search_space[variable_index].size(); j++)
			search_space_indexes[variable_index].push_back(j);
	std::vector<std::vector<unsigned>> tasks_vec;
	
	tasks_vec.resize(N_total);
	unsigned long long task_index = 0;
	search_space_point cur_point;
	while (SSPEMDD_utils::next_cartesian(search_space_indexes, index_arr, cur_point_indexes)) {
		cur_point = fromPointIndexesToPoint(cur_point_indexes);
		tasks_vec[task_index][0] = ;
		task_index++;
	}
	sstream_out << "tasks_vec.size() " << tasks_vec.size() << std::endl;
	
	unsigned send_task_count = 0;
	unsigned processed_task_count = 0;
	sstream_out << "task_array_len " << task_array_len << std::endl;
	sstream_out << "tasks_vec[0].size() " << tasks_vec[0].size() << std::endl;
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		for (int j = 0; j < (int)task_array_len; j++)
			task[j] = tasks_vec[send_task_count][j];
		MPI_Send(task_array, task_array_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << std::endl;
	std::ofstream ofile("mpi_out");
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();

	// get results and send new tasks on idle computing processes
	while (processed_task_count < cws_all_cartesians.size()) {
		MPI_Recv(task, task_array_len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		sstream_out << "processed_task_count " << processed_task_count << std::endl;
		
		cur_point.cb   = task[0];
		cur_point.rhob = task[1];
		cur_point.R    = task[2];
		cur_point.tau  = task[3];
		cur_point.cws.clear();
		for (unsigned i = 4; i < search_space.size(); i++)
			cur_point.cws.push_back(task[i]);
		
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
			//sstream_out << "send_task_count " << send_task_count << std::endl;
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
	task_len = record_point.cws.size() + 4;
	task = new double[task_len];
	result_len = task_len + 1;
	result = new double[result_len];
}

void sspemdd_parallel::deallocateArrays()
{
	delete[] task;
	delete[] result;
}