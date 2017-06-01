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
	corecount ( 0 )
{}

sspemdd_parallel::~sspemdd_parallel()
{}

void sspemdd_parallel::MPI_main()
{
	if (rank == 0)
		control_process();
	else if (rank > 0)
		computing_process();
}

void sspemdd_parallel::control_process()
{
#ifdef _MPI
	MPI_Status status;
	std::stringstream sstream_out;
	mpi_start_time = MPI_Wtime();
	
	std::cout << "control_process() started" << std::endl;

	sstream_out << "MPI control process" << std::endl;
	sstream_out << "Start residual is: " << record_point.residual << std::endl;
	sstream_out << "Search space:" << std::endl;
	sstream_out << cb1 << " < c_b < " << cb2 << std::endl;
	sstream_out << R1 << " < Range < " << R2 << std::endl;
	sstream_out << rhob1 << " < rho_b < " << rhob2 << std::endl;
	sstream_out << tau1 << " < tau < " << tau2 << std::endl;
	sstream_out << "cw1_arr :" << std::endl;
	for (auto &x : cw1_arr)
		sstream_out << x << " ";
	sstream_out << std::endl;
	sstream_out << "cw2_arr :" << std::endl;
	for (auto &x : cw2_arr)
		sstream_out << x << " ";
	sstream_out << std::endl;
	sstream_out << "ncpl_arr :" << std::endl;
	for (auto &x : ncpl_arr)
		sstream_out << x << " ";
	sstream_out << std::endl;
	
	std::vector<std::vector<double>> point_values_vec;
	point_values_vec.resize(N_total);
	sstream_out << "point_values_vec as vector of tasks" << std::endl;
	sstream_out << "point_values_vec.size() " << point_values_vec.size() << std::endl;
	unsigned long long index = 0;
	std::vector<double> cur_point_values;
	std::vector<int> index_arr;
	while (SSPEMDD_utils::next_cartesian(search_space, index_arr, cur_point_values))
		point_values_vec[index++] = cur_point_values;
	std::cout << "next_cartesian() finished" << std::endl;

	sstream_out << "point_values_vec[0].size() " << point_values_vec[0].size() << std::endl;
	task_len = point_values_vec[0].size() + 1; // point data + task index
	task = new double[task_len];
	result_len = 2; // calculated residual + index of a task
	result = new double[result_len];
	sstream_out << "task_len " << task_len << std::endl;
	
	unsigned send_task_count = 0;
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		//sstream_out << "before filling task" << std::endl;
		//std::cout << sstream_out.str();
		for (unsigned j = 0; j < task_len - 1; j++)
			task[j] = point_values_vec[send_task_count][j];
		task[task_len - 1] = (double)send_task_count;
		//sstream_out << "sending task" << std::endl;
		//for (unsigned j = 0; j < task_len; j++)
		//	sstream_out << task[j] << " ";
		//sstream_out << std::endl;
		//std::cout << sstream_out.str();
		MPI_Send(task, task_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << std::endl;
	std::cout << "first send_task_count " << send_task_count << std::endl;

	std::ofstream ofile("mpi_out");
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	unsigned processed_task_count = 0;
	unsigned received_task_index;
	double received_residual;
	
	// get results and send new tasks on idle computing processes
	while (processed_task_count < point_values_vec.size()) {
		MPI_Recv(result, result_len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		if (processed_task_count % 10000 == 0)
			sstream_out << std::endl << "processed_task_count " << processed_task_count << std::endl;
		
		received_residual = result[0];
		received_task_index = (unsigned)result[1];
		//sstream_out << "received_residual " << received_residual << std::endl;
		//sstream_out << "received_task_index  " << received_task_index << std::endl;
		
		if (received_residual < record_point.residual) {
			record_point = fromDoubleVecToPoint(point_values_vec[received_task_index]);
			record_point.residual = received_residual;
			
			std::cout << "Control process, new residual minimum : "  << received_residual << std::endl;
			sstream_out << std::endl << "Control process, new residual minimum:" << std::endl;
			sstream_out << "err = " << record_point.residual << ", parameters:" << std::endl;
			sstream_out << "c_b = " << record_point.cb << 
				           ", rho_b= " << record_point.rhob << 
						   ", tau = " << record_point.tau <<
				           ", R = " << record_point.R << std::endl;
			sstream_out << "cws_min :" << std::endl;
			for (auto &x : record_point.cws)
				sstream_out << x << " ";
			sstream_out << std::endl;
			sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << std::endl;
		}
		// if free tasks for sending
		if (send_task_count < point_values_vec.size()) {
			for (unsigned j = 0; j < task_len - 1; j++)
				task[j] = point_values_vec[send_task_count][j];
			task[task_len - 1] = send_task_count;
			MPI_Send(task, task_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			send_task_count++;
		}
		else {
			// send stop-messages
			for (unsigned j = 0; j < task_len; j++)
				task[j] = STOP_MESSAGE;
			MPI_Send(task, task_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		}
		
		ofile.open("mpi_out", std::ios_base::app);
		ofile << sstream_out.rdbuf();
		sstream_out.clear(); sstream_out.str("");
		ofile.close(); ofile.clear();
	}
	
	sstream_out << std::endl << "SEARCH ENDED!" << std::endl;
	sstream_out << "err = " << record_point.residual << ", parameters:" << std::endl;
	sstream_out << "c_b = " << record_point.cb <<
		", rho_b = " << record_point.rhob <<
		", R = " << record_point.R <<
		", tau = " << record_point.tau << std::endl;
	sstream_out << "cws :" << std::endl;
	for (auto &x : record_point.cws)
		sstream_out << x << " ";
	sstream_out << std::endl;
	sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << std::endl;
	
	ofile.open("mpi_out", std::ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	delete[] task;
	delete[] result;
	
	std::cout << "finilizing process " << rank << std::endl;
	MPI_Finalize();
#endif
}

void sspemdd_parallel::computing_process()
{
#ifdef _MPI
	MPI_Status status;
	task_len = search_space.size() + 1;
	task = new double[task_len];
	result_len = 2; // index of a task + calculated residual
	result = new double[result_len];
	std::vector<double> cur_point_values_vec;
	search_space_point cur_point;
	cur_point_values_vec.resize(task_len-1);
	double task_index;
	
	if (rank == 1) {
		std::cout << "computing process, rank " << rank << std::endl;
		std::cout << "task_len " << task_len << std::endl;
		std::cout << "result_len " << result_len << std::endl;
	}

	std::stringstream cur_process_points_sstream;

	cur_process_points_sstream << "cb rhob R tau ";
	for (unsigned i = 0; i < cw1_arr.size(); i++) {
		cur_process_points_sstream << "cw" << i << " ";
	}
	cur_process_points_sstream << "residual";
	cur_process_points_sstream << std::endl;
	
	for (;;) {
		MPI_Recv(task, task_len, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (rank == 1) {
			std::cout << "received task " << std::endl;
			for ( unsigned i=0; i < task_len; i++)
				std::cout << task[i] << " " << std::endl;
			std::cout << std::endl;
		}
		
		// if stop-message then finalize
		if (task[task_len-1] == -1) {
			std::cout << "rank " << rank << " received stop-message" << std::endl;
			break;
		}
		//std::cout << "1" << std::endl;
		for ( unsigned i=0; i < task_len - 1; i++ )
			cur_point_values_vec[i] = task[i];
		//std::cout << "2" << std::endl;
		task_index = task[task_len - 1];
		//std::cout << "3" << std::endl;
		cur_point = fromDoubleVecToPoint(cur_point_values_vec);
		//std::cout << "4" << std::endl;
		
		if (rank == 1)
			std::cout << "received task_index " << task_index << std::endl;
		
		fillDataComputeResidual(cur_point); // calculated residual is written to cur_point

		cur_process_points_sstream << cur_point.cb << " " 
			                       << cur_point.rhob << " "
								   << cur_point.R << " "
								   << cur_point.tau << " ";
		for (unsigned i = 0; i < cur_point.cws.size(); i++)
			cur_process_points_sstream << cur_point.cws[i] << " ";
		cur_process_points_sstream << cur_point.residual;
		cur_process_points_sstream << std::endl;
		
		result[0] = cur_point.residual;
		result[1] = task_index;
		MPI_Send(result, result_len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	delete[] task;
	delete[] result;

	std::stringstream cur_process_file_name_sstream;
	cur_process_file_name_sstream << "points_process" << rank;
	std::ofstream cur_process_file(cur_process_file_name_sstream.str().c_str());
	cur_process_file << cur_process_points_sstream.rdbuf();
	cur_process_file.close();

	MPI_Finalize();

#endif
}