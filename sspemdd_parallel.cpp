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
	stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	cout << endl << "control_process() started" << endl;

	vector<vector<double>> depths_vec;
	createDepthsArray(depths_vec);
	cout << "depths_vec.size() " << depths_vec.size() << endl;
	
	double *task = new double[TASK_LEN];
	double *result = new double[RESULT_LEN];
	
	unsigned send_task_count = 0;
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		//sstream_out << "before filling task" << std::endl;
		//std::cout << sstream_out.str();
		unsigned elements_to_send = depths_vec[send_task_count].size() + 1;
		cout << "elements_to_send " << elements_to_send << endl;
		for (unsigned j = 0; j < elements_to_send - 1; j++)
			task[j] = depths_vec[send_task_count][j];
		task[elements_to_send] = (double)send_task_count;
		for (unsigned j = elements_to_send + 1; j < TASK_LEN; j++)
			task[j] = -1;
		MPI_Send(task, TASK_LEN, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << std::endl;

	std::ofstream ofile("mpi_out");
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	unsigned processed_task_count = 0;
	unsigned received_task_index;
	double received_residual;

	// get results and send new tasks on idle computing processes
	while (processed_task_count < depths_vec.size()) {
		MPI_Recv(result, RESULT_LEN, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		if (processed_task_count % 10 == 0)
			sstream_out << std::endl << "processed_task_count " << processed_task_count << std::endl;
		
		received_task_index = (unsigned)result[0];
		received_residual = result[1];
		//sstream_out << "received_residual " << received_residual << std::endl;
		//sstream_out << "received_task_index  " << received_task_index << std::endl;
		
		if (received_residual < record_point.residual) {
			record_point.residual = received_residual;
			record_point.cb   = result[2];
			record_point.rhob = result[3];
			record_point.R    = result[4];
			record_point.tau  = result[5];
			unsigned cur_depths_number = depths_vec[received_task_index].size();
			record_point.cws.resize(cur_depths_number);
			record_point.depths.resize(cur_depths_number);
			for (unsigned j = 0; j < cur_depths_number; j++)
				record_point.cws[j] = result[6 + j];
			for (unsigned j = 0; j < cur_depths_number; j++)
				record_point.depths[j] = result[6 + cur_depths_number + j];

			if (record_point.depths != depths_vec[received_task_index]) {
				cerr << "record_point.depths != depths_vec[received_task_index]" << endl;
				cerr << "received_task_index " << received_task_index << endl;
				MPI_Abort(MPI_COMM_WORLD, 0);
				return;
			}
			
			std::cout << "Control process, new residual minimum : "  << received_residual << std::endl;
			sstream_out << std::endl << "Control process, new residual minimum:" << std::endl;
			sstream_out << "err = " << record_point.residual << ", parameters:" << std::endl;
			sstream_out << "c_b = " << record_point.cb << 
				           ", rho_b= " << record_point.rhob << 
						   ", tau = " << record_point.tau <<
				           ", R = " << record_point.R << std::endl;
			sstream_out << "cws :" << std::endl;
			for (auto &x : record_point.cws)
				sstream_out << x << " ";
			sstream_out << std::endl;
			sstream_out << "depths :" << std::endl;
			for (auto &x : record_point.depths)
				sstream_out << x << " ";
			sstream_out << std::endl;
			sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << std::endl;
		}
		// if free tasks for sending
		if (send_task_count < depths_vec.size()) {
			unsigned elements_to_send = depths_vec[send_task_count].size() + 1;
			cout << "elements_to_send " << elements_to_send << endl;
			for (unsigned j = 0; j < elements_to_send - 1; j++)
				task[j] = depths_vec[send_task_count][j];
			task[elements_to_send] = (double)send_task_count;
			for (unsigned j = elements_to_send + 1; j < TASK_LEN; j++)
				task[j] = -1;
			MPI_Send(task, TASK_LEN, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			send_task_count++;
		}
		else {
			// send stop-messages
			for (unsigned j = 0; j < TASK_LEN; j++)
				task[j] = STOP_MESSAGE;
			MPI_Send(task, TASK_LEN, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
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
	double *task = new double[TASK_LEN];
	double *result = new double[RESULT_LEN];
	double task_index;
	vector<double> depths;

	stringstream cur_process_points_sstream;
	search_space_point local_point_record;
	
	for (;;) {
		MPI_Recv(task, TASK_LEN, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (rank == 1) {
			cout << "computing process, rank " << rank << endl;
			cout << "received task " << endl;
			for ( unsigned j=0; j < TASK_LEN; j++)
				cout << task[j] << " " << endl;
			cout << endl;
		}
		for (unsigned j = 0; j < TASK_LEN; j++) {
			if (task[j] == -1)
				break;
			depths[j] = task[j];
		}
		task_index = task[depths.size()];

		if (rank == 1) {
			cout << "received depths " << endl;
			for (auto &x : depths)
				cout << x << " " << endl;
			cout << endl;
			cout << "received task_index " << task_index << endl;
		}
		
		// if stop-message then finalize
		if (task_index == -1) {
			cout << "rank " << rank << " received stop-message" << std::endl;
			break;
		}

		init(depths);
		local_point_record = findLocalMinHillClimbing(depths);
		
		result[0] = task_index;
		result[1] = local_point_record.residual;
		result[2] = local_point_record.cb;
		result[3] = local_point_record.rhob;
		result[4] = local_point_record.R;
		result[5] = local_point_record.tau;
		unsigned cur_depths_number = depths.size();
		for (unsigned j=0; j < cur_depths_number; j++) {
			result[6 + j] = local_point_record.cws[j];
			result[6 + cur_depths_number + j] = local_point_record.depths[j];
		}
		
		MPI_Send(result, RESULT_LEN, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

		/*cur_process_points_sstream << cur_point.cb << " " 
			                       << cur_point.rhob << " "
								   << cur_point.R << " "
								   << cur_point.tau << " ";
		for (unsigned i = 0; i < cur_point.cws.size(); i++)
			cur_process_points_sstream << cur_point.cws[i] << " ";
		cur_process_points_sstream << cur_point.residual;
		cur_process_points_sstream << std::endl;*/
	}
	delete[] task;
	delete[] result;

	/*std::stringstream cur_process_file_name_sstream;
	cur_process_file_name_sstream << "points_process" << rank;
	std::ofstream cur_process_file(cur_process_file_name_sstream.str().c_str());
	cur_process_file << cur_process_points_sstream.rdbuf();
	cur_process_file.close();*/

	MPI_Finalize();

#endif
}