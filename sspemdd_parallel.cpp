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
	MPI_Status status, cur_status;
	stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	cout << endl << "control_process() started" << endl;

	vector<vector<double>> depths_vec;
	createDepthsArray(depths_vec);
	cout << "depths_vec.size() " << depths_vec.size() << endl;
	
	double *task = new double[TASK_LEN];
	double *result = new double[RESULT_LEN];
	
	int send_task_count = 0;
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		sendTask(task, send_task_count, computing_process_index, depths_vec[send_task_count]);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << endl;

	std::ofstream ofile("mpi_out");
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	unsigned processed_task_count = 0;
	int received_task_index;
	double received_residual;
	double task_processing_time;
	int stop_message = -1;
	
	// get results and send new tasks on idle computing processes
	while (processed_task_count < depths_vec.size()) {
		MPI_Recv( &received_task_index,  1, MPI_INT, MPI_ANY_SOURCE,    MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		cur_status = status;
		MPI_Recv( &task_processing_time, 1, MPI_DOUBLE, cur_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		MPI_Recv( result, RESULT_LEN, MPI_DOUBLE, cur_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		processed_task_count++;
		if (processed_task_count % 100 == 0) {
			sstream_out << endl << "processed_task_count " << processed_task_count;
			sstream_out << " , time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
		}
		
		received_residual = result[0];
		
		if (received_residual < record_point.residual) {
			record_point.residual = received_residual;
			record_point.cb   = result[1];
			record_point.rhob = result[2];
			record_point.R    = result[3];
			record_point.tau  = result[4];
			unsigned cur_depths_size = depths_vec[received_task_index].size();
			record_point.depths.resize(cur_depths_size);
			unsigned cur_cws_size = cur_depths_size - 1;
			record_point.cws.resize(cur_cws_size);
			for (unsigned j = 0; j < cur_cws_size; j++)
				record_point.cws[j] = result[5 + j];
			for (unsigned j = 0; j < cur_depths_size; j++)
				record_point.depths[j] = result[5 + cur_cws_size + j];
			
			if (record_point.depths != depths_vec[received_task_index]) {
				cerr << "record_point.depths != depths_vec[received_task_index]" << endl;
				cerr << "received_task_index " << received_task_index << endl;
				MPI_Abort(MPI_COMM_WORLD, 0);
				return;
			}
			
			//cout << "Control process, new residual minimum : "  << received_residual << endl;
			sstream_out << endl << "Control process, new residual minimum:" << endl;
			sstream_out << "task processing time " << task_processing_time << " s" << endl;
			sstream_out << "err = " << record_point.residual << ", parameters:" << endl;
			sstream_out << "c_b = " << record_point.cb << 
				           ", rho_b= " << record_point.rhob << 
						   ", tau = " << record_point.tau <<
				           ", R = " << record_point.R << endl;
			sstream_out << "cws :" << endl;
			for (auto &x : record_point.cws)
				sstream_out << x << " ";
			sstream_out << endl;
			sstream_out << "depths :" << endl;
			for (auto &x : record_point.depths)
				sstream_out << x << " ";
			sstream_out << endl;
			sstream_out << "time from start " << MPI_Wtime() - mpi_start_time << " s" << endl;
		}
		// if free tasks for sending
		if (send_task_count < depths_vec.size()) {
			sendTask(task, send_task_count, status.MPI_SOURCE, depths_vec[send_task_count]);
			send_task_count++;
		}
		else
			MPI_Send( &stop_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		
		ofile.open("mpi_out", std::ios_base::app);
		ofile << sstream_out.rdbuf();
		sstream_out.clear(); sstream_out.str("");
		ofile.close(); ofile.clear();
	}
	
	sstream_out << endl << "SEARCH ENDED!" << endl;
	sstream_out << "err = " << record_point.residual << ", parameters:" << endl;
	sstream_out << "c_b = " << record_point.cb <<
		", rho_b = " << record_point.rhob <<
		", R = " << record_point.R <<
		", tau = " << record_point.tau << endl;
	sstream_out << "cws :" << endl;
	for (auto &x : record_point.cws)
		sstream_out << x << " ";
	sstream_out << endl;
	sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << endl;
	
	ofile.open("mpi_out", std::ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	delete[] task;
	delete[] result;
	
	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}

void sspemdd_parallel::sendTask(double *task, int task_index, unsigned process_index, vector<double> depths)
{
	unsigned cur_depths_size = depths.size();
	for (unsigned j = 0; j < cur_depths_size; j++)
		task[j] = depths[j];
	for (unsigned j = cur_depths_size; j < TASK_LEN; j++)
		task[j] = -1;
#ifdef _MPI
	MPI_Send( &task_index, 1,        MPI_INT, process_index, 0, MPI_COMM_WORLD);
	MPI_Send( task,        TASK_LEN, MPI_DOUBLE,   process_index, 0, MPI_COMM_WORLD);
#endif
}

void sspemdd_parallel::computing_process()
{
#ifdef _MPI
	MPI_Status status;
	double *task = new double[TASK_LEN];
	double *result = new double[RESULT_LEN];
	int task_index;
	vector<double> depths;

	stringstream cur_process_points_sstream;
	search_space_point local_point_record;
	
	for (;;) {
		MPI_Recv( &task_index, 1,        MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		// if stop-message then finalize
		if (task_index == -1) {
			cout << "rank " << rank << " received stop-message" << endl;
			break;
		}
		MPI_Recv( task,        TASK_LEN, MPI_DOUBLE,   0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

		if (rank == 1) {
			cout << "computing process, rank " << rank << endl;
			cout << "received task " << endl;
			for ( unsigned j = 0; j < TASK_LEN; j++)
				cout << task[j] << " " << endl;
			cout << endl;
		}
		depths.clear();
		for (unsigned j = 0; j < TASK_LEN; j++) {
			if (task[j] == -1)
				break;
			depths.push_back( task[j] );
		}
		
		if (rank == 1) {
			cout << "received depths " << endl;
			for (auto &x : depths)
				cout << x << " " << endl;
			cout << endl;
			cout << "received task_index " << task_index << endl;
		}

		init(depths);
		double processing_time = MPI_Wtime();
		local_point_record = findLocalMinHillClimbing(depths);
		processing_time = MPI_Wtime() - processing_time;
		if ( (rank == 1) && (verbosity > 0) )
			cout << "task " << task_index << " completed" << endl;
		
		result[0] = local_point_record.residual;
		result[1] = local_point_record.cb;
		result[2] = local_point_record.rhob;
		result[3] = local_point_record.R;
		result[4] = local_point_record.tau;
		unsigned cur_cws_size = local_point_record.cws.size();
		for (unsigned j = 0; j < cur_cws_size; j++)
			result[5 + j] = local_point_record.cws[j];
		unsigned cur_depths_size = depths.size();
		for (unsigned j = 0; j < cur_depths_size; j++)
			result[5 + cur_cws_size + j] = local_point_record.depths[j];
		for (unsigned j = 5 + cur_cws_size + cur_depths_size; j < RESULT_LEN; j++)
			result[j] = -1;
		
		MPI_Send( &task_index,      1,          MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD );
		MPI_Send( &processing_time, 1,          MPI_DOUBLE,   0, 0, MPI_COMM_WORLD );
		MPI_Send( result,          RESULT_LEN, MPI_DOUBLE,   0, 0, MPI_COMM_WORLD );
		
		/*cur_process_points_sstream << cur_point.cb << " " 
			                       << cur_point.rhob << " "
								   << cur_point.R << " "
								   << cur_point.tau << " ";
		for (unsigned i = 0; i < cur_point.cws.size(); i++)
			cur_process_points_sstream << cur_point.cws[i] << " ";
		cur_process_points_sstream << cur_point.residual;
		cur_process_points_sstream << endl;*/
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