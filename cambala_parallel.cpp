#include "cambala_parallel.h"
#include "cambala_sequential.h"
#include "cambala_utils.h"

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

CAMBALA_parallel::CAMBALA_parallel() :
	corecount ( 0 )
{}

CAMBALA_parallel::~CAMBALA_parallel()
{}

void CAMBALA_parallel::MPI_main()
{
	cout << "corecount " << corecount << endl;
	if (launch_type == "ils") {
		if (rank == 0)
			controlProcessIls();
		else if (rank > 0)
			computingProcessIls();
	}
	else if (launch_type == "bruteforce") {
		if (rank == 0)
			controlProcessBruteforce();
		else if (rank > 0)
			computingProcessBruteforce();
	}
}

void CAMBALA_parallel::controlProcessIls()
{
/*
#ifdef _MPI
	MPI_Status status, cur_status;
	stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	cout << endl << "control_process() started" << endl;

	vector<vector<double>> depths_vec;
	createDepthsArray(depths_vec);
	cout << "depths_vec.size() " << depths_vec.size() << endl;
	
	double *task = new double[ILS_TASK_LEN];
	double *result = new double[ILS_RESULT_LEN];
	
	int send_task_count = 0;
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		sendTaskIls(task, send_task_count, computing_process_index, depths_vec[send_task_count]);
		send_task_count++;
	}
	sstream_out << "send_task_count " << send_task_count << endl;
	sstream_out << "tasks number " << depths_vec.size() << endl;
	
	ofstream ofile(output_filename, ios_base::app);
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
		MPI_Recv( result, ILS_RESULT_LEN, MPI_DOUBLE, cur_status.MPI_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status );
		processed_task_count++;
		if (processed_task_count % 1000 == 0) {
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
			sstream_out << "task index " << received_task_index << endl;
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
			sendTaskIls(task, send_task_count, status.MPI_SOURCE, depths_vec[send_task_count]);
			send_task_count++;
		}
		else
			MPI_Send( &stop_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
		
		ofile.open(output_filename, ios_base::app);
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
	sstream_out << "depths :" << endl;
	for (auto &x : record_point.depths)
		sstream_out << x << " ";
	sstream_out << endl;
	sstream_out << "final time " << MPI_Wtime() - mpi_start_time << " s" << endl;
	
	ofile.open(output_filename, ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	delete[] task;
	delete[] result;
	
	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
*/
}

void CAMBALA_parallel::sendTaskIls(double *task, int task_index, unsigned process_index, vector<double> depths)
{
	unsigned cur_depths_size = depths.size();
	for (unsigned j = 0; j < cur_depths_size; j++)
		task[j] = depths[j];
	for (unsigned j = cur_depths_size; j < ILS_TASK_LEN; j++)
		task[j] = -1;
#ifdef _MPI
	MPI_Send( &task_index, 1,        MPI_INT, process_index, 0, MPI_COMM_WORLD);
	MPI_Send( task, ILS_TASK_LEN, MPI_DOUBLE,   process_index, 0, MPI_COMM_WORLD);
#endif
}

void CAMBALA_parallel::computingProcessIls()
{
/*
#ifdef _MPI
	MPI_Status status;
	double *task = new double[ILS_TASK_LEN];
	double *result = new double[ILS_RESULT_LEN];
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
		MPI_Recv( task, ILS_TASK_LEN, MPI_DOUBLE,   0, MPI_ANY_TAG, MPI_COMM_WORLD, &status );

		if (rank == 1) {
			cout << "computing process, rank " << rank << endl;
			cout << "received task " << endl;
			for ( unsigned j = 0; j < ILS_TASK_LEN; j++)
				cout << task[j] << " " << endl;
			cout << endl;
		}
		depths.clear();
		for (unsigned j = 0; j < ILS_TASK_LEN; j++) {
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
		for (unsigned j = 5 + cur_cws_size + cur_depths_size; j < ILS_RESULT_LEN; j++)
			result[j] = -1;
		
		MPI_Send( &task_index,      1,          MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD );
		MPI_Send( &processing_time, 1,          MPI_DOUBLE,   0, 0, MPI_COMM_WORLD );
		MPI_Send( result,          ILS_RESULT_LEN, MPI_DOUBLE,   0, 0, MPI_COMM_WORLD );
	}
	delete[] task;
	delete[] result;

	MPI_Finalize();

#endif
	*/
}

void CAMBALA_parallel::controlProcessBruteforce()
{
	stringstream sstream_out;
	
	vector<double> h_vec;
	for (unsigned i = 0; i < nh; i++)
		h_vec.push_back(h1 + (nh == 1 ? 0 : i*(h2 - h1) / (nh - 1)));

	sstream_out << "h_vec.size() " << h_vec.size() << endl;
	for (auto &x : h_vec)
		sstream_out << x << " ";
	sstream_out << endl;

	vector<vector<double>> total_depths_vec;
	
	for (auto &h : h_vec) {
		vector<vector<double>> depths_vec;
		createDepthsArray(h, depths_vec);
		for (auto &depths : depths_vec)
			total_depths_vec.push_back(depths);
		if (total_depths_vec.size() >= MAX_DEPTHS_VECTORS) {
			sstream_out << "WARNING. total_depths_vec.size() >= MAX_DEPTHS_VECTORS" << endl;
			sstream_out << "total_depths_vec.size() changed to " << MAX_DEPTHS_VECTORS << endl;
			total_depths_vec.resize(MAX_DEPTHS_VECTORS);
			break;
		}
	}
	sstream_out << "total_depths_vec.size() " << total_depths_vec.size() << endl;
	
	ofstream ofile(output_filename, ios_base::app);;
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
	
	for (unsigned i = 0; i < total_depths_vec.size(); i++) {
		controlProcessFixedDepths(total_depths_vec[i], i);
		ofstream ofile(output_filename, ios_base::app);;
		ofile << endl << "***" << endl << 
			"depths # " << i << " out of " << total_depths_vec.size() << " processed" << endl;
		ofile.close(); ofile.clear();
	}
#ifdef _MPI
	// send stop-messages to all computing processes
	int stop_message = -1;
	for (unsigned i=1; i< corecount; i++)
		MPI_Send(&stop_message, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	cout << "stop-messages were sent" << endl;
	MPI_Finalize();
#endif
}

void CAMBALA_parallel::controlProcessFixedDepths(const vector<double> depths, const unsigned depths_index)
{
#ifdef _MPI
	MPI_Status status;
	stringstream sstream_out;
	mpi_start_time = MPI_Wtime();

	cout << "control_process() started of depths # " << depths_index << endl;
	
	init(depths);

	sstream_out << endl;
	sstream_out << "MPI control process" << endl;
	sstream_out << "Start residual is: " << record_point.residual << endl;
	sstream_out << "Search space:" << endl;
	sstream_out << cb1 << " < c_b < " << cb2 << endl;
	sstream_out << R1 << " < Range < " << R2 << endl;
	sstream_out << rhob1 << " < rho_b < " << rhob2 << endl;
	sstream_out << tau1 << " < tau < " << tau2 << endl;
	sstream_out << "cw1_arr :" << endl;
	for (auto &x : cw1_arr)
		sstream_out << x << " ";
	sstream_out << endl;
	sstream_out << "cw2_arr :" << endl;
	for (auto &x : cw2_arr)
		sstream_out << x << " ";
	sstream_out << endl;
	sstream_out << "ncpl_arr :" << endl;
	for (auto &x : ncpl_arr)
		sstream_out << x << " ";
	sstream_out << endl;

	vector<vector<double>> point_values_vec;
	point_values_vec.resize(N_total);
	sstream_out << "point_values_vec as vector of tasks" << endl;
	sstream_out << "point_values_vec.size() " << point_values_vec.size() << endl;
	unsigned long long index = 0;
	vector<double> cur_point_values;
	vector<int> index_arr;
	while (CAMBALA_utils::next_cartesian(search_space, index_arr, cur_point_values))
		point_values_vec[index++] = cur_point_values;
	cout << "next_cartesian() finished" << endl;

	sstream_out << "point_values_vec[0].size() " << point_values_vec[0].size() << endl;
	unsigned task_len = point_values_vec[0].size() + 1; // point data + task index
	double *task = new double[task_len];
	unsigned result_len = 2; // calculated residual + index of a task
	double *result = new double[result_len];
	sstream_out << "task_len " << task_len << endl;
	sstream_out << "result_len " << result_len << endl;

	unsigned send_task_count = 0;
	int depths_p_len = depths.size();
	double *depths_p = new double[depths_p_len];
	cout << "depths_p" << endl;
	for (unsigned i = 0; i < depths_p_len; i++)
		depths_p[i] = depths[i];
	
	// sending first part of tasks
	for (int computing_process_index = 1; computing_process_index < corecount; computing_process_index++) {
		for (unsigned j = 0; j < task_len - 1; j++)
			task[j] = point_values_vec[send_task_count][j];
		task[task_len - 1] = (double)send_task_count;
		if (verbosity > 0) {
			cout << "Sending from control process" << endl;
			cout << "depths_p_len " << depths_p_len << endl;
			cout << "depths_p" << endl;
			for (unsigned i = 0; i < depths_p_len; i++)
				cout << depths_p[i] << " ";
			cout << endl;
			cout << "task_len " << task_len << endl;
			cout << "task" << endl;
			for (unsigned i = 0; i < task_len; i++)
				cout << task[i] << " ";
			cout << endl;
		}
		MPI_Send(&depths_p_len, 1, MPI_INT, computing_process_index, 0, MPI_COMM_WORLD);
		MPI_Send(depths_p, depths_p_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		MPI_Send(&task_len, 1, MPI_INT, computing_process_index, 0, MPI_COMM_WORLD);
		MPI_Send(task, task_len, MPI_DOUBLE, computing_process_index, 0, MPI_COMM_WORLD);
		send_task_count++;
	}
	delete[] depths_p;
	sstream_out << "send_task_count " << send_task_count << endl;
	cout << "first send_task_count " << send_task_count << endl;

	ofstream ofile(output_filename, ios_base::app);;
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();

	unsigned processed_task_count = 0;
	unsigned received_task_index;
	double received_residual;
	double previous_record_time = MPI_Wtime();
	unsigned long long record_count = 0;
	int zero_message = 0;
	
	// get results and send new tasks on idle computing processes
	while (processed_task_count < point_values_vec.size()) {
		MPI_Recv(result, result_len, MPI_DOUBLE, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		processed_task_count++;
		if (processed_task_count % 100000 == 0) {
			ofstream ofile(output_filename, ios_base::app);
			ofile << endl << "processed_task_count " << processed_task_count << " out of " << point_values_vec.size() << endl;
			ofile << "time from start " << MPI_Wtime() - mpi_start_time << endl;
			ofile.close(); ofile.clear();
		}

		received_residual = result[0];
		received_task_index = (unsigned)result[1];
		//sstream_out << "received_residual " << received_residual << endl;
		//sstream_out << "received_task_index  " << received_task_index << endl;

		if (received_residual < record_point.residual) {
			record_count++;
			record_point = fromDoubleVecToPoint(point_values_vec[received_task_index]);
			record_point.depths = depths;
			record_point.residual = received_residual;

			double current_time = MPI_Wtime();
			if (current_time - previous_record_time >= 1000) {
				reportRecordPoint(record_point, record_count);
				previous_record_time = MPI_Wtime();
			}
		}
		// if free tasks for sending
		if (send_task_count < point_values_vec.size()) {
			for (unsigned j = 0; j < task_len - 1; j++)
				task[j] = point_values_vec[send_task_count][j];
			task[task_len - 1] = send_task_count;
			MPI_Send(&zero_message, 1, MPI_INT, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			MPI_Send(task, task_len, MPI_DOUBLE, status.MPI_SOURCE, 0, MPI_COMM_WORLD);
			send_task_count++;
		}
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

	ofile.open(output_filename, ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();

	delete[] task;
	delete[] result;

	cout << "finilizing process " << rank << endl;
	MPI_Finalize();
#endif
}

void CAMBALA_parallel::computingProcessBruteforce()
{
#ifdef _MPI
	MPI_Status status;
	unsigned result_len = 2; // index of a task + calculated residual
	double *result = new double[result_len];
	vector<double> cur_point_values_vec;
	search_space_point cur_point;
	double task_index;
	double *task;

	if (rank == 1) {
		cout << "computing process, rank " << rank << endl;
		cout << "result_len " << result_len << endl;
	}

	stringstream cur_process_points_sstream;

	int message;
	double d_message;
	vector<double> depths;
	unsigned task_len;
	bool isTaskReceived = false;

	for (;;) {
		MPI_Recv(&message, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
		if (verbosity > 1)
		if (message < 0) {
			// if stop-message then finalize
			cout << "rank " << rank << " received stop-message" << endl;
			break;
		}
		else if (message > 0) {
			// get new depths array
			unsigned depths_len = message;
			double *depths_p = new double[depths_len];
			cout << "depths_p_len " << depths_len << endl;
			MPI_Recv(depths_p, depths_len, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			depths.resize(depths_len);
			for (unsigned j = 0; j < depths_len; j++)
				depths[j] = depths_p[j];
			delete[] depths_p;
			if (rank == 1) {
				cout << "received new depths" << endl;
				for (auto &x : depths)
					cout << x << " ";
				cout << endl;
			}
			MPI_Recv(&task_len, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
			if (isTaskReceived)
				delete[] task;
			task = new double[task_len];
			isTaskReceived = true;

			init(depths);

			cur_process_points_sstream << "cb rhob R tau ";
			for (unsigned i = 0; i < cw1_arr.size(); i++)
				cur_process_points_sstream << "cw" << i << " ";
			cur_process_points_sstream << "residual";
		}
		MPI_Recv(task, task_len, MPI_DOUBLE, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

		if (rank == 1) {
			cout << "received task " << endl;
			for (unsigned i = 0; i < task_len; i++)
				cout << task[i] << " " << endl;
			cout << endl;
		}

		cur_point_values_vec.resize(task_len - 1);
		for (unsigned i = 0; i < task_len - 1; i++)
			cur_point_values_vec[i] = task[i];
		task_index = task[task_len - 1];
		cur_point = fromDoubleVecToPoint(cur_point_values_vec);
		cur_point.depths = depths;
		
		if ( (rank == 1) && ( verbosity > 0 ) )
			cout << "received task_index " << task_index << endl;

		fillDataComputeResidual(cur_point); // calculated residual is written to cur_point

		cur_process_points_sstream << cur_point.cb << " "
			<< cur_point.rhob << " "
			<< cur_point.R << " "
			<< cur_point.tau << " ";
		for (unsigned i = 0; i < cur_point.cws.size(); i++)
			cur_process_points_sstream << cur_point.cws[i] << " ";
		for (unsigned i = 0; i < cur_point.depths.size(); i++)
			cur_process_points_sstream << cur_point.depths[i] << " ";
		cur_process_points_sstream << cur_point.residual;
		cur_process_points_sstream << endl;

		result[0] = cur_point.residual;
		result[1] = task_index;
		MPI_Send(result, result_len, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	delete[] task;
	delete[] result;
	
	stringstream cur_process_file_name_sstream;
	cur_process_file_name_sstream << "points_process" << rank;
	ofstream cur_process_file;
	cur_process_file.open(cur_process_file_name_sstream.str().c_str(),ios_base::app);
	cur_process_file << cur_process_points_sstream.rdbuf();
	cur_process_file.close();
	cur_process_points_sstream.clear(); cur_process_points_sstream.str("");

	MPI_Finalize();
#endif
}

void CAMBALA_parallel::reportRecordPoint( search_space_point record_point, unsigned long long record_count )
{
#ifdef _MPI
	stringstream sstream_out;
	sstream_out << endl << "record_count " << record_count << endl;
	sstream_out << "control process, new residual minimum" << endl;
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

	ofstream ofile(output_filename, ios_base::app);
	ofile << sstream_out.rdbuf();
	sstream_out.clear(); sstream_out.str("");
	ofile.close(); ofile.clear();
#endif
}