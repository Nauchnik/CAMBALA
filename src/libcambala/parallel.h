#ifndef CAMBALA_PARALLEL_H
#define CAMBALA_PARALLEL_H

#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include "cambala_sequential.h"

const double STOP_MESSAGE = -1;
const unsigned ILS_TASK_LEN = 10;
const unsigned ILS_RESULT_LEN = 20;
const unsigned MAX_DEPTHS_VECTORS = 1000000;

class CAMBALA_parallel : public CAMBALA_sequential
{
public:
	CAMBALA_parallel();
	~CAMBALA_parallel();
	int corecount;
	double mpi_start_time;
	void MPI_main();
private:
	void controlProcessIls();
	void computingProcessIls();
	void sendTaskIls(double *task, int task_index, unsigned process_index, vector<double> depths);

	void controlProcessBruteforce();
	void controlProcessFixedDepths(const vector<double> depths, const unsigned depths_index);
	void computingProcessBruteforce();
	
	void reportRecordPoint( Point record_point, unsigned long long record_count );
};

#endif