#ifndef SSPEMDD_PARALLEL_H
#define SSPEMDD_PARALLEL_H

#ifdef _MPI
#include <mpi.h>
#endif

#include <iostream>
#include <vector>
#include <algorithm>
#include "sspemdd_sequential.h"

const double STOP_MESSAGE = -1;

class sspemdd_parallel : public sspemdd_sequential
{
public:
	sspemdd_parallel();
	~sspemdd_parallel();
	int rank;
	int corecount;
	double mpi_start_time;
	unsigned task_len;
	unsigned result_len;
	double *task;
	double *result;
	void MPI_main();
private:
	void control_process();
	void computing_process();
};

#endif