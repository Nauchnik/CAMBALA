#ifndef SOLVERS_INTERFACE_H_
#define SOLVERS_INTERFACE_H_

#include "types.h"
#include <vector>

class Solver
{
public:
	virtual void Solve () = 0;
	virtual void LoadData(std::vector<double> depths) = 0;
	virtual Point getBestPoint() = 0
	virtual ~Solver(){};
};

#endif
