#ifndef SOLVERS_INTERFACE_H_
#define SOLVERS_INTERFACE_H_

#include <vector>
#include "types.h"
#include "residual/selector.h"

class Solver
{
public:
	virtual void Solve () = 0;
	//virtual void LoadData(std::vector<double> depths) = 0;
	virtual Point getBestPoint() = 0;
	virtual void SetResidualCalculatorSelector(ResCalcSelector* p_sel) = 0;
	virtual void LoadSearchSpaceDims(SearchSpaceDims ssd) = 0;
	//virtual ~Solver(){};
};

#endif
