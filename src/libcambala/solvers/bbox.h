#ifndef SOLVERS_BBOX_H_
#define SOLVERS_BBOX_H_

#include "solvers/interface.h"
#include <vector>
#include <mpproblem.hpp>
#include <box/box.hpp>

class BlackBox : public Solver
{
public:
	BlackBox (std::string name); 
	BlackBox (); 
	~BlackBox();
private:
	void DoSolve();
	void DoLoadSearchSpaceDims(SearchSpaceDims ssd);
	SearchSpaceDims ssd_;
};


class AcousticsObjective : public COMPI::Functor <double>
{
public:
	AcousticsObjective(ResCalc* c, std::vector <double> mask) : rc_(c), ssdv_mask(mask) { }
	double func(const double* x);
private:
	ResCalc* rc_;
	std::vector <double> ssdv_mask;
};

#endif 
