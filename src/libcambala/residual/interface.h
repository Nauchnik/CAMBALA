#ifndef RESIDUE_CALCULATOR_H_
#define RESIDUE_CALCULATOR_H_
#include "types.h"
#include "model.h"

class ResidualCalculator
{
public:
	virtual void LoadImmutableData (const Model& m) = 0;
	virtual double CalculatePointResidual(Point p) = 0;
	virtual double CalculatePointResidual(const Model& m, Point p) = 0;
	//virtual ~ResidualCalculator(){};
};

using ResCals = ResidualCalculator;

#endif
