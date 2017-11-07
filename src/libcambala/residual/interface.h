#ifndef RESIDUAL_INTERFACE_H_
#define RESIDUAL_INTERFACE_H_
#include <string>
#include "types.h"
// Using "Template Method" design pattern!

class ResidualCalculator
{
public:
	void LoadModel (const Model& m);
	double ComputeResidual(Point& p);
	double ComputeResidual(const Model& m, Point& p);
	std::string getName();

	virtual ~ResidualCalculator();

protected:
	std::string name_;
private:
	virtual void DoLoadModel (const Model& m) = 0;
	virtual double DoComputeResidual(Point& p) = 0;
};

using ResCalc = ResidualCalculator;

#endif
