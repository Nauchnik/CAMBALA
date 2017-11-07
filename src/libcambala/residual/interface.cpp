#include "residual/interface.h"
#include "easylogging++.h"




double ResidualCalculator::ComputeResidual(const Model& m, Point& p)
{
	LoadModel(m);
	return ComputeResidual(p);
}

void ResidualCalculator::LoadModel(const Model& m)
{
	LOG(DEBUG) << "ResCalc " << getName() << " loading model data";
	DoLoadModel(m);
}

double ResidualCalculator::ComputeResidual(Point& p)
{
	LOG(DEBUG) << "ResCalc " << getName() << " starting residual calculation";
	double result = DoComputeResidual(p);
	LOG(DEBUG) << "Residual: " << result;
	return result;
}

/*
ResidualCalculator::ResidualCalculator(std::string name)
{
	name_ = name;
	LOG(INFO) << "Created residual calculator \""<< name_ << "\"." ;
}

ResidualCalculator::ResidualCalculator() { }
*/

std::string ResidualCalculator::getName() { return name_; }
ResidualCalculator::~ResidualCalculator(){}
