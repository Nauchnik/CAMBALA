#include "cambala.h"
#include "residual/cpu.h"
#include "rescalc_data.cpp"

#include "easylogging++.h"

INITIALIZE_EASYLOGGINGPP

int main(int argc, char *argv[])
{
	ResCalc* cpu64 = new ResCalcCPU64<double>("cpu64");
	Model m;
	Point p;
	SetTestData(m,p);
	LOG(INFO) << "Original residual: " << p.residual;
	cpu64->CalculatePointResidual(m,p);
	LOG(INFO) << "Calculated residual: " << p.residual;
	return(0);
}

	
