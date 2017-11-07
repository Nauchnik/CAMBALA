#ifndef CAMBALA_H_
#define CAMBALA_H_
#include <map>
#include "residual/interface.h"
#include "scenario.h"


//using namespace std;


class CAMBALA
{
public:
	CAMBALA();
	void Solve(const Scenario& c);
	double directPointCalc(Point p);
	void reportFinalResult();
	void AddResidualCalculator(std::string name, ResCalc* rc);
	std::map <std::string, ResCalc*> calcs_;
};

vector<double> makeDepths(double h, double H, const vector <Dim>& d, const vector<Dim>& cw);

#endif
