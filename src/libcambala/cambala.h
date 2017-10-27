#ifndef CAMBALA_H_
#define CAMBALA_H_
#include "scenario.h"
#include "residual/selector.h"


//using namespace std;


class CAMBALA
{
public:
	void Solve(const Scenario& c);
	double directPointCalc(Point p);
	void reportFinalResult();
	ResCalcSelector res_calc_sel_;
};

vector<double> makeDepths(double h, double H, const vector <Dim>& d, const vector<Dim>& cw);

#endif
