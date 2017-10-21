#ifndef CAMBALA_H_
#define CAMBALA_H_
#include "scenario.h"


using namespace std;


class CAMBALA
{
public:
	void Solve(const Scenario& c);
	void directPointCalc(Point p);
	void reportFinalResult();
};

#endif
