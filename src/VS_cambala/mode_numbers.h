#ifndef MODE_NUMBERS_H
#define MODE_NUMBERS_H

#include <vector>

using namespace std;

class mode_numbers 
{
public:
	mode_numbers();
	double freq;
	double iModesSubset;
	vector<double> c1s;
	vector<double> c2s;
	vector<double> rhos;
	vector<unsigned> Ns_points;
	vector<double> depths;
	vector<double> zr;
};

#endif
