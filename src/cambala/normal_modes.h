#ifndef NORMAL_MODES_H
#define NORMAL_MODES_H

#include <vector>

using namespace std;

class normal_modes
{
public:
	normal_modes();
	double freq;
	double iModesSubset;
	unsigned ordRich;
	vector<double> c1s;
	vector<double> c2s;
	vector<double> rhos;
	vector<unsigned> Ns_points;
	vector<double> depths;
	vector<double> zr;

	vector<double> comp_wnumb_extr_lin_dz();
};

#endif
