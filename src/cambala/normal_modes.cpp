#include "normal_modes.h"
//#include "compute.h"

//using namespace CAMBALA_compute;

normal_modes::normal_modes():
	freq(0.0),
	iModesSubset(-1.0)
{}

vector<double> normal_modes::comp_wnumb_extr_lin_dz()
{
	vector<double> out_wnum;
	//double omeg1 = 2 * LOCAL_M_PI*freq;
	//out_wnum = compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, iModesSubset, ordRich);
	return out_wnum;
}

