#ifndef TYPES_H_
#define TYPES_H_
#include <vector>
#include "constants.h"
struct Point
{
	double R = START_HUGE_VALUE;
	double rhob = START_HUGE_VALUE;
	double cb = START_HUGE_VALUE;
	double tau = START_HUGE_VALUE;
	std::vector<double> cws;
	double residual = START_HUGE_VALUE;

	bool operator==(const Point& a) const
	{
		return (R == a.R && tau == a.tau && rhob == a.rhob && cb == a.cb && cws == a.cws);
	}
	
};

Point CAMBALA_sequential::fromDoubleVecToPoint(vector<double> double_vec)
{
	Point point;
	point.cb   = double_vec[0];
	point.rhob = double_vec[1];
	point.R    = double_vec[2];
	point.tau  = double_vec[3];
	for (unsigned i = 4; i < double_vec.size(); i++)
		point.cws.push_back(double_vec[i]);
	point.residual = START_HUGE_VALUE;
	return point;
}


void CAMBALA_sequential::fromPointToFile(const Point &point, ofstream &ofile)
{
	ofile << point.residual << " " << point.cb << " " << point.rhob << " "
		<< point.R << " " << point.tau << " ";
	for (unsigned i = 0; i < point.cws.size(); i++)
		ofile << point.cws[i] << " ";
}
// function for BOINC client application
Point CAMBALA_sequential::fromStrToPoint(string str)
{
	Point point;
	stringstream sstream;
	sstream << str;
	sstream >> point.residual >> point.cb >> point.rhob >> point.R >> point.tau;
	double val;
	while (sstream >> val)
		point.cws.push_back(val);
	return point;
}
void CAMBALA_sequential::PrintPoint(const Point& point)
{
	cout << endl;
	cout << endl << "New residual minimum:" << endl;
	cout << "err = " << point.residual << ", parameters:" << endl;
	cout << "c_b = " << point.cb
		<< ", rho_b = " << point.rhob
		<< ", tau = " << point.tau
		<< ", R = " << point.R << endl;
	cout << "cws_min :" << endl;
	for (auto x : point.cws)
		cout << x << " ";
	cout << endl;
	cout << endl;
	cout << "Ns_points " << endl;
	for (auto x : Ns_points)
		cout << x << " ";
	cout << endl;
	cout << endl;
}

struct Dim
{
	double l = 0;
	double s = 0;
	double r = 0;
};

struct SearchSpaceDims
{
	Dim cb;
	Dim rhob;
	Dim R;
	Dim tau;
	vector <Dim> cw;
};



/*
struct reduced_search_space_attribute
{
	bool R;
	bool rhob;
	bool cb;
	bool tau;
	std::vector<bool> cws;
};
*/
#endif
