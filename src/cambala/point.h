#ifndef POINT_H
#define POINT_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>

const double START_HUGE_VALUE = 1e100;

using namespace std;

namespace CAMBALA_point {
	struct search_space_point
	{
		double R;
		double rhob;
		double cb;
		double tau;
		vector<double> cws;
		double residual;
		vector<double> depths;

		bool operator==(const search_space_point& a) const
		{
			return (R == a.R && tau == a.tau && rhob == a.rhob && cb == a.cb && cws == a.cws && depths == a.depths);
		}
		friend bool operator<(const search_space_point& a, const search_space_point& b)
		{
			return (a.residual < b.residual);
		}
	};

	struct reduced_search_space_attribute
	{
		bool R;
		bool rhob;
		bool cb;
		bool tau;
		vector<bool> cws;
	};

	inline search_space_point fromDoubleVecToPoint(const vector<double> double_vec)
	{
		search_space_point point;
		point.cb = double_vec[0];
		point.rhob = double_vec[1];
		point.R = double_vec[2];
		point.tau = double_vec[3];
		for (unsigned i = 4; i < double_vec.size(); i++)
			point.cws.push_back(double_vec[i]);
		point.residual = START_HUGE_VALUE;
		return point;
	}

	// function for BOINC client application
	inline search_space_point fromStrToPoint(const string str, const unsigned long long cws_count )
	{
		search_space_point point;
		stringstream sstream;
		sstream << str;
		sstream >> point.residual >> point.cb >> point.rhob >> point.R >> point.tau;
		double val = 0;
		for (unsigned long long i = 0; i < cws_count; i++) {
			sstream >> val;
			point.cws.push_back(val);
		}
		while (sstream >> val)
			point.depths.push_back(val);
		
		return point;
	}

	inline void fromPointToFile(const search_space_point point, ofstream &ofile)
	{
		ofile << point.residual << " " << point.cb << " " << point.rhob << " "
			<< point.R << " " << point.tau << " ";
		for (unsigned i = 0; i < point.cws.size(); i++)
			ofile << point.cws[i] << " ";
		for (unsigned i = 0; i < point.depths.size(); i++)
			ofile << point.depths[i] << " ";
	}

	inline bool isNan(const double x) { // check if a given valune is not a number (NaN)
		return x != x;
	}
	
	inline bool isCorrectCalculatedPoint(const search_space_point point) 
	{
		if ((point.residual <= 0) || (point.cb <= 0) || (point.R <= 0) || (point.rhob <= 0))
			return false;
		if (isNan(point.residual) || isNan(point.cb) || isNan(point.R) || isNan(point.rhob))
			return false;
		
		for (auto &x : point.cws)
			if (x <= 0) return false;
		for (auto &x : point.cws)
			if (isNan(x)) return false;
		
		for (auto &x : point.depths)
			if (x < 0) return false;
		for (auto &x : point.depths)
			if (isNan(x)) return false;
		
		return true;
	}

	inline string strPointData(const search_space_point point)
	{
		stringstream sstream;
		sstream << "err = " << point.residual << ", parameters :\n";
		sstream << "c_b = " << point.cb << endl
			<< "tau = " << point.tau << endl
			<< "rho_b = " << point.rhob << endl
			<< "R = " << point.R << endl;
		sstream << "cws :\n";
		for (unsigned j=0; j < point.cws.size()-1; j++)
			sstream << point.cws[j] << " ";
		sstream << point.cws[point.cws.size() - 1] << endl;
		sstream << "depths :\n";
		for (unsigned j = 0; j < point.depths.size() - 1; j++)
			sstream << point.depths[j] << " ";
		sstream << point.depths[point.depths.size() - 1] << endl;
		sstream << endl;
		return sstream.str();
	}
}

#endif
