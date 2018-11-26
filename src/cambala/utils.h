#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

namespace CAMBALA_utils{
	inline bool is_pos_number(const string& s)
	{
		if (s.empty())
			return false;
		for (auto x : s)
			if (!isdigit(x))
				return false;
		return true;
	}

	inline void getThreeValuesFromStr(string str, double &val1, double &val2, double &val3)
	{
		val1 = val3 = -1;
		val2 = 1;
		string word1, word2, word3;
		for (auto &x : str)
			if (x == ':')
				x = ' ';
		stringstream sstream;
		sstream << str;
		sstream >> word1 >> word2 >> word3;
		istringstream(word1) >> val1;
		istringstream(word2) >> val2;
		istringstream(word3) >> val3;
		if (val3 == -1)
			val3 = val1;
	}

	inline string doubleToStr(const double d) { 
		stringstream sstream; 
		sstream << d;
		return sstream.str();
	}

	inline string doubleVecToStr(vector<double> vec) {
		if (vec.empty())
			return "";
		string str;
		for (unsigned i=0; i<vec.size()-1; i++)
			str += doubleToStr(vec[i]) + " ";
		str += doubleToStr(vec[vec.size()-1]);
		return str;
	}

	inline double getMidVecValue(vector<double> vec) 
	{
		return (vec.size() == 1) ? vec[0] : (vec[vec.size() / 2]);
	}
	
	// construct all combinations of some parameters
	template< typename T >
	bool next_cartesian(vector<T> &vii, vector<int> &index_arr, T &cur_vi)
	{
		if (index_arr.size() == 0) { // init
			index_arr.resize(vii.size());
			//for( auto &x : index_arr )
			//	x = 0;
			for (vector<int> ::iterator it = index_arr.begin(); it != index_arr.end(); ++it)
				*it = 0;
		}
		if (index_arr[0] == -1)
			return false;
		// get current value
		cur_vi.resize(vii.size());
		for (unsigned i = 0; i < index_arr.size(); ++i)
			cur_vi[i] = vii[i][index_arr[i]];
		// check if last iteration
		bool IsLastValue = true;
		for (unsigned i = 0; i < index_arr.size(); ++i) {
			if (index_arr[i] != vii[i].size() - 1) {
				IsLastValue = false;
				break;
			}
		}
		if (IsLastValue)
			index_arr[0] = -1; // condition of stopping
		else {
			// find last changable row to increase its value
			int last_changable = index_arr.size() - 1;
			while (last_changable != -1){
				if (index_arr[last_changable] < (int)(vii[last_changable].size() - 1))
					break;
				--last_changable;
			}
			index_arr[last_changable]++;
			for (unsigned i = last_changable + 1; i < index_arr.size(); ++i)
				index_arr[i] = 0;
		}

		return true;
	}
};

#endif