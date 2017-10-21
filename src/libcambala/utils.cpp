#include "utils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm> // for correct std::remove 


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

//TODO:rewrite me as separate template class!
vector <vector <double>> DoubleVecVecRead(string filename)
{
	vector <vector <double>> vecvec;
	ifstream file(filename.c_str());
	if (!file.is_open())
	{
		cerr << "file " << filename << " wasn't opened" << endl;
		exit(1);
	}
	stringstream myLineStream;
	string myLine;
	while (getline(file, myLine))
	{ 
		// delete windows endline symbol for correct reading
		myLine.erase(remove(myLine.begin(), myLine.end(), '\r'), myLine.end());
		myLineStream << myLine;

		vector<double> buffvect;
		while (!myLineStream.eof())
		{
			double buff;
			myLineStream >> buff;
			buffvect.push_back(buff);
		}

		vecvec.push_back(buffvect);
		myLineStream.str(""); myLineStream.clear();
	}
	file.close();
	return vecvec;
}


vector <double> DoubleVecVecGetFirstColumn(vector <vector <double>> vecvec)
{
	vector <double> first_column;
	for (auto fl: vecvec)
		first_column.push_back(fl[0]);
	return first_column;
}

vector <vector <double>> DoubleVecVecGetOtherColumns(vector <vector <double>> vecvec)
{
	vector <vector <double>> out;
	for (auto fl: vecvec)
		out.push_back(vector<double>(fl.begin()+1, fl.end()));
	return out;
}
