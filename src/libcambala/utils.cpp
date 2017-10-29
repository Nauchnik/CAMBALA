#include "utils.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <algorithm> // for correct std::remove 

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
