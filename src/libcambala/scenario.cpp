#include "scenario.h"
#include "iostream"
#include "fstream"
#include "utils.h"

//extern string output_filename;
//extern int verbosity;

void getThreeValuesFromStr(string str, double &val1, double &val2, double &val3)
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

Dim getDimFromStr(stringstream& sstream)
{
	string word;
	sstream >> word;
	double l,s,r;
	getThreeValuesFromStr(word, l, s, r);
	return Dim (l,s,r);
}


void outDim(stringstream& ss, Dim x) { 
	//TODO: move me to Dim class and make me produce string
	ss << x.l << " " << x.s << " " << x.r << endl;
}

void Scenario::print()
{
	//TODO: make me produce output in actual script file format!
	stringstream ss;
	ss << "Parameters :" << endl;
	ss << "object_function_type_ " << object_function_type_ << endl;
	ss << "ppm_ " << ppm_ << endl;
	ss << "iterated_local_search_runs_ " << iterated_local_search_runs_ << endl;
	ss << "H_ " << H_ << endl;
	ss << "h_ " << hDim_.l << endl;
	ss << dtimesFileName_ << endl;
	ss << spmagFileName_  << endl;
	ss << "cw_ :" << endl;
	for (auto &x : ssd_.cw)
		outDim(ss, x);
	ss << "depths_ :" << endl;
	for (auto &x : depthsDim_)
		outDim(ss, x);
	ss << "R_ :" << endl;
	outDim(ss, ssd_.R);
	ss << "cb_ :" << endl;
	outDim(ss, ssd_.cb);
	ss << "tau_ :" << endl;
	outDim(ss, ssd_.tau);

	/*
	ofstream ofile(output_filename, ios_base::out);
	ofile << ss.str();
	ofile.close();
	*/
}

int Scenario::readFile(string scenarioFileName)
{
	/*
	if (verbosity > 0)
		cout << "scenarioFileName " << scenarioFileName << endl;
	*/
	ifstream scenarioFile(scenarioFileName.c_str());

	if (!scenarioFile.is_open())
	{
		cerr << "scenarioFile with the name " << scenarioFileName << " wasn't openend" << endl;
		return -1;
	}

	//TODO: replace with something more readable, like case, pattern
	//matching lib or PEGTL
	string str, word, tmp_word;
	stringstream sstream;
	while (getline(scenarioFile, str))
	{
		if ((str == "") || (str[0] == '%'))
			continue;
		sstream << str;
		sstream >> word;
		if (word.find("dtimes_file") != string::npos)
			sstream >> dtimesFileName_;
		else if (word.find("spmag_file") != string::npos)
			sstream >> spmagFileName_;
		else if (word == "H")
			sstream >> H_;
		else if (word == "h")
			sstream >> hDim_.l;
		else if ((word.size() >= 2) && (word[0] == 'c') && (word[1] == 'w'))
		{
			word = word.substr(2, word.size() - 2);
			unsigned cw_index = 0;
			istringstream(word) >> cw_index; // FIXME: we don't use the index actually!
			ssd_.cw.push_back(getDimFromStr(sstream));
		}
		else if ((word.size() == 2) && (word[0] == 'd') && (isdigit(word[1])))
		{
			int d_index;
			word = word.substr(1, word.size() - 1);
			istringstream(word) >> d_index;// FIXME: we don't use the index actually!
			depthsDim_.push_back (getDimFromStr(sstream));
		}
		else if (word == "R")
			ssd_.R = getDimFromStr(sstream);
		else if (word == "rhob")
			ssd_.rhob = getDimFromStr(sstream);
		else if (word == "cb")
			ssd_.cb = getDimFromStr(sstream);
		else if (word == "tau")
			ssd_.tau = getDimFromStr(sstream);
		else if (word == "ils_iterations")
			sstream >> iterated_local_search_runs_;
		else if (word == "function_type")
			sstream >> object_function_type_;
		else if (word == "ppm_")
			sstream >> ppm_;
		sstream.str(""); sstream.clear();
	}
	
	//TODO: replace this checks with except!
	if (ppm_ == 0) 
	{ // if ppm_ wasn't set directly
		if ((object_function_type_ == "uniform") || (object_function_type_ == "weighted"))
			ppm_ = 2;
		else if ((object_function_type_ == "uniform2") || (object_function_type_ == "weighted2"))
			ppm_ = 1;
		else {
			cerr << "unknown object_function_type_ " << object_function_type_ << endl;
			exit(1);
		}
	}

	if (!ssd_.cw.size())
	{
		cerr << "!cw1_.size()" << endl;
		return -1;
	}
	if (!H_)
	{
		cerr << "!H" << endl;
		return -1;
	}

	print();

	/*
	if (verbosity > 0)
		cout << "readScenario() finished" << endl;
		*/

	return 0;
}

void Scenario::readInputDataFromFiles()
{
	auto expdelays_vv = DoubleVecVecRead(dtimesFileName_);
	freqs_ = DoubleVecVecGetFirstColumn(expdelays_vv);
	modal_delays_ = DoubleVecVecGetOtherColumns(expdelays_vv);
	if ("weighted")
	{
		auto spmag_vv = DoubleVecVecRead(spmagFileName_);
		spmag_ = DoubleVecVecGetOtherColumns(spmag_vv);
		//TODO: add verification of freqs between modal_delays and
		//spmag files
	}
}

Scenario::Scenario(string scenarioFileName)
{
	readFile(scenarioFileName);
	print();
}
