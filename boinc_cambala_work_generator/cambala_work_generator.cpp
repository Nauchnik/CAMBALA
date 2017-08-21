/*****************************************************************************************
Work generator for Acoustics@home. 2017 (c).
Format of WUs: scenarios of Acoustics@home.
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/
#include "cambala_sequential.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

int main( int argc, char **argv )
{
#ifdef _DEBUG
	argc = 2;
	argv[1] = "./43_hydro_R_uniform260.txt";
#endif

	if (argc < 2) {
		cerr << "Usage: scenario_file" << endl;
		return 1;
	}
	
	string scenario_file_name = argv[1];
	
	CAMBALA_sequential cambala_seq;
	int retval;
	retval = cambala_seq.readScenario(scenario_file_name);
	if (retval) {
		cerr << "cambala_seq.readScenario() falied " << endl;
		exit(retval);
	}
	vector<vector<double>> depths_vec;
	cambala_seq.createDepthsArray(cambala_seq.h1, depths_vec);
	retval = cambala_seq.init(depths_vec[0]);
	if (retval) {
		cerr << "cambala_seq.init() falied " << endl;
		exit(retval);
	}
	
	stringstream sstream_static_strings, sstream;
	ifstream ifile(scenario_file_name.c_str());
	string str, word, word1, word2, word3;
	reduced_search_space_attribute reduced_s_s_a;
	reduced_s_s_a.cb = reduced_s_s_a.R = reduced_s_s_a.tau = reduced_s_s_a.rhob = false;
	reduced_s_s_a.cws.resize(cambala_seq.cw1_arr.size());
	for (unsigned i=0; i<reduced_s_s_a.cws.size(); i++)
		reduced_s_s_a.cws[i] = false;
	
	unsigned cw_index;
	while (getline(ifile, str)) {
		sstream << str;
		word1 = word2 = word3 = "";
		sstream >> word1 >> word2 >> word3;
		if (word3 == "server") {
			if ( word1 == "R" )
				reduced_s_s_a.R = true;
			else if (word1 == "cb")
				reduced_s_s_a.cb = true;
			else if (word1 == "rhob")
				reduced_s_s_a.rhob = true;
			else if (word1 == "tau")
				reduced_s_s_a.tau = true;
			else if ((word1.size() >= 2) && (word1[0] == 'c') && (word1[1] == 'w')) {
				word = word1.substr(2, word1.size() - 2);
				istringstream(word) >> cw_index;
				reduced_s_s_a.cws[cw_index] = true;
			}
		}
		else
			sstream_static_strings << str << endl;
		sstream.str("");
		sstream.clear();
	}
	ifile.close();

	cout << "reduced_s_s_a.R " << reduced_s_s_a.R << endl;
	cout << "reduced_s_s_a.cb " << reduced_s_s_a.cb << endl;
	cout << "reduced_s_s_a.rhob " << reduced_s_s_a.rhob << endl;
	cout << "reduced_s_s_a.tau " << reduced_s_s_a.tau << endl;
	for (unsigned i=0; i < reduced_s_s_a.cws.size(); i++)
		cout << "reduced_s_s_a.cws index " << i << " " << reduced_s_s_a.cws[i] << endl;
	cout << endl;
	
	cambala_seq.reduceSearchSpace(reduced_s_s_a);
	vector<search_space_point> reduced_points_vec = cambala_seq.getSearchSpacePointsVec(depths_vec[0]);
	cout << "reduced_points_vec.size() " << reduced_points_vec.size() << endl;

	string str_to_remove = "./";
	unsigned pos = scenario_file_name.find(str_to_remove);
	if (pos != string::npos)
		scenario_file_name.erase(pos, str_to_remove.length());
	fstream temp_wu_file;
	for (unsigned long long i = 0; i < reduced_points_vec.size(); i++) {
		sstream << scenario_file_name << "-wu" << i+1;
		string wu_name = sstream.str();
		sstream.str(""); sstream.clear();
		string cur_wu_input_file_name = "input_" + wu_name;

		temp_wu_file.open("tmp_wu_file", ios_base::out);
		if (!temp_wu_file.is_open()) {
			cerr << "Failed to create tmp_wu_file" << endl;
			exit(1);
		}

		// write input data to WU file
		temp_wu_file << sstream_static_strings.str();
		if (reduced_s_s_a.cb == true)
			temp_wu_file << "cb " << reduced_points_vec[i].cb << endl;
		if (reduced_s_s_a.rhob == true)
			temp_wu_file << "rhob " << reduced_points_vec[i].rhob << endl;
		if (reduced_s_s_a.R == true)
			temp_wu_file << "R " << reduced_points_vec[i].R << endl;
		if (reduced_s_s_a.tau == true)
			temp_wu_file << "tau " << reduced_points_vec[i].tau << endl;
		for (unsigned j=0; j<reduced_s_s_a.cws.size(); j++) {
			if (reduced_s_s_a.cws[j] == true)
				temp_wu_file << "cw" << j << " " << reduced_points_vec[i].cws[j] << endl;
		}
		temp_wu_file.close();
		temp_wu_file.clear();
		
		string system_str = "cp tmp_wu_file `./bin/dir_hier_path " + cur_wu_input_file_name + "`";
		//cout << "before system command : " << system_str << endl; 
		system(system_str.c_str());
		//cout << "done" << endl;
		system_str = "./bin/create_work -appname sspemdd -wu_name " + wu_name +
			" -wu_template ./templates/workunit_cambala.xml" +
			" -result_template ./templates/result_cambala.xml " + cur_wu_input_file_name;
		cout << system_str << endl;
		system(system_str.c_str());
		//cout << "done" << endl;
	}
	
	return 0;
}