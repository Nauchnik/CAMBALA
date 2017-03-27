/*****************************************************************************************
Work generator for Acoustics@home. 2017 (c).
Format of WUs: scenarios of Acoustics@home.
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/
#include "sspemdd_sequential.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

int main( int argc, char **argv )
{
#ifdef _DEBUG
	argc = 2;
	argv[1] = "310_hydro_R_uniform260_work_generator.txt";
#endif

	if (argc < 2) {
		std::cerr << "Usage: scenario_file" << std::endl;
		return 1;
	}
	
	std::string scenario_file_name = argv[1];
	
	sspemdd_sequential sspemdd_seq;
	sspemdd_seq.readScenario(scenario_file_name);
	sspemdd_seq.init();

	std::stringstream sstream_static_strings, sstream;
	std::ifstream ifile(scenario_file_name.c_str());
	std::string str, word, word1, word2, word3;
	reduced_search_space_attribute reduced_s_s_a;
	reduced_s_s_a.cb = reduced_s_s_a.R = reduced_s_s_a.tau = reduced_s_s_a.rhob = false;
	reduced_s_s_a.cws.resize(sspemdd_seq.cw1_arr.size());
	for (unsigned i=0; i<reduced_s_s_a.cws.size(); i++)
		reduced_s_s_a.cws[i] = false;
	
	unsigned cw_index;
	while (getline(ifile, str)) {
		sstream << str;
		sstream >> word1 >> word2 >> word3;
		if (word3 == "server") {
			if ( word1 == "R" )
				reduced_s_s_a.R = true;
			else if (word1 == "tau")
				reduced_s_s_a.tau = true;
			else if (word1 == "cb")
				reduced_s_s_a.tau = true;
			else if (word1 == "rhob")
				reduced_s_s_a.tau = true;
			else if ((word1.size() >= 2) && (word1[0] == 'c') && (word1[1] == 'w')) {
				word = word1.substr(2, word1.size() - 2);
				std::istringstream(word) >> cw_index;
				reduced_s_s_a.cws[cw_index] = true;
			}
		}
		else
			sstream_static_strings << str << std::endl;
		sstream.str("");
		sstream.clear();
	}
	ifile.close();

	sspemdd_seq.reduceSearchSpace(reduced_s_s_a);
	std::vector<search_space_point> points_vec = sspemdd_seq.getSearchSpacePointsVec();
	std::cout << "points_vec.size() " << points_vec.size() << std::endl;

	std::fstream temp_wu_file;
	for (unsigned long long i = 0; i < points_vec.size(); i++) {
		sstream << "310_hydro_R_uniform260" << "-wu" << i+1;
		std::string wu_name = sstream.str();
		sstream.str(""); sstream.clear();
		std::string cur_wu_input_file_name = "input_" + wu_name;

		temp_wu_file.open("tmp_wu_file", std::ios_base::out);
		if (!temp_wu_file.is_open()) {
			std::cerr << "Failed to create tmp_wu_file" << std::endl;
			exit(1);
		}

		// write input data to WU file
		temp_wu_file << sstream_static_strings.str();
		if (reduced_s_s_a.cb == true)
			temp_wu_file << "cb " << points_vec[i].cb << std::endl;
		if (reduced_s_s_a.rhob == true)
			temp_wu_file << "rhob " << points_vec[i].rhob << std::endl;
		if (reduced_s_s_a.R == true)
			temp_wu_file << "R " << points_vec[i].R << std::endl;
		if (reduced_s_s_a.tau == true)
			temp_wu_file << "tau " << points_vec[i].tau << std::endl;
		for (unsigned j=0; j<reduced_s_s_a.cws.size(); j++) {
			if (reduced_s_s_a.cws[j] == true)
				temp_wu_file << "cw" << j << " " << points_vec[i].cws[j] << std::endl;
		}
		temp_wu_file.close();
		temp_wu_file.clear();
		
		std::string system_str = "cp tmp_wu_file `./bin/dir_hier_path " + cur_wu_input_file_name + "`";
		//std::cout << "before system command : " << system_str << std::endl; 
		system(system_str.c_str());
		//std::cout << "done" << std::endl;
		system_str = "./bin/create_work -appname sspemdd -wu_name " + wu_name +
			" -wu_template ./templates/workunit_sspemdd.xml" +
			" -result_template ./templates/result_sspemdd.xml " + cur_wu_input_file_name;
		std::cout << system_str << std::endl;
		system(system_str.c_str());
		//std::cout << "done" << std::endl;
	}
	
	return 0;
}