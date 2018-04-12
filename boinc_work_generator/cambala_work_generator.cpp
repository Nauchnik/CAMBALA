/*****************************************************************************************
Work generator for Acoustics@home. 2017 (c).
Format of WUs: scenarios of Acoustics@home.
Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/
#ifndef _WIN32
#include <mysql.h>
#endif

#include "sequential.h"

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

const long long MIN_WUS_FOR_CREATION = 100;
const long long MAX_WUS_FOR_CREATION = 2500;

struct cambala_boinc_config
{
	long long unsent_wus_required;
	long long seconds_between_launches;
	long long total_wus;
	long long created_wus;
};

bool isTest = false;
stringstream sstream_wu_common;
stringstream sstream_scenario_wout_created_wus;
CAMBALA_sequential cambala_seq;
cambala_boinc_config c_b_config;
reduced_search_space_attribute reduced_s_s_a;

#ifndef _WIN32
long long getCountOfUnsentWUs(string pass_file_name);
int processQuery( MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec );
#endif
int readConfig( const string config_file_name );
int generateWUs();
void getReducedSpaceAttrib();
vector<search_space_point> generateWUsFixedDepths( const vector<double> depths );
void runSystemWUsGeneration( vector<search_space_point> &points_to_generate );

string scenario_file_name = "";
string pass_file_name = "";
long long unsent_wus_count = -1;
long long wus_for_creation = 0;

int main( int argc, char **argv )
{
#ifdef _DEBUG
	argc = 3;
	scenario_file_name = "..\\boinc_work_generator\\906_uniform_dynamic_d_2n_boinc";
	//isTest = true;
	unsent_wus_count = 2400;
#else
	if ( argc < 3 ) {
		cerr << "Usage : program scenario_file_name pass_file_name [-test]" << endl;
		return -1;
	}
	scenario_file_name = argv[1];
	string pass_file_name = argv[2]; // read password for database
	if (argc >3) {
		string tmp_str = argv[3];
		if ( (tmp_str == "-test") || (tmp_str == "--test") )
			isTest = true;
	}
#endif

	cout << "scenario_file_name " << scenario_file_name << endl;
	cout << "pass_file_name " << pass_file_name << endl;
	cout << "isTest " << isTest << endl;
	
	int retval = readConfig( scenario_file_name );
	if (retval) {
		cerr << "readConfig() falied \n";
		exit(retval);
	}
	
#ifndef _WIN32
	unsent_wus_count = getCountOfUnsentWUs(pass_file_name);
#endif
	if (unsent_wus_count < 0) {
		cerr << "unsent_wus_count " << unsent_wus_count << endl;
		exit(-1);
	}
	cout << "unsent_wus_count " << unsent_wus_count << endl;
	
	wus_for_creation = MAX_WUS_FOR_CREATION - unsent_wus_count;
	
	cout << "wus_for_creation " << wus_for_creation << endl;
	
	if ( (wus_for_creation >= MIN_WUS_FOR_CREATION) && 
	          (wus_for_creation <= MAX_WUS_FOR_CREATION) && 
              (c_b_config.created_wus < c_b_config.total_wus) )
		generateWUs();
	else {
		cout << "generation is not required recently \n";
		cout << "MIN_WUS_FOR_CREATION " << MIN_WUS_FOR_CREATION << endl;
	}
	
	return 0;
}

int readConfig(const string config_file_name)
{
	cout << "readConfig()" << endl;
	ifstream config_file(config_file_name.c_str());
	if (!config_file.is_open()) {
		cerr << "config_file " << config_file_name << " wasn't opened" << endl;
		return -1;
	}
	c_b_config.unsent_wus_required = -1;
	c_b_config.seconds_between_launches = -1;
	c_b_config.total_wus = -1;
	c_b_config.created_wus = -1;

	string str;
	while (getline(config_file, str)) {
		stringstream sstream;
		string word1 = "", word2 = "";
		sstream << str;
		sstream >> word1 >> word2;
		if (word1 == "unsent_wus_required")
			istringstream(word2) >> c_b_config.unsent_wus_required;
		if (word1 == "seconds_between_launches")
			istringstream(word2) >> c_b_config.seconds_between_launches;
		if (word1 == "total_wus")
			istringstream(word2) >> c_b_config.total_wus;
		if (word1 == "created_wus")
			istringstream(word2) >> c_b_config.created_wus;
	}
	config_file.close();

	if (c_b_config.unsent_wus_required < MIN_WUS_FOR_CREATION) {
		cerr << "c_b_config.unsent_wus_required " << c_b_config.unsent_wus_required << endl;
		return -1;
	}
	if (c_b_config.seconds_between_launches <= 0) {
		cerr << "c_b_config.seconds_between_launches " << c_b_config.seconds_between_launches << endl;
		return -1;
	}
	if (c_b_config.total_wus < MIN_WUS_FOR_CREATION) {
		cerr << "c_b_config.total_wus " << c_b_config.total_wus << endl;
		return -1;
	}
	if (c_b_config.created_wus < 0) {
		cerr << "c_b_config.created_wus " << c_b_config.created_wus << endl;
		return -1;
	}

	cout << "c_b_config.unsent_wus_required " << c_b_config.unsent_wus_required << endl;
	cout << "c_b_config.seconds_between_launches " << c_b_config.seconds_between_launches << endl;
	cout << "c_b_config.total_wus " << c_b_config.total_wus << endl;
	cout << "c_b_config.created_wus " << c_b_config.created_wus << endl;

	return 0;
}

#ifndef _WIN32
long long getCountOfUnsentWUs(string pass_file_name)
{
	long long unsent_count;
	char *host = "localhost";
    char *db;
	char *user;
    char *pass;
	MYSQL *conn;
	
	ifstream pass_file;
	pass_file.open( pass_file_name.c_str() );
	if ( !pass_file.is_open() ) {
		cerr << "pass_file " << pass_file_name << " wasn't opened" << endl;
		exit(1);
	}
	string str;
	getline( pass_file, str );
	db = new char[str.length() + 1];
	strcpy( db, str.c_str() );
	db[str.length()] = NULL;
	getline( pass_file, str );
	user = new char[str.length() + 1];
	strcpy( user, str.c_str() );
	user[str.length()] = NULL;
	getline( pass_file, str );
	pass = new char[str.length() + 1];
	strcpy( pass, str.c_str() );
	pass[str.length()] = NULL;
	cout << "db "   << db   << endl;
	cout << "user " << user << endl;
	cout << "pass " << pass << endl;
	
	conn = mysql_init(NULL);
	if(conn == NULL)
		cerr << "Error: can't create MySQL-descriptor\n" << endl;
	
	if(!mysql_real_connect(conn, host, user, pass, db, 0, NULL, 0)) {
		cerr << "Error: can't connect to MySQL server" << endl;
		exit(1);
	}
	delete[] db;
	delete[] user;
	delete[] pass;

	vector< vector<stringstream *> > result_vec;
	str = "SELECT COUNT(*) FROM workunit WHERE id IN(SELECT workunitid FROM result WHERE server_state = 2)";
	cout << str << endl;

	processQuery( conn, str, result_vec );
	*result_vec[0][0] >> unsent_count;
	result_vec.clear();
	mysql_close(conn);
	
	return unsent_count;
}

int processQuery( MYSQL *conn, string str, vector< vector<stringstream *> > &result_vec )
{
	MYSQL_RES *res;
	MYSQL_ROW row;
	int num_fields;
	
	if ( mysql_query(conn, str.c_str()) != 0 ) {
		cerr << "Error: can't execute SQL-query\n";
		return -1;
	}
	
	res = mysql_store_result( conn );

	if( res == NULL ) {
		cerr << "Error: can't get the result description\n";
		return -1;
	}

	num_fields = mysql_num_fields(res);
	stringstream *sstream_p;
	vector<stringstream *> result_data;

	if ( mysql_num_rows( res ) > 0 ) {
		while((row = mysql_fetch_row(res)) != NULL) {
			for( int i = 0; i < num_fields; ++i ) {
				sstream_p = new stringstream();
				*sstream_p << row[i]; // get value
				result_data.push_back( sstream_p );
			}
			result_vec.push_back( result_data );
			result_data.clear();
		}
	}

	mysql_free_result(res);

	return 0;
}
#endif

int generateWUs()
{
	cout << "generateWUs()" << endl;
	
	int retval = cambala_seq.readScenario(scenario_file_name);
	if (retval) {
		cerr << "cambala_seq.readScenario() falied " << endl;
		exit(retval);
	}
	
	cambala_seq.createDepthsArray();
	cout << "depths_vec size " << cambala_seq.depths_vec.size() << endl;

	getReducedSpaceAttrib();

	long long wus_sum = 0;
	long long skipped_wus = 0;
	long long depths_from = 0, depths_to = 0;
	vector<search_space_point> points;
	for (unsigned i = 0; i < cambala_seq.depths_vec.size(); i++) {
		points = generateWUsFixedDepths(cambala_seq.depths_vec[i]);
		wus_sum += points.size();
		if (wus_sum <= c_b_config.created_wus) {
			depths_from = i + 1;
			skipped_wus += points.size();
		}
		else if ( (!depths_to) && (wus_sum > c_b_config.created_wus + wus_for_creation) )
			depths_to = i;
	}
	
	cout << "depths_from " << depths_from << endl;
	cout << "depths_to " << depths_to << endl;
	cout << "calculated wus_sum " << wus_sum << endl;
	if ((c_b_config.total_wus <= 0) || (c_b_config.total_wus > wus_sum)) {
		c_b_config.total_wus = wus_sum;
		cout << "total_wus was changed to " << c_b_config.total_wus << endl;
	}

	if (isTest) {
		cout << "Test launch, don't generate WUs \n";
		return 0;
	}

	vector<search_space_point> tmp_wus, wus_to_generate;
	for (unsigned i=depths_from; i<=depths_to; i++) {
		points = generateWUsFixedDepths(cambala_seq.depths_vec[i]);
		for (auto &x : points)
			tmp_wus.push_back(x);
	}

	long long wu_from = c_b_config.created_wus - skipped_wus;
	long long wu_to = c_b_config.created_wus - skipped_wus + wus_for_creation;
	if ( (wu_from < 0) || (wu_to >= tmp_wus.size()) ) {
		cerr << "incorrect wu_from or wu_to " << wu_from << " " << wu_to << endl;
		exit;
	}
	for (unsigned i = wu_from; i < wu_to; i++)
		wus_to_generate.push_back(tmp_wus[i]);
	
	runSystemWUsGeneration(wus_to_generate);
	
	return 0;
}

void getReducedSpaceAttrib()
{
	stringstream sstream;
	ifstream ifile(scenario_file_name.c_str());
	string str, word, word1, word2, word3;
	reduced_s_s_a.cb = reduced_s_s_a.R = reduced_s_s_a.tau = reduced_s_s_a.rhob = false;
	reduced_s_s_a.cws.resize(cambala_seq.cw1_arr.size());
	for (unsigned i = 0; i < reduced_s_s_a.cws.size(); i++)
		reduced_s_s_a.cws[i] = false;

	unsigned cw_index;
	while (getline(ifile, str)) {
		if (str.find("created_wus") == string::npos)
			sstream_scenario_wout_created_wus << str << endl;
		sstream << str;
		word1 = word2 = word3 = "";
		sstream >> word1 >> word2 >> word3;
		if (word3 == "server") {
			if (word1 == "R")
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
		else if ((word1[0] != 'd') || (!isdigit(word1[1])))
			sstream_wu_common << str << endl;
		sstream.str("");
		sstream.clear();
	}
	ifile.close();

	cout << "reduced_s_s_a.R " << reduced_s_s_a.R << endl;
	cout << "reduced_s_s_a.cb " << reduced_s_s_a.cb << endl;
	cout << "reduced_s_s_a.rhob " << reduced_s_s_a.rhob << endl;
	cout << "reduced_s_s_a.tau " << reduced_s_s_a.tau << endl;
	for (unsigned i = 0; i < reduced_s_s_a.cws.size(); i++)
		cout << "reduced_s_s_a.cws index " << i << " " << reduced_s_s_a.cws[i] << endl;
	cout << endl;
}

vector<search_space_point> generateWUsFixedDepths(const vector<double> depths)
{
	int retval = cambala_seq.init(depths);
	if (retval) {
		cerr << "cambala_seq.init() falied " << endl;
		exit(retval);
	}

	cambala_seq.reduceSearchSpace(reduced_s_s_a);
	vector<search_space_point> reduced_points_vec = cambala_seq.getSearchSpacePointsVec(depths);

	return reduced_points_vec;
}

void runSystemWUsGeneration( vector<search_space_point> &wus_to_generate )
{
	string str_to_remove = "./";
	size_t pos = scenario_file_name.find(str_to_remove);
	if (pos != string::npos)
		scenario_file_name.erase(pos, str_to_remove.length());
	fstream temp_wu_file;
	stringstream sstream;
	
	for (unsigned i = 0; i < wus_to_generate.size(); i++) {
		sstream << scenario_file_name << "-wu" << i + 1;
		string wu_name = sstream.str();
		sstream.str(""); sstream.clear();
		string cur_wu_input_file_name = "input_" + wu_name;

		temp_wu_file.open("tmp_wu_file", ios_base::out);
		if (!temp_wu_file.is_open()) {
			cerr << "Failed to create tmp_wu_file" << endl;
			exit(1);
		}

		// write input data to WU file
		temp_wu_file << sstream_wu_common.str();
		if (reduced_s_s_a.cb == true)
		temp_wu_file << "cb " << wus_to_generate[i].cb << endl;
		if (reduced_s_s_a.rhob == true)
		temp_wu_file << "rhob " << wus_to_generate[i].rhob << endl;
		if (reduced_s_s_a.R == true)
		temp_wu_file << "R " << wus_to_generate[i].R << endl;
		if (reduced_s_s_a.tau == true)
		temp_wu_file << "tau " << wus_to_generate[i].tau << endl;
		for (unsigned j = 0; j<reduced_s_s_a.cws.size(); j++) {
			if (reduced_s_s_a.cws[j] == true)
				temp_wu_file << "cw" << j << " " << wus_to_generate[i].cws[j] << endl;
		}
		for (unsigned j = 1; j < reduced_s_s_a.cws.size(); j++) {
			unsigned depths_index = j - 1;
			if (depths_index >= wus_to_generate[i].depths.size()) {
				cerr << "depths_index >= depths.size() \n";
				cerr << depths_index << " >= " << wus_to_generate[i].depths.size() << "\n";
			}
			temp_wu_file << "d" << j << " " << wus_to_generate[i].depths[depths_index] << endl;
		}
		temp_wu_file.close();
		temp_wu_file.clear();
		
		string system_str = "cp tmp_wu_file `./bin/dir_hier_path " + cur_wu_input_file_name + "`";
		//cout << "before system command : " << system_str << endl;
		if (!isTest)
			system(system_str.c_str());
		//cout << "done" << endl;
		system_str = "./bin/create_work -appname sspemdd -wu_name " + wu_name +
		" -wu_template ./templates/workunit_cambala.xml" +
		" -result_template ./templates/result_cambala.xml " + cur_wu_input_file_name;
		cout << system_str << endl;
		if (!isTest)
			system(system_str.c_str());
	}
	c_b_config.created_wus += wus_to_generate.size();
	cout << "new c_b_config.created_wus " << c_b_config.created_wus << endl;
	
	ofstream ofile(scenario_file_name.c_str());
	ofile << sstream_scenario_wout_created_wus.rdbuf();
	ofile << "created_wus " << c_b_config.created_wus;
	ofile.close();
}