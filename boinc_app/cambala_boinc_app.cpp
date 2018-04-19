/*****************************************************************************************
// CAMBALA-based BOINC client application for Acoustics@home -- Copyright(c) 2017
// CAMBALA stands for Coupled Acoustic Modes
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS)
*****************************************************************************************/

#ifdef _WIN32
#include "boinc_win.h"
#else
#include "config.h"
#include <cstdio>
#include <cctype>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <csignal>
#include <unistd.h>
#endif

#include "str_util.h"
#include "util.h"
#include "filesys.h"
#include "boinc_api.h"
#include "mfile.h"

#include <sstream>
#include <string>
#include <fstream>

#include "sequential.h"
#include "utils.h"
#include "point.h"

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

using namespace std;

int do_work(search_space_point &cur_record_point);
int do_checkpoint( const unsigned long long &total_points, 
	               const unsigned long long &processed_points, 
	               const search_space_point &cur_record_point );

int main(int argc, char **argv)
{
    char buf[256];
	int retval = boinc_init();
    if ( retval ) {
        fprintf(stderr, "%s APP: boinc_init() returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit( retval );
    }
	
	search_space_point cur_record_point;
	cur_record_point.residual = START_HUGE_VALUE;

	retval = do_work(cur_record_point);
	if ( retval ) {
		fprintf( stderr, "%s APP: do_work() returned \n" );
		exit(retval);
	}

	// resolve, open and write answer to output file
	string output_file_name;
    boinc_resolve_filename_s( OUTPUT_FILENAME, output_file_name);
	ofstream output_file(output_file_name.c_str());
    if ( !output_file.is_open() ) {
        fprintf(stderr, "%s APP: app output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
		exit(retval);
    }
	
	fromPointToFile(cur_record_point, output_file);
	
	output_file.close();
	
    boinc_finish(0);
}

int do_work (search_space_point &cur_record_point)
{
	string input_file_name = "";
	string str = "';

	// open the input file (resolve logical name first)
	boinc_resolve_filename_s(INPUT_FILENAME, input_file_name);

	if (input_file_name == "")
		return false;

	CAMBALA_sequential cambala_seq;
	int retval;
	retval = cambala_seq.readScenario(input_file_name);
	if (retval) {
		fprintf(stderr, "APP: readScenario() failed %d\n", retval);
		exit(retval);
	}
	retval = cambala_seq.readInputDataFromFiles();
	if (retval) {
		fprintf(stderr, "APP: readInputDataFromFiles() failed %d\n", retval);
		exit(retval);
	}
	vector<vector<double>> depths_vec = cambala_seq.createDepthsArray();
	if (depths_vec.size() > 1) {
		fprintf(stderr, "APP: depths_vec.size() %zd\n", depths_vec.size());
		exit(1);
	}
	vector<double> depths = depths_vec[0];
	depths_vec.clear();
	cout << "depths :";
	for (unsigned i = 0; i < depths.size(); i++)
		cout << " " << depths[i];
	cout << endl;
	retval = cambala_seq.init(depths);
	if (retval) {
		fprintf(stderr, "APP: init() failed %d\n", retval);
		exit(retval);
	}
	vector<search_space_point> points_vec = cambala_seq.getSearchSpacePointsVec(depths);
	
	unsigned long long total_points = points_vec.size();

	// read data from the checkpoint file if such exists
	string chpt_file_name;
	fstream chpt_file;
	unsigned long long processed_points = 0;
	boinc_resolve_filename_s(CHECKPOINT_FILE, chpt_file_name);
	chpt_file.open(chpt_file_name.c_str(), ios_base::in);
	if (chpt_file.is_open()) {
		getline(chpt_file, str);
		istringstream(str) >> processed_points;
		getline(chpt_file, str);
		unsigned cws_count = cambala_seq.cw1_arr.size();
		cur_record_point = fromStrToPoint(str, cws_count);
		chpt_file.close();
		cout << "point from chpt file" << endl;
	}
	
	if (processed_points == total_points) // exit if all points already processed
		return true;

	if (!total_points) {
		retval = -1;
		fprintf(stderr, "APP: total_points == 0");
		exit(retval);
	}
	
	double dval;
	unsigned long long checkpoint_every_point;
	if ((cambala_seq.object_function_type == "uniform") || (cambala_seq.object_function_type == "weighted"))
		checkpoint_every_point = 10;
	else 
		checkpoint_every_point = 1;
	
	for ( unsigned long long i = processed_points; i < total_points; i++) {
		dval = cambala_seq.fillDataComputeResidual(points_vec[i]);
		if (dval < cur_record_point.residual)
			cur_record_point = points_vec[i];
		
		// checkpoint current results
		//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
		if ((i+1) % checkpoint_every_point == 0) {
			retval = do_checkpoint(total_points, i + 1, cur_record_point);
			if (retval) {
				fprintf(stderr, "APP: checkpoint failed %d\n", retval);
				exit(retval);
			}
			boinc_checkpoint_completed();
		}
		cout << "processed " << i+1 << " out of " << total_points << endl;
        //}
	}
	
	return 0;
}

int do_checkpoint( const unsigned long long &total_points, 
	               const unsigned long long &processed_points, 
	               const search_space_point &cur_record_point )
{
    string resolved_name = "";
	
	ofstream temp_ofile( "temp" );
	if (!temp_ofile.is_open()) {
		fprintf(stderr, "APP: do_checkpoint() temp ofile wasn't opened");
		return -1;
	}
	
	temp_ofile << processed_points << endl;
	fromPointToFile(cur_record_point, temp_ofile);
    temp_ofile.close();
	
    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
	if (resolved_name == "") {
		fprintf(stderr, "APP: resolved_name is empty");
		return -1;
	}
	
	int retval = 0;
    retval = boinc_rename( "temp", resolved_name.c_str() );
	if ( retval ) {
		fprintf(stderr, "APP: do_checkpoint() boinc_rename() returned %d", retval);
		return -1;
	}

	if (!processed_points) {
		fprintf(stderr, "APP: do_checkpoint() processed_points == 0");
		return -1;
	}

	if (!total_points) {
		fprintf(stderr, "APP: do_checkpoint() total_points == 0");
		return -1;
	}

	boinc_fraction_done( (double)processed_points / (double)total_points );

    return 0;
}

#ifdef _WIN32
int WINAPI WinMain(HINSTANCE hInst, HINSTANCE hPrevInst, LPSTR Args, int WinMode) {
    LPSTR command_line;
    char* argv[100];
    int argc;

    command_line = GetCommandLine();
    argc = parse_command_line( command_line, argv );
    return main(argc, argv);
}
#endif