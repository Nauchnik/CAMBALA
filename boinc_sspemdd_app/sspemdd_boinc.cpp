/*****************************************************************************************
// SSPEMDD-based BOINC client application for Acoustics@home -- Copyright(c) 2017
// SSPEMDD stands for Sound Speed Profile Estimator from Modal Delay Data
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS),
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

#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"

#define CHECKPOINT_FILE "chpt"
#define INPUT_FILENAME "in"
#define OUTPUT_FILENAME "out"

search_space_point fromStrToPoint(std::string str);
bool do_work( const std::string &input_file_name, 
	          const unsigned long long &processed_points, 
	          search_space_point &current_record_point );
int do_checkpoint( const unsigned long long &total_points, 
	               const unsigned long long &processed_points, 
	               const search_space_point &current_record_point );
void fromPointToFile(const search_space_point &point, std::ofstream &ofile);

int main(int argc, char **argv) 
{
    char buf[256];
	int retval = boinc_init();
    if ( retval ) {
        fprintf(stderr, "%s boinc_init returned %d\n",
            boinc_msg_prefix(buf, sizeof(buf)), retval
        );
        exit( retval );
    }
	
	std::string input_file_name, output_file_name, chpt_file_name;
	std::string str;
	std::ofstream output_file;
	std::fstream chpt_file;
	unsigned long long processed_points = 0;
	
	// open the input file (resolve logical name first)
	boinc_resolve_filename_s( INPUT_FILENAME, input_file_name);
	
	// read data from the checkpoint file if such exists
	search_space_point cur_record_point;
	cur_record_point.residual = START_HUGE_VALUE;
    boinc_resolve_filename_s( CHECKPOINT_FILE, chpt_file_name);
	chpt_file.open( chpt_file_name.c_str(), std::ios_base::in );
	if ( chpt_file.is_open() ) {
		getline(chpt_file, str);
		std::istringstream(str) >> processed_points;
		getline(chpt_file, str);
		cur_record_point = fromStrToPoint(str);
		chpt_file.close();
		std::cout << "point from chpt file" << std::endl;
	}
	
	if ( !do_work( input_file_name, processed_points, cur_record_point ) ) {
		fprintf( stderr, "APP: do_work() failed:\n" );
		perror("do_work");
        exit(1);
	}

	// resolve, open and write answer to output file
    boinc_resolve_filename_s( OUTPUT_FILENAME, output_file_name);
	output_file.open( output_file_name.c_str() );
    if ( !output_file.is_open() ) {
        fprintf(stderr, "%s APP: app output open failed:\n",
            boinc_msg_prefix(buf, sizeof(buf))
        );
		exit(-1);
    }
	
	fromPointToFile(cur_record_point, output_file);
	
	output_file.close();
	
    boinc_finish(0);
}

search_space_point fromStrToPoint(std::string str)
{
	search_space_point point;
	std::stringstream sstream;
	sstream << str;
	sstream >> point.residual >> point.cb >> point.rhob >> point.R >> point.tau;
	double val;
	while (sstream >> val)
		point.cws.push_back(val);
	return point;
}

bool do_work( const std::string &input_file_name, 
	          const unsigned long long &processed_points,
	          search_space_point &current_record_point)
{
	sspemdd_sequential sspemdd_seq;

	int retval;
	retval = sspemdd_seq.readScenario(input_file_name);
	if (retval) {
		fprintf(stderr, "APP: readScenario() failed %d\n", retval);
		exit(retval);
	}
	retval = sspemdd_seq.readInputDataFromFiles();
	if (retval) {
		fprintf(stderr, "APP: readInputDataFromFiles() failed %d\n", retval);
		exit(retval);
	}
	retval = sspemdd_seq.init();
	if (retval) {
		fprintf(stderr, "APP: init() failed %d\n", retval);
		exit(retval);
	}
	std::vector<search_space_point> points_vec = sspemdd_seq.getSearchSpacePointsVec();
	
	unsigned long long total_points = points_vec.size();
	
	if (processed_points == total_points) // exit if all points already processed
		return true;

	if (!total_points) {
		retval = -1;
		fprintf(stderr, "APP: total_points == 0", retval);
		exit(retval);
	}
	
	double dval;
	for ( unsigned long long i = processed_points; i < total_points; i++) {
		dval = sspemdd_seq.fillDataComputeResidual(points_vec[i]);
		if (dval < current_record_point.residual)
			current_record_point = points_vec[i];
		
		// checkpoint current results
		//if ( ( boinc_is_standalone() ) || ( boinc_time_to_checkpoint() ) ) {
		retval = do_checkpoint(total_points, i+1, current_record_point );
        if (retval) {
			fprintf(stderr, "APP: checkpoint failed %d\n", retval );
            exit(retval);
        }
		boinc_checkpoint_completed();
		std::cout << "processed " << i+1 << " out from " << total_points << std::endl;
        //}
	}

	return true;
}

int do_checkpoint( const unsigned long long &total_points, 
	               const unsigned long long &processed_points, 
	               const search_space_point &current_record_point )
{
	int retval = 0;
    std::string resolved_name;
	
	std::ofstream temp_ofile( "temp" );
	if ( !temp_ofile.is_open() ) return 1;

	temp_ofile << processed_points << std::endl;
	fromPointToFile(current_record_point, temp_ofile);
    temp_ofile.close();
	
    boinc_resolve_filename_s(CHECKPOINT_FILE, resolved_name);
    retval = boinc_rename( "temp", resolved_name.c_str() );

	boinc_fraction_done( (double)processed_points / (double)total_points );

    return retval;
}

void fromPointToFile(const search_space_point &point, std::ofstream &ofile)
{
	ofile << point.residual << " " << point.cb << " " << point.rhob << " "
		       << point.R << " " << point.tau << " ";
	for (unsigned i = 0; i < point.cws.size(); i++)
		ofile << point.cws[i] << " ";
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