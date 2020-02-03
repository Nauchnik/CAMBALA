#include <iostream>
#include <string>
#include <chrono>
#include "normal_modes.h"

int main(int argc, char ** argv)
{
	string scenarioFileName = "";
#ifdef _DEBUG
	//scenarioFileName = "../../scenarios/CambalaNormalModes/test1_scenario.txt";
	scenarioFileName = "../../scenarios/CambalaNormalModes_ac_modes/case_1/case1_scenario.txt";
#else
	if (argc < 2) {
		cout << "Usage: program normal_modes_scenario_name\n";
		exit(-1);
	}
	scenarioFileName = argv[1];
#endif
	cout << "scenarioFileName " << scenarioFileName << endl;
	chrono::high_resolution_clock::time_point start_t = chrono::high_resolution_clock::now();

	NormalModes n_m;
	n_m.read_scenario(scenarioFileName);
	n_m.compute_for_all_depths();

	chrono::high_resolution_clock::time_point cur_t = chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = chrono::duration_cast<chrono::duration<double>>(cur_t - start_t);
	cout  << "elapsed " << time_span.count() << " seconds\n";
	//n_m.print_wnumbers();
	//n_m.print_mfunctions_zr();
	//n_m.print_modal_group_velocities();

	return 0;
}