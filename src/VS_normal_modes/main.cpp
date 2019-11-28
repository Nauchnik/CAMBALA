#include <iostream>
#include <string>
#include "normal_modes.h"

int main(int argc, char ** argv)
{
	string scenarioFileName = "";
#ifdef _DEBUG
	//scenarioFileName = "normal_modes_scenario.txt";
	scenarioFileName = "../../scenarios/CambalaNormalModes/test2_scenario.txt";
#else
	if (argc < 2) {
		cout << "Usage: program normal_modes_scenario_name\n";
		exit(-1);
	}
	scenarioFileName = argv[1];
#endif
	cout << "scenarioFileName " << scenarioFileName << endl;

	NormalModes n_m;
	n_m.read_data(scenarioFileName);
	n_m.compute_khs();
	n_m.print_khs();
	n_m.compute_mfunctions_zr();
	n_m.print_mfunctions_zr();

	return 0;
}