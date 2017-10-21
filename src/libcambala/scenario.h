
class Scenario
{
public:
	double H_;
	double h_;
	string dtimesFileName_;
	string spmagFileName_;
	SearchSpaceDims ssd_;
	vector <Dim> depthsDim_;
	size_t iterated_local_search_runs_;
	string object_function_type_;
	size_t ppm_;
	vector <double> freqs_; // loaded from 1st column of dtimesFileName_
	vector <vector<double>>modal_delays_; // loaded from dtimesFilename_
	vector <vector<double>> spmag_; // loaded from smpagFilename_
	int Scenario::readFile(string scenarioFileName);
	void Scenario::print();
};
