#include "cambala.h"
#include "cambala_utils.h"
#include <iostream>
#include <complex>
#include <time.h>
#include <stdexcept>

/*
CAMBALA::CAMBALA() :
	launch_type("bruteforce"),
	object_function_type_("uniform"),
	output_filename("cambala_out"),
	depths_filename("cambala_depths_out"),
	H_(0),
	nh_(1),
	ncb(1),
	nrhob(1),
	nR(1),
	ntau(1),
	cb1(2000.0),
	cb2(2000.0),
	R1(3400.0),
	R2(3600.0),
	tau1(0.0),
	tau2(0.0),
	rhob1(2.0),
	rhob2(2.0),
	n_layers_w(1),
	iterated_local_search_runs_(10),
	verbosity(1),
	isTimeDelayPrinting(false),
	ppm_(0),
	rank(0)
{
	srand((unsigned)time(NULL));
	start_chrono_time = chrono::high_resolution_clock::now();
}

*/

void CAMBALA::reportFinalResult()
{
	/*
	// fix final time
	chrono::high_resolution_clock::time_point t2;
	chrono::duration<double> time_span;
	t2 = chrono::high_resolution_clock::now();
	time_span = chrono::duration_cast<chrono::duration<double>>(t2 - start_chrono_time);

	ofstream ofile(output_filename, ios_base::app);

	ofile << endl;
	ofile << "total solving time (chrono) " << time_span.count() << endl;
	ofile << "SEARCH ENDED!" << endl;
	ofile << "RESULTING VALUE:" << endl;
	ofile << "err = " << record_point_.residual << ", parameters:" << endl;
	ofile << "c_b = " << record_point_.cb << endl
			  << "tau = " << record_point_.tau << endl
			  << "rho_b = " << record_point_.rhob << endl
			  << "R = " << record_point_.R << endl;
	ofile << "cws :" << endl;
	for (auto &x : record_point_.cws)
		ofile << x << " ";
	ofile << endl;
	ofile << "total solving time " << time_span.count() << endl;

	ofile.close();
	*/
}


/*
void CAMBALA::printDelayTime(double R, vector<unsigned> mode_numbers, vector<vector<double>> modal_group_velocities)
{
	string ofile_name = "delayTimeOutput_" + object_function_type_ + ".txt";
	ofstream ofile(ofile_name);
	for (unsigned ii = 0; ii < freqs.size(); ii++) {
		unsigned mnumb = mode_numbers.at(ii);
		ofile << freqs.at(ii) << "\t";
		for (unsigned jj = 0; jj < mnumb; jj++)
			ofile << R / modal_group_velocities[ii][jj] << "\t";
		ofile << endl;
	}
	ofile.close();
}
*/

double CAMBALA::directPointCalc( Point point )
{
	isTimeDelayPrinting = true;
	return fillDataComputeResidual(point);
}

void CAMBALA::Solve(const Scenario& c)
{
	for (auto cur_h: getDimGrid(c.h_))
	{

		Solver* s;
		/*
		if (launch_type == "hillclimbing")
			s = new HillClimbing;
		else if (launch_type == "bbox")
			s = new BBox;
		else if (launch_type == "bruteforce")
		*/
			s = new BruteForce;

		Model m;
		m.depths = makeDepths(cur_h, c.H_, c.depthsDim_, c.cwDim);
		m.freqs = c.freqs_;
		m.weight_coeffs = c.spmag_;
		res_calc_sel_.LoadModel(m);

		s->SetResidualCalculator(res_calc_sel_);
		s->LoadSearchSpaceDims(c.ssd_);
		s->Solve();
		Point point = s->getBestPoint();
		delete s;

		if (point.residual < record_point_.residual)
		{
			record_point_ = point;
			if (verbosity > 0)
				PrintPoint(record_point_);
		}
		cout << "Processed " << (cur_h-h.l)/h.s << " out of " << h.r/h.s << " h (max depths)" << endl;
	}
}


vector<double> makeDepths(double h, double H, const vector <Dim>& d, const vector<Dim>& cw)
{
	vector <double> depths;
	size_t n_layers_w = cw.size();
	double layer_thickness_w = h / n_layers_w;
	for (unsigned jj = 1; jj <= n_layers_w; jj+)
		depths.push_back(layer_thickness_w*jj);
	depths.push_back(H);
}

