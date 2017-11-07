//#include <iostream>
#include "cambala.h"
#include "utils.h"
#include "types.h"
#include "solvers/interface.h"
#include "solvers/discrete.h"
#include "solvers/bruteforce.h"
#include "solvers/hillclimbing.h"
#define ELPP_STL_LOGGING
#include "easylogging++.h"


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
	/*
	//isTimeDelayPrinting = true;
	return fillDataComputeResidual(point);
	*/
}

void CAMBALA::Solve(const Scenario& c)
{
	Point record_point;
	for (auto cur_h: getDimGrid(c.hDim_))
	{
		Solver* s;
		/*
		if (launch_type == "hillclimbing")
			s = new HillClimbing;
		else if (launch_type == "bbox")
			s = new BBox;
		else if (launch_type == "bruteforce")
			s = new BruteForce;
		*/
			s = new HillClimbing();
			//s = new BruteForce;

		//TODO: rewrite model init as a separate function/method
		Model m;
		m.depths = makeDepths(cur_h, c.H_, c.depthsDim_, c.ssd_.cw);
		m.freqs = c.freqs_;
		m.weight_coeffs = c.spmag_;
		m.Ns_points.push_back((unsigned)round(c.ppm_*m.depths[0]));
		m.exp_delays = c.modal_delays_;

		for (size_t i=1; i<m.depths.size(); ++i)
			m.Ns_points.push_back((unsigned)round(c.ppm_*(m.depths[i] - m.depths[i-1])));

		//LOG(DEBUG) << "c.ppm , m.depths[0] "<< c.ppm_ << " " << m.depths[0] ;
		LOG(DEBUG) << "depths freqs weight_coeffs Ns_points "
			<< m.depths.size() << " " 
			<< m.freqs.size()<< " " 
			<< m.weight_coeffs.size() << " "
			<< m.depths.size() << " "
		        << m.Ns_points.size();


		ResCalc* rc = calcs_["fast"];
		rc->LoadModel(m);

		s->SetResidualCalculator(rc);
		s->LoadSearchSpaceDims(c.ssd_);
		s->Solve();
		Point point = s->getBestPoint();
		delete s;

		if (point.residual < record_point.residual)
		{
			record_point = point;
			//if (verbosity > 0) PrintPoint(record_point_);
		}
		// TODO: make a general output method for Point!
		LOG(INFO) << "Record point found: " 
			<< record_point.residual << " Coords: "
			<< record_point.R << " "
			<< record_point.rhob << " "
			<< record_point.cb << " "
			<< record_point.tau <<  " "
			<< record_point.cws;
		//cout << "Processed " << (cur_h-h.l)/h.s << " out of " << h.r/h.s << " h (max depths)" << endl;
	}
}

vector<double> makeDepths(double h, double H, const vector <Dim>& d, const vector<Dim>& cw)
{
	LOG(DEBUG) << "makeDepths: "<< h << " " << H << " ";
	vector <double> depths;
	size_t n_layers_w = cw.size();
	double layer_thickness_w = h / n_layers_w;
	for (unsigned jj = 1; jj <= n_layers_w; ++jj) 
		depths.push_back(layer_thickness_w*jj);
	depths.push_back(H);
	return depths;
}


void CAMBALA::AddResidualCalculator(std::string name, ResCalc* rc)
{
	LOG(DEBUG) << "Adding residual calculator \""<< rc->getName() << "\"" << " as \"" << name << "\".";
	//TODO: add more semantic checks
	calcs_[name] = rc;
	//add default calculators
	if( calcs_["default"]==nullptr)
	{
		calcs_["default"] = rc;
		calcs_["precise"] = rc;
		calcs_["fast"] = rc;
		LOG(DEBUG) << " Default calc now is \"" << calcs_["default"]->getName() << "\"." ;
	}
}


CAMBALA::CAMBALA()
{
	// residual calculators names
	calcs_["default"] = nullptr;
	calcs_["precise"] = nullptr;
	calcs_["fast"] = nullptr;
}

