/*****************************************************************************************
// SSPEMDD: Sound Speed Profile Estimator from Modal Delay Data -- Copyright(c) 2015-2017
// Pavel Petrov (Il'ichev Pacific Oceanological Institute of FEB RAS),
// Oleg Zaikin (Matrosov Institute for System Dynamics and Control Theory of SB RAS)
*****************************************************************************************/

#include "sspemdd_sequential.h"

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <cmath>

using namespace std;

int main()
{
	string scenarioFileName = "test_scenario.txt";

    //объявить экземпляр
	sspemdd_sequential sspemdd;
	sspemdd.verbosity = 0;
    //читать сценарий
	//sspemdd.readScenario(scenarioFileName);
    //создание массивов c1s, c2s,...
    //    sspemdd.init();

	vector<double> depths, c1s, c2s, rhos, freqs, zr, Rr;
	vector<unsigned> Ns_points;
        double d0, fc, zc, rc;


    unsigned ppm = 2;

    //string ProfileFName = "ssp6m.txt";
    //sspemdd.load_profile_deep_water(ProfileFName,ppm,depths,c1s,c2s,rhos, Ns_points);

    //
    string ProfileFName = "test_pek.txt";

    sspemdd.load_layers_data(ProfileFName,depths,c1s,c2s,rhos, Ns_points);

    cout << " MEDIA PARAMS " << endl;
    for (unsigned jj=0; jj<depths.size(); jj++){
        if (jj==0) { d0 = 0; } else { d0 = depths.at(jj-1); }
        cout << d0 << "-" << depths.at(jj) << "\t" << c1s.at(jj) << "\t" << c2s.at(jj) << "\t" << rhos.at(jj) << "\t" << Ns_points.at(jj) << endl;
    }
    cout << endl;
    cout << " ---- " << endl;

    freqs.clear();

    ifstream Myfile;
    Myfile.open("freqs.txt");

	if (!Myfile.is_open()) {
		cerr << "freqs.txt wasn't opened" << endl;
		return 1;
	}

    while (!(Myfile.eof()))
	{
		Myfile >> fc;
        freqs.push_back(fc);
		//cout << cc << "--" << dc << endl;
	}
	Myfile.close();


	cout << " FREQUENCIES " << endl;
    cout << freqs.at(0) << "--" << freqs.back() << endl;

    Myfile.open("srcrec.txt");

	if (!Myfile.is_open()) {
		cerr << "srcrec.txt wasn't opened" << endl;
		return 1;
	}

    while (!(Myfile.eof()))
	{
		Myfile >> rc;
		Myfile >> zc;
        zr.push_back(zc);
        Rr.push_back(rc);
		//cout << cc << "--" << dc << endl;
	}
	Myfile.close();

	cout << " SOURCE " << endl;
    cout << "R=" << Rr.at(0) << "; zs=" << zr.at(0) << endl;
    cout << " RECEIVERS " << endl;
    for (unsigned jj=1; jj<zr.size(); jj++){
        cout << "Rr" << jj << "=" << Rr.at(jj) << "\t"  << "zr" << jj << "=" << zr.at(jj) << endl;
    }
    cout << " ---- " << endl;

    // ### INITIALIZATION OF VARIABLES

    unsigned int ordRich = 3;
	double omeg1;
	cout.precision(12);

	// ### WAVENUMBERS BLOCK ###


	vector<double> result_vec;
	ofstream ofile("output.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {
		omeg1 = 2 * LOCAL_M_PI*freqs[i];
		result_vec = sspemdd.compute_wnumbers_extrap2(omeg1, depths, c1s, c2s, rhos, Ns_points, 0.6, ordRich);
		for (unsigned j = 0; j < result_vec.size(); j++)
			ofile << result_vec.at(j) << "\t";
		ofile << endl;


		cout << " f= " << omeg1/(2 * LOCAL_M_PI) << endl;
        cout << " ---- " << endl;
		cout << " k_j " << endl;
        for (unsigned jj=0; jj<result_vec.size(); jj++){
            cout << result_vec.at(jj) << "\t";
        }
        cout << endl;
        cout << " ---- " << endl;

	}
	ofile.close();


    // ### MODE FUNCTIONS BLOCK ###
    sspemdd.compute_all_mfunctions(omeg1, depths, c1s, c2s, rhos, Ns_points, result_vec );



    // ### GROUP VELOCITIES BLOCK ###

    vector<vector<double>> mgvs;
    vector<unsigned> mnumbs;
    sspemdd.compute_modal_grop_velocities2(freqs, 0.005, depths, c1s, c2s, rhos, Ns_points,-1.0,ordRich, mgvs,mnumbs);

    ofile.open("mgvs2.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {

		for (unsigned j = 0; j < mnumbs.at(i); j++)
			ofile << mgvs.at(i).at(j) << "\t";
		ofile << endl;


		cout << " f= " << freqs.at(i) << endl;
        cout << " ---- " << endl;
		cout << " vg_j " << endl;
        for (unsigned jj=0; jj<mnumbs.at(i); jj++){
            cout << mgvs.at(i).at(jj) << "\t";
        }
        cout << endl;
        cout << " ---- " << endl;

	}
	ofile.close();

    sspemdd.compute_modal_grop_velocities(freqs, 0.005, depths, c1s, c2s, rhos, Ns_points,-1.0,ordRich, mgvs,mnumbs);

    ofile.open("mgvs.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {

		for (unsigned j = 0; j < mnumbs.at(i); j++)
			ofile << mgvs.at(i).at(j) << "\t";
		ofile << endl;


		cout << " f= " << freqs.at(i) << endl;
        cout << " ---- " << endl;
		cout << " vg_j " << endl;
        for (unsigned jj=0; jj<mnumbs.at(i); jj++){
            cout << mgvs.at(i).at(jj) << "\t";
        }
        cout << endl;
        cout << " ---- " << endl;

	}
	ofile.close();



    // ### SOUND PRESSURE FIELD COMPUTATION BLOCK ###
/*

    double freq_threshold = 0.5;

    ofstream Pfile("PHelm.txt");
    vector<complex<double>> PHelm;
    for (unsigned ff=0; ff<freqs.size(); ff++ ){
        cout << "f=" << freqs.at(ff) << ": " << endl;

        if (freqs.at(ff)>freq_threshold) {
            PHelm = sspemdd.compute_cpl_pressure(freqs.at(ff),depths,c1s,c2s,rhos,Ns_points,Rr,zr,1/sqrt(2),ordRich);

            Pfile << freqs.at(ff);

            for (unsigned jj=0; jj<PHelm.size(); jj++ ){

                Pfile << "\t" << PHelm.at(jj).real() << "  " << PHelm.at(jj).imag();

            }
            Pfile << endl;
        } else {

            for (unsigned jj=0; jj<Rr.size(); jj++ ){
                Pfile << "\t" << 0 << "  " << 0;
            }
            Pfile << endl;
        }

    }
    Pfile.close();
*/




	cout << "test passed" << endl;

    complex<double> U(0.0,LOCAL_M_PI/4);
    cout << sqrt(2)*exp(U) << endl;
    cout << sqrt(2)*exp(U) + Iu << endl;
    U = (complex<double>(4.5,2.7));


	return 0;
}
