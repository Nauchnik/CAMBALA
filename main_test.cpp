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

int main()
{
	std::string scenarioFileName = "test_scenario.txt";

    //объявить экземпляр
	sspemdd_sequential sspemdd;
	sspemdd.verbosity = 0;
    //читать сценарий
	//sspemdd.readScenario(scenarioFileName);
    //создание массивов c1s, c2s,...
    //    sspemdd.init();

	std::vector<double> depths, c1s, c2s, rhos, freqs, zr, Rr;
	std::vector<unsigned> Ns_points;
        double d0, fc, zc, rc;


    unsigned ppm = 2;

    //std::string ProfileFName = "ssp6m.txt";
    //sspemdd.load_profile_deep_water(ProfileFName,ppm,depths,c1s,c2s,rhos, Ns_points);

    //
    std::string ProfileFName = "test_pek.txt";

    sspemdd.load_layers_data(ProfileFName,depths,c1s,c2s,rhos, Ns_points);

    std::cout << " MEDIA PARAMS " << std::endl;
    for (unsigned jj=0; jj<depths.size(); jj++){
        if (jj==0) { d0 = 0; } else { d0 = depths.at(jj-1); }
        std::cout << d0 << "-" << depths.at(jj) << "\t" << c1s.at(jj) << "\t" << c2s.at(jj) << "\t" << rhos.at(jj) << "\t" << Ns_points.at(jj) << std::endl;
    }
    std::cout << std::endl;
    std::cout << " ---- " << std::endl;

    freqs.clear();

    std::ifstream Myfile;
    Myfile.open("freqs.txt");

    while (!(Myfile.eof()))
	{
		Myfile >> fc;
        freqs.push_back(fc);
		//std::cout << cc << "--" << dc << std::endl;
	}
	Myfile.close();


	std::cout << " FREQUENCIES " << std::endl;
    std::cout << freqs.at(0) << "--" << freqs.back() << std::endl;

    Myfile.open("srcrec.txt");

    while (!(Myfile.eof()))
	{
		Myfile >> rc;
		Myfile >> zc;
        zr.push_back(zc);
        Rr.push_back(rc);
		//std::cout << cc << "--" << dc << std::endl;
	}
	Myfile.close();

	std::cout << " SOURCE " << std::endl;
    std::cout << "R=" << Rr.at(0) << "; zs=" << zr.at(0) << std::endl;
    std::cout << " RECEIVERS " << std::endl;
    for (unsigned jj=1; jj<zr.size(); jj++){
        std::cout << "Rr" << jj << "=" << Rr.at(jj) << "\t"  << "zr" << jj << "=" << zr.at(jj) << std::endl;
    }
    std::cout << " ---- " << std::endl;

    // ### INITIALIZATION OF VARIABLES

    unsigned int ordRich = 3;
	double omeg1;
	std::cout.precision(12);

	// ### WAVENUMBERS BLOCK ###


	std::vector<double> result_vec;
	std::ofstream ofile("output.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {
		omeg1 = 2 * LOCAL_M_PI*freqs[i];
		result_vec = sspemdd.compute_wnumbers_extrap2(omeg1, depths, c1s, c2s, rhos, Ns_points, 0.6, ordRich);
		for (unsigned j = 0; j < result_vec.size(); j++)
			ofile << result_vec.at(j) << "\t";
		ofile << std::endl;


		std::cout << " f= " << omeg1/(2 * LOCAL_M_PI) << std::endl;
        std::cout << " ---- " << std::endl;
		std::cout << " k_j " << std::endl;
        for (unsigned jj=0; jj<result_vec.size(); jj++){
            std::cout << result_vec.at(jj) << "\t";
        }
        std::cout << std::endl;
        std::cout << " ---- " << std::endl;

	}
	ofile.close();


    // ### MODE FUNCTIONS BLOCK ###
    sspemdd.compute_all_mfunctions(omeg1, depths, c1s, c2s, rhos, Ns_points, result_vec );



    // ### GROUP VELOCITIES BLOCK ###

    std::vector<std::vector<double>> mgvs;
    std::vector<unsigned> mnumbs;
    sspemdd.compute_modal_grop_velocities2(freqs, 0.005, depths, c1s, c2s, rhos, Ns_points,-1.0,ordRich, mgvs,mnumbs);

    ofile.open("mgvs2.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {

		for (unsigned j = 0; j < mnumbs.at(i); j++)
			ofile << mgvs.at(i).at(j) << "\t";
		ofile << std::endl;


		std::cout << " f= " << freqs.at(i) << std::endl;
        std::cout << " ---- " << std::endl;
		std::cout << " vg_j " << std::endl;
        for (unsigned jj=0; jj<mnumbs.at(i); jj++){
            std::cout << mgvs.at(i).at(jj) << "\t";
        }
        std::cout << std::endl;
        std::cout << " ---- " << std::endl;

	}
	ofile.close();

    sspemdd.compute_modal_grop_velocities(freqs, 0.005, depths, c1s, c2s, rhos, Ns_points,-1.0,ordRich, mgvs,mnumbs);

    ofile.open("mgvs.txt");
	ofile.precision(12);

	for (unsigned i = 0; i < freqs.size(); i++) {

		for (unsigned j = 0; j < mnumbs.at(i); j++)
			ofile << mgvs.at(i).at(j) << "\t";
		ofile << std::endl;


		std::cout << " f= " << freqs.at(i) << std::endl;
        std::cout << " ---- " << std::endl;
		std::cout << " vg_j " << std::endl;
        for (unsigned jj=0; jj<mnumbs.at(i); jj++){
            std::cout << mgvs.at(i).at(jj) << "\t";
        }
        std::cout << std::endl;
        std::cout << " ---- " << std::endl;

	}
	ofile.close();



    // ### SOUND PRESSURE FIELD COMPUTATION BLOCK ###
/*

    double freq_threshold = 0.5;

    std::ofstream Pfile("PHelm.txt");
    std::vector<std::complex<double>> PHelm;
    for (unsigned ff=0; ff<freqs.size(); ff++ ){
        std::cout << "f=" << freqs.at(ff) << ": " << std::endl;

        if (freqs.at(ff)>freq_threshold) {
            PHelm = sspemdd.compute_cpl_pressure(freqs.at(ff),depths,c1s,c2s,rhos,Ns_points,Rr,zr,1/sqrt(2),ordRich);

            Pfile << freqs.at(ff);

            for (unsigned jj=0; jj<PHelm.size(); jj++ ){

                Pfile << "\t" << PHelm.at(jj).real() << "  " << PHelm.at(jj).imag();

            }
            Pfile << std::endl;
        } else {

            for (unsigned jj=0; jj<Rr.size(); jj++ ){
                Pfile << "\t" << 0 << "  " << 0;
            }
            Pfile << std::endl;
        }

    }
    Pfile.close();
*/




	std::cout << "test passed" << std::endl;

    std::complex<double> U(0.0,LOCAL_M_PI/4);
    std::cout << sqrt(2)*exp(U) << std::endl;
    std::cout << sqrt(2)*exp(U) + Iu << std::endl;
    U = (std::complex<double>(4.5,2.7));


	return 0;
}
