/* 
 * File:   test_hydro_R_weighted.cpp
 * Author: Oleg Zaikin
 *
 * Created on March 1, 2017, 11:27 PM
 */

#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include "acoustics_hydro_R_weighted.hpp"

int main(int argc, char *argv[])
{
    const int n = 5;
    const double eps = 1e-10;
    const double vref = 3.15003e-07;
	sspemdd_sequential sspemdd_seq;
	sspemdd_seq.verbosity = 0;
	sspemdd_seq.readScenario("39_hydro_R_weighted260.txt");
	std::vector<std::pair<double, double>> vPair;
	vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
	vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
	vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHydroRWeightedProblemFactory ahrwpf(vPair);
    COMPI::MPProblem<double>* prob = ahrwpf.getProblem();
    double x[n] = { 6981, 1480, 1485, 1450, 1470 };
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = "    << v    << "\n";
	std::cout << "vref = " << vref << "\n";
	std::cout << "eps = "  << eps  << "\n";
    SG_ASSERT(SGABS(v - vref) < eps);
    return 0;
}