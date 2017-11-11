#include <iostream>
#include <utility>
#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <sspemdd_sequential.h>
#include "acoustics_bottom_R_uniform.hpp"

int main(int argc, char *argv[]) {
    const int n = 3;
    const double eps = 1e-5;
    const double vref = 0.0130504;

    sspemdd_sequential sspemdd_seq;
    sspemdd_seq.verbosity = 0;
    sspemdd_seq.readScenario("311_bottom_R_uniform260.txt");
    std::vector<std::pair<double, double>> vPair;
    vPair.push_back(std::make_pair(sspemdd_seq.R1, sspemdd_seq.R2));
    vPair.push_back(std::make_pair(sspemdd_seq.rhob1, sspemdd_seq.rhob2));
    vPair.push_back(std::make_pair(sspemdd_seq.cb1, sspemdd_seq.cb2));
    ACOUSTIC::AcousticsHomogWaterUniformProblemFactory ahwpf(vPair);
    COMPI::MPProblem<double>* prob = ahwpf.getProblem();
    double x[n] = {7016, 1.7, 1735};
    double v = prob->mObjectives.at(0)->func(x);
    std::cout << "v = " << v << "\n";
    std::cout << "vref = " << vref << "\n";
    std::cout << "eps = " << eps << "\n";
    SG_ASSERT(SGABS(v - vref) < eps);
    return 0;
}