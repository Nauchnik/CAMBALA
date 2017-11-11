/* 
 * File:   acoustics_homog_water.hpp
 * Author: Oleg Zaikin
 *
 * Created on Devember 9, 2016, 8:12 PM
 */

#ifndef ACOUSTICS_HOMOG_WATER_UNIFORM_HPP
#define ACOUSTICS_HOMOG_WATER_UNIFORM_HPP

#include <mpproblem.hpp>
#include <box/box.hpp>
#include <sspemdd_sequential.h>
#include <vector>
#include <utility>

namespace ACOUSTIC {

    class AcousticsHomogWaterUniformObjective : public COMPI::Functor <double> {
    public:

        AcousticsHomogWaterUniformObjective(int n) : mN(n) {
        }

        double func(const double* x) {
			sspemdd_sequential sspemdd_seq;
			sspemdd_seq.verbosity = 0;
			sspemdd_seq.readScenario("311_bottom_R_uniform260.txt");
			sspemdd_seq.readInputDataFromFiles();
            sspemdd_seq.init();
            search_space_point cur_point;
			// const values
			cur_point.tau = sspemdd_seq.tau1;
			cur_point.cws = sspemdd_seq.cw1_arr;
			// variable values
			cur_point.R = x[0];
            cur_point.rhob = x[1];
            cur_point.cb = x[2];
            return sspemdd_seq.fill_data_compute_residual(cur_point);
        }

    private:
        int mN;
    };

    /**
     * Factory to produce instances of DeJong optimization problem
     */
    class AcousticsHomogWaterUniformProblemFactory {
    public:

        /**
         * n = 3
         * x0 = radius (distance)
         * x1 = bottom density 
         * x2 = bottom sound speed 
         */
        AcousticsHomogWaterUniformProblemFactory(const std::vector<std::pair<double, double>> &vPair) : mVPair(vPair) {
			mN = vPair.size();
        }

        COMPI::MPProblem<double>* getProblem() const {
            COMPI::MPProblem<double>* prob = new COMPI::MPProblem<double>();
            prob->mVarTypes.assign(mN, COMPI::MPProblem<double>::VariableTypes::GENERIC);
            prob->mObjectives.push_back(std::make_shared<AcousticsHomogWaterUniformObjective>(mN));
            prob->mBox = new snowgoose::Box<double>(mN);
			for (int i = 0; i < mN; i++) {
				prob->mBox->mA[i] = mVPair[i].first;
				prob->mBox->mB[i] = mVPair[i].second;
			}
            return prob;
        }

    private:
        int mN;
		std::vector<std::pair<double, double>> mVPair;
    };
}

#endif /* ACOUSTICS_HOMOG_WATER_UNIFORM_HPP */

