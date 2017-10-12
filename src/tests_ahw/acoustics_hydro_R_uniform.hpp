/* 
 * File:   acoustics_hydro_R_unifrom.hpp
 * Author: Oleg Zaikin
 *
 * Created on March 1, 2017, 11:57 PM
 */

#ifndef ACOUSTICS_HYDRO_R_UNIFORM_HPP
#define ACOUSTICS_HYDRO_R_UNIFORM_HPP

#include <mpproblem.hpp>
#include <box/box.hpp>
#include <sspemdd_sequential.h>
#include <vector>
#include <utility>

namespace ACOUSTIC {

    class AcousticsHydroRUniformObjective : public COMPI::Functor <double> {
    public:

		AcousticsHydroRUniformObjective(int n) : mN(n) {
        }

        double func(const double* x) {
			sspemdd_sequential sspemdd_seq;
			sspemdd_seq.verbosity = 0;
			sspemdd_seq.readScenario("310_hydro_R_uniform260.txt");
			sspemdd_seq.readInputDataFromFiles();
            sspemdd_seq.init();
            search_space_point cur_point;
			// const values
			cur_point.tau  = sspemdd_seq.tau1;
			cur_point.cws  = sspemdd_seq.cw1_arr;
			cur_point.cb   = sspemdd_seq.cb1;
			cur_point.rhob = sspemdd_seq.rhob1;
			// variable values
			cur_point.R      = x[0];
			cur_point.cws[1] = x[1];
			cur_point.cws[2] = x[2];
			cur_point.cws[3] = x[3];
			cur_point.cws[4] = x[4];
			
            return sspemdd_seq.fill_data_compute_residual(cur_point);
        }

    private:
        int mN;
    };

    /**
     * Factory to produce instances of DeJong optimization problem
     */
    class AcousticsHydroRUniformProblemFactory {
    public:

        /**
         * n = 5
         * x0 = radius (distance)
         * x1-x4 = hydro sound speed 
         */
		AcousticsHydroRUniformProblemFactory(const std::vector<std::pair<double, double>> &vPair) : mVPair(vPair) {
			mN = vPair.size();
        }

        COMPI::MPProblem<double>* getProblem() const {
            COMPI::MPProblem<double>* prob = new COMPI::MPProblem<double>();
            prob->mVarTypes.assign(mN, COMPI::MPProblem<double>::VariableTypes::GENERIC);
            prob->mObjectives.push_back(std::make_shared<AcousticsHydroRUniformObjective>(mN));
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

#endif /* ACOUSTICS_HYDRO_R_UNIFORM_HPP */

