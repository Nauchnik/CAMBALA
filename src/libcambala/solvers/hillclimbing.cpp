#include "solvers/hillclimbing.h"
#include <iostream>

Point CAMBALA_sequential::getNonRandomStartPoint()
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws

	Point point;

	point.cb = (cb1 == cb2) ? cb1 : (cb2 - cb1) / 2;
	point.rhob = (rhob1 == rhob2) ? rhob1 : (rhob2 - rhob1) / 2;
	point.R = (R1 == R2) ? R1 : (R2 - R1) / 2;
	point.tau = (tau1 == tau2) ? tau1 : (tau2 - tau1) / 2;

	point.cws.resize(cw1_arr.size());
	for ( unsigned i = 0; i < point.cws.size(); i++ ) {
		if (cw1_arr[i] == cw2_arr[i])
			point.cws[i] = cw1_arr[i];
		else if ((START_CW_VALUE >= cw1_arr[i]) && (START_CW_VALUE <= cw2_arr[i])) // if cws should be modified
			point.cws[i] = START_CW_VALUE;
		else
			point.cws[i] = (cw2_arr[i] - cw1_arr[i]) / 2;
	}

	return point;
}
HillClimbing::Solve()
{
	if (verbosity > 1)
		cout << "findLocalMinHillClimbing" << endl;
	unsigned u_val = 1;
	for (unsigned i = 0; i < depths_.size() - 2; i++)
		u_val *= (unsigned)depths_[i];
	srand(u_val);

	// choose random point in the search space
	/*for (unsigned i = 0; i < search_space.size(); i++) // i stands for variable_index
	local_record_point_indexes[i] = rand() % search_space[i].size(); // get random index
	cur_point_indexes = local_record_point_indexes;*/

	Point local_record_point = getNonRandomStartPoint(depths_);
	vector<unsigned> local_record_point_indexes = fromPointToPointIndexes( local_record_point );
	fillDataComputeResidual(local_record_point); // calculated residual is written to cur_point

	if (point.residual < record_point.residual)
	{
		record_point = point;
		if (verbosity > 0)
			PrintPoint(record_point);
	}

	bool isCheckRequired = false;
	for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
		if (search_space[i].size() > 1) {
			isCheckRequired = true;
			break;
		}
	}
	if ( (!isCheckRequired) && (verbosity > 0) ) {
		cout << "1 element in a search space, fast exit" << endl;
		return local_record_point;
	}

	checked_points.reserve(N_total);
	checked_points.push_back(local_record_point);
	unsigned skipped_points = 0;
	bool isContinueDimension;
	vector<unsigned> cur_point_indexes;
	Point cur_point;
	// launch iterations of hill climbing
	for (unsigned run_index = 0; run_index < iterated_local_search_runs; run_index++) {
		bool isLocalMin;
		cout << endl;
		cout << "iteration " << run_index << " of ILS" << endl;
		do { // do while local min not reached
			isLocalMin = true; // if changing of every variable will not lead to a record updata, then a local min reached
			for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
				if (search_space[i].size() == 1) {
					//cout << "one value of a variable, skip it" << endl;
					continue;
				}
				//cout << "variable_index " << variable_index << endl;
				cur_point_indexes = local_record_point_indexes;
				vector<unsigned> point_indexes_before_increase = cur_point_indexes;
				unsigned index_from = cur_point_indexes[i]; // don't check index twice
				if (verbosity > 0)
					cout << "index_from " << index_from << endl;
				bool isDecreaseTurn = false;
				bool isTriggerIncToDec = false;
				do { // change value of a variabble while it leads to updating of a record
					double old_record_residual = local_record_point.residual;
					if (isDecreaseTurn) {
						if (isTriggerIncToDec) { // move to a point before increase
							cur_point_indexes = point_indexes_before_increase;
							isTriggerIncToDec = false;
						}
						if (cur_point_indexes[i] == 0)
							cur_point_indexes[i] = search_space[i].size() - 1;
						else
							cur_point_indexes[i]--;
					}
					else {
						cur_point_indexes[i]++;
						if (cur_point_indexes[i] == search_space[i].size())
							cur_point_indexes[i] = 0;
					}
					if (cur_point_indexes[i] == index_from) {
						if (verbosity > 0)
							cout << "cur_point_indexes[i] == index_from. Break iteration." << endl;
						break;
					}
					if (verbosity > 0) {
						cout << "checking index " << cur_point_indexes[i] <<
							", max index " << search_space[i].size() - 1 << endl;
						cout << "cur_point_indexes" << endl;
						for (unsigned j = 0; j < cur_point_indexes.size(); j++)
							cout << cur_point_indexes[j] << " ";
						cout << endl;
					}
					cur_point = fromPointIndexesToPoint(cur_point_indexes, depths_);
					if (find(checked_points.begin(), checked_points.end(), cur_point) != checked_points.end()) {
						skipped_points++;
						continue;
					}
					double d_val = fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
					checked_points.push_back(cur_point);
					isContinueDimension = false;
					if (d_val < old_record_residual) { // new record was found
						local_record_point.residual = d_val;
						local_record_point_indexes = cur_point_indexes;
						isLocalMin = false;
						isContinueDimension = true;
					}
					if ((!isContinueDimension) && (!isDecreaseTurn)) {
						isDecreaseTurn = true; // try to decrease current value
						if (verbosity > 0)
							cout << "isDecreaseTurn " << isDecreaseTurn << endl;
						isContinueDimension = true;
						isTriggerIncToDec = true;
						continue;
					}
				} while (isContinueDimension);
			}
		} while (!isLocalMin);

		//cout << endl << "*** local minimum in hill climbing" << endl;
		//cout << "local record of residual " << record_point.residual << endl;
		//cout << "-----" << endl;
		if (verbosity > 0)
			cout << "new random cur_point_indexes : " << endl;

		for(;;) {
			// prmutate current global minimum point to obtain a new start point
			for (unsigned i = 0; i < search_space.size(); i++) {
				if (search_space[i].size() == 1) {
					cur_point_indexes[i] = 0;
					if (verbosity > 0)
						cout << cur_point_indexes[i] << " ";
					continue;
				}
				unsigned rand_numb = rand();
				if (rand_numb % 3 == 0)
					cur_point_indexes[i] = local_record_point_indexes[i];
				else
					cur_point_indexes[i] = (rand_numb % search_space[i].size());
			}
			cur_point = fromPointIndexesToPoint(cur_point_indexes, depths_);
			if (find(checked_points.begin(), checked_points.end(), cur_point) == checked_points.end())
				break;
		}
		if (verbosity > 0) {
			cout << "new random point" << endl;
			for (unsigned j = 0; j < cur_point_indexes.size(); j++)
				cout << cur_point_indexes[j] << " ";
			cout << endl;
		}

		fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
		checked_points.push_back(cur_point);
		// new start point
		local_record_point = cur_point;
		local_record_point_indexes = cur_point_indexes;

		cout << "checked_points size " << checked_points.size() << endl;
		cout << "skipped_points " << skipped_points << endl;
		cout << "---" << endl;
	}

	recordPoint_ = local_record_point;
}
