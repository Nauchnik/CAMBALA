#include "solvers/ils.h"
#include <iostream>

void ILS::Solve()
{
	for (unsigned i = 0; i < depths_vec_.size(); i++)
	{
		init(depths_vec_[i]);
		SearchSpace search_space = loadValuesToSearchSpaceVariables(search_space_dims);
		//CAMBALA_seq.findGlobalMinBruteForce(depths_vec[i]);
		findLocalMinHillClimbing(depths_vec_[i]);
		cout << "Processed " << i + 1 << " out of " << depths_vec_.size() << " depths" << endl;
	}
}
