#ifndef SOLVERS_ILS_H_
#define SOLVERS_ILS_H_
#include "solvers/interface.h"
#include <vector>

class ILS : public Solver
{
public:
	void Solve ();
	void LoadData (std::vector<std::vector<double>> depths_vec); { depths_vec_ = depths_vec; }
	Point getBestPoint() { return recordPoint_ };
private:
	std::vector<std::vector<double>> depths_vec_;
	Point recordPoint_;
}
#endif
