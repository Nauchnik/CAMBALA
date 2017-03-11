#include <vector>
#include <cstdlib>
//#include <cmath>
//#include <cfloat>
#include <algorithm>
#include <limits>
#include "linalg.h"
using EvVec = std::vector <float>;

struct BidiagMat
{
	std::vector <float> md;
	std::vector <float> sd;
	BidiagMat (int sz = 0)
	{
		md.resize(sz,0);
		sd.resize(sz,0);
	}
};


struct Interval
{
	float l;
	float u;
};


Interval computeGerschgorin(const BidiagMat& mat)
{
	float lg = std::numeric_limits<float>::max();
	float ug = std::numeric_limits<float>::lowest();

	auto& d = mat.md;
	auto& s = mat.sd;
	size_t mat_size = mat.md.size();

	// compute bounds
	for( unsigned int i = 1; i < (mat_size - 1); ++i)
	{
		// sum over the absolute values of all elements of row i
		float sum_abs_ni = std::fabs( s[i-1]) + std::fabs( s[i]);

		lg = std::min( lg, d[i] - sum_abs_ni);
		ug = std::max( ug, d[i] + sum_abs_ni);
	}  

	// first and last row, only one superdiagonal element

	// first row
	lg = std::min( lg, d[0] - std::fabs( s[0]));
	ug = std::max( ug, d[0] + std::fabs( s[0]));

	// last row
	lg = std::min( lg, d[mat_size-1] - std::fabs( s[mat_size-2]));
	ug = std::max( ug, d[mat_size-1] + std::fabs( s[mat_size-2]));

	// increase interval to avoid side effects of fp arithmetic
	float bnorm = std::max( std::fabs( ug), std::fabs( lg));

	// these values depend on the implementation of floating count that is
	// employed in the following
	const float EPSILON = std::numeric_limits<float>::epsilon();
	float psi_0 = 11 * EPSILON * bnorm;
	float psi_n = 11 * EPSILON * bnorm;

	lg = lg - bnorm * 2 * mat_size * EPSILON - psi_0;
	ug = ug + bnorm * 2 * mat_size * EPSILON + psi_n;

	ug = std::max( lg, ug);
	return (Interval{lg,ug});
}

EvVec alglibEv(const BidiagMat& mat, Interval eiv)
{
	const int mat_size = mat.md.size();

	alglib::real_1d_array md, sd;
	md.setlength (mat_size);
	sd.setlength (mat_size - 1);

	for (int i = 0; i < mat_size; ++i)
		md[i] = mat.md[i];

	for (int i = 0; i < mat_size-1; ++i)
		sd[i] = mat.sd[i];

	alglib::ae_int_t eigen_count = 0;
	alglib::real_2d_array eigenvectors;
	alglib::smatrixtdevdr (md, sd, mat_size, 0, (double)eiv.l, (double)eiv.u, 
			eigen_count, eigenvectors);

	EvVec out;
	for (int i = 0; i<eigen_count; ++i)
		out.push_back(float(md[i]));
	return std::move(out);
}

BidiagMat generateTestMatrix(const int mat_size)
{
	BidiagMat mat(mat_size);
	for (size_t i = 0; i < mat.md.size(); ++i)
	{
		mat.md[i] = i*i*i;
		mat.sd[i] = i*i;
	}
	return std::move(mat);
}

int main (const int argc, const char** argv)
{

	const int mat_size = std::atoi(argv[1]);

	const BidiagMat mat = generateTestMatrix(mat_size);

	auto gg = computeGerschgorin(mat);
	EvVec ev_control = alglibEv(mat, gg);
	//EvVec ev_cpu = computeEv(mat, gg);

	//CompareEvs(ev_cpu, ev_control);
}

