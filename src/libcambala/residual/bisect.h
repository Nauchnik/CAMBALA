#include "assert.h"
#include "stdlib.h"
#include <iostream>
#include <math.h>

template <typename ftype> 
struct Interval
{
	ftype ll;
	ftype rl;
};

struct Interval_e
{
	int lc;
	int rc;
};

template<class T>
#ifdef __CUDACC__
__host__  __device__
#endif 
T min( const T& lhs, const T& rhs) { return (lhs < rhs) ? lhs : rhs; }

template<class T>
#ifdef __CUDACC__
__host__  __device__
#endif
T max( const T& lhs, const T& rhs) { return (lhs < rhs) ? rhs : lhs; }

template <typename ftype> 
inline ftype sign_f( const ftype& val) { return (val < 0.0f) ? -1.0f : 1.0f; }

template <typename ftype> 
inline void StoreEv (ftype eiv, ftype* eigenvalues, int& eigenvalues_c) { eigenvalues[eigenvalues_c++] = eiv; }

template <typename ftype> 
inline void PushInterval (const Interval<ftype> ival, const Interval_e& ival_e, Interval<ftype>* intervals, Interval_e* intervals_e, int& intervals_c)
{
	intervals[intervals_c] = ival;
	intervals_e[intervals_c] = ival_e;
	++intervals_c;	
}



#if FTYPE==float
  static const float precision = 0.000001;
#else
  #if FTYPE==double
    static const double precision = 0.00000001;
  #else
    #error Unknown primary storage type: FTYPE
  #endif
#endif


#define MIN_ABS_INTERVAL (ftype)5.0e-37

template <typename ftype> 
int computeSmallerEigenvals (const ftype x, const ftype* md, const ftype* sd, const int mat_size)
{
	//assert (sd[0]==0.0f);
	// perform (optimized) Gaussian elimination to determine the number
	// of eigenvalues that are smaller than x
	const ftype* sd_shifted = sd-1;
	ftype  delta = 1.0f;
	int count = 0;
	delta = md[0] - x;
	count += (delta < 0);
	for (int i = 1; i < mat_size; ++i)
	{
		delta = md[i] - x - (sd_shifted[i] * sd_shifted[i]) / delta;
		count += (delta < 0);
	}
	assert (count>=0);
	assert (count<=mat_size);
	return count;

}

template <typename ftype> 
ftype midpoint(const ftype& left, const ftype& right)
{
	// try to avoid overflow
	ftype mid;
	if (sign_f (left) == sign_f (right))
		mid = left + (right - left) * 0.5f;
	else
		mid = (left + right) * 0.5f;

	return mid;
}
template <typename ftype> 
inline ftype midpoint(const Interval<ftype>& curint) { return midpoint (curint.ll, curint.rl);};

template <typename ftype> 
bool intervalConverged (const Interval<ftype>& c)
{
	ftype t0 = fabs (c.rl - c.ll);
	ftype t1 = max (fabs(c.ll), fabs(c.rl)) * precision;
	return (t0 <= max (MIN_ABS_INTERVAL, t1));

}


template <typename ftype> 
void PushOrStoreInterval (const Interval<ftype> curint, const Interval_e curint_e,
		Interval<ftype>* intervals, Interval_e* intervals_e, int& intervals_c,
		ftype* eigenvalues, int& eigenvalues_c)
{
	if (intervalConverged (curint))
		for (int i=0; i<(curint_e.rc-curint_e.lc); ++i)
			StoreEv (midpoint (curint), eigenvalues, eigenvalues_c);
	else
		PushInterval (curint, curint_e, intervals, intervals_e, intervals_c);
	//std::cout <<  "interval  " <<  curint.ll << " " << curint.rl << " " << intervalConverged(curint) << "\n" << std::flush;
}

template <typename ftype> 
int bisectGPU (const ftype* md, const ftype* sd, const int mat_size,
		const ftype ll, const ftype rl,
		ftype* eigenvalues)
{
	assert(ll <= rl);
	int eigenvalues_c = 0; // Number of eigenvalues found

	// We can't put these into a single structure, since CUDA will
	// drop it from registers into local memory
	Interval<ftype> intervals[MAX_WNUMS]; // Intervals borders
	Interval_e intervals_e[MAX_WNUMS]; // Eigenvalue counts for intervals borders
	int intervals_c = 0; // Intervals stack height

	Interval_e toplevel_e = Interval_e {
				computeSmallerEigenvals (ll, md, sd, mat_size),
				computeSmallerEigenvals (rl, md, sd, mat_size)};
	PushOrStoreInterval (Interval<ftype> {ll, rl}, toplevel_e,  intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);

	while (intervals_c > 0)
	{
		// Pop intervals stack
		--intervals_c;
		Interval_e curint_e = intervals_e[intervals_c];
		Interval<ftype> curint = intervals[intervals_c];
		int lc = curint_e.lc;
		int rc = curint_e.rc;
		//std::cout << mat_size << " " << intervals_c << " Lim: " << curint.ll << " " << curint.rl << "lrc " << lc << " " << rc << "\n" << std::flush;

		ftype m = midpoint (curint);
		int mc = computeSmallerEigenvals (m, md, sd, mat_size);
		if (mc-lc > 0)
			PushOrStoreInterval (Interval<ftype> {curint.ll, m}, Interval_e {lc, mc}, intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);
		if (rc-mc > 0)
			PushOrStoreInterval (Interval<ftype> {m, curint.rl}, Interval_e {mc, rc}, intervals, intervals_e, intervals_c, eigenvalues, eigenvalues_c);
	}
	return eigenvalues_c;
}
