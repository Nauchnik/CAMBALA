#ifndef BISECT_CPU_H_
#define BISECT_CPU_H_
struct Interval
{
	float ll;
	float rl;
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

inline float sign_f( const float& val) { return (val < 0.0f) ? -1.0f : 1.0f; }

int computeSmallerEigenvals (const float x, const float* md, const float* sd, const int mat_size);
float midpoint(const float& left, const float& right);
inline float midpoint(const Interval& curint) { return midpoint (curint.ll, curint.rl);};
bool intervalConverged (const Interval& c);
void PushOrStoreInterval (const Interval curint, Interval* intervals, int& intervals_c, float* eigenvalues, int& eigenvalues_c);
int bisectCpu (const float* md, const float* sd, const int mat_size, const float ll, const float rl, float* eigenvalues);
inline void StoreEv (float eiv, float* eigenvalues, int& eigenvalues_c) { eigenvalues[eigenvalues_c++] = eiv; }
inline void PushInterval (const Interval ival, const Interval_e& ival_e, Interval* intervals, Interval_e* intervals_e, int& intervals_c)
{
	intervals[intervals_c] = ival;
	intervals_e[intervals_c] = ival_e;
	++intervals_c;	
}
#endif
