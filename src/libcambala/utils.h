#ifndef UTILS_H_
#define UTILS_H_
#include <vector>
#include <string>

using namespace std;

template< typename T > bool next_cartesian(vector<T> &vii, vector<int> &index_arr, T &cur_vi);
vector <vector <double>> DoubleVecVecRead(string filename);
vector <double> DoubleVecVecGetFirstColumn(vector <vector <double>> vecvec);
vector <vector <double>> DoubleVecVecGetOtherColumns(vector <vector <double>> vecvec);




#endif
