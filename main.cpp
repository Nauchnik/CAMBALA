// +-----------------------------------------------------------------------------------+
// | Client application for the volunteer computing project Acoustics@home             |                       
// +-----------------------------------------------------------------------------------+
// | Pacific Oceanological Institute, Institute for System Dynamics and Control Theory |   
// +-----------------------------------------------------------------------------------+
// | Authors: Pavel Petrov, Oleg Zaikin                                                |
// +-----------------------------------------------------------------------------------+

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include "linalg.h"

using namespace std;

// layer data 
struct layer
{
	vector<double> zend; // vector of finite depths
	double cbeg;         // velocity of a sound at the beginning of a layer
	double cend;         // velocity of a sound at the end of a layer
	double dbeg;         // density at the beginning of a layer
	double dend;         // density at the end of a layer
};

// data of corresponding function output 
struct eigenfunctions_velocities
{
	vector<double> eigenfunctions;
	vector<double> velocities;
};

// prototypes of functions
vector<double> calc_chanel_mods_wave_numbers( double &depth_step, double &freq, vector<double> &sound_velocity,
											  vector<double> &density, vector<double> &point_indexes );
vector<double> calc_wave_numbers_richardson( layer  &cur_layer, double &initial_depth_step, double &freq );
eigenfunctions_velocities calc_eigenfunctions_velocities( layer  &cur_layer, double &initial_depth_step,   
														  double &freq, double &wave_number );
int main( int argc, char **argv )
{
	layer  cur_layer;
	double initial_depth_step;
	double freq;
	double wave_number;

	// test filling
	initial_depth_step = 2.0;
	freq = 0.1;
	cur_layer.cbeg = 0.2;
	cur_layer.cend = 0.8;
	cur_layer.dbeg = 0.5;
	cur_layer.dend = 0.9;
	for ( unsigned i=0; i < 5; i++ )
		cur_layer.zend.push_back( i );
	//
	
	vector<double> cur_wave_numbers;
	eigenfunctions_velocities cur_eigenfunc_veloc;
	cur_wave_numbers = calc_wave_numbers_richardson( cur_layer, initial_depth_step, freq );
	cout << "calc_wave_numbers_richardson() done" << endl;
	cur_eigenfunc_veloc = calc_eigenfunctions_velocities( cur_layer, initial_depth_step, freq, wave_number );
	cout << "calc_eigenfunctions_velocities() done" << endl;
	
	return 0;
}

/* function for calculating wave number of channel mods
Функция для расчета волновых чисел канальных мод
Вход:
dz -- шаг по глубине;
f -- частота звука;
вектор c = (ci); -- скорости звука с шагом dz; Длина = N. Значения N -- "большие".
вектор d = (di); -- плотности; Длина N.
вектор m = (mj); -- индексы точек, где параметры среды терпят разрыв 
                     (точки, где заканчивается один слой и начинается другой). Длина M <=10.
Работа: сформировать трехдиагональные матрицы, запустить солвер собственных значений для отрезка [omega/cmax omega/cmin].
Выход: собственные значения kj^2 акустической спектральной задачи */
vector<double> calc_chanel_mods_wave_numbers( double &depth_step, 
										      double &freq, // sound frequency
											  vector<double> &sound_velocity,
											  vector<double> &density,
											  vector<double> &point_indexes )
{
	vector<double> spectr_problem_eigenvalues;
	
	// make tridiagonal matrix and find its eigenvalues in given interval
	// ...
	stringstream sstream;
	// calculate eigenvalues and eigenvectors
	int n=2000;
	double from = 0, to = 0.01; // interval for eigenvalues
	
	sstream << "n : " << n << endl;
	alglib::real_2d_array A, eigenvectors; // V - собств вектор
	alglib::real_1d_array eigenvalues; // Lm -собств знач
	A.setlength(n,n);
	eigenvectors.setlength(n,n);
	eigenvalues.setlength(n);
	
	// fill matrix by zeros
	for ( int i=0; i < n; i++ )
		for ( int j=0; j < n; j++ )
			A[i][j] = 0.0;
	
	// fill tridiagonal matrix
	// make Dirichlet case matrix
	// На главной диагонали стоят числа -2/(h^2), на под- и над- диагоналях стоят числа 1/(h^2).
	double h = 2.0;
	for ( int i=0; i < n; i++ ) {
		A[i][i] = -2/pow(h,2);
		if ( i != n-1 ) {
			A[i+1][i] = 1/pow(h,2);
			A[i][i+1] = 1/pow(h,2);
		}
	}
	
	sstream << "A :" << endl;
	sstream << "first diagonal above main :" << endl;
	for ( int i=0; i < n-1; i++ )
		sstream << A[i][i+1] << " ";
	sstream << endl;
	
	sstream << "main diagonal :" << endl;
	for ( int i=0; i < n; i++ )
		sstream << A[i][i] << " ";
	sstream << endl;
	
	sstream << "first diagonal below main :" << endl;
	for ( int i=0; i < n-1; i++ )
		sstream << A[i+1][i] << " ";
	sstream << endl;
	
	sstream << "interval : (" << from << ", " << to << "]" << endl;
	alglib::ae_int_t eigen_count = 0;
	//alglib::smatrixevd(A,n,1,0,eigenvalues,eigenvectors);
	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	alglib::smatrixevdr( A, n, 0, 0, from, to, eigen_count, eigenvalues, eigenvectors); // bisection method
	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	
	sstream << "solving time : " << time_span.count() << endl;
	sstream << "eigenvalues count in interval : " << eigen_count << endl; 
	//cout << "eigenvalues count : " << eigenvalues.length() << endl;

	for ( int i=0; i < eigenvalues.length(); i++ ) {
		sstream << "eigenvalue # " << i << " : " << eigenvalues[i] << endl;
		/*sstream << "eigenvector # " << i << " : " << endl;
		for ( int j=0; j < n; j++ )
			sstream << eigenvectors[i][j] << " ";
		if ( i < eigenvalues.length() - 1 ) 
			sstream << endl;*/
	}
	
	ofstream ofile( "out" );
	ofile << sstream.str();
	ofile.close();
	sstream.clear(); sstream.str("");

	cout << "finding eigenvalues done" << endl;

	// fill vector wave_numbers
	// test filling
	spectr_problem_eigenvalues.resize(10);
	for( unsigned i=0; i < spectr_problem_eigenvalues.size(); i++ )
		spectr_problem_eigenvalues[i] = i;
	// 
	
	return spectr_problem_eigenvalues;
}

// function for calculating wave numbers with Richardson extrapolation
vector<double> calc_wave_numbers_richardson( layer  &cur_layer, // info about current layer
									         double &initial_depth_step, // initial step of depth
									         double &freq ) // sound frequency
{
	vector<double> wave_numbers;
	
	// call function
	vector<double> channel_mods_wave_numbers, sound_velocity, density, point_indexes;
	// fill vectors sound_velocity, density, point_indexes
	unsigned channel_mods_iterations = 1; // test filling
	double cur_depth_step;
	for ( unsigned i = 0; i < channel_mods_iterations; i++ ) {
		cur_depth_step = initial_depth_step + i*0.1;
		channel_mods_wave_numbers = calc_chanel_mods_wave_numbers( cur_depth_step, freq, sound_velocity, density, point_indexes );
	}

	// fill vector wave_number
	// ...
	
	return wave_numbers;
}

/* Function for calculating eigenfunctions and group velocities
Функция для расчета собственных фуннкций и групповых скоростей:
Вход:
Как в функции (2) + значение волнового числа;
Работа:
-- Найти собственные функции методом Рунге-Кутта или явно с помощью функций Эйри.
-- Проинтегрировать их и нормировать на 1.
-- Найти групповую скорость. 
Выход:
собственная функции psij, нормированные на 1.
групповая скорость моды.*/
eigenfunctions_velocities calc_eigenfunctions_velocities( layer  &cur_layer,          // info about current layer
									                      double &initial_depth_step, // initial step of depth
									                      double &freq,               // sound frequency
											              double &wave_number )       // value of a wave number
{
	eigenfunctions_velocities eigenfunc_veloc; // for output data
	
	// fill eigenfunctions and velocities
	// test filling
	eigenfunc_veloc.eigenfunctions.resize(10);
	for( unsigned i=0; i < eigenfunc_veloc.eigenfunctions.size(); i++ )
		eigenfunc_veloc.eigenfunctions[i] = i;
	eigenfunc_veloc.velocities.resize(10);
	for( unsigned i=0; i < eigenfunc_veloc.velocities.size(); i++ )
		eigenfunc_veloc.velocities[i] = i;
	
	return eigenfunc_veloc;
}

