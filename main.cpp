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




vector<double> compute_trapped_modes(double &omeg, vector<double> &c, vector<double> &rho, vector<int> &interface_idcs, vector<double> &meshsizes);


int main( int argc, char **argv )
{

	double freq = 50;
    double c_w = 1500;
    double c_b = 2000;
    double rho_w = 1;
    double rho_b = 2;
    double dz = 1;
    int nz = 201;
    int ib = 90;    //at POINT = 89 we have water, at POINT = 90 we have bottom
                    //   ii = 89,                  at ii = 90

    double omeg = 2*3.141592653589793*freq;

    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh { dz,dz };
    vector<int> input_interf_idcs { ib };
    vector<double> out_wnum;

    for ( int ii = 0; ii<nz; ii++ ){
        if (ii<ib){                   //
            input_c.push_back(c_w);
            input_rho.push_back(rho_w);
        }
        else {
            input_c.push_back(c_b);
            input_rho.push_back(rho_b);
        }

    }

    cout << "freq = " << freq << endl;
    cout << "omega = " << omeg << endl;

    out_wnum = compute_trapped_modes(omeg, input_c, input_rho,input_interf_idcs, input_mesh);

   for (int ii=0; ii<out_wnum.size();  ii++) {

            cout << ii << "->" << out_wnum.at(ii) << endl;


    }

	return 0;
}




vector<double> compute_trapped_modes(           double &omeg, // sound frequency
											  vector<double> &c,
											  vector<double> &rho,
											  vector<int> &interface_idcs,
											  vector<double> &meshsizes)
{
	// prepare the three diagonals of the matrix to be passed to the EIG function
    // for the c = c_j, j=0... N_points
    // interfaces are at z = z_m,  interface_idcs = {m}, if empty then we have NO interfaces
    // mesh size in the j-th layer is meshsizes.at(j); this vector has at least one element,
    // for the k-th interface interface_idcs.at(k-1) we have meshsizes meshsizes.at(k-1) and meshsizes.at(k) below and above the interface respectively
    // for c(interface_idcs.at(k-1)) the value of c is the one BELOW the k-th interface
    //(i.e. for the water-bottom interface at the boundary we take the value from the bottom)


	vector<double> md;
	vector<double> ud;
	vector<double> ld;

    int N_points = c.size();
    int layer_number = 0;
    double dz = meshsizes.at(layer_number);
    double dz_next = dz;
    double q = 0;
    double cp, cm, dp, dm, cmin, cmax, kappamax, kappamin;
    int next_interface_idx;

    if ( interface_idcs.size() > 0 )
    {
        next_interface_idx = interface_idcs.at(0)-1;
    }
    else
    {
        next_interface_idx = N_points;
    }

    cmin = c.at(0);
    cmax = c.at(0);

    for( int ii=0; ii < N_points; ii++ ){
        if (c.at(ii) < cmin) { cmin = c.at(ii); }
        if (c.at(ii) > cmax) { cmax = c.at(ii); }

    }

    kappamax = omeg/cmin;
    kappamin = omeg/cmax;


    for(int ii=0; ii < N_points-2; ii++ ){

        // ordinary point


        ud.push_back( 1/(dz*dz) );
        ld.push_back( 1/(dz*dz) );
        md.push_back(-2/(dz*dz)  + omeg*omeg/(c.at(ii+1)*c.at(ii+1)) );


        // special case of the point at the interface

        if (ii == next_interface_idx) {         //ii -- z(ii+1), z(0) = 0
            layer_number = layer_number + 1;    // вообще ii=89 -- вода, в ii=90 -дно,
                                                // здесь ii = 89 -- интерфейс, уже дно

            cp = c.at(ii+1);
            dp = rho.at(ii+1);
            cm = c.at(ii);
            dm = rho.at(ii);
            q = 1/( dz_next*dm + dz*dp );

            dz_next = meshsizes.at(layer_number);


            ld.at(ii) = 2*q*dp/dz;
            md.at(ii) = -2*q*( dz_next*dp + dz*dm )/(dz*dz_next) + omeg*omeg*q*( dz*dp*cp*cp + dz_next*dm*cm*cm )/( cp*cp*cm*cm ) ;
            ud.at(ii) = 2*q*dm/dz_next;

            if ( interface_idcs.size() > layer_number )
            {
                next_interface_idx = interface_idcs.at(layer_number) - 1;
            }
            else
            {
                next_interface_idx = N_points;
            }

            dz = dz_next;
        }

    }

    // HERE WE CALL THE EIG ROUTINE!!!
    //input: diagonals ld, md, ud + interval [0 k_max]
    //output: wnumbers2 = wave numbers squared

    alglib::real_2d_array A, eigenvectors; // V - собств вектор
	alglib::real_1d_array eigenvalues; // Lm -собств знач
	A.setlength(N_points-2,N_points-2);
	eigenvectors.setlength(N_points-2,N_points-2);
	eigenvalues.setlength(N_points-2);
    alglib::ae_int_t eigen_count = 0;




	// fill matrix by zeros
	for (int ii=0; ii < N_points-2; ii++ )
		for ( int jj=0; jj < N_points-2; jj++ )
			A[ii][jj] = 0.0;


	// fill tridiagonal matrix
	// make Dirichlet case matrix
	// Ќа главной диагонали сто€т числа -2/(h^2), на под- и над- диагонал€х сто€т числа 1/(h^2).
	for ( int ii=0; ii < N_points-2; ii++ ) {
		A[ii][ii] = md.at(ii);


		if ( ii >0 ) {
            A[ii][ii-1] = ld.at(ii);
		}

		if ( ii < N_points-3 ) {
            A[ii][ii+1] = ud.at(ii);
		}

	}


    ofstream myFile("thematrixdiags.txt");
    for (int ii=0; ii<N_points-2; ii++){
        myFile << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
    }


	vector<double> wnumbers;




    //alglib::smatrixevdr( A, N_points-2, 0, 0, kappamin, kappamax, eigen_count, eigenvalues, eigenvectors); // CHECK CALL ARGUMENTS!
    alglib::smatrixevdr( A, N_points-2, 0, 0, kappamin*kappamin, kappamax*kappamax, eigen_count, eigenvalues, eigenvectors); // CHECK CALL ARGUMENTS!

    for (int ii=0; ii<eigen_count; ii++) {
            wnumbers.push_back( sqrt(eigenvalues[eigen_count-ii-1]) );

    }


    return wnumbers;
}


/* function for calculating wave number of channel mods
‘ункци€ дл€ расчета волновых чисел канальных мод
¬ход:
dz -- шаг по глубине;
f -- частота звука;
вектор c = (ci); -- скорости звука с шагом dz; ƒлина = N. «начени€ N -- "большие".
вектор d = (di); -- плотности; ƒлина N.
вектор m = (mj); -- индексы точек, где параметры среды терп€т разрыв
                     (точки, где заканчиваетс€ один слой и начинаетс€ другой). ƒлина M <=10.
–абота: сформировать трехдиагональные матрицы, запустить солвер собственных значений дл€ отрезка [omega/cmax omega/cmin].
¬ыход: собственные значени€ kj^2 акустической спектральной задачи */

/*

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
	double from = -0.001, to = 0.001; // interval for eigenvalues

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
	// Ќа главной диагонали сто€т числа -2/(h^2), на под- и над- диагонал€х сто€т числа 1/(h^2).
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
	alglib::smatrixevdr( A, n, 1, 0, from, to, eigen_count, eigenvalues, eigenvectors); // bisection method
	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);

	sstream << "solving time : " << time_span.count() << endl;
	sstream << "eigenvalues count in interval : " << eigen_count << endl;
	//cout << "eigenvalues count : " << eigenvalues.length() << endl;

	for ( int i=0; i < eigenvalues.length(); i++ ) {
		sstream << "eigenvalue # " << i << " : " << eigenvalues[i] << endl;
		sstream << "eigenvector # " << i << " : " << endl;
		for ( int j=0; j < n; j++ )
			sstream << eigenvectors[i][j] << " ";
		if ( i < eigenvalues.length() - 1 )
			sstream << endl;
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


*/
