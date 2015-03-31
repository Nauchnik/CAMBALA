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
#include <iomanip>
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




vector<double> compute_wnumbers(double &omeg, vector<double> &c, vector<double> &rho, vector<unsigned> &interface_idcs, vector<double> &meshsizes,unsigned flOnlyTrapped);

vector<double> compute_wnumbers_extrap(double &omeg, vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, unsigned flOnlyTrapped,unsigned &ordRich);

int main( int argc, char **argv )
{

	double freq = 50;
    double c_w = 1500;
    double c_b = 2000;
    double rho_w = 1;
    double rho_b = 2;
    double dz = 1;
    unsigned nz = 501;
    unsigned ib = 90;    //at POINT = 89 we have water, at POINT = 90 we have bottom
                    //   ii = 89,                  at ii = 90

    cout.precision(15);

    double omeg = 2*3.141592653589793*freq;

    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh { dz,dz };
    vector<unsigned> input_interf_idcs { ib };
    vector<double> out_wnum;

    for ( unsigned ii = 0; ii<nz; ii++ ){
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

/*    out_wnum = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,1);

   for (unsigned ii=0; ii<out_wnum.size();  ii++) {

            cout << ii << "->" << out_wnum.at(ii) << endl;


    }
*/
    cout << "NEW: extrapolation" << endl;

    vector<double> depths {90,1500};
    vector<double> c1s  {1500,2000};
    vector<double> c2s  {1500,2000};
    vector<double> rhos  {1,2};
    vector<unsigned> Ns_points  {60,940};
    unsigned rord = 3;

    out_wnum = compute_wnumbers_extrap(omeg,depths,c1s,c2s,rhos,Ns_points,1,rord);

    cout << "Extrapolated ev:" << endl;
    for (unsigned ii=0; ii<out_wnum.size();  ii++) {

            cout << ii << "->" << out_wnum.at(ii) << endl;


    }

	return 0;


}




vector<double> compute_wnumbers_extrap(       double &omeg, // sound frequency
											  vector<double> &depths,
											  vector<double> &c1s,
											  vector<double> &c2s,
											  vector<double> &rhos,
											  vector<unsigned> &Ns_points,
											  unsigned flOnlyTrapped,
											  unsigned &ordRich)
/*  subroutine for computing wavenumbers for a given waveguide structure
    the computation is performed by the FD method for certain meshsize,
    Richardson extrapolation is used to improve the


Layer structure:

depths_{i-1}----c=c1s_i--
*
*
... <-Ns_points_i       rho = rhos_i
*
*
depths_{i}------c=c2s_i--

Other parameters:
omeg -- cyclic frequency, omeg = 2*Pi*f;
flOnlyTrapped -- flag to determine the mode subset: set to 0 to compute all propagating modes, i.e. such that k^2>=0, otherwise only trapped modes are computed
ordRich -- order of the Richardson extrapolation;

the top of the first layer is z=0
	*/


{
    vector<double> coeff_extrap;
    if (ordRich == 1) {
        coeff_extrap.assign({1});
    }
    else if (ordRich == 2){
        coeff_extrap.assign({-1,2});
    }
    else if (ordRich == 3){
        coeff_extrap.assign({0.5, -4, 4.5});
        //coeff_extrap.assign({0.1, -0.6, 1.5});
    }
    else if (ordRich == 4){
        coeff_extrap.assign({-1/double(6),4,-13.5,32/double(3)});
    }
    else {
        ordRich = 3;
        coeff_extrap.assign({0.5, -4, 4.5});
    }

    cout << "Richardson coeffs" << endl;
    for (int ii=0; ii<coeff_extrap.size() ; ii++ ){
        cout << coeff_extrap.at(ii) << endl;
    }


    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh;
    vector<unsigned> input_interf_idcs;
    vector<double> out_wnum;
    vector<double> wnum_extrapR;
    double zc = 0;
    double zp = 0;
    double dz = 0;
    unsigned m_wnum = 1000;

    unsigned n_layers = depths.size();
    unsigned n_points_total = 0;
    unsigned n_points_layer = 0;
// outer loop for Richardson coefficient rr
    for (unsigned rr = 1; rr <= ordRich; rr++){

        input_c.clear();
        input_rho.clear();
        input_interf_idcs.clear();
        input_mesh.clear();
        out_wnum.clear();

        input_c.push_back(0);
        input_rho.push_back(0);
        n_points_total = 1;
        zp = 0;

        for (unsigned ll = 0; ll<n_layers; ll++){
            zc = depths.at(ll);
            n_points_layer = Ns_points.at(ll)*rr;
            dz = (zc - zp)/(  n_points_layer  );
            input_mesh.push_back(  dz  );
            input_c.at(n_points_total-1) = c1s.at(ll) ;
            input_rho.at(n_points_total-1) = rhos.at(ll) ;

            n_points_total = n_points_total + n_points_layer;

            for (unsigned kk=1; kk<= n_points_layer; kk++) {
                input_rho.push_back(rhos.at(ll));
                input_c.push_back( c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk/n_points_layer );
            }

            if (ll < n_layers - 1) {
                input_interf_idcs.push_back(n_points_total-1);
            }
            zp = zc;
        }

        cout << "rr=" << rr << endl;

        out_wnum = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,flOnlyTrapped);
        m_wnum = min(m_wnum, out_wnum.size() );

        if (rr == 1) { wnum_extrapR.assign(m_wnum,0);}


        for (unsigned mm=0; mm<m_wnum; mm++ ) {
            wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum.at(mm)*coeff_extrap.at(rr-1);
        }



            for (unsigned ii=0; ii<out_wnum.size();  ii++) {

                cout << ii << "->" << out_wnum.at(ii) << endl;
            }

    }


return wnum_extrapR;
}



vector<double> compute_wnumbers(           double &omeg, // sound frequency
											  vector<double> &c,
											  vector<double> &rho,
											  vector<unsigned> &interface_idcs,
											  vector<double> &meshsizes,
											  unsigned flOnlyTrapped                 // set flOnlyTrapped = 0 to compute all propagating modes, i.e. such that k^2>=0
											  )
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
    vector<double> wnumbers;

    int N_points = c.size();
    unsigned layer_number = 0;
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

    if (flOnlyTrapped == 0 )
        kappamin = 0;


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

    alglib::real_2d_array eigenvectors; // V - собств вектор
	alglib::real_1d_array eigenvalues; // Lm -собств знач
	//eigenvectors.setlength(N_points-2,N_points-2);
	eigenvalues.setlength(N_points-2);
    alglib::ae_int_t eigen_count = 0;

    //new function !!!
    alglib::real_1d_array main_diag, second_diag;
    main_diag.setlength(N_points-2);
    second_diag.setlength(N_points-3);





    for ( int ii=0; ii < N_points-3; ii++ ) {
        second_diag[ii] = sqrt(ud.at(ii)*ld.at(ii+1));
        main_diag[ii] = md.at(ii);
	}
    main_diag[N_points-3] = md.at(N_points-3);


    ofstream myFile("thematrixdiags.txt");
    for (int ii=0; ii<N_points-2; ii++){
        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
    }
    myFile.close();

    ofstream myFile1("thematrixdiags_sym.txt");
    for (int ii=0; ii<N_points-3; ii++){
        myFile1 << std::fixed << std::setprecision(16) << main_diag[ii] << "  " << second_diag[ii] << endl;
    }
    myFile1.close();


    alglib::smatrixtdevdr(main_diag,  second_diag,  N_points-2,  0, kappamin*kappamin,  kappamax*kappamax, eigen_count,eigenvectors);

    for (int ii=0; ii<eigen_count; ii++) {
        wnumbers.push_back( sqrt(main_diag[eigen_count-ii-1]) );
    }






    return wnumbers;
}

