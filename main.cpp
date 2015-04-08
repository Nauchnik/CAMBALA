// +-----------------------------------------------------------------------------------+
// | Client application for the volunteer computing project Acoustics@home             |
// +-----------------------------------------------------------------------------------+
// | Pacific Oceanological Institute, Institute for System Dynamics and Control Theory |
// +-----------------------------------------------------------------------------------+
// | Authors: Pavel Petrov, Oleg Zaikin                                                |
// +-----------------------------------------------------------------------------------+

#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include "linalg.h"
#include <cmath>

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
vector<double> compute_wnumbers_extrap_lin_dz(double &omeg, vector<double> &depths,vector<double> &c1s,vector<double> &c2s,vector<double> &rhos,vector<unsigned> &Ns_points, unsigned flOnlyTrapped,unsigned &ordRich);
int compute_modal_grop_velocities( vector<double> &freqs, double deltaf, vector<double> &depths, vector<double> &c1s, vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, unsigned flOnlyTrapped, unsigned &ordRich, vector<vector<double>> &modal_group_velocities, vector<unsigned> &mode_numbers );
double compute_modal_delays_residual_uniform( vector<double> &freqs, vector<double> &depths, vector<double> &c1s, vector<double> &c2s, vector<double> &rhos, vector<unsigned> &Ns_points, double R, vector<vector<double>> &experimental_delays, vector<unsigned> &experimental_mode_numbers);

int main( int argc, char **argv )
{
/*
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

    double omeg = 2*M_PI*freq;

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

    out_wnum = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,1);

   for (unsigned ii=0; ii<out_wnum.size();  ii++) {

            cout << ii << "->" << out_wnum.at(ii) << endl;


    }

    cout << "NEW: extrapolation" << endl;

    vector<double> depths {90,600};
    vector<double> c1s  {1500,2000};
    vector<double> c2s  {1500,2000};
    vector<double> rhos  {1,2};
    vector<unsigned> Ns_points  {180,1020};
    unsigned rord = 3;
    vector<double> freqs {20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100};
    vector<vector<double>> modal_group_velocities;
    vector<unsigned> mode_numbers;
    double deltaf = 0.5;

    out_wnum = compute_wnumbers_extrap_lin_dz(omeg,depths,c1s,c2s,rhos,Ns_points,1,rord);

    cout << "Extrapolated ev:" << endl;
    for (unsigned ii=0; ii<out_wnum.size();  ii++) {

            cout << ii << "->" << out_wnum.at(ii) << endl;


    }

    compute_modal_grop_velocities( freqs, deltaf, depths,c1s, c2s, rhos, Ns_points, 1, rord, modal_group_velocities, mode_numbers );

    ofstream myFile("mgv.txt");
//    for (int ii=0; ii<N_points-2; ii++){
//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
//    }
//

    for (unsigned ii=0; ii<freqs.size();  ii++) {
            cout << "f=" << freqs.at(ii) << endl;

            for (unsigned jj=0; jj<mode_numbers.at(ii);  jj++) {

            cout << modal_group_velocities[ii][jj] << endl;
            myFile << std::fixed << std::setprecision(16) << modal_group_velocities[ii][jj] << "\t";

            }
            myFile << endl;
    }

    myFile.close();

*/
	vector<double> depths{90,600};
    vector<double> c1s{1500,2000};
    vector<double> c2s{1500,2000};
    vector<double> rhos{1,2};
    vector<unsigned> Ns_points{180,1020};
    unsigned rord = 3;
    vector<vector<double>> modal_group_velocities;
    vector<unsigned> mode_numbers;
    vector<vector<double>> modal_delays;
    vector<double> freqs;

    freqs.clear();
    double buff;
    vector<double> buffvect;
    mode_numbers.clear();
    string myLine;
    double residual;
    double R = 3500;

    //SEARCH SPACE
    double cb1 = 1600;
    double cb2 = 2600;
    double cbmin, cb_cur;
    unsigned ncb = 10;

    double rhob1 = 1;
    double rhob2 = 3;
    double rhobmin, rhob_cur;
    unsigned nrhob = 2;

    double R1 = 3000;
    double R2 = 4000;
    double Rmin, R_cur;
    unsigned nR = 10;

    double resmin = 10000;

	// fix start time 
	chrono::high_resolution_clock::time_point t1 = std::chrono::high_resolution_clock::now();
	
	ifstream myFileSynth("dtimes_synth_1.txt"); // delays for R = 3500, cw = 1500, cb = 2000, rhow = 1, rhob = 2;

    while ( getline(myFileSynth,myLine) ){
        stringstream myLineStream(myLine);
        myLineStream >> buff;
        freqs.push_back(buff);
        buffvect.clear();

        while (!myLineStream.eof()) {
            myLineStream >> buff;
            buffvect.push_back(buff);
			mode_numbers.push_back((unsigned)buffvect.size());
        }

        modal_delays.push_back(buffvect);
    }
    myFileSynth.close();

    // sample call of the residual computation routine
    residual = compute_modal_delays_residual_uniform( freqs, depths,c1s, c2s, rhos, Ns_points, R, modal_delays, mode_numbers);
    cout << "Residual is: " << residual << endl;

    // brute force minimum search
    cout << "BRUTE FORCE MINIMUM SEARCH" << endl;
    cout << "Search space:" << endl;
    cout << cb1 << "< c_b <" << cb2 << ", step" << (cb2 - cb1)/ncb << endl;
    cout << rhob1 << "< rho_b <" << rhob2 << ", step" << (rhob2 - rhob1)/nrhob << endl;
    cout << R1 << "< Range <" << R2 << ", step" << (R2 - R1)/nR << endl;

    for (unsigned ii=0; ii<=ncb; ii++) {
        cb_cur = cb1 + ii*(cb2 - cb1)/ncb;
        c1s.at(1) = cb_cur;
        c2s.at(1) = cb_cur;

        for (unsigned jj=0; jj<=nrhob; jj++) {
            rhob_cur = rhob1 + jj*(rhob2 - rhob1)/nrhob;
            rhos.at(1) = rhob_cur;

            for (unsigned kk=0; kk<=nR; kk++) {
                R_cur = R1 + kk*(R2 - R1)/nR;

                residual = compute_modal_delays_residual_uniform( freqs, depths,c1s, c2s, rhos, Ns_points, R_cur, modal_delays, mode_numbers);

                if (residual < resmin) {
                    resmin = residual;
                    cbmin = cb_cur;
                    rhobmin = rhob_cur;
                    Rmin = R_cur;
                    cout << "New residual minimum:" << endl;
                    cout << "err=" << resmin << ", parameters:" << endl;
                    cout << "c_b=" << cbmin << ", rho_b=" << rhobmin << ", R=" << Rmin <<  endl;

                }
            }
        }
    }

	// fix final time
	chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
	chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);
	
    cout << "SEARCH ENDED!" << endl;
    cout << "RESULTING VALUE:" << endl;
    cout << "err=" << resmin << ", parameters:" << endl;
    cout << "c_b=" << cbmin << ", rho_b=" << rhobmin << ", R=" << Rmin <<  endl;
	cout << "time " << time_span.count() << endl;

	ofstream out("out_time");
	out << time_span.count() << " s";
	out.close();
	
	return 0;
}

/*
A routine for computing delay residual.
Arguments:
1) Environment: five arrays of the same length: depth, c1s, c2s, rhos, Ns_points;
(each entry describes one layer as described in the comments to compute_wnumbers_extrap() )

2) Source-receive distance: R -- distance from the source to the receiver;

3) Experimental data: modal delays:
    -- experimental_mode_numbers: number of modes for each frequency in the recorded signal
    -- experimental_delays: experimental_delays[ii][jj] is the delay of jj+1-th mode for the frequency freqs[ii]

The routine computes the (uniform) residual (misfit) of experimental data and the "theoretical" delays for a given environment model.

It should be used as follows: for a set of environment models the residual should be computed. The minimal value of the residual indicates
the most "adequate" model.
*/
double compute_modal_delays_residual_uniform( vector<double> &freqs,
                                     vector<double> &depths,
                                     vector<double> &c1s,
                                     vector<double> &c2s,
                                     vector<double> &rhos,
                                     vector<unsigned> &Ns_points,
                                     double R,
                                     vector<vector<double>> &experimental_delays,
                                     vector<unsigned> &experimental_mode_numbers
                                     )
{
    unsigned rord = 3;
    unsigned flTrappedOnly = 1;
    double deltaf = 0.5;
    double residual = 0;
    unsigned mnumb;
    double mdelay;

    vector<vector<double>> modal_group_velocities;
    vector<unsigned> mode_numbers;

    compute_modal_grop_velocities( freqs, deltaf, depths,c1s, c2s, rhos, Ns_points, flTrappedOnly, rord, modal_group_velocities, mode_numbers );

    for (unsigned ii=0; ii<freqs.size();  ii++) {
		mnumb = min(mode_numbers.at(ii), experimental_mode_numbers.at(ii) );
        for (unsigned jj=0; jj<mnumb;  jj++) {
			if (experimental_delays[ii][jj]>0) {
				mdelay =  R/modal_group_velocities[ii][jj];
                residual = residual + pow(experimental_delays[ii][jj]-mdelay,2);
            }
        }
    }

	residual = sqrt(residual);

	return residual;
}


int compute_modal_grop_velocities(      vector<double> &freqs,
                                        double deltaf,
                                        vector<double> &depths,
                                        vector<double> &c1s,
                                        vector<double> &c2s,
                                        vector<double> &rhos,
                                        vector<unsigned> &Ns_points,
                                        unsigned flOnlyTrapped,
                                        unsigned &ordRich,
                                        vector<vector<double>> &modal_group_velocities,
                                        vector<unsigned> &mode_numbers
                                        )
{
    mode_numbers.clear();
    modal_group_velocities.clear();

    vector<double> out_wnum1;
    vector<double> out_wnum2;
    vector<double> mgv_ii;
    unsigned nwnum;
	unsigned nfr = (unsigned)freqs.size();
    double omeg1, omeg2;

    for (unsigned ii=0; ii<nfr; ii++) {
        out_wnum1.clear();
        out_wnum2.clear();
        mgv_ii.clear();
        omeg1 = 2*M_PI*(freqs.at(ii) + deltaf/2);
        out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1,depths,c1s,c2s,rhos,Ns_points,1,ordRich);
		nwnum = (unsigned)out_wnum1.size();

        /*
        cout << "f=" << freqs.at(ii) << "Hz" << endl;

        for (unsigned jj=0; jj < nwnum; jj++ )
        {
            cout << "k_" << jj+1 << "=" << out_wnum1.at(jj) << endl;
        }
        */

		omeg2 = 2 *M_PI*(freqs.at(ii) - deltaf / 2);
        out_wnum2 = compute_wnumbers_extrap_lin_dz(omeg2,depths,c1s,c2s,rhos,Ns_points,1,ordRich);
        nwnum = min(nwnum, (unsigned)out_wnum2.size());

        for (unsigned jj=0; jj < nwnum; jj++ )
        {
            mgv_ii.push_back(  (omeg1 - omeg2)/(out_wnum1.at(jj) - out_wnum2.at(jj) )  );
        }

        modal_group_velocities.push_back( mgv_ii );
        mode_numbers.push_back( nwnum );
    }

	return 0;
}


/*
General considerations:
1) It is better to choose mesh in such a way that mesh size is the same for all z. This gives better accuracy!
(the reason: if the diagonal elements in the upper and lower diagonals vary with row index j, then accuracy is crippled, the
 closer are u_j to each other -- the better, ideally they should be equal)
2) compute_wnumbers_extrap_lin_dz() makes the number of points within each layer a multiple of 12;
it is better to set all the numbers to multiple of 12 in advance
3) Richardson extrapolation of the order 3 gives reasonable accuracy
*/
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
	switch (ordRich) {
	case 1 : 
		coeff_extrap.assign({ 1 });
		break;
	case 2 :
		coeff_extrap.assign({ -1, 2 });
		break;
	case 4 :
		coeff_extrap.assign({ -1 / double(6), 4, -13.5, 32 / double(3) });
	default :
		ordRich = 3;
		coeff_extrap.assign({ 0.5, -4, 4.5 });
		//coeff_extrap.assign({0.1, -0.6, 1.5});
		break;
	}
//
//    cout << "Richardson coeffs" << endl;
//    for (unsigned ii=0; ii<coeff_extrap.size() ; ii++ ){
//        cout << coeff_extrap.at(ii) << endl;
//    }

    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh;
    vector<unsigned> input_interf_idcs;
    vector<double> out_wnum2;
    vector<double> wnum_extrapR;
    double zc = 0;
    double zp = 0;
    double dz = 0;
    unsigned m_wnum = 1000;

	unsigned n_layers = (unsigned)depths.size();
    unsigned n_points_total = 0;
    unsigned n_points_layer = 0;

// outer loop for Richardson coefficient rr
    for (unsigned rr = 1; rr <= ordRich; rr++){
        input_c.clear();
        input_rho.clear();
        input_interf_idcs.clear();
        input_mesh.clear();
        out_wnum2.clear();

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

        //cout << "rr=" << rr << endl;

        out_wnum2 = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,flOnlyTrapped);
		m_wnum = min(m_wnum, (unsigned)out_wnum2.size());

        if (rr == 1) { wnum_extrapR.assign(m_wnum,0);}

        for (unsigned mm=0; mm<m_wnum; mm++ ) {
            wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr-1);
        }

        /*
            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
            }
        */
    }

	for (unsigned mm=0; mm<m_wnum; mm++ ) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}

	return wnum_extrapR;
}



vector<double> compute_wnumbers_extrap_lin_dz( double &omeg, // sound frequency
											   vector<double> &depths,
											   vector<double> &c1s,
											   vector<double> &c2s,
											   vector<double> &rhos,
											   vector<unsigned> &Ns_points,
											   unsigned flOnlyTrapped,
											   unsigned &ordRich )
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
    vector<double> input_c;
    vector<double> input_rho;
    vector<double> input_mesh;
    vector<unsigned> input_interf_idcs;
    vector<double> out_wnum2;
    vector<double> wnum_extrapR;
    double zc = 0;
    double zp = 0;
    double dz = 0;
    unsigned m_wnum = 1000;

	unsigned n_layers = (unsigned)depths.size();
    unsigned n_points_total = 0;
    unsigned n_points_layer = 0;

    vector<double> coeff_extrap;
	switch (ordRich) {
	case 1 :
		coeff_extrap.assign({1});
		break;
	case 2 :
		coeff_extrap.assign({1.333333333333333, -0.333333333333333});
		break;
	case 4 :
		// coeff_extrap.assign({1.595325630252102, -0.788449059052564, 0.216346153846154, -0.023222725045691});
		coeff_extrap.assign({1.6, -0.8, 0.228571428571429, -0.028571428571429});
		break;
	default:
		ordRich = 3;
		coeff_extrap.assign({1.5, -0.6, 0.1});
		//coeff_extrap.assign({0.1, -0.6, 1.5});
		break;
	}
	
    // number of points in each layer is multiple of 12
    // this allows us to use nz_ii = nz/ii, ii = 1,2,3,4
    for (unsigned ii=0; ii < n_layers; ii++ ){
        Ns_points.at(ii) = 12*(Ns_points.at(ii)/12);
    }

//    cout << "Richardson coeffs" << endl;
//    for (int ii=0; ii<coeff_extrap.size() ; ii++ ){
//        cout << coeff_extrap.at(ii) << endl;
//    }

// outer loop for Richardson coefficient rr
    for (unsigned rr = 1; rr <= ordRich; rr++){
        input_c.clear();
        input_rho.clear();
        input_interf_idcs.clear();
        input_mesh.clear();
        out_wnum2.clear();

        input_c.push_back(0);
        input_rho.push_back(0);
        n_points_total = 1;
        zp = 0;

        for (unsigned ll = 0; ll<n_layers; ll++){
            zc = depths.at(ll);
            n_points_layer = Ns_points.at(ll)/rr;
            dz = (zc - zp)/(  n_points_layer  );

//            cout << "np=" << n_points_layer << "  " << "dz=" << dz <<endl;

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

//        cout << "rr=" << rr << endl;

        out_wnum2 = compute_wnumbers(omeg, input_c, input_rho,input_interf_idcs, input_mesh,flOnlyTrapped);
		m_wnum = min(m_wnum, (unsigned)out_wnum2.size());

        if (rr == 1) { wnum_extrapR.assign(m_wnum,0);}

        for (unsigned mm=0; mm<m_wnum; mm++ ) {
            wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr-1);
        }

//            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
//
//                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
//            }
//            cout << endl;
    }

	for (unsigned mm=0; mm<m_wnum; mm++ ) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}
	return wnum_extrapR;
}


vector<double> compute_wnumbers( double &omeg, // sound frequency
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
    vector<double> wnumbers2;

	int N_points = (unsigned)c.size();
    unsigned layer_number = 0;
    double dz = meshsizes.at(layer_number);
    double dz_next = dz;//    ofstream myFile("thematrixdiags.txt");
//    for (int ii=0; ii<N_points-2; ii++){
//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
//    }
//    myFile.close();
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
	eigenvalues.setlength(N_points-2);
    alglib::ae_int_t eigen_count = 0;

    alglib::real_1d_array main_diag, second_diag;
    main_diag.setlength(N_points-2);
    second_diag.setlength(N_points-3);

    for ( int ii=0; ii < N_points-3; ii++ ) {
        second_diag[ii] = sqrt(ud.at(ii)*ld.at(ii+1));
        main_diag[ii] = md.at(ii);
	}
    main_diag[N_points-3] = md.at(N_points-3);

//    ofstream myFile("thematrixdiags.txt");
//    for (int ii=0; ii<N_points-2; ii++){
//        myFile << std::fixed << std::setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
//    }
//    myFile.close();

//    ofstream myFile1("thematrixdiags_sym.txt");
//    for (int ii=0; ii<N_points-3; ii++){
//        myFile1 << std::fixed << std::setprecision(16) << main_diag[ii] << "  " << second_diag[ii] << endl;
//    }
//    myFile1.close();

    alglib::smatrixtdevdr(main_diag,  second_diag,  N_points-2,  0, kappamin*kappamin,  kappamax*kappamax, eigen_count,eigenvectors);

    for (int ii=0; ii<eigen_count; ii++) {
        wnumbers2.push_back( main_diag[eigen_count-ii-1] );
    }

    return wnumbers2;
}

