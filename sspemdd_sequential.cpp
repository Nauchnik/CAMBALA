#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"
#include <iostream>
#include <complex>
#include <time.h>
#include <stdexcept>


double RK4(double omeg2, // sound frequency
    double kh2,
	double deltah,
	double c1,
	double c2,
	unsigned Np,
	vector<double> &phi0,
	vector<double> &dphi0)
{

    double f11,f12,f21,f22,f31,f32,f41,f42, cc;
    double h = deltah/Np;
    double layer_int = 0.0;


    for (unsigned kk = 0; kk <Np; kk++) {

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

        f11 = dphi0.back();
        cc = c1 + (c2-c1)*kk*h/deltah;
        cc = cc*cc;
        f12 = (kh2 - omeg2/cc  )*phi0.back();

        f21 = dphi0.back() + 0.5*f12*h;
        cc = c1 + (c2-c1)*(kk + 0.5)*h/deltah;
        cc = cc*cc;
        f22 = (kh2 - omeg2/cc  )*(phi0.back() + 0.5*f11*h) ;

        f31 = dphi0.back() + 0.5*f22*h;
        cc = c1 + (c2-c1)*(kk + 0.5)*h/deltah;
        cc = cc*cc;
        f32 = (kh2 - omeg2/cc  )*(phi0.back() + 0.5*f21*h) ;

        f41 = dphi0.back() + f32*h;
        cc = c1 + (c2-c1)*(kk + 1)*h/deltah;
        cc = cc*cc;
        f42 = (kh2 - omeg2/cc  )*(phi0.back() + f31*h) ;

        phi0.push_back( phi0.back() + h*( f11 + 2*f21 + 2*f31 + f41 )/6 );
        dphi0.push_back( dphi0.back() + h*( f12 + 2*f22 + 2*f32 + f42 )/6 );

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

    }

    return layer_int;
}


double Layer_an_exp(double omeg2, // sound frequency
    double kh2,
	double deltah,
	double c,
	unsigned Np,
	vector<double> &phi0,
	vector<double> &dphi0)
{

    double c1, c2, kv;
    double h = deltah/Np;
    double layer_int = 0.0;

    kv = sqrt( kh2 - omeg2/(c*c) );
    c1 = 0.5*(phi0.back() - dphi0.back()/kv);
    c2 = 0.5*(phi0.back() + dphi0.back()/kv);
    c2 = 0;


/*
    cout << "kv c1 c2" << endl;
    cout << kv << endl;
    cout << c1 << endl;
    cout << c2 << endl;
    cout << h << endl;
    cout << "-----" << endl;
    cout << phi0.back() << endl;
    cout << dphi0.back()/kv << endl;
*/
    for (unsigned kk = 0; kk <Np; kk++) {

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();



        phi0.push_back(     c1*exp( -kv*((kk+1)*h) )    +   c2*exp( kv*((kk+1)*h) )      ) ;
        dphi0.push_back(    -c1*kv*exp( -kv*((kk+1)*h) )    +   c2*kv*exp( kv*((kk+1)*h) )  ) ;

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();
        //cout << layer_int << endl;
    }

    return layer_int;
}





/*
double Euler(double omeg2, // sound frequency
    double kh2,
	double deltah,
	double c1,
	double c2,
	unsigned Np,
	vector<double> &phi0,
	vector<double> &dphi0)
{

    double f11,f12,f21,f22,f31,f32,f41,f42, cc;
    double h = deltah/Np;
    double layer_int = 0.0;

    cout << "h value" << endl;
    cout << h << endl;

    for (unsigned kk = 0; kk <Np; kk++) {

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

        f11 = dphi0.back();
        cc = c1 + (c2-c1)*kk*h/deltah;
        cc = cc*cc;
        f12 = (kh2 - omeg2/cc  )*phi0.back();

        if (kk<2){
            cout << "f11 f12 value" << endl;
            cout << phi0.back() << endl;
            cout << dphi0.back() << endl;
            cout << f11 << endl;
            cout << f12 << endl;
            cout << phi0.back() + h*f11 << endl;
            cout << dphi0.back() + h*f12 << endl;
        }

        phi0.push_back( phi0.back() + h*f11 );
        dphi0.push_back( dphi0.back() + h*f12 );

        layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

    }

    return layer_int;
}
*/

sspemdd_sequential::sspemdd_sequential() :
	object_function_type("weighted"),
	h(0),
	H(0),
	ncb(1),
	nrhob(1),
	nR(1),
	ntau(1),
	cb1(2000.0),
	cb2(2000.0),
	R1(3400.0),
	R2(3600.0),
	tau1(0.0),
	tau2(0.0),
	rhob1(2.0),
	rhob2(2.0),
	n_layers_w(1),
	iterated_local_search_runs(10),
	verbosity(1),
	N_total(1),
	isTimeDelayPrinting(false),
	rank(0)
{
	record_point.cb       = START_HUGE_VALUE;
	record_point.rhob     = START_HUGE_VALUE;
	record_point.R        = START_HUGE_VALUE;
	record_point.tau      = START_HUGE_VALUE;
	record_point.residual = START_HUGE_VALUE;
	srand((unsigned)time(NULL));
	start_chrono_time = std::chrono::high_resolution_clock::now();
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

//tau_comment: tau is appended here as a new argument for residual computation

void sspemdd_sequential::load_layers_data(
    std::string LayersFName,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points)
{
    std::ifstream Myfile(LayersFName);
    double d, c1,c2,rho;
    unsigned nsp;

    depths.clear();
    c1s.clear();
    c2s.clear();
    rhos.clear();
    Ns_points.clear();

	while (!(Myfile.eof()))
	{
		Myfile >> d;
		Myfile >> c1;
		Myfile >> c2;
		Myfile >> rho;
		Myfile >> nsp;

        depths.push_back(d);
        c1s.push_back(c1);
        c2s.push_back(c2);
        rhos.push_back(rho);
        Ns_points.push_back(nsp);

	}
	Myfile.close();

}

void sspemdd_sequential::load_profile_deep_water(
    std::string ProfileFName,
    unsigned ppm,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points)
{
    std::ifstream Myfile(ProfileFName);
    double cp, cc, dc, dp;

    Myfile >> dp;
    Myfile >> cp;

    unsigned npc;

	while (!(Myfile.eof()))
	{
		Myfile >> dc;
		Myfile >> cc;

        depths.push_back(dc);
        c1s.push_back(cp);
        c2s.push_back(cc);
        rhos.push_back(1);

        npc = (unsigned)ppm*(dc - dp);
        Ns_points.push_back(npc);

        cp = cc;
        dp = dc;

	}
	Myfile.close();


}

double sspemdd_sequential::compute_modal_delays_residual_uniform(vector<double> &freqs,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double R,
	double tau,
	vector<vector<double>> &experimental_delays,
	vector<unsigned> &experimental_mode_numbers
	)
{
	unsigned rord = 3;

	double iModesSubset = -1.0;
	double deltaf = 0.05;

	double residual = 0;
	unsigned mnumb;
	double mdelay;
	//2016.04.27:Pavel: now we use RMS as the residual
    unsigned nRes = 0;

	vector<vector<double>> modal_group_velocities;
	vector<unsigned> mode_numbers;

	compute_modal_grop_velocities(freqs, deltaf, depths, c1s, c2s, rhos, Ns_points, iModesSubset, rord, modal_group_velocities, mode_numbers);

	for (unsigned ii = 0; ii<freqs.size(); ii++) {
		//2016.04.27:Pavel: mnumb = std::min(mode_numbers.at(ii), experimental_mode_numbers.at(ii));
		mnumb = experimental_mode_numbers.at(ii);
		for (unsigned jj = 0; jj<mnumb; jj++) {
			if (experimental_delays[ii][jj]>0) {
                nRes = nRes + 1;
                //2016.04.27:Pavel:
                if ( jj<mode_numbers.at(ii) ) {
                    mdelay = R / modal_group_velocities[ii][jj];
                }
                else if ( (ii+1<freqs.size()) && (jj<mode_numbers.at(ii+1))  ) {
                    mdelay = R / modal_group_velocities[ii+1][jj];
                }
                else {
                    mdelay = 0;
                }
				//tau_comment: this is the very place where it comes into play in the computation
				//please check the search block!
				residual = residual + pow(experimental_delays[ii][jj] + tau - mdelay, 2);
			}
		}
	}
    //2016.04.27:Pavel: RMS
	residual = sqrt(residual/nRes);

	return residual;
}

// New version from 17.05.2017, group velocities computed using perturbative approach

double sspemdd_sequential::compute_modal_delays_residual_uniform2(vector<double> &freqs,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double R,
	double tau,
	vector<vector<double>> &experimental_delays,
	vector<unsigned> &experimental_mode_numbers
	)
{
	unsigned rord = 3;
	double iModesSubset = 1/sqrt(2);
	double deltaf = 0.05;
	double residual = 0;
	unsigned mnumb;
	double mdelay;
	//2016.04.27:Pavel: now we use RMS as the residual
    unsigned nRes = 0;

	vector<vector<double>> modal_group_velocities;
	vector<unsigned> mode_numbers;

	compute_modal_grop_velocities2(freqs, deltaf, depths, c1s, c2s, rhos, Ns_points, iModesSubset, rord, modal_group_velocities, mode_numbers);
	/*cout << "iModesSubset " << iModesSubset << endl;
	cout << "freqs" << endl;
	for (unsigned i = 0; i < 100; i++)
		cout << freqs[i] << " ";
	cout << endl;
	cout << "mode_numbers" << endl;
	for (unsigned i=0; i < 100; i++)
		cout << mode_numbers[i] << " ";
	cout << endl;*/
	
	for (unsigned ii = 0; ii<freqs.size(); ii++) {
		//2016.04.27:Pavel: mnumb = std::min(mode_numbers.at(ii), experimental_mode_numbers.at(ii));
		mnumb = experimental_mode_numbers.at(ii);
		for (unsigned jj = 0; jj<mnumb; jj++) {
			if (experimental_delays[ii][jj]>0) {
                nRes = nRes + 1;
                //2016.04.27:Pavel:
                if ( jj<mode_numbers.at(ii) ) {
                    mdelay = R / modal_group_velocities[ii][jj];
                }
                else {
                    mdelay = R / modal_group_velocities[ii][mode_numbers.at(ii)-1];
                }
				//tau_comment: this is the very place where it comes into play in the computation
				//please check the search block!
				residual = residual + pow(experimental_delays[ii][jj] + tau - mdelay, 2);
			}
		}
	}
    //2016.04.27:Pavel: RMS
	residual = sqrt(residual/nRes);

	if (isTimeDelayPrinting)
		printDelayTime(R, mode_numbers, modal_group_velocities);

	return residual;
}

//2016.12.31:Pavel: a residual functions where the "experimental" spectrogram modulud is taken as the weight coefficients
//this is a simplest nonuniform residual function

double sspemdd_sequential::compute_modal_delays_residual_weighted(vector<double> &freqs,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double R,
	double tau,
	vector<vector<double>> &experimental_delays,
    vector<vector<double>> &weight_coeffs,   //2016.12.31:Pavel: this is a key parameter controlling the weights
	vector<unsigned> &experimental_mode_numbers
	)
{
	unsigned rord = 3;
	double iModesSubset = -1.0;
	double deltaf = 0.05;
	double residual = 0;
	unsigned mnumb;
	double mdelay;
	//2016.04.27:Pavel: now we use RMS as the residual
    unsigned nRes = 0;

	vector<vector<double>> modal_group_velocities;
	vector<unsigned> mode_numbers;

	compute_modal_grop_velocities(freqs, deltaf, depths, c1s, c2s, rhos, Ns_points, iModesSubset, rord, modal_group_velocities, mode_numbers);

	for (unsigned ii = 0; ii<freqs.size(); ii++) {
		//2016.04.27:Pavel: mnumb = std::min(mode_numbers.at(ii), experimental_mode_numbers.at(ii));
		mnumb = experimental_mode_numbers.at(ii);
		for (unsigned jj = 0; jj<mnumb; jj++) {
			if (experimental_delays[ii][jj]>0) {
                nRes = nRes + 1;
                //2016.04.27:Pavel:
                if ( jj<mode_numbers.at(ii) ) {
                    mdelay = R / modal_group_velocities[ii][jj];
                }
                else if ( (ii+1<freqs.size()) && (jj<mode_numbers.at(ii+1))  ) {
                    mdelay = R / modal_group_velocities[ii+1][jj];
                }
                else {
                    mdelay = 0;
                }
				//tau_comment: this is the very place where it comes into play in the computation
				//please check the search block!
                                //2016.12.31:Pavel: weight coefficients are included
				residual = residual + weight_coeffs[ii][jj]*pow(experimental_delays[ii][jj] + tau - mdelay, 2);
			}
		}
	}
    //2016.04.27:Pavel: RMS
	double d = (double)(residual / (double)nRes);
	residual = sqrt(d);

	/*std::ofstream ofile("R_mgv");
	for (unsigned ii = 0; ii < freqs.size(); ii++) {
		mnumb = mode_numbers.at(ii);
		ofile << freqs.at(ii) << "\t";
		for (unsigned jj = 0; jj < mnumb; jj++)
			ofile << R / modal_group_velocities[ii][jj] << "\t";
		ofile << endl;
	}
	ofile.close();*/

	return residual;
}

int sspemdd_sequential::compute_wnumbers_bb(vector<double> &freqs,
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
	vector<double> mgv_ii;
	unsigned nwnum;
	unsigned nfr = (unsigned)freqs.size();
	double omeg1;

	for (unsigned ii = 0; ii<nfr; ii++) {
		out_wnum1.clear();
		mgv_ii.clear();
		omeg1 = 2 * LOCAL_M_PI*(freqs.at(ii) + deltaf / 2);
		out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, -1.0, ordRich);
		nwnum = (unsigned)out_wnum1.size();

		for (unsigned jj = 0; jj < nwnum; jj++)
		{
			mgv_ii.push_back(out_wnum1.at(jj));
		}

		modal_group_velocities.push_back(mgv_ii);
		mode_numbers.push_back(nwnum);
	}

	return 0;
}

int sspemdd_sequential::compute_modal_grop_velocities(vector<double> &freqs,
	double deltaf,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double iModesSubset,
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

	for (unsigned ii = 0; ii<nfr; ii++) {
		out_wnum1.clear();
		out_wnum2.clear();
		mgv_ii.clear();
		omeg1 = 2 * LOCAL_M_PI*(freqs.at(ii) + deltaf);
		out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, -1.0, ordRich);
		nwnum = (unsigned)out_wnum1.size();

		/*
		cout << "f=" << freqs.at(ii) << "Hz" << endl;

		for (unsigned jj=0; jj < nwnum; jj++ )
		{
		cout << "k_" << jj+1 << "=" << out_wnum1.at(jj) << endl;
		}
		*/

		omeg2 = 2 * LOCAL_M_PI*(freqs.at(ii));
		out_wnum2 = compute_wnumbers_extrap_lin_dz(omeg2, depths, c1s, c2s, rhos, Ns_points, -1.0, ordRich);
		//nwnum = std::min(nwnum, (unsigned)out_wnum2.size());
        nwnum = (unsigned)out_wnum2.size();

		for (unsigned jj = 0; jj < nwnum; jj++)
		{
			mgv_ii.push_back((omeg1 - omeg2) / (out_wnum1.at(jj) - out_wnum2.at(jj)));
		}

		modal_group_velocities.push_back(mgv_ii);
		mode_numbers.push_back(nwnum);
	}

	return 0;
}


int sspemdd_sequential::compute_modal_grop_velocities2(vector<double> &freqs,
	double deltaf,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double iModesSubset,
	unsigned &ordRich,
	vector<vector<double>> &modal_group_velocities,
	vector<unsigned> &mode_numbers
)
{
	mode_numbers.clear();
	modal_group_velocities.clear();

	vector<double> out_wnum;

	vector<double> mgv_ii, phi, dphi;
	unsigned nwnum;
	unsigned nfr = (unsigned)freqs.size();
	double omeg, vgc;

	for (unsigned ii = 0; ii<nfr; ii++) {
		out_wnum.clear();

		mgv_ii.clear();
		omeg = 2 * LOCAL_M_PI*freqs.at(ii);
		out_wnum = compute_wnumbers_extrap2(omeg, depths, c1s, c2s, rhos, Ns_points, iModesSubset , ordRich);
		nwnum = (unsigned)out_wnum.size();

		for (unsigned jj = 0; jj < nwnum; jj++)
		{

			    compute_wmode1(omeg, depths, c1s, c2s, rhos, Ns_points, out_wnum.at(jj), phi, dphi);

                vgc = compute_wmode_vg(omeg, depths, c1s, c2s, rhos, Ns_points, out_wnum.at(jj), phi);
                mgv_ii.push_back( vgc );
		}

		modal_group_velocities.push_back(mgv_ii);
		mode_numbers.push_back(nwnum);
	}

	return 0;
}

/*

compute_mfunctions_zr() computs the mode functions corresponding to the media parameters described by
the arrays [depths,c1s,c2s,rhos,Ns_points] for a given set of the wavenumbers (e.g. obtained from the function
compute_wnumbers_extrap() ). The functions are computed at the receiver depths from the array "zr". The values in zr are assumed to be
sorted in ascending order

*/
void sspemdd_sequential::compute_mfunctions_zr(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	vector<double> &khs,
	vector<double> &zr,
	vector<vector<double>> &mfunctions_zr)
	{
        vector<double> phi;
        vector<double> dphi;

        vector<unsigned> i_zr;
        vector<double> t_zr;


        vector<double> z, phim_zr;

        double zp;
        unsigned cur_layer = 0;
        unsigned cur_points = 0;
        unsigned nzr = zr.size();
        unsigned i_inside_l = 0;



    mfunctions_zr.clear();

    zp = 0;

    for (unsigned jj=0; jj<nzr; jj++ ){

        while ( depths.at(cur_layer)< zr.at(jj) ) {
            cur_points = cur_points + Ns_points.at(cur_layer);
            zp = depths.at(cur_layer);
            cur_layer = cur_layer + 1;
        }

        i_inside_l = (unsigned) ( (zr.at(jj) - zp)*Ns_points.at(cur_layer)/( depths.at(cur_layer) - zp ) );
        i_inside_l = std::min(i_inside_l , Ns_points.at(cur_layer) );
        i_zr.push_back( cur_points + i_inside_l );
        t_zr.push_back( (zr.at(jj) - zp)*Ns_points.at(cur_layer)/( depths.at(cur_layer) - zp ) - i_inside_l );



        // For non-ordered set of zr!!!! Slows the interpolation down!
        cur_layer = 0;
        cur_points = 0;
    }


    for (unsigned ii = 0; ii < khs.size(); ii++) {
		double kh = khs.at(ii);
		compute_wmode1(omeg, depths, c1s, c2s, rhos, Ns_points, kh, phi, dphi);

        phim_zr.clear();
		for (unsigned jj=0; jj<nzr; jj++ ) {

            phim_zr.push_back(  (1-t_zr.at(jj))*phi.at( i_zr.at(jj)) + t_zr.at(jj)*phi.at(i_zr.at(jj)+1)   );

		}

        mfunctions_zr.push_back( phim_zr );
	}

//    //mfunctions output to a file
//    std::ofstream ofile("mfunctionszr.txt");
//
//    for (unsigned jj = 0; jj < nzr; jj++) {
//			ofile << zr.at(jj) << " ";
//	}
//	ofile << endl;
//
//	for (unsigned ii = 0; ii < khs.size(); ii++) {
//		for (unsigned jj = 0; jj < nzr; jj++)
//			ofile << mfunctions_zr[ii][jj] << " ";
//		ofile << endl;
//	}
//	ofile.close();
}

vector<std::complex<double>> sspemdd_sequential::compute_cpl_pressure(double f,
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	vector<double> &Rr,
	vector<double> &zr,
    double iModesSubset,
	unsigned ordRich)
	{

        vector<vector<double>> modefunctions;
        vector<double> khs;
        vector<std::complex<double>> PHelm;


        double omeg = 2 * LOCAL_M_PI*f;
        double R;
        unsigned nzr = zr.size();
        std::complex<double> Prc;


        khs = compute_wnumbers_extrap_lin_dz(omeg, depths, c1s, c2s, rhos, Ns_points, iModesSubset, ordRich);

        unsigned nmod = khs.size();


//        // TEST
//        cout << nmod << " modes " << endl;
//
//        // TEST
//        //for (unsigned qq=0; qq<nmod; qq++){
//        for (unsigned qq=0; qq<20; qq++){
//
//            cout  << " kh" << qq <<" = " << khs.at(qq) << endl;
//
//        }


        if (nmod>0) {



            compute_mfunctions_zr(omeg, depths, c1s, c2s, rhos, Ns_points, khs,zr, modefunctions);

//            // TEST
//            for (unsigned ss=0; ss<20; ss++){
//            //for (unsigned ss=0; ss<nmod; ss++){
//
//                    cout  << " mf" << ss << " = " << modefunctions.at(ss).at(0) << endl;
//
//            }


            //PHelm(1:nr) = PHelm(1:nr) + psiz*psizs*(1i*exp(-1i*pi/4)./sqrt(8*pi*R) ).*exp(1i*wnum(mm)*R)/sqrt(wnum(mm));

            //modefunctions -- vector of vectors, that represent the values of certain mode at all zr
            for (unsigned ii = 1; ii < nzr; ii++) {
                Prc = std::complex<double>(0.0,0.0);
                R = Rr.at(ii);

                for (unsigned jj = 0; jj < nmod; jj++) {
                    //PHelm.back() = PHelm.back() +

                    /*
                    //TEST
                    cout << "phi1zr=" << modefunctions.at(jj).at(ii) << endl;
                    cout << "phi1zs=" << modefunctions.at(jj).at(0) << endl;
                    cout << "R=" << R << endl;
                    cout << "k=" << khs.at(jj) << endl;
                    cout << "exp=" << exp( Iu*khs.at(jj)*R ) << endl;
                    */

                    Prc = Prc + exp( Iu*khs.at(jj)*R )*modefunctions.at(jj).at(ii)*modefunctions.at(jj).at(0)/sqrt(khs.at(jj));
                }
                Prc = Iu*exp(-Iu*LOCAL_M_PI/4.0)*Prc/sqrt(8*LOCAL_M_PI*R);
                PHelm.push_back( Prc );
            }

        } else {
            // TEST
            cout << 0 << " modes " << endl;

            Prc = std::complex<double>(0.0,0.0);
            for (unsigned ii = 1; ii < nzr; ii++) {
                PHelm.push_back( Prc );
            }
        }


        return PHelm;

	}



/*

compute_all_mfunctions() computs the mode functions corresponding to the media parameters described by
the arrays [depths,c1s,c2s,rhos,Ns_points] for a given set of the wavenumbers (e.g. obtained from the function
compute_wnumbers_extrap() ). The functions are written to the file "mfunctions.txt" line-by-line, the first line contains
the values of depth

*/
void sspemdd_sequential::compute_all_mfunctions(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	vector<double> &khs)
	{
        vector<double> phi;
        vector<double> dphi;

        vector<double> z;
        double h,z0;


    std::ofstream ofile("mfunctions.txt");



	for (unsigned ii = 0; ii < depths.size(); ii++) {

		if (ii>0){
            h = (depths.at(ii) - depths.at(ii-1))/Ns_points.at(ii);
            z0 = depths.at(ii-1);
		}
		else {
		    h = (depths.at(ii))/Ns_points.at(ii);
		    z0 = 0;
		}

		for (unsigned jj = 0; jj < Ns_points.at(ii); jj++){
            z.push_back( z0 + h*jj );
            ofile << z.back() << " ";
		}
	}
	z.push_back( depths.back() );
	ofile << z.back() << endl;
    ofile << endl;

	for (unsigned i = 0; i < khs.size(); i++) {
		double kh = khs[i];
		compute_wmode1(omeg, depths, c1s, c2s, rhos, Ns_points, kh, phi, dphi);
		for (unsigned j = 0; j < phi.size(); j++)
			ofile << phi[j] << " ";
		ofile << endl;
	}
	ofile.close();


	}


/*

compute_wmode() computs the mode function "phi" and its derivative "dphi" for media parameters described by
the arrays [depths,c1s,c2s,rhos,Ns_points] and for a given horizontal wavenumber "kh". The functions are normalized in
the standard way (using inverse density as a weight function). The Runge-Kutta (4th order) scheme is used for solving the
ODE in the matrix formulation. Due to the round-off errors it is unstable in the bottom for the discrete spectrum modes, as
the ODE is stiff there (solution involves a decaying exponential).

*/

void sspemdd_sequential::compute_wmode(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double kh,
	vector<double> &phi,
	vector<double> &dphi)
{
    double deltah;
    double phi2int =0;
    double layer_int =0;

    phi.clear();
    dphi.clear();

    phi.push_back(0.0);
    dphi.push_back(1.0);


    unsigned n_layers = (unsigned)depths.size();

    for (unsigned ll=0; ll<n_layers; ll++) {

        if (ll == 0 ){ deltah = depths.at(0);}
        else {
            deltah = depths.at(ll) - depths.at(ll-1);
            dphi.back() = rhos.at(ll)*dphi.back()/rhos.at(ll-1);
        }

        if ((ll==n_layers-1) && ( omeg/c1s.at(ll) < kh  )) {
            // in the bottom layer where the mode functions exhibits decay use
            layer_int = Layer_an_exp(omeg*omeg, kh*kh, deltah, c1s.at(ll), Ns_points.at(ll) , phi, dphi);
        } else {
            // use the Runge-Kutta 4th order scheme in regular layers
            layer_int = RK4(omeg*omeg, kh*kh, deltah, c1s.at(ll), c2s.at(ll), Ns_points.at(ll) , phi, dphi);
        }

        phi2int = phi2int + layer_int/rhos.at(ll);

        // TEST
        //cout << layer_int << endl;
        // TEST
    }

    double phiNorm = sqrt(phi2int);
    unsigned nz = (unsigned)phi.size();

    // TEST

    //cout << phi2int << endl;
    //cout << phiNorm << endl;
    // TEST

    for (unsigned ll = 0; ll < nz; ll++ ) {
        phi.at(ll) = phi.at(ll)/phiNorm;
        dphi.at(ll) = dphi.at(ll)/phiNorm;
    }

}

/*

compute_wmode1() computs the mode function "phi" and its derivative "dphi" for media parameters described by
the arrays [depths,c1s,c2s,rhos,Ns_points] and for a given horizontal wavenumber "kh". The functions are normalized in
the standard way (using inverse density as a weight function). The Runge-Kutta (4th order) scheme is used for solving the
ODE in the matrix formulation.

The parallel shooting is used in this version. For the discrete spectrum modes (refracted-refracted modes). We find the last layer L
(from the bottom) where k(z) <  kh  (hence we have the decaying exponent). Starting from the bottom layer we solve the ODE to L-1 th layer
going in the negative direction of z axis. Then we match the solutions coming from the top and the bottom.

*/

void sspemdd_sequential::compute_wmode1(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double kh,
	vector<double> &phi,
	vector<double> &dphi)
{
    double deltah,phiNorm;
    double phi2int =0.0;
    double bphi2int =0.0;

    double layer_int =0.0;

    vector<double> bphi, bdphi;

    phi.clear();
    dphi.clear();

    phi.push_back(0.0);
    dphi.push_back(1.0);




    unsigned n_layers = (unsigned)depths.size();
    unsigned L = n_layers;

    while ((  kh > omeg/( std::min( c1s.at(L-1), c2s.at(L-1) ) )  ) && (L>1) ) {
        L = L - 1;
    }



    // shooting from the surface
    for (unsigned ll=0; ll<L; ll++) {

        if (ll == 0 ){ deltah = depths.at(0);}
        else {
            deltah = depths.at(ll) - depths.at(ll-1);
            dphi.back() = rhos.at(ll)*dphi.back()/rhos.at(ll-1);
        }

        layer_int = RK4(omeg*omeg, kh*kh, deltah, c1s.at(ll), c2s.at(ll), Ns_points.at(ll) , phi, dphi);

        phi2int = phi2int + layer_int/rhos.at(ll);

        //TEST --> overflow problem
        if (phi2int >  1.0e+50) {
            for (unsigned qq=0; qq<phi.size(); qq++) {
                phi.at(qq) = phi.at(qq)/(1.0e+20);
                dphi.at(qq) = dphi.at(qq)/(1.0e+20);
            }
            phi2int = phi2int/(1.0e+40);
        }

    }

    // shooting from the bottom

    if (L < n_layers){

        bphi.push_back(0.0);
        bdphi.push_back(1.0);

        for (unsigned ll=n_layers-1; ll>=L; ll--) {

            deltah = depths.at(ll) - depths.at(ll-1);

            layer_int = RK4(omeg*omeg, kh*kh, deltah, c2s.at(ll), c1s.at(ll), Ns_points.at(ll) , bphi, bdphi);
            bphi2int = bphi2int + layer_int/rhos.at(ll);


//            //TEST
//            cout  << "ll="<< ll <<" bphi2int " << bphi2int  << endl;
//            cout  << "ll="<< ll <<" layer_int " << layer_int  << endl;
//            cout  << "ll="<< ll <<" deltah " << deltah  << endl;
//            cout  << "ll="<< ll <<" kh " << kh  << endl;
//            cout  << "ll="<< ll <<" c2 " << c2s.at(ll)  << endl;
//            cout  << "ll="<< ll <<" c1 " << c1s.at(ll)  << endl;
//            cout  << "ll="<< ll <<" nsp " << Ns_points.at(ll)  << endl;
//            cout  << "ll="<< ll <<" nsp-1 " << Ns_points.at(ll-1)  << endl;
//            cout  << "ll="<< ll <<" nsp-2 " << Ns_points.at(ll-2)  << endl;
//            for (unsigned qq=0; qq<bphi.size(); qq++  ){
//                cout  << "qq="<< qq <<" bphi " << bphi.at(qq) <<" bdphi " << bdphi.at(qq)  << endl;
//            }

            bdphi.back() = rhos.at(ll-1)*bdphi.back()/rhos.at(ll);

            //TEST --> overflow problem
            if (bphi2int >  1.0e+50) {

                for (unsigned qq=0; qq<bphi.size(); qq++) {
                    bphi.at(qq) = bphi.at(qq)/(1.0e+20);
                    bdphi.at(qq) = bdphi.at(qq)/(1.0e+20);
                }
                bphi2int = bphi2int/(1.0e+40);
            }

        }

        // matching the shooting solutions

        double cmatching = phi.back()/bphi.back();



        for (int ll = bphi.size()-2; ll>=0; ll-- ){

            phi.push_back( cmatching*bphi.at(ll) );
            dphi.push_back( cmatching*bdphi.at(ll) );


        }

        cmatching = cmatching*cmatching;
        phiNorm = sqrt(phi2int + bphi2int*cmatching);

    } else {

        phiNorm = sqrt(phi2int);

    }
//    //TEST
//    cout  << " phiNorm " << phiNorm  << endl;
//    cout  << " phi2int " << phi2int  << endl;
//    cout  << " bphi2int " << bphi2int  << endl;
//    cout  << " cmatching " << phi.back()/bphi.back()  << endl;
//    cout  << " L " << L  << endl;
//    cout  << " n_layers " << n_layers  << endl;
//
//    throw std::invalid_argument("Ururu");

//        // TEST
//        for (unsigned qq=0; qq<nmod; qq++){
//
//            cout  << " kh" << qq <<" = " << khs.at(qq) << endl;
//
//        }


    unsigned nz = (unsigned)phi.size();



    for (unsigned ll = 0; ll < nz; ll++ ) {
        phi.at(ll) = phi.at(ll)/phiNorm;
        dphi.at(ll) = dphi.at(ll)/phiNorm;
    }

}



double sspemdd_sequential::compute_wmode_vg(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double kh,
	vector<double> &phi)
{
    double deltah, h, cc, termc, vg;
    double vg_int =0.0;
    double layer_int;
    unsigned Np, nc;
    
    unsigned n_layers = (unsigned)depths.size();

    nc = 0;

    for (unsigned ll=0; ll<n_layers; ll++) {

        if (ll == 0 ){ deltah = depths.at(0);}
        else {
            deltah = depths.at(ll) - depths.at(ll-1);
        }

        Np = Ns_points.at(ll);
        h = deltah/Np;

        cc = c1s.at(ll);
        layer_int = 0.0;

        termc = phi.at(nc)*phi.at(nc)/( cc*cc);

        for (unsigned jj = 0; jj < Np; jj++){
            layer_int = layer_int + termc;
            nc = nc + 1;
            cc = c1s.at(ll) + ( c2s.at(ll)-c1s.at(ll) )*(jj+1)*h/deltah;
            termc = phi.at(nc)*phi.at(nc)/( cc*cc );
            layer_int = layer_int + termc;
        }

        layer_int = 0.5*layer_int*h/rhos.at(ll);

        vg_int = vg_int + layer_int;


    }

    vg = 1/( omeg*vg_int/kh );

    return vg;
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
vector<double> sspemdd_sequential::compute_wnumbers_extrap(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double iModesSubset,
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
	iModesSubset -- controls the subset of the modes computed: iModesSubset <0 -> trapped modes only; 0<=iModesSubset<1 -> a subset of modes with wavenumbers within [iModesSubset*kmax kmax], in particular iModesSubset = 0 computes all propagating modes
	ordRich -- order of the Richardson extrapolation;

	the top of the first layer is z=0
	*/
{
	vector<double> coeff_extrap;
	switch (ordRich) {
	case 1:
		coeff_extrap.assign({ 1 });
		break;
	case 2:
		coeff_extrap.assign({ -1, 2 });
		break;
	case 4:
		coeff_extrap.assign({ -1 / double(6), 4, -13.5, 32 / double(3) });
		break;
	default:
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

	unsigned m_wnum = 100000;

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
			dz = (zc - zp) / (n_points_layer);
			input_mesh.push_back(dz);
			input_c.at(n_points_total - 1) = c1s.at(ll);
			input_rho.at(n_points_total - 1) = rhos.at(ll);

			n_points_total = n_points_total + n_points_layer;

			for (unsigned kk = 1; kk <= n_points_layer; kk++) {
				input_rho.push_back(rhos.at(ll));
				input_c.push_back(c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk / n_points_layer);
			}

			if (ll < n_layers - 1) {
				input_interf_idcs.push_back(n_points_total - 1);
			}
			zp = zc;
		}

		//cout << "rr=" << rr << endl;

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh, iModesSubset);
		m_wnum = std::min(m_wnum, (unsigned)out_wnum2.size());

		if (rr == 1) { wnum_extrapR.assign(m_wnum, 0); }

		for (unsigned mm = 0; mm<m_wnum; mm++) {
			wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr - 1);
		}

		/*
		for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
		cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
		}
		*/
	}

	for (unsigned mm = 0; mm<m_wnum; mm++) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}

	return wnum_extrapR;
}



vector<double> sspemdd_sequential::compute_wnumbers_extrap2(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double iModesSubset,
	unsigned &ordRich)
	/*
	14.05.17 Pavel
	subroutine for computing wavenumbers for a given waveguide structure
	the computation is performed by the FD method for certain meshsize,
	Richardson extrapolation is used to improve the accuracy

    This version features the meshes with 1 x, 2 x, 4 x, 8 x Ns_points numbers of points

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
	iModesSubset -- controls the subset of the modes computed: iModesSubset <0 -> trapped modes only; 0<=iModesSubset<1 -> a subset of modes with wavenumbers within [iModesSubset*kmax kmax], in particular iModesSubset = 0 computes all propagating modes
	ordRich -- order of the Richardson extrapolation;

	the top of the first layer is z=0
	*/
{
	vector<double> coeff_extrap;
	switch (ordRich) {
	case 1:
		coeff_extrap.assign({ 1 });
		break;
	case 2:
		coeff_extrap.assign({ -0.333333333333333, 1.333333333333333});
		break;
	case 4:
		coeff_extrap.assign({ -0.000352733686067, 0.029629629629630, -0.474074074074074, 1.444797178130511 });
		break;
	default:
		ordRich = 3;
		coeff_extrap.assign({ 0.022222222222222, -0.444444444444444, 1.422222222222222 });
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
	unsigned m_wnum = 100000;

	unsigned n_layers = (unsigned)depths.size();
	unsigned n_points_total = 0;
	unsigned n_points_layer = 0;

    unsigned r_mult = 1; //Richardson factor for the number of points
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
			n_points_layer = Ns_points.at(ll)*r_mult;
			dz = (zc - zp) / (n_points_layer);
			input_mesh.push_back(dz);
			input_c.at(n_points_total - 1) = c1s.at(ll);
			input_rho.at(n_points_total - 1) = rhos.at(ll);

			n_points_total = n_points_total + n_points_layer;

			for (unsigned kk = 1; kk <= n_points_layer; kk++) {
				input_rho.push_back(rhos.at(ll));
				input_c.push_back(c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk / n_points_layer);
			}

			if (ll < n_layers - 1) {
				input_interf_idcs.push_back(n_points_total - 1);
			}
			zp = zc;
		}

		//cout << "rr=" << rr << endl;

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh, iModesSubset);
		m_wnum = std::min(m_wnum, (unsigned)out_wnum2.size());

		if (rr == 1) { wnum_extrapR.assign(m_wnum, 0); }

		for (unsigned mm = 0; mm<m_wnum; mm++)
			wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr - 1);

		r_mult = r_mult*2;

		/*
		for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
		cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
		}
		*/
	}

	for (unsigned mm = 0; mm<m_wnum; mm++)
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));

	return wnum_extrapR;
}




vector<double> sspemdd_sequential::compute_wnumbers_extrap_lin_dz(double &omeg, // sound frequency
	vector<double> &depths,
	vector<double> &c1s,
	vector<double> &c2s,
	vector<double> &rhos,
	vector<unsigned> &Ns_points,
	double iModesSubset,
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
	iModesSubset -- controls the subset of the modes computed: iModesSubset <0 -> trapped modes only; 0<=iModesSubset<1 -> a subset of modes with wavenumbers within [iModesSubset*kmax kmax], in particular iModesSubset = 0 computes all propagating modes
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

	unsigned m_wnum = 100000;

	unsigned n_layers = (unsigned)depths.size();
	unsigned n_points_total = 0;
	unsigned n_points_layer = 0;

	vector<double> coeff_extrap;
	switch (ordRich) {
	case 1:
		coeff_extrap.assign({ 1 });
		break;
	case 2:
		coeff_extrap.assign({ 1.333333333333333, -0.333333333333333 });
		break;
	case 4:
		// coeff_extrap.assign({1.595325630252102, -0.788449059052564, 0.216346153846154, -0.023222725045691});
		coeff_extrap.assign({ 1.6, -0.8, 0.228571428571429, -0.028571428571429 });
		break;
	default:
		ordRich = 3;
		coeff_extrap.assign({ 1.5, -0.6, 0.1 });
		//coeff_extrap.assign({0.1, -0.6, 1.5});
		break;
	}

	// number of points in each layer is multiple of 12
	// this allows us to use nz_ii = nz/ii, ii = 1,2,3,4
	for (unsigned ii = 0; ii < n_layers; ii++){
        if ( Ns_points.at(ii) % 12 > 0 ){
            Ns_points.at(ii) = 12 * (1+(Ns_points.at(ii)/12));
        }
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
			n_points_layer = Ns_points.at(ll) / rr;
			dz = (zc - zp) / (n_points_layer);

			//            cout << "np=" << n_points_layer << "  " << "dz=" << dz <<endl;

			input_mesh.push_back(dz);
			input_c.at(n_points_total - 1) = c1s.at(ll);
			input_rho.at(n_points_total - 1) = rhos.at(ll);

			n_points_total = n_points_total + n_points_layer;

			for (unsigned kk = 1; kk <= n_points_layer; kk++) {
				input_rho.push_back(rhos.at(ll));
				input_c.push_back(c1s.at(ll) + (c2s.at(ll) - c1s.at(ll))*kk / n_points_layer);
			}

			if (ll < n_layers - 1) {
				input_interf_idcs.push_back(n_points_total - 1);
			}
			zp = zc;
		}

		//        cout << "rr=" << rr << endl;

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh, iModesSubset);
		m_wnum = std::min(m_wnum, (unsigned)out_wnum2.size());

		if (rr == 1) { wnum_extrapR.assign(m_wnum, 0); }

		for (unsigned mm = 0; mm<m_wnum; mm++) {
			wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr - 1);
		}

		//            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
		//
		//                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
		//            }
		//            cout << endl;
	}

	for (unsigned mm = 0; mm<m_wnum; mm++) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}
	return wnum_extrapR;
}

vector<double> sspemdd_sequential::compute_wnumbers(double &omeg, // sound frequency
	vector<double> &c,
	vector<double> &rho,
	vector<unsigned> &interface_idcs,
	vector<double> &meshsizes,
	double iModesSubset                 // set iModesSubset <0 -> trapped modes only; 0<=iModesSubset<1 -> a subset of modes with wavenumbers
	)                                   // within [iModesSubset*kmax kmax], in particular iModesSubset = 0 computes all propagating modes
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

	if (interface_idcs.size() > 0)
	{
		next_interface_idx = interface_idcs.at(0) - 1;
	}
	else
	{
		next_interface_idx = N_points;
	}

	cmin = c.at(0);
	cmax = c.at(0);

	for (int ii = 0; ii < N_points; ii++){
		if (c.at(ii) < cmin) { cmin = c.at(ii); }
		if (c.at(ii) > cmax) { cmax = c.at(ii); }
	}

	kappamax = omeg / cmin;
	kappamin = omeg / cmax;

	if (iModesSubset >= 0) {

        if (iModesSubset < 1) {
            kappamin = iModesSubset*kappamax;
        }
        else {
            throw std::invalid_argument("Invalid iModeSubset: use either -1 or a value from [0 1)");
        }

	}




	for (int ii = 0; ii < N_points - 2; ii++){
		// ordinary point
		ud.push_back(1 / (dz*dz));
		ld.push_back(1 / (dz*dz));
		md.push_back(-2 / (dz*dz) + omeg*omeg / (c.at(ii + 1)*c.at(ii + 1)));

		// special case of the point at the interface

		if (ii == next_interface_idx) {         //ii -- z(ii+1), z(0) = 0
			layer_number = layer_number + 1;    // âîîáùå ii=89 -- âîäà, â ii=90 -äíî,
			// çäåñü ii = 89 -- èíòåðôåéñ, óæå äíî
			cp = c.at(ii + 1);
			dp = rho.at(ii + 1);
			cm = c.at(ii);
			dm = rho.at(ii);


			dz_next = meshsizes.at(layer_number);
            q = 1 / (dz_next*dm + dz*dp);  // Ïîìåíÿòü ìåñòàìè ñ ïðåäûäóùåé ñòðîêîé??

			ld.at(ii) = 2 * q*dp / dz;
			md.at(ii) = -2 * q*(dz_next*dp + dz*dm) / (dz*dz_next) + omeg*omeg*q*(dz*dp*cp*cp + dz_next*dm*cm*cm) / (cp*cp*cm*cm);
			ud.at(ii) = 2 * q*dm / dz_next;

			if (interface_idcs.size() > layer_number)
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

	alglib::real_2d_array eigenvectors; // V - ñîáñòâ âåêòîð
	alglib::real_1d_array eigenvalues; // Lm -ñîáñòâ çíà÷
	eigenvalues.setlength(N_points - 2);
	alglib::ae_int_t eigen_count = 0;

	alglib::real_1d_array main_diag, second_diag;
	main_diag.setlength(N_points - 2);
	second_diag.setlength(N_points - 3);

	for (int ii = 0; ii < N_points - 3; ii++) {
		second_diag[ii] = sqrt(ud.at(ii)*ld.at(ii + 1));
		main_diag[ii] = md.at(ii);
	}
	main_diag[N_points - 3] = md.at(N_points - 3);

/* TEST: the sparse matrix diagonals output
	    std::ofstream myFile("thematrixdiags.txt");
	    for (int ii=0; ii<N_points-2; ii++){
	        myFile << std::fixed  << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
	    }
	    myFile.close();

	    std::ofstream myFile1("thematrixdiags_sym.txt");
	    for (int ii=0; ii<N_points-3; ii++){
	        myFile1 << std::fixed  << main_diag[ii] << "  " << second_diag[ii] << endl;
	    }
	    myFile1.close();
*/
	alglib::smatrixtdevdr(main_diag, second_diag, N_points - 2, 0, kappamin*kappamin, kappamax*kappamax, eigen_count, eigenvectors);

	for (int ii = 0; ii<eigen_count; ii++) {
		wnumbers2.push_back(main_diag[eigen_count - ii - 1]);
	}

	return wnumbers2;
}

//tau_comment: search and output,
//check for tau!
int sspemdd_sequential::init(vector<double> depths)
{
	if ((!rank) && (verbosity > 0)) {
		cout << "init() started" << endl;
		cout << "depths: ";
		for (auto &x : depths)
			cout << x << " ";
		cout << endl;
	}
	
	n_layers_w = depths.size() - 1;
	
	if (!n_layers_w) {
		cerr << "n_layers_w == 0" << endl;
		return -1;
	}

	ncpl_arr = ncpl_init_arr;
	ncpl_arr.resize(n_layers_w);
	cw1_arr = cw1_init_arr;
	cw1_arr.resize(n_layers_w);
	cw2_arr = cw2_init_arr;
	cw2_arr.resize(n_layers_w);

	c1s.resize(n_layers_w + 1);
	for (auto &x : c1s)
		x = 1500;
	c2s.resize(n_layers_w + 1);
	for (auto &x : c2s)
		x = 1500;
	rhos.resize(n_layers_w + 1);
	for (auto &x : rhos)
		x = 1;
	Ns_points.resize(n_layers_w + 1);
	unsigned ppm = 1;
	Ns_points[0] = (unsigned)round(ppm*depths[0]);
	for (unsigned i=1; i < depths.size(); i++ )
		Ns_points[i] = (unsigned)round(ppm*(depths[i] - depths[i-1]));

	N_total = nR*nrhob*ncb*ntau;
	for (auto &x : ncpl_arr)
		N_total *= (unsigned long long)x;

	if (!N_total) {
		cerr << "N_total == 0" << endl;
		cerr << "nR nrhob ncb ntau : " << endl;
		cerr << nR << " " << nrhob << " " << ncb << " " << ntau << endl;
		cerr << "ncpl_arr : " << endl;
		for (auto &x : ncpl_arr)
			cerr << x << " ";
		cerr << endl;
		cerr << "depths : " << endl;
		for (auto &x : depths)
			cerr << x << " ";
		cerr << endl;
		return -1;
	}
	
	if ( (!rank) && (verbosity > 0) )
		cout << "N_total " << N_total << endl;
	
	if (cw1_arr.size() != (depths.size() - 1)) {
		cerr << "cw1_arr.size() != (depths.size() - 1)";
		cerr << endl;
		exit(1);
	}
	loadValuesToSearchSpaceVariables();
	
	if ( (!rank) && (verbosity > 0) )
		cout << "init() finished" << endl;

	return 0;
}

int sspemdd_sequential::createDepthsArray(vector<vector<double>> &depths_vec)
{
	if (d1_arr.size() == 0) {
		n_layers_w = cw1_arr.size();
		double layer_thickness_w = h / n_layers_w;
		vector<double> depths;
		for (unsigned jj = 1; jj <= n_layers_w; jj++)
			depths.push_back(layer_thickness_w*jj);
		depths.push_back(H);
		depths_vec.push_back(depths);
	}
	else {
		vector<vector<double>> search_space_depths;
		search_space_depths.resize(d1_arr.size());
		for (unsigned i = 0; i < d2_arr.size(); i++) {
			double cur_val = d2_arr[i];
			for (;;) {
				search_space_depths[i].push_back(cur_val);
				cur_val -= d_step[i];
				if (cur_val < d1_arr[i])
					break;
			}
		}

		vector<int> index_arr;
		vector<double> tmp_depths;
		vector<vector<double>> ::iterator it;
		double cur_treshold;
		while (SSPEMDD_utils::next_cartesian(search_space_depths, index_arr, tmp_depths)) {
			vector<double> depths;
			cur_treshold = tmp_depths[0] + 3;
			depths.push_back(tmp_depths[0]); // at least 1 water layer must exist
			for (unsigned i = 1; i < tmp_depths.size(); i++) {
				if (tmp_depths[i] >= cur_treshold) {
					depths.push_back(tmp_depths[i]);
					cur_treshold = tmp_depths[i] + 2;
				}
			}
			it = find(depths_vec.begin(), depths_vec.end(), depths);
			if (it == depths_vec.end())
				depths_vec.push_back(depths);
		}

		for (auto &x : depths_vec) {
			x.push_back(h);
			x.push_back(H);
		}
	}
	cout << "depths_vec.size() " << depths_vec.size() << endl;
	
	return 0;
}

void sspemdd_sequential::reportFinalResult()
{
	// fix final time
	std::chrono::high_resolution_clock::time_point t2;
	std::chrono::duration<double> time_span;
	t2 = std::chrono::high_resolution_clock::now();
	time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - start_chrono_time);

	cout << endl;
	cout << "total solving time " << time_span.count() << endl;
	cout << "SEARCH ENDED!" << endl;
	cout << "RESULTING VALUE:" << endl;
	cout << "err = " << record_point.residual << ", parameters:" << endl;
	cout << "c_b = " << record_point.cb << endl
			  << "tau = " << record_point.tau << endl
			  << "rho_b = " << record_point.rhob << endl
			  << "R = " << record_point.R << endl;
	cout << "cws :" << endl;
	for (auto &x : record_point.cws)
		cout << x << " ";
	cout << endl;
	cout << "depths " << endl;
	for (auto &x : record_point.depths)
		cout << x << " ";
	cout << endl;
	cout << "total solving time " << time_span.count() << endl;
}

void sspemdd_sequential::findGlobalMinBruteForce(vector<double> depths)
{
	cout << "findGlobalMinBruteForce()" << endl;

	vector<search_space_point> search_space_points_vec = getSearchSpacePointsVec(depths);
	cout << "search_space_points_vec.size() " << search_space_points_vec.size() << endl;
	
	for (auto &x : search_space_points_vec)
		fillDataComputeResidual(x); // calculated residual is written to cur_point
}

vector<search_space_point> sspemdd_sequential::getSearchSpacePointsVec(vector<double> depths)
{
	vector<int> index_arr;
	vector<unsigned> cur_point_indexes;
	vector<vector<unsigned>> search_space_indexes;
	search_space_indexes.resize(search_space.size());
	for (unsigned i = 0; i < search_space.size(); i++)
		for (unsigned j = 0; j < search_space[i].size(); j++)
			search_space_indexes[i].push_back(j);

	vector<search_space_point> points_vec;
	while (SSPEMDD_utils::next_cartesian(search_space_indexes, index_arr, cur_point_indexes))
		points_vec.push_back(fromPointIndexesToPoint(cur_point_indexes, depths));

	return points_vec;
}

void sspemdd_sequential::reduceSearchSpace(reduced_search_space_attribute &reduced_s_s_a)
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4, ...] - cws
	// search_space_variables[..., ...] - depths
	if (reduced_s_s_a.cb == false) {
		search_space[0].resize(1);
		search_space[0][0] = cb1;
	}
	if (reduced_s_s_a.rhob == false) {
		search_space[1].resize(1);
		search_space[1][0] = rhob1;
	}
	if (reduced_s_s_a.R == false) {
		search_space[2].resize(1);
		search_space[2][0] = R1;
	}
	if (reduced_s_s_a.tau == false) {
		search_space[3].resize(1);
		search_space[3][0] = tau1;
	}
	for (unsigned i=0; i < reduced_s_s_a.cws.size(); i++) {
		if (reduced_s_s_a.cws[i] == false) {
			search_space[4 + i].resize(1);
			search_space[4 + i][0] = cw1_arr[i];
		}
	}
}

double sspemdd_sequential::fillDataComputeResidual( search_space_point &point )
{ // finally specify sound speed in water
  // the parameters are transformed into the arrays c1s, c2s, rhos
	for (unsigned jj = 0; jj < n_layers_w - 1; jj++) {
		c1s.at(jj) = point.cws.at(jj);
		c2s.at(jj) = point.cws.at(jj + 1);
		rhos.at(jj) = 1;
	}
	c1s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	c2s.at(n_layers_w - 1) = point.cws.at(n_layers_w - 1);
	rhos.at(n_layers_w - 1) = 1;
	c1s.at(n_layers_w) = point.cb;
	c2s.at(n_layers_w) = point.cb;
	rhos.at(n_layers_w) = point.rhob;
	vector<double> depths = point.depths;
	if (depths.size() == 0) {
		std::cerr << "depths.size() == 0" << endl;
		exit(-1);
	}
	
	if (verbosity > 1) {
		/*for (unsigned jj = 0; jj <= n_layers_w; jj++)
			cout << "Layer #" << jj + 1 << ": c=" << c1s.at(jj) << "..." << c2s.at(jj) << "; rho=" << rhos.at(jj) << "; np=" << Ns_points.at(jj) << endl;
		cout << residual << endl << endl;*/
		cout << "depths : ";
		for (auto &x : depths)
			cout << x << " ";
		cout << endl;
		cout << "Ns_points : ";
		for (auto &x : Ns_points)
			cout << x << " ";
		cout << endl;
	}

	if (object_function_type == "uniform") {
		point.residual = compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, mode_numbers);
	}
	else if (object_function_type == "uniform2") {
		point.residual = compute_modal_delays_residual_uniform2(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, mode_numbers);
	}
	else if (object_function_type == "weighted") {
		point.residual = compute_modal_delays_residual_weighted(freqs, depths, c1s, c2s, rhos, Ns_points,
			point.R, point.tau, modal_delays, weight_coeffs, mode_numbers);
	}
	else {
		std::cerr << "unknown object_function_type " << object_function_type << endl;
		exit(1);
	}

	if ( verbosity > 1 )
		cout << "point.residual " << point.residual << endl;
	
	if (point.residual < record_point.residual) {
		record_point = point;
		if (verbosity > 0) {
			cout << endl;
			cout << endl << "New residual minimum:" << endl;
			cout << "err = " << record_point.residual << ", parameters:" << endl;
			cout << "c_b = " << record_point.cb
				<< ", rho_b = " << record_point.rhob
				<< ", tau = " << record_point.tau
				<< ", R = " << record_point.R << endl;
			cout << "cws_min :" << endl;
			for (auto &x : record_point.cws)
				cout << x << " ";
			cout << endl;
			cout << "depths " << endl;
			for (auto &x : record_point.depths)
				cout << x << " ";
			cout << endl;
			cout << "Ns_points " << endl;
			for (auto &x : Ns_points)
				cout << x << " ";
			cout << endl;
			cout << endl;
		}
	}

	return point.residual;
}

void sspemdd_sequential::loadValuesToSearchSpaceVariables()
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws
	vector<double> tmp_vec;
	search_space.clear();

	// fill search_space_variables[0] with cb
	tmp_vec.resize(ncb);
	for (unsigned long long i = 0; i < ncb; i++)
		tmp_vec[i] = cb1 + (ncb == 1 ? 0 : i*(cb2 - cb1) / (ncb - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[1] with rhob
	tmp_vec.resize(nrhob);
	for (unsigned long long i = 0; i < nrhob; i++)
		tmp_vec[i] = rhob1 + (nrhob == 1 ? 0 : i*(rhob2 - rhob1) / (nrhob - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[2] with R
	tmp_vec.resize(nR);
	for (unsigned long long i = 0; i < nR; i++)
		tmp_vec[i] = R1 + (nR == 1 ? 0 : i*(R2 - R1) / (nR - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[3] with tau
	tmp_vec.resize(ntau);
	for (unsigned long long i = 0; i < ntau; i++)
		tmp_vec[i] = tau1 + (ntau == 1 ? 0 : i*(tau2 - tau1) / (ntau - 1));
	search_space.push_back(tmp_vec);

	// fill search_space_variables[4-...] with cws
	for (unsigned long long i = 0; i < cw1_arr.size(); i++) {
		tmp_vec.resize(ncpl_arr[i]);
		for (unsigned long long j = 0; j < ncpl_arr[i]; j++)
			tmp_vec[j] = cw1_arr[i] + (ncpl_arr[i] == 1 ? 0 : j*(cw2_arr[i] - cw1_arr[i]) / (ncpl_arr[i] - 1));
		search_space.push_back(tmp_vec);
	}

	if (!rank)
		cout << "loadValuesToSearchSpaceVariables() finished" << endl;
}

search_space_point sspemdd_sequential::findLocalMinHillClimbing(vector<double> depths)
{
	cout << "findLocalMinHillClimbing" << endl;

	// choose random point in the search space
	/*for (unsigned i = 0; i < search_space.size(); i++) // i stands for variable_index
	local_record_point_indexes[i] = rand() % search_space[i].size(); // get random index
	cur_point_indexes = local_record_point_indexes;*/
	
	search_space_point local_record_point = getNonRandomStartPoint(depths);
	vector<unsigned> local_record_point_indexes = fromPointToPointIndexes( local_record_point );
	fillDataComputeResidual(local_record_point); // calculated residual is written to cur_point

	bool isCheckRequired = false;
	for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
		if (search_space[i].size() > 1) {
			isCheckRequired = true;
			break;
		}
	}
	if ( (!isCheckRequired) && (verbosity > 0) ) {
		cout << "1 element in a search space, fast exit" << endl;
		return local_record_point;
	}

	checked_points.reserve(N_total);
	checked_points.push_back(local_record_point);
	unsigned skipped_points = 0;
	bool isContinueDimension;
	vector<unsigned> cur_point_indexes;
	search_space_point cur_point;
	// launch iterations of hill climbing
	for (unsigned run_index = 0; run_index < iterated_local_search_runs; run_index++) {
		bool isLocalMin;
		cout << endl;
		cout << "iteration " << run_index << " of ILS" << endl;
		do { // do while local min not reached
			isLocalMin = true; // if changing of every variable will not lead to a record updata, then a local min reached
			for (unsigned i = 0; i < search_space.size(); i++) { // i stands for variable_index
				if (search_space[i].size() == 1) {
					//cout << "one value of a variable, skip it" << endl;
					continue;
				}
				//cout << "variable_index " << variable_index << endl;
				cur_point_indexes = local_record_point_indexes;
				vector<unsigned> point_indexes_before_increase = cur_point_indexes;
				unsigned index_from = cur_point_indexes[i]; // don't check index twice
				if (verbosity > 0)
					cout << "index_from " << index_from << endl;
				bool isDecreaseTurn = false;
				bool isTriggerIncToDec = false;
				do { // change value of a variabble while it leads to updating of a record
					double old_record_residual = local_record_point.residual;
					if (isDecreaseTurn) {
						if (isTriggerIncToDec) { // move to a point before increase
							cur_point_indexes = point_indexes_before_increase;
							isTriggerIncToDec = false;
						}
						if (cur_point_indexes[i] == 0)
							cur_point_indexes[i] = search_space[i].size() - 1;
						else
							cur_point_indexes[i]--;
					}
					else {
						cur_point_indexes[i]++;
						if (cur_point_indexes[i] == search_space[i].size())
							cur_point_indexes[i] = 0;
					}
					if (cur_point_indexes[i] == index_from) {
						if (verbosity > 0)
							cout << "cur_point_indexes[i] == index_from. Break iteration." << endl;
						break;
					}
					if (verbosity > 0) {
						cout << "checking index " << cur_point_indexes[i] <<
							", max index " << search_space[i].size() - 1 << endl;
						cout << "cur_point_indexes" << endl;
						for (unsigned j = 0; j < cur_point_indexes.size(); j++)
							cout << cur_point_indexes[j] << " ";
						cout << endl;
					}
					cur_point = fromPointIndexesToPoint(cur_point_indexes, depths);
					if (std::find(checked_points.begin(), checked_points.end(), cur_point) != checked_points.end()) {
						skipped_points++;
						continue;
					}
					double d_val = fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
					checked_points.push_back(cur_point);
					isContinueDimension = false;
					if (d_val < old_record_residual) { // new record was found
						local_record_point.residual = d_val;
						local_record_point_indexes = cur_point_indexes;
						isLocalMin = false;
						isContinueDimension = true;
					}
					if ((!isContinueDimension) && (!isDecreaseTurn)) {
						isDecreaseTurn = true; // try to decrease current value
						if (verbosity > 0)
							cout << "isDecreaseTurn " << isDecreaseTurn << endl;
						isContinueDimension = true;
						isTriggerIncToDec = true;
						continue;
					}
				} while (isContinueDimension);
			}
		} while (!isLocalMin);

		//cout << endl << "*** local minimum in hill climbing" << endl;
		//cout << "local record of residual " << record_point.residual << endl;
		//cout << "-----" << endl;
		if (verbosity > 0)
			cout << "new random cur_point_indexes : " << endl;
		
		for(;;) {
			// prmutate current global minimum point to obtain a new start point
			for (unsigned i = 0; i < search_space.size(); i++) {
				if (search_space[i].size() == 1) {
					cur_point_indexes[i] = 0;
					if (verbosity > 0)
						cout << cur_point_indexes[i] << " ";
					continue;
				}
				unsigned rand_numb = rand();
				if (rand_numb % 3 == 0)
					cur_point_indexes[i] = local_record_point_indexes[i];
				else
					cur_point_indexes[i] = (rand_numb % search_space[i].size());
			}
			cur_point = fromPointIndexesToPoint(cur_point_indexes, depths);
			if (std::find(checked_points.begin(), checked_points.end(), cur_point) == checked_points.end())
				break;
		}
		if (verbosity > 0) {
			cout << "new random point" << endl;
			for (unsigned j = 0; j < cur_point_indexes.size(); j++)
				cout << cur_point_indexes[j] << " ";
			cout << endl;
		}
		
		fillDataComputeResidual(cur_point); // calculated residual is written to cur_point
		checked_points.push_back(cur_point);
		// new start point
		local_record_point = cur_point;
		local_record_point_indexes = cur_point_indexes;

		cout << "checked_points size " << checked_points.size() << endl;
		cout << "skipped_points " << skipped_points << endl;
		cout << "---" << endl;
	}

	return local_record_point;
}

search_space_point sspemdd_sequential::fromPointIndexesToPoint( vector<unsigned> cur_point_indexes, 
	                                                            vector<double> depths)
{
	search_space_point point;
	point.cb   = search_space[0][cur_point_indexes[0]];
	point.rhob = search_space[1][cur_point_indexes[1]];
	point.R    = search_space[2][cur_point_indexes[2]];
	point.tau  = search_space[3][cur_point_indexes[3]];
	for (unsigned i = 4; i < search_space.size(); i++)
		point.cws.push_back(search_space[i][cur_point_indexes[i]]);
	point.depths = depths;
	point.residual = START_HUGE_VALUE;
	return point;
}

vector<unsigned> sspemdd_sequential::fromPointToPointIndexes( search_space_point point )
{
	vector<unsigned> cur_point_indexes;
	cur_point_indexes.resize(search_space.size());

	for (unsigned j = 0; j<search_space[0].size(); j++)
		if ( search_space[0][j] == point.cb ) {
			cur_point_indexes[0] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[1].size(); j++)
		if (search_space[1][j] == point.rhob) {
			cur_point_indexes[1] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[2].size(); j++)
		if (search_space[2][j] == point.R) {
			cur_point_indexes[2] = j;
			break;
		}
	for (unsigned j = 0; j<search_space[3].size(); j++)
		if (search_space[3][j] == point.tau) {
			cur_point_indexes[3] = j;
			break;
		}
	for (unsigned i = 4; i < search_space.size(); i++) {
		for (unsigned j = 0; j < search_space[i].size(); j++)
			if (search_space[i][j] == point.cws[i - 4]) {
				cur_point_indexes[i] = j;
				break;
			}
	}

	return cur_point_indexes;
}

search_space_point sspemdd_sequential::fromDoubleVecToPoint(vector<double> double_vec)
{
	search_space_point point;
	point.cb   = double_vec[0];
	point.rhob = double_vec[1];
	point.R    = double_vec[2];
	point.tau  = double_vec[3];
	unsigned depths_number = (double_vec.size() - 4) / 2;
	for (unsigned i = 4; i < 4 + depths_number; i++)
		point.cws.push_back(double_vec[i]);
	for (unsigned i = 4 + depths_number; i < double_vec.size(); i++)
		point.depths.push_back(double_vec[i]);
	point.residual = START_HUGE_VALUE;
	return point;
}

double sspemdd_sequential::getRecordResidual()
{
	return record_point.residual;
}

void sspemdd_sequential::getThreeValuesFromStr(std::string str, double &val1, double &val2, double &val3)
{
	val1 = val3 = -1;
	val2 = 1;
	std::string word1, word2, word3;
	for (auto &x : str)
		if (x == ':')
			x = ' ';
	std::stringstream sstream;
	sstream << str;
	sstream >> word1 >> word2 >> word3;
	std::istringstream(word1) >> val1;
	std::istringstream(word2) >> val2;
	std::istringstream(word3) >> val3;
	if (val3 == -1)
		val3 = val1;
}

int sspemdd_sequential::readScenario(std::string scenarioFileName)
{
// read constant and variable values from a scenario file
	if ( (!rank) && (verbosity > 0) )
		cout << "scenarioFileName " << scenarioFileName << endl;
	std::ifstream scenarioFile(scenarioFileName.c_str());

	if (!scenarioFile.is_open()) {
		std::cerr << "scenarioFile with the name " << scenarioFileName << " wasn't openend" << endl;
		return -1;
	}

	std::string str, word, tmp_word;
	std::stringstream sstream;
	unsigned cw_index = 0, d_index = 0;
	double cur_val_step = 0, cur_val1 = 0, cur_val2 = 0;
	while (getline(scenarioFile, str)) {
		if ((str == "") || (str[0] == '%'))
			continue;
		sstream << str;
		sstream >> word;
		if (word.find("dtimes_file") != std::string::npos)
			sstream >> dtimesFileName;
		else if (word.find("spmag_file") != std::string::npos)
			sstream >> spmagFileName;
		else if (word == "h")
			sstream >> h;
		else if (word == "H")
			sstream >> H;
		else if ((word.size() >= 2) && (word[0] == 'c') && (word[1] == 'w')) {
			word = word.substr(2, word.size()-2);
			std::istringstream(word) >> cw_index;
			if (cw1_init_arr.size() < cw_index + 1)
				cw1_init_arr.resize(cw_index + 1);
			if (cw2_init_arr.size() < cw_index + 1)
				cw2_init_arr.resize(cw_index + 1);
			if (ncpl_init_arr.size() < cw_index + 1)
				ncpl_init_arr.resize(cw_index + 1);
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			cw1_init_arr[cw_index] = cur_val1;
			cw2_init_arr[cw_index] = cur_val2;
			if (cur_val1 == cur_val2)
				ncpl_init_arr[cw_index] = 1;
			else
				ncpl_init_arr[cw_index] = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if ((word.size() == 2) && (word[0] == 'd') && (isdigit(word[1]))) {
			word = word.substr(1, word.size() - 1);
			std::istringstream(word) >> d_index;
			d_index--;
			if (d1_arr.size() < d_index + 1) {
				d1_arr.resize(d_index + 1);
				d2_arr.resize(d_index + 1);
				d_step.resize(d_index + 1);
			}
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			d1_arr[d_index] = cur_val1;
			d2_arr[d_index] = cur_val2;
			d_step[d_index] = cur_val_step;
		}
		else if (word == "R") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			R1 = cur_val1;
			R2 = cur_val2;
			if (R1 == R2)
				nR = 1;
			else
				nR = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "rhob") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			rhob1 = cur_val1;
			rhob2 = cur_val2;
			if (rhob1 == rhob2)
				nrhob = 1;
			else
				nrhob = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "cb") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			cb1 = cur_val1;
			cb2 = cur_val2;
			if (cb1 == cb2)
				ncb = 1;
			else
				ncb = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		else if (word == "tau") {
			sstream >> word;
			getThreeValuesFromStr(word, cur_val1, cur_val_step, cur_val2);
			tau1 = cur_val1;
			tau2 = cur_val2;
			if (tau1 == tau2)
				ntau = 1;
			else
				ntau = (unsigned)(ceil((cur_val2 - cur_val1) / cur_val_step)) + 1;
		}
		sstream.str(""); sstream.clear();
	}
	
	if (!cw1_init_arr.size()) {
		std::cerr << "!cw1_init_arr.size()" << endl;
		return -1;
	}
	if (!h || !H) {
		std::cerr << "!h || !H" << endl;
		return -1;
	}

	if ( ((rank == 0) || (rank == 1)) && (verbosity > 0) ) {
		cout << "Parameters :" << endl;
		cout << "cw1_init_arr :" << endl;
		for (auto &x : cw1_init_arr)
			cout << x << " ";
		cout << endl;
		cout << "cw2_init_arr :" << endl;
		for (auto &x : cw2_init_arr)
			cout << x << " ";
		cout << endl;
		cout << "ncpl_init_arr :" << endl;
		for (auto &x : ncpl_init_arr)
			cout << x << " ";
		cout << endl;
		cout << "n_layers_w " << n_layers_w << endl;
		cout << "nR " << nR << endl;
		cout << "R1 " << R1 << endl;
		cout << "R2 " << R2 << endl;
		cout << "ntau " << ntau << endl;
		cout << "tau1 " << tau1 << endl;
		cout << "tau2 " << tau2 << endl;
		cout << "nrhob " << nrhob << endl;
		cout << "rhob1 " << rhob1 << endl;
		cout << "rhob2 " << rhob2 << endl;
		cout << "ncb " << ncb << endl;
		cout << "cb1 " << cb1 << endl;
		cout << "cb2 " << cb2 << endl;
		cout << "dtimes_file " << dtimesFileName << endl;
		cout << "spmag_file " << spmagFileName << endl;

		cout << "readScenario() finished" << endl;
	}

	return 0;
}

int sspemdd_sequential::readInputDataFromFiles()
{	
	std::ifstream dtimesFile(dtimesFileName.c_str());
	if (!dtimesFile.is_open()) {
		std::cerr << "dtimesFile " << dtimesFileName << " wasn't opened" << endl;
		return -1;
	}
	std::stringstream myLineStream;
	std::string myLine;
	double buff;
	vector<double> buffvect;
	// reading the "experimental" delay time data from a file
	while (std::getline(dtimesFile, myLine)) {
		myLine.erase(std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete windows endline symbol for correct reading
		myLineStream << myLine;
		myLineStream >> buff;
		freqs.push_back(buff);

		buffvect.clear();
		while (!myLineStream.eof()) {
			myLineStream >> buff;
			buffvect.push_back(buff);
			mode_numbers.push_back((unsigned)buffvect.size());
		}

		modal_delays.push_back(buffvect);
		myLineStream.str(""); myLineStream.clear();
	}
	dtimesFile.close();

	weight_coeffs.clear();
	std::ifstream spmagFile(spmagFileName.c_str());
	if (spmagFile.is_open()) {
		buffvect.clear();
		while (std::getline(spmagFile, myLine)) {
			myLine.erase(std::remove(myLine.begin(), myLine.end(), '\r'), myLine.end()); // delete windows endline symbol for correct reading
			myLineStream << myLine;
			myLineStream >> buff;

			buffvect.clear();
			while (!myLineStream.eof()) {
				myLineStream >> buff;
				buffvect.push_back(buff);
			}

			weight_coeffs.push_back(buffvect);
			myLineStream.str(""); myLineStream.clear();
		}
		spmagFile.close();

		if ( (!rank) && (verbosity > 0) ) {
			cout << "weight_coeffs.size() " << weight_coeffs.size() << endl;
			cout << "weight_coeffs first 10 lines : " << endl;
			for (unsigned i = 0; i < 10; i++) {
				for (auto &x : weight_coeffs[i])
					cout << x << " ";
				cout << endl;
			}
		}
	}
	else {
		object_function_type = "uniform2";
		if ( (!rank) && (verbosity > 0) )
			cout << "spmagFile " << spmagFileName << " wasn't opened" << endl;
	}

	if ( (!rank) && (verbosity > 0) ){
		cout << "object_function_type changed to " << object_function_type << endl;
		cout << "readInputDataFromFiles() finished " << endl;
	}
	return 0;
}

search_space_point sspemdd_sequential::getNonRandomStartPoint( vector<double> depths )
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3] - tau
	// search_space_variables[4...] - cws

	search_space_point point;

	point.cb = (cb1 == cb2) ? cb1 : (cb2 - cb1) / 2;
	point.rhob = (rhob1 == rhob2) ? rhob1 : (rhob2 - rhob1) / 2;
	point.R = (R1 == R2) ? R1 : (R2 - R1) / 2;
	point.tau = (tau1 == tau2) ? tau1 : (tau2 - tau1) / 2;
	
	point.cws.resize(cw1_arr.size());
	for ( unsigned i = 0; i < point.cws.size(); i++ ) {
		if (cw1_arr[i] == cw2_arr[i])
			point.cws[i] = cw1_arr[i];
		else if ((START_CW_VALUE >= cw1_arr[i]) && (START_CW_VALUE <= cw2_arr[i])) // if cws should be modified
			point.cws[i] = START_CW_VALUE;
		else
			point.cws[i] = (cw2_arr[i] - cw1_arr[i]) / 2;
	}

	point.depths = depths;

	return point;
}

void sspemdd_sequential::printDelayTime(double R, vector<unsigned> mode_numbers, vector<vector<double>> modal_group_velocities)
{
	ofstream ofile("delayTimeOutput.txt");
	for (unsigned ii = 0; ii < freqs.size(); ii++) {
		unsigned mnumb = mode_numbers.at(ii);
		ofile << freqs.at(ii) << "\t";
		for (unsigned jj = 0; jj < mnumb; jj++)
			ofile << R / modal_group_velocities[ii][jj] << "\t";
		ofile << endl;
	}
	ofile.close();
}

void sspemdd_sequential::directPointCalc( search_space_point point )
{
	isTimeDelayPrinting = true;
	fillDataComputeResidual(point);
}