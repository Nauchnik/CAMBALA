#include "normal_modes.h"
#include "linalg.h"

//using namespace CAMBALA_compute;

NormalModes::NormalModes():
	iModesSubset(-1.0)
{}

vector<double> NormalModes::compute_wnumbers_extrap_lin_dz(double &omeg)
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

	unsigned n_layers = (unsigned)M_depths.size();
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
	for (unsigned ii = 0; ii < n_layers; ii++) {
		if (M_Ns_points.at(ii) % 12 > 0) {
			M_Ns_points.at(ii) = 12 * (1 + (M_Ns_points.at(ii) / 12));
		}
	}

	//    cout << "Richardson coeffs" << endl;
	//    for (int ii=0; ii<coeff_extrap.size() ; ii++ ){
	//        cout << coeff_extrap.at(ii) << endl;
	//    }

	// outer loop for Richardson coefficient rr
	for (unsigned rr = 1; rr <= ordRich; rr++) {
		input_c.clear();
		input_rho.clear();
		input_interf_idcs.clear();
		input_mesh.clear();
		out_wnum2.clear();

		input_c.push_back(0);
		input_rho.push_back(0);
		n_points_total = 1;
		zp = 0;

		for (unsigned ll = 0; ll < n_layers; ll++) {
			zc = M_depths.at(ll);
			n_points_layer = M_Ns_points.at(ll) / rr;
			dz = (zc - zp) / (n_points_layer);

			//            cout << "np=" << n_points_layer << "  " << "dz=" << dz <<endl;

			input_mesh.push_back(dz);
			input_c.at(n_points_total - 1) = M_c1s.at(ll);
			input_rho.at(n_points_total - 1) = M_rhos.at(ll);

			n_points_total = n_points_total + n_points_layer;

			for (unsigned kk = 1; kk <= n_points_layer; kk++) {
				input_rho.push_back(M_rhos.at(ll));
				input_c.push_back(M_c1s.at(ll) + (M_c2s.at(ll) - M_c1s.at(ll))*kk / n_points_layer);
			}

			if (ll < n_layers - 1) {
				input_interf_idcs.push_back(n_points_total - 1);
			}
			zp = zc;
		}

		//        cout << "rr=" << rr << endl;

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh);
		m_wnum = min(m_wnum, (unsigned)out_wnum2.size());

		if (rr == 1) { wnum_extrapR.assign(m_wnum, 0); }

		for (unsigned mm = 0; mm < m_wnum; mm++) {
			wnum_extrapR.at(mm) = wnum_extrapR.at(mm) + out_wnum2.at(mm)*coeff_extrap.at(rr - 1);
		}

		//            for (unsigned ii=0; ii<out_wnum2.size();  ii++) {
		//
		//                cout << ii << "->" << sqrt(out_wnum2.at(ii)) << endl;
		//            }
		//            cout << endl;
	}

	for (unsigned mm = 0; mm < m_wnum; mm++) {
		wnum_extrapR.at(mm) = sqrt(wnum_extrapR.at(mm));
	}
	return wnum_extrapR;
}


vector<double> NormalModes::compute_wnumbers(
	double &omeg,
	vector<double> &c,
	vector<double> &rho,
	vector<unsigned> &interface_idcs,
	vector<double> &meshsizes
)                                   
{
	// set iModesSubset < 0->trapped modes only; 0 <= iModesSubset < 1->a subset of modes with wavenumbers
	// within [iModesSubset*kmax kmax], in particular iModesSubset = 0 computes all propagating modes
	//
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
	//        myFile << fixed << setprecision(16) << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
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

	for (int ii = 0; ii < N_points; ii++) {
		if (c.at(ii) < cmin) { cmin = c.at(ii); }
		if (c.at(ii) > cmax) { cmax = c.at(ii); }
	}

	kappamax = omeg / cmin;
	kappamin = omeg / cmax;

	if (iModesSubset >= 0) {

		if (iModesSubset < 1) {
			kappamin = iModesSubset * kappamax;
		}
		else {
			throw invalid_argument("Invalid iModeSubset: use either -1 or a value from [0 1)");
		}

	}




	for (int ii = 0; ii < N_points - 2; ii++) {
		// ordinary point
		ud.push_back(1 / (dz*dz));
		ld.push_back(1 / (dz*dz));
		md.push_back(-2 / (dz*dz) + omeg * omeg / (c.at(ii + 1)*c.at(ii + 1)));

		// special case of the point at the interface

		if (ii == next_interface_idx) {         //ii -- z(ii+1), z(0) = 0
			layer_number = layer_number + 1;    // âîîáùå ii=89 -- âîäà, â ii=90 -äíî,
			// çäåñü ii = 89 -- èíòåðôåéñ, óæå äíî
			cp = c.at(ii + 1);
			dp = rho.at(ii + 1);
			cm = c.at(ii);
			dm = rho.at(ii);


			dz_next = meshsizes.at(layer_number);
			q = 1 / (dz_next*dm + dz * dp);  // Ïîìåíÿòü ìåñòàìè ñ ïðåäûäóùåé ñòðîêîé??

			ld.at(ii) = 2 * q*dp / dz;
			md.at(ii) = -2 * q*(dz_next*dp + dz * dm) / (dz*dz_next) + omeg * omeg*q*(dz*dp*cp*cp + dz_next * dm*cm*cm) / (cp*cp*cm*cm);
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
			ofstream myFile("thematrixdiags.txt");
			for (int ii=0; ii<N_points-2; ii++){
				myFile << fixed  << ld.at(ii) << "  " << md.at(ii) << "  " << ud.at(ii) << endl;
			}
			myFile.close();

			ofstream myFile1("thematrixdiags_sym.txt");
			for (int ii=0; ii<N_points-3; ii++){
				myFile1 << fixed  << main_diag[ii] << "  " << second_diag[ii] << endl;
			}
			myFile1.close();
	*/
	alglib::smatrixtdevdr(main_diag, second_diag, N_points - 2, 0, kappamin*kappamin, kappamax*kappamax, eigen_count, eigenvectors);

	for (int ii = 0; ii < eigen_count; ii++) {
		wnumbers2.push_back(main_diag[eigen_count - ii - 1]);
	}

	return wnumbers2;
}

double NormalModes::RK4(double omeg2, // sound frequency
	double kh2,
	double deltah,
	double c1,
	double c2,
	unsigned Np,
	vector<double> &phi0,
	vector<double> &dphi0)
{

	double f11, f12, f21, f22, f31, f32, f41, f42, cc;
	double h = deltah / Np;
	double layer_int = 0.0;


	for (unsigned kk = 0; kk < Np; kk++) {

		layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

		f11 = dphi0.back();
		cc = c1 + (c2 - c1)*kk*h / deltah;
		cc = cc * cc;
		f12 = (kh2 - omeg2 / cc)*phi0.back();

		f21 = dphi0.back() + 0.5*f12*h;
		cc = c1 + (c2 - c1)*(kk + 0.5)*h / deltah;
		cc = cc * cc;
		f22 = (kh2 - omeg2 / cc)*(phi0.back() + 0.5*f11*h);

		f31 = dphi0.back() + 0.5*f22*h;
		cc = c1 + (c2 - c1)*(kk + 0.5)*h / deltah;
		cc = cc * cc;
		f32 = (kh2 - omeg2 / cc)*(phi0.back() + 0.5*f21*h);

		f41 = dphi0.back() + f32 * h;
		cc = c1 + (c2 - c1)*(kk + 1)*h / deltah;
		cc = cc * cc;
		f42 = (kh2 - omeg2 / cc)*(phi0.back() + f31 * h);

		phi0.push_back(phi0.back() + h * (f11 + 2 * f21 + 2 * f31 + f41) / 6);
		dphi0.push_back(dphi0.back() + h * (f12 + 2 * f22 + 2 * f32 + f42) / 6);

		layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();

	}

	return layer_int;
}

double NormalModes::Layer_an_exp(double omeg2, // sound frequency
	double kh2,
	double deltah,
	double c,
	unsigned Np,
	vector<double> &phi0,
	vector<double> &dphi0)
{

	double c1, c2, kv;
	double h = deltah / Np;
	double layer_int = 0.0;

	kv = sqrt(kh2 - omeg2 / (c*c));
	c1 = 0.5*(phi0.back() - dphi0.back() / kv);
	c2 = 0.5*(phi0.back() + dphi0.back() / kv);
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
	for (unsigned kk = 0; kk < Np; kk++) {

		layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();



		phi0.push_back(c1*exp(-kv * ((kk + 1)*h)) + c2 * exp(kv*((kk + 1)*h)));
		dphi0.push_back(-c1 * kv*exp(-kv * ((kk + 1)*h)) + c2 * kv*exp(kv*((kk + 1)*h)));

		layer_int = layer_int + 0.5*h*phi0.back()*phi0.back();
		//cout << layer_int << endl;
	}

	return layer_int;
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

void NormalModes::load_layers_data(string LayersFName)
{
	ifstream Myfile(LayersFName);
	double d, c1, c2, rho;
	unsigned nsp;

	M_depths.clear();
	M_c1s.clear();
	M_c2s.clear();
	M_rhos.clear();
	M_Ns_points.clear();

	while (!(Myfile.eof()))
	{
		Myfile >> d;
		Myfile >> c1;
		Myfile >> c2;
		Myfile >> rho;
		Myfile >> nsp;

		M_depths.push_back(d);
		M_c1s.push_back(c1);
		M_c2s.push_back(c2);
		M_rhos.push_back(rho);
		M_Ns_points.push_back(nsp);

	}
	Myfile.close();
}

void NormalModes::load_profile_deep_water(string ProfileFName, const unsigned ppm)
{
	ifstream Myfile(ProfileFName);
	double cp, cc, dc, dp;

	Myfile >> dp;
	Myfile >> cp;

	unsigned npc;

	while (!(Myfile.eof()))
	{
		Myfile >> dc;
		Myfile >> cc;

		M_depths.push_back(dc);
		M_c1s.push_back(cp);
		M_c2s.push_back(cc);
		M_rhos.push_back(1);

		npc = (unsigned)abs(ppm*(dc - dp));
		M_Ns_points.push_back(npc);

		cp = cc;
		dp = dc;

	}
	Myfile.close();
}

double NormalModes::compute_modal_delays_residual_uniform(
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

	compute_modal_grop_velocities(deltaf, modal_group_velocities, mode_numbers);

	for (unsigned ii = 0; ii < freqs.size(); ii++) {
		//2016.04.27:Pavel: mnumb = min(mode_numbers.at(ii), experimental_mode_numbers.at(ii));
		mnumb = experimental_mode_numbers.at(ii);
		for (unsigned jj = 0; jj < mnumb; jj++) {
			if (experimental_delays[ii][jj] > 0) {
				nRes = nRes + 1;
				//2016.04.27:Pavel:
				if (jj < mode_numbers.at(ii)) {
					mdelay = R / modal_group_velocities[ii][jj];
				}
				else if ((ii + 1 < freqs.size()) && (jj < mode_numbers.at(ii + 1))) {
					mdelay = R / modal_group_velocities[ii + 1][jj];
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
	residual = sqrt(residual / nRes);

	//if (isTimeDelayPrinting)
	//	printDelayTime(R, freqs, mode_numbers, modal_group_velocities);

	return residual;
}

int NormalModes::compute_modal_grop_velocities(
	double deltaf,
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

	for (unsigned ii = 0; ii < nfr; ii++) {
		out_wnum1.clear();
		out_wnum2.clear();
		mgv_ii.clear();
		omeg1 = 2 * M_PI*(freqs.at(ii) + deltaf);
		out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1);
		nwnum = (unsigned)out_wnum1.size();

		/*
		cout << "f=" << freqs.at(ii) << "Hz" << endl;

		for (unsigned jj=0; jj < nwnum; jj++ )
		{
		cout << "k_" << jj+1 << "=" << out_wnum1.at(jj) << endl;
		}
		*/

		omeg2 = 2 * M_PI*(freqs.at(ii));
		out_wnum2 = compute_wnumbers_extrap_lin_dz(omeg2);
		//nwnum = min(nwnum, (unsigned)out_wnum2.size());
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