#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"
#include <iostream>
#include <time.h>

sspemdd_sequential::sspemdd_sequential() :
	ncb(0),
	nrhob(0),
	nR(0),
	cb1(0.0),
	cb2(0.0),
	cw1(0.0),
	cw2(0.0),
	R1(0.0),
	R2(0.0),
	rhob1(0.0),
	rhob2(0.0),
	n_layers_w(0),
	launchType(0),
	iterated_local_search_runs(10)
{
	record_point.cb = 1e50;
	record_point.rhob = 1e50;
	record_point.R = 1e50;
	record_point.residual = 1e100;
	srand(time(NULL));
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
double sspemdd_sequential::compute_modal_delays_residual_uniform(std::vector<double> &freqs,
	std::vector<double> &depths,
	std::vector<double> &c1s,
	std::vector<double> &c2s,
	std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points,
	double R,
	std::vector<std::vector<double>> &experimental_delays,
	std::vector<unsigned> &experimental_mode_numbers
	)
{
	unsigned rord = 3;
	unsigned flTrappedOnly = 1;
	double deltaf = 0.5;
	double residual = 0;
	unsigned mnumb;
	double mdelay;

	std::vector<std::vector<double>> modal_group_velocities;
	std::vector<unsigned> mode_numbers;

	compute_modal_grop_velocities(freqs, deltaf, depths, c1s, c2s, rhos, Ns_points, flTrappedOnly, rord, modal_group_velocities, mode_numbers);

	for (unsigned ii = 0; ii<freqs.size(); ii++) {
		mnumb = std::min(mode_numbers.at(ii), experimental_mode_numbers.at(ii));
		for (unsigned jj = 0; jj<mnumb; jj++) {
			if (experimental_delays[ii][jj]>0) {
				mdelay = R / modal_group_velocities[ii][jj];
				residual = residual + pow(experimental_delays[ii][jj] - mdelay, 2);
			}
		}
	}

	residual = sqrt(residual);

	return residual;
}

int sspemdd_sequential::compute_modal_grop_velocities(std::vector<double> &freqs,
	double deltaf,
	std::vector<double> &depths,
	std::vector<double> &c1s,
	std::vector<double> &c2s,
	std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points,
	unsigned flOnlyTrapped,
	unsigned &ordRich,
	std::vector<std::vector<double>> &modal_group_velocities,
	std::vector<unsigned> &mode_numbers
	)
{
	mode_numbers.clear();
	modal_group_velocities.clear();

	std::vector<double> out_wnum1;
	std::vector<double> out_wnum2;
	std::vector<double> mgv_ii;
	unsigned nwnum;
	unsigned nfr = (unsigned)freqs.size();
	double omeg1, omeg2;

	for (unsigned ii = 0; ii<nfr; ii++) {
		out_wnum1.clear();
		out_wnum2.clear();
		mgv_ii.clear();
		omeg1 = 2 * LOCAL_M_PI*(freqs.at(ii) + deltaf / 2);
		out_wnum1 = compute_wnumbers_extrap_lin_dz(omeg1, depths, c1s, c2s, rhos, Ns_points, 1, ordRich);
		nwnum = (unsigned)out_wnum1.size();

		/*
		cout << "f=" << freqs.at(ii) << "Hz" << endl;

		for (unsigned jj=0; jj < nwnum; jj++ )
		{
		cout << "k_" << jj+1 << "=" << out_wnum1.at(jj) << endl;
		}
		*/

		omeg2 = 2 * LOCAL_M_PI*(freqs.at(ii) - deltaf / 2);
		out_wnum2 = compute_wnumbers_extrap_lin_dz(omeg2, depths, c1s, c2s, rhos, Ns_points, 1, ordRich);
		nwnum = std::min(nwnum, (unsigned)out_wnum2.size());

		for (unsigned jj = 0; jj < nwnum; jj++)
		{
			mgv_ii.push_back((omeg1 - omeg2) / (out_wnum1.at(jj) - out_wnum2.at(jj)));
		}

		modal_group_velocities.push_back(mgv_ii);
		mode_numbers.push_back(nwnum);
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
std::vector<double> sspemdd_sequential::compute_wnumbers_extrap(double &omeg, // sound frequency
	std::vector<double> &depths,
	std::vector<double> &c1s,
	std::vector<double> &c2s,
	std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points,
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
	std::vector<double> coeff_extrap;
	switch (ordRich) {
	case 1:
		coeff_extrap.assign({ 1 });
		break;
	case 2:
		coeff_extrap.assign({ -1, 2 });
		break;
	case 4:
		coeff_extrap.assign({ -1 / double(6), 4, -13.5, 32 / double(3) });
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

	std::vector<double> input_c;
	std::vector<double> input_rho;
	std::vector<double> input_mesh;
	std::vector<unsigned> input_interf_idcs;
	std::vector<double> out_wnum2;
	std::vector<double> wnum_extrapR;
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

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh, flOnlyTrapped);
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

std::vector<double> sspemdd_sequential::compute_wnumbers_extrap_lin_dz(double &omeg, // sound frequency
	std::vector<double> &depths,
	std::vector<double> &c1s,
	std::vector<double> &c2s,
	std::vector<double> &rhos,
	std::vector<unsigned> &Ns_points,
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
	std::vector<double> input_c;
	std::vector<double> input_rho;
	std::vector<double> input_mesh;
	std::vector<unsigned> input_interf_idcs;
	std::vector<double> out_wnum2;
	std::vector<double> wnum_extrapR;
	double zc = 0;
	double zp = 0;
	double dz = 0;
	unsigned m_wnum = 1000;

	unsigned n_layers = (unsigned)depths.size();
	unsigned n_points_total = 0;
	unsigned n_points_layer = 0;

	std::vector<double> coeff_extrap;
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
		Ns_points.at(ii) = 12 * (Ns_points.at(ii) / 12);
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

		out_wnum2 = compute_wnumbers(omeg, input_c, input_rho, input_interf_idcs, input_mesh, flOnlyTrapped);
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

std::vector<double> sspemdd_sequential::compute_wnumbers(double &omeg, // sound frequency
	std::vector<double> &c,
	std::vector<double> &rho,
	std::vector<unsigned> &interface_idcs,
	std::vector<double> &meshsizes,
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

	std::vector<double> md;
	std::vector<double> ud;
	std::vector<double> ld;
	std::vector<double> wnumbers2;

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

	if (flOnlyTrapped == 0)
		kappamin = 0;

	for (int ii = 0; ii < N_points - 2; ii++){
		// ordinary point
		ud.push_back(1 / (dz*dz));
		ld.push_back(1 / (dz*dz));
		md.push_back(-2 / (dz*dz) + omeg*omeg / (c.at(ii + 1)*c.at(ii + 1)));

		// special case of the point at the interface

		if (ii == next_interface_idx) {         //ii -- z(ii+1), z(0) = 0
			layer_number = layer_number + 1;    // вообще ii=89 -- вода, в ii=90 -дно,
			// здесь ii = 89 -- интерфейс, уже дно
			cp = c.at(ii + 1);
			dp = rho.at(ii + 1);
			cm = c.at(ii);
			dm = rho.at(ii);
			q = 1 / (dz_next*dm + dz*dp);

			dz_next = meshsizes.at(layer_number);

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

	alglib::real_2d_array eigenvectors; // V - собств вектор
	alglib::real_1d_array eigenvalues; // Lm -собств знач
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

	alglib::smatrixtdevdr(main_diag, second_diag, N_points - 2, 0, kappamin*kappamin, kappamax*kappamax, eigen_count, eigenvectors);

	for (int ii = 0; ii<eigen_count; ii++) {
		wnumbers2.push_back(main_diag[eigen_count - ii - 1]);
	}

	return wnumbers2;
}

void sspemdd_sequential::init()
{
	// TODO move to constructor
	std::vector<double> tmp_vec{ 1490, 1490, 1480, 1465, 1460 }; // use when ncpl = 1
	cws_fixed = tmp_vec;
	record_point.cws.resize(n_layers_w);
	for (unsigned i = 0; i < record_point.cws.size(); i++)
		record_point.cws[i] = cw2;
}

void sspemdd_sequential::findGlobalMinBruteForce()
{
	// make cws_all_cartesians - all cartesians of all speeds in water
	std::vector<std::vector<double>> cws_vii; // all variants for every depth
	std::vector<int> index_arr;
	std::vector<double> cws_vi;
	std::vector<std::vector<double>> cws_all_cartesians;
	bool isAdding;
	if (ncpl == 1)
		cws_all_cartesians.push_back(cws_fixed);
	else {
		for (unsigned ncpl_cur = 0; ncpl_cur < ncpl; ncpl_cur++)
			cws_vi.push_back(cw1 + ncpl_cur*(cw2 - cw1) / (ncpl - 1));
		for (unsigned i = 0; i < n_layers_w; i++)
			cws_vii.push_back(cws_vi);
		while (SSPEMDD_utils::next_cartesian(cws_vii, index_arr, cws_vi)) {
			switch (launchType) {
			case 1:
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 2:
				isAdding = true;
				for (unsigned cws_vi_index = 0; cws_vi_index < cws_vi.size() - 1; cws_vi_index++)
					if (cws_vi[cws_vi_index] <= cws_vi[cws_vi_index + 1])
						isAdding = false;
				if (isAdding)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 3:
				if (cws_vi[0] != 1490)
					isAdding = false;
				else {
					isAdding = true;
					for (unsigned cws_vi_index = 0; cws_vi_index < cws_vi.size() - 1; cws_vi_index++)
						if (cws_vi[cws_vi_index] <= cws_vi[cws_vi_index + 1])
							isAdding = false;
				}
				if (isAdding)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 7: // like case 1
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			case 8: // like case 1
				if (cws_vi[0] == 1490)
					cws_all_cartesians.push_back(cws_vi);
				break;
			default:
				cws_all_cartesians.push_back(cws_vi);
				break;
			}
		}
	}

	std::cout << "cws_all_cartesians.size() " << cws_all_cartesians.size() << std::endl;
	double cb_cur, rhob_cur, R_cur;
	search_space_point cur_point;
	
	// inverting for bottom halfspace parameters + sound speed in water
	for (unsigned cur_ncb = 0; cur_ncb < ncb; cur_ncb++) {
		for (unsigned cur_nrhob = 0; cur_nrhob < nrhob; cur_nrhob++) {
			for (unsigned cur_nR = 0; cur_nR < nR; cur_nR++) {
				// specify bottom parameters;
				if (ncb > 1) { cb_cur = cb1 + cur_ncb*(cb2 - cb1) / (ncb - 1); }
				else { cb_cur = cb1; }
				if (nrhob > 1) { rhob_cur = rhob1 + cur_nrhob * (rhob2 - rhob1) / (nrhob - 1); }
				else { rhob_cur = rhob1; }
				// specify range
				if (nR > 1) { R_cur = R1 + cur_nR*(R2 - R1) / (nR - 1); }
				else { R_cur = R1; }
				cur_point.cb = cb_cur;
				cur_point.rhob = rhob_cur;
				cur_point.R = R_cur;
				for (auto &cws_cur : cws_all_cartesians) {
					cur_point.cws = cws_cur;
					fill_data_compute_residual(cur_point);
				}
			}
		}
	}
}

void sspemdd_sequential::fill_data_compute_residual( search_space_point &point)
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
	
	//for (unsigned jj = 0; jj <= n_layers_w; jj++)
	//	std::cout << "Layer #" << jj + 1 << ": c=" << c1s.at(jj) << "..." << c2s.at(jj) << "; rho=" << rhos.at(jj) << "; np=" << Ns_points.at(jj) << std::endl;
	point.residual = compute_modal_delays_residual_uniform(freqs, depths, c1s, c2s, rhos, Ns_points, point.R, modal_delays, mode_numbers);
	//std::cout << residual << std::endl << std::endl;

	if (point.residual < record_point.residual) {
		record_point.residual = point.residual;
		record_point.cb       = point.cb;
		record_point.rhob     = point.rhob;
		record_point.R        = point.R;
		record_point.cws      = point.cws;
		std::cout << std::endl;
		std::cout << std::endl << "New residual minimum:" << std::endl;
		std::cout << "err = " << record_point.residual << ", parameters:" << std::endl;
		std::cout << "c_b = " << record_point.cb <<
			", rho_b = " << record_point.rhob <<
			", R = " << record_point.R << std::endl;
		std::cout << "cws_min :" << std::endl;
		for (auto &x : record_point.cws)
			std::cout << x << " ";
		std::cout << std::endl;
	}
}

void sspemdd_sequential::loadValuesToSearchSpaceVariables()
{
	// search_space_variables[0] - cb
	// search_space_variables[1] - rhob
	// search_space_variables[2] - R
	// search_space_variables[3...] - cws
	std::vector<double> tmp_vec;
	
	// fill search_space_variables[0] with cb
	tmp_vec.resize(ncb);
	for (unsigned i = 0; i < ncb; i++)
		tmp_vec[i] = cb1 + (ncb == 1 ? 0 : i*(cb2 - cb1) / (ncb - 1) );
	search_space.push_back(tmp_vec);

	// fill search_space_variables[1] with rhob
	tmp_vec.resize(nrhob);
	for (unsigned i = 0; i < nrhob; i++)
		tmp_vec[i] = rhob1 + (nrhob == 1 ? 0 : i*(rhob2 - rhob1) / (nrhob - 1));
	search_space.push_back(tmp_vec);
	
	// fill search_space_variables[2] with cb
	tmp_vec.resize(nR);
	for (unsigned i = 0; i < nR; i++)
		tmp_vec[i] = R1 + (nR == 1 ? 0 : i*(R2 - R1) / (nR - 1));
	search_space.push_back(tmp_vec);
	
	// fill search_space_variables[3-...] with cws
	tmp_vec.resize(ncpl);
	for (unsigned i = 0; i < ncpl; i++)
		tmp_vec[i] = cw1 + (ncpl == 1 ? 0 : i*(cw2 - cw1) / (ncpl - 1));
	for (unsigned i = 0; i < n_layers_w; i++)
		search_space.push_back(tmp_vec);
}

void sspemdd_sequential::findLocalMinHillClimbing()
{
	loadValuesToSearchSpaceVariables();
	// choose random point in the search space
	std::vector<unsigned> cur_point_indexes, record_point_indexes;
	search_space_point cur_point, global_record_point;
	cur_point_indexes.resize(search_space.size());
	for (unsigned variable_index = 0; variable_index < search_space.size(); variable_index++)
		cur_point_indexes[variable_index] = rand() % search_space[variable_index].size(); // get random index
	record_point_indexes = cur_point_indexes;
	
	// calculate residual in the start point
	cur_point = fromPointIndexesToPoint( cur_point_indexes );
	fill_data_compute_residual( cur_point ); // calculated residual is written to cur_point
	global_record_point = record_point;
	
	// launch hill climbing for variables
	bool isLocalMin;
	bool isRecordUpdateInDimension;
	unsigned checked_points_number = 0;
	unsigned index_from;
	double old_record_residual;
	
	for (unsigned run_index = 0; run_index < iterated_local_search_runs; run_index++) {
		std::cout << "run " << run_index << " of ILS" << std::endl;
		do { // do while local min not reached
			isLocalMin = true; // if changing of every variable will not lead to a record updata, then a local min reached
			for (unsigned variable_index = 0; variable_index < search_space.size(); variable_index++) {
				if (search_space[variable_index].size() == 1) {
					//std::cout << "one value of a variable, skip it" << std::endl;
					continue;
				}
				std::cout << "variable_index " << variable_index << std::endl;
				cur_point_indexes = record_point_indexes;
				index_from = cur_point_indexes[variable_index]; // don't check index twice
				std::cout << "index_from " << index_from << std::endl;
				do { // change value of a variabble while it leads to updating of a record
					old_record_residual = record_point.residual;
					cur_point_indexes[variable_index]++;
					if (cur_point_indexes[variable_index] == search_space[variable_index].size())
						cur_point_indexes[variable_index] = 0;
					if (cur_point_indexes[variable_index] == index_from) {
						std::cout << "cur_point_indexes[variable_index] == index_from. Break iteration." << std::endl;
						break;
					}
					std::cout << "checking index " << cur_point_indexes[variable_index] <<
						", max index " << search_space[variable_index].size() - 1 << std::endl;
					cur_point = fromPointIndexesToPoint(cur_point_indexes);
					fill_data_compute_residual(cur_point); // calculated residual is written to cur_point
					checked_points_number++;
					std::cout << "checked_points_number " << checked_points_number << std::endl;
					std::cout << "-----" << std::endl;
					if (record_point.residual < old_record_residual) { // new record was found
						record_point_indexes = cur_point_indexes;
						isLocalMin = false;
						isRecordUpdateInDimension = true;
					}
					else
						isRecordUpdateInDimension = false;
				} while (isRecordUpdateInDimension);
			}
		} while (!isLocalMin);

		if ( record_point.residual < global_record_point.residual ) {
			global_record_point = record_point;
			std::cout << "on run_index " << run_index << " new global minimum with residual "
				      << global_record_point.residual << std::endl;
		}
		
		std::cout << std::endl << "*** local minimum in hill climbing" << std::endl;
		std::cout << "local record of residual " << record_point.residual << std::endl;
		std::cout << "-----" << std::endl;
		std::cout << "new random cur_point_indexes : " << std::endl;
		for (unsigned variable_index = 0; variable_index < search_space.size(); variable_index++) {
			cur_point_indexes[variable_index] = (rand() % search_space[variable_index].size());
			std::cout << cur_point_indexes[variable_index] << " ";
		}
		
		cur_point = fromPointIndexesToPoint(cur_point_indexes);
		fill_data_compute_residual(cur_point); // calculated residual is written to cur_point
		record_point = cur_point;
		record_point_indexes = cur_point_indexes;
	}
	
	// during optimization record_point is a local minimum, finaly it's the global minimum
	record_point = global_record_point;
}

search_space_point sspemdd_sequential::fromPointIndexesToPoint(std::vector<unsigned> cur_point_indexes)
{
	search_space_point point;
	point.cb   = search_space[0][cur_point_indexes[0]];
	point.rhob = search_space[1][cur_point_indexes[1]];
	point.R    = search_space[2][cur_point_indexes[2]];
	for (unsigned i = 3; i < search_space.size(); i++)
		point.cws.push_back(search_space[i][cur_point_indexes[i]]);
	point.residual = 1e100;
	return point;
}