
void EvalPointCPU( Point &point,
		const std::vector<double> &freqs_d,
		const std::vector<unsigned> &Ns_points_d,
		const std::vector<double> &depths_d,
		const std::vector<std::vector<double>> &modal_delays);

double CAMBALA_sequential::computeResidual(
		Point &point,
		const std::vector<double> &freqs_d,
		const std::vector<unsigned> &Ns_points_d,
		const std::vector<double> &depths_d,
		const std::vector<std::vector<double>> &modal_delays)
{ 
	vector<double> c1s;
	vector<double> c2s;
	vector<double> rhos;
	vector<unsigned> Ns_points;
	
	// finally specify sound speed in water
	// the parameters are transformed into the arrays c1s, c2s, rhos
	if (verbosity > 1)
		cout << "fillDataComputeResidual()" << endl;

	if (point.cws.size() != depths.size() - 1)
	{
		cerr << "point.cws.size() != depths.size() - 1" << endl;
		cerr << point.cws.size() << " " << depths.size() - 1 << endl;
		exit(1);
	}
	if (depths.size() == 0)
	{
		cerr << "depths.size() == 0" << endl;
		exit(-1);
	}

	if (verbosity > 1)
	{
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

		point.residual = compute_modal_delays_residual_weighted2(freqs, depths, c1s, c2s, rhos, Ns_points, point.R, point.tau, modal_delays, weight_coeffs, mode_numbers);


	if ( verbosity > 1 )
		cout << "point.residual " << point.residual << endl;


	return point.residual;
}
