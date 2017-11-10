#include "solvers/bbox.h"
#define ELPP_STL_LOGGING
#include "easylogging++.h"

#include <common/sgerrcheck.hpp>
#include <common/utilmacro.hpp>
#include <pointgen/randpointgen.hpp>
#include "mcplusvcd.hpp"
void BlackBox::DoLoadSearchSpaceDims(SearchSpaceDims ssd)
{
	ssd_ = ssd;
}

// Hujak-hujak -> production! 
Point PointReduced2Point (std::vector <double> mask, const double* x)
{
	//TODO: safety checks!
	//TODO: rewrite me proper!
	std::vector <double> ptmp(mask.size());
	for (size_t i=0,j=0; i<mask.size(); ++i)
	{
		ptmp[i] = mask[i];
		if (mask[i] == INFINITY)
		{
			ptmp[i] = x[j++];
		}
	}
	Point p;
	p.R = ptmp[0];
	p.rhob = ptmp[1];
	p.cb = ptmp[2];
	p.tau = ptmp[3];
	for (int i = 4; i<ptmp.size(); ++i)
	{
		p.cws.push_back(ptmp[i]);
	}
	return p;
}

// even more hujak-hujak!!
void Point2Reduced (std::vector <double> mask, Point p, double* x)
{
	std::vector <double> ptmp(mask.size());
	ptmp[0] = p.R;
	ptmp[1] = p.rhob;
	ptmp[2] = p.cb;
	ptmp[3] = p.tau;
	for (int i = 4; i<ptmp.size(); ++i)
		ptmp[i] = p.cws[i];

	for (size_t i=0, j=0; i<mask.size(); ++i)
	{
		if (mask[i]==INFINITY)
			x[j++]=x[i];
	}
}
void BlackBox::DoSolve()
{
	//TODO: safety checks!
	// Copy dims to temporary vector
	std::vector <Dim> ssdv;
	ssdv.push_back(ssd_.R);
	ssdv.push_back(ssd_.rhob);
	ssdv.push_back(ssd_.cb);
	ssdv.push_back(ssd_.tau);
	for (auto& d: ssd_.cw)
		ssdv.push_back(d);

	// Mask of inactive dimensions
	std::vector <double> mask(ssdv.size());
	std::vector <Dim> ssdv_reduced;
	for (size_t i = 0; i < ssdv.size(); ++i)
	{
		Dim& d = ssdv[i];
		mask[i] = INFINITY;
		if (d.l == d.r)
		{
			mask[i] = d.l;
		}
		else
		{
			ssdv_reduced.push_back(d);
		}
	}

	int mN = ssdv_reduced.size();
	snowgoose::Box<double> mbox(mN);
	for (size_t i = 0; i < ssdv_reduced.size(); ++i)
	{
		mbox.mA[i] = ssdv_reduced[i].l;
		mbox.mB[i] = ssdv_reduced[i].r;
	}

	COMPI::MPProblem<double> prob = COMPI::MPProblem<double>();

	prob.mVarTypes.assign(mN, COMPI::MPProblem<double>::VariableTypes::GENERIC);
	prob.mObjectives.push_back(std::make_shared<AcousticsObjective>(rc_, mask));

	prob.mBox = &mbox;

	// Setup solver
	const int n = prob.mVarTypes.size();
	const int numberOfPoints = 4;
	const double minimalGranularity = 1e-4;

	MCplusVCD mcsearch(prob, minimalGranularity, numberOfPoints);

	// Setup initial point
	double x[n];

	if (startingPoint_ == Point())
	{
		snowgoose::BoxUtils::getCenter(*(prob.mBox), (double*) x);
	}
	else
	{
		Point2Reduced(mask, startingPoint_, x);
	}
	

	double v = prob.mObjectives.at(0)->func(x);

	// Run solver
	LOG(DEBUG) << "Searching with " << mcsearch.about() << "\n";
	mcsearch.search(x, v);

	// Print results
	LOG(DEBUG) << "Found x = " << snowgoose::VecUtils::vecPrint(n, x) << "\n";
	LOG(DEBUG) << "Objective value = " << v << "\n";

	Point p;
	p.R    = x[0];
	p.rhob = x[1];
	p.cb   = x[2];
	p.tau  = x[3];
	for (size_t i = 0; i < mN; ++i)
		p.cws.push_back(x[i+4]);
	p.residual = v;

	UpdateRecord(p);
}

BlackBox::BlackBox (std::string name)
{
	name_ = name;
}

BlackBox::BlackBox ()
{
	name_ = "BlackBox";
}

BlackBox::~BlackBox(){}


double AcousticsObjective::func(const double* x)
{
	Point p = PointReduced2Point(ssdv_mask, x);
	return rc_->ComputeResidual(p);
}

