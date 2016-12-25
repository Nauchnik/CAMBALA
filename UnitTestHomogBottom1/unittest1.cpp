#include "stdafx.h"
#include "CppUnitTest.h"
#include "sspemdd_sequential.h"
#include "sspemdd_utils.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;

namespace UnitTestHomogBottom1
{		
	TEST_CLASS(UnitTest1)
	{
	public:
		
		TEST_METHOD(TestMethod1)
		{
			/*
			RESULTING VALUE:
			err = 0.00267443, parameters:
			c_b = 1830
			tau = -0.005
			rho_b = 1.45
			R = 8000
			cws : 1500 
			*/
			sspemdd_sequential sspemdd_seq;
			sspemdd_seq.readDataFromFile("../8000_extracted.txt", 0);
			sspemdd_seq.init();
			search_space_point cur_point;
			cur_point.cb = 1830;
			cur_point.rhob = 1.45;
			cur_point.tau = -0.005;
			cur_point.R = 8000;
			cur_point.cws.resize(1);
			cur_point.cws[0] = 1500;
			double d = sspemdd_seq.fill_data_compute_residual(cur_point);
			Assert::AreEqual(sspemdd_seq.fill_data_compute_residual(cur_point), 0.0026744323772541706);
		}
	};
}