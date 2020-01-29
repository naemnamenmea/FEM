#include "stdafx.hpp"

#include "test_runner.hpp"
#include "tests.hpp"

void RunAllTests()
{
	TestRunner tr;
	// Numerical integration tests
	 RUN_TEST(tr, TestBasic1d);
	 RUN_TEST(tr, TestBasic2d);
	 RUN_TEST(tr, TestBasic3d);

	// Input tests
	RUN_TEST(tr, TestFEMChekhovInputNodes);
	RUN_TEST(tr, TestFEMChekhovInputMaterials);
	RUN_TEST(tr, TestFEMChekhovInputFiniteElements);
}
