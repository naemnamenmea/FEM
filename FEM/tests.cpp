#include "stdafx.hpp"

#include "test_runner.h"
#include "tests.hpp"

// Numerical integration tests
void TestGaussIntegrBasic();
void TestGaussIntegr1d();
void TestGaussIntegr2d();
void TestGaussIntegr3d();

// Input tests
void TestFEMChekhovInputNodes();
void TestFEMChekhovInputMaterials();
void TestFEMChekhovInputFiniteElements();

namespace TEST
{
void RunAllTests()
{
	TestRunner tr;

	// Input tests
	RUN_TEST(tr, TestFEMChekhovInputNodes);
	RUN_TEST(tr, TestFEMChekhovInputMaterials);
	RUN_TEST(tr, TestFEMChekhovInputFiniteElements);

	// Numerical integration tests
	//RUN_TEST(tr, TestGaussIntegrBasic);
	//RUN_TEST(tr, TestGaussIntegr1d);
	//RUN_TEST(tr, TestGaussIntegr2d);
	//RUN_TEST(tr, TestGaussIntegr3d);
}
}  // namespace TEST