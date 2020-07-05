#include "test_runner.hpp"

const std::string TEST_DATA_DIR = "../data/";

const std::string& GetTestDataDir()
{
	return TEST_DATA_DIR;
}

// Numerical integration tests
void TestGaussIntegrBasic();
void TestGaussIntegr1dOrder();
void TestGaussIntegr2dOrder();
void TestGaussIntegr3dOrder();
void TestGaussIntegr1d();
void TestGaussIntegr2d();
void TestGaussIntegr3d();

// Input tests
// void TestFEMChekhovInputNodes();
// void TestFEMChekhovInputMaterials();
// void TestFEMChekhovInputFiniteElements();

void Test3dFiniteElement();

void TestFEVolumeCalc();
void TestFEVolumeCalc3dCube();
void TestFEVolumeCalc2dRectangle();
void TestFEVolumeCalc1dRod();

void RunAllTests()
{
	TestRunner tr;

	// Input tests
	// RUN_TEST(tr, TestFEMChekhovInputNodes);
	// RUN_TEST(tr, TestFEMChekhovInputMaterials);
	// RUN_TEST(tr, TestFEMChekhovInputFiniteElements);

	// Numerical integration tests
	RUN_TEST(tr, TestGaussIntegrBasic);
	RUN_TEST(tr, TestGaussIntegr1dOrder);
	RUN_TEST(tr, TestGaussIntegr2dOrder);
	RUN_TEST(tr, TestGaussIntegr3dOrder);
	RUN_TEST(tr, TestGaussIntegr1d);
	RUN_TEST(tr, TestGaussIntegr2d);
	RUN_TEST(tr, TestGaussIntegr3d);
	// return;
	RUN_TEST(tr, Test3dFiniteElement);

	RUN_TEST(tr, TestFEVolumeCalc);
	RUN_TEST(tr, TestFEVolumeCalc1dRod);
	RUN_TEST(tr, TestFEVolumeCalc2dRectangle);
	RUN_TEST(tr, TestFEVolumeCalc3dCube);
}

int main()
{
	RunAllTests();

	return 0;
}