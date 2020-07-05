#include "FiniteElementPool.hpp"
#include "FiniteElementModel.hpp"
#include "profile.hpp"
#include "tchar.h"
#include <fstream>

const std::string& GetTestDataDir();

static const std::string GetInData()
{
	return GetTestDataDir() + _T("/Folder/");
}

#define PATH GetInData()

namespace
{
void RunTest(const std::string& filename)
{
	std::ostringstream message;
	message << "filename: \"" << filename << "\"";
	tFE_model fem;
	fem.ImportFromFile((PATH + filename).c_str());
	double femVolume = 0.0;

	std::ostringstream os;
	{
		LOG_DURATION_STREAM(os, "calc_time");
		femVolume = fem.Volume();
	}
	message << " | fem (" << fem.Name() << ") Volume: " << femVolume << " | " << os.str();
	std::cerr << message.str();
}
}  // namespace

void TestFEVolumeCalc()
{
	tFE_model fem;
	fem.ImportFromFile("../data/3dcube1000el.xml");
	double feVolume = 0.0;
	double femVolume = 0.0;
	{
		LOG_DURATION("Single FE volume");
		feVolume = fem.GetFE(1).Volume();
	}
	std::cout << "fe #1 Volume: " << feVolume << std::endl;

	{
		LOG_DURATION("Total Volume");
		femVolume = fem.Volume();
	}
	std::cout << "fem Volume: " << femVolume << std::endl;
}

void TestFEVolumeCalc3dCube()
{
	RunTest("3dcube1000el.xml");
	RunTest("3dcube2000el.xml");
	RunTest("3dcube4000el.xml");
}

void TestFEVolumeCalc2dRectangle()
{
	RunTest("2drectangle1000el.xml");
	RunTest("2drectangle2000el.xml");
	RunTest("2drectangle4000el.xml");
}

void TestFEVolumeCalc1dRod()
{
	RunTest("1drod1000el.xml");
	RunTest("1drod2000el.xml");
	RunTest("1drod4000el.xml");
}
