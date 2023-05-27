#include "stdafx.hpp"

#include <fstream>
#include <iomanip>
#include <functional>
#include <iostream>
#include <memory>
#include <type_traits>
#include <sstream>
#include <string>
#include <vector>
#include "app_constants.hpp"
#include "auxiliary.hpp"
#include "fem.hpp"
#include "numeric_math.hpp"
#include "point3d.hpp"
#include "cmdLineParams.hpp"
#include "test_runner.h"
#include "tests.hpp"

int main(int argc, const char* argv[])
{
	TEST::RunAllTests();

	try
	{
		FiniteElementModel<double, 3> FEModel;
		CmdLineParams cmdParams(argc, argv);
		FEModel.ReadData(cmdParams.GetInputFilePath());
		FEModel.WriteData(cmdParams.GetOutputFilePath());

		std::cout << "FEM executed successfully." << std::endl;
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown error occured." << std::endl;
	}

	return 0;
}
