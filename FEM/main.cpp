#include "stdafx.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
#include "app_constants.hpp"
#include "auxiliary.hpp"
#include "fem.hpp"
#include "numeric_math.hpp"
#include "point3d.hpp"
#include "test_runner.hpp"
#include "tests.hpp"

int main(unsigned int argc, const char* argv[])
{
	RunAllTests();

	try
	{
		FiniteElementModel<3> FEModel;
		try
		{
			FEModel.SetParams(argc, argv);
		}
		catch (std::exception& e)
		{
			std::cerr << e.what() << std::endl;
			exit(1);
		}

		try
		{
			FEModel.ReadData("data/in/2dim_model_001.txt", FEM::IO_FORMAT::CHEKHOV);
		}
		catch (std::exception& e)
		{
			std::cerr << e.what() << std::endl;
		}

		FEModel.WriteData("data/out/2dim_model_001.txt", FEM::IO_FORMAT::CHEKHOV);
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << std::endl;
	}
	catch (...)
	{
		std::cerr << "Unknown error occured" << std::endl;
	}

	// std::system("pause");

	return 0;
}
