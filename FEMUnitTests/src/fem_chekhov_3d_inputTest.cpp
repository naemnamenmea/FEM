#include "stdafx.hpp"

#include <filesystem>
#include "auxiliary.hpp"
#include "fem.hpp"
#include "math_constants.hpp"
#include "tests_helper.hpp"

namespace
{
namespace fs = std::filesystem;

fs::path inputDir = "data";
fs::path outputDir = inputDir;
std::string fileName1 = "2dim_model_001.txt";
std::string fileName2 = "no_end_tag.txt";
}  // namespace

void TestFEMChekhovInputNodes()
{
	fs::path inputFile = inputDir / fileName1;

	{
		// test node:
		// 15 36.3 8.185 0. free free fixed

		FiniteElementModel<double, 3> feModel;
		feModel.ReadData(inputFile);
		const auto& nodes = feModel.GetNodes();

		ASSERT_EQUAL(nodes.size(), 45);

		auto it = nodes.find(15);

		ASSERT(it != nodes.end());

		const auto& node = it->second;
		const auto& coord = node.GetCoord();
		const auto& freedomDegrees = node.GetDegreesOfFreedom();

		{
			const Tensor1s<double> expected(36.3, 8.185, 0.);

			ASSERT_EQUAL(coord, expected);
		}

		{
			using FD = FEM::DEGREES_OF_FREEDOM;
			decltype(freedomDegrees) expected = {FD::FREE, FD::FREE, FD::FIXED};

			ASSERT_EQUAL(freedomDegrees, expected);
		}
	}
}

void TestFEMChekhovInputMaterials()
{
	fs::path inputFile = inputDir / fileName1;

	{
		// test material:
		// Al_Alloy 0.0028

		FiniteElementModel<double, 3> feModel;
		feModel.ReadData(inputFile);
		const auto& materials = feModel.GetMaterials();

		{
			auto it = materials.find({"Al_Alloy"});
			ASSERT(it != materials.end());

			const auto& material = *it;
			const auto& materialProp = material.GetProperties();

			ASSERT_DOUBLE_EQUAL(materialProp.GetDensity(), 0.0028);
		}

		{
			auto it = materials.find({"Steel"});
			ASSERT(it != materials.end());

			const auto& material = *it;
			const auto& materialProp = material.GetProperties();

			ASSERT_DOUBLE_EQUAL(materialProp.GetDensity(), 0.0078);
		}
	}
}

void TestFEMChekhovInputFiniteElements()
{
	fs::path inputFile = inputDir / fileName1;

	{
		// test finite element:

		FiniteElementModel<double, 3> feModel;
		feModel.ReadData(inputFile);
		const auto& finiteElements = feModel.GetFiniteElements();

		ASSERT_EQUAL(finiteElements.size(), 56);

		// 17 Quad4 25 26 29 28 Al_Alloy .1
		{
			auto it = finiteElements.find(17);

			ASSERT(it != finiteElements.end());

			const auto& fe = it->second;

			ASSERT_EQUAL(fe->SpaceDimension(), 2);

			ASSERT_EQUAL(fe->GetMaterial().GetName(), "Al_Alloy");

			if (fe->SpaceDimension() < 3)
			{
				ASSERT_DOUBLE_EQUAL(fe->GetParameter(), 0.1);
			}
		}

		// 37 Rod2 25 28 Steel 1.
		{
			auto it = finiteElements.find(37);

			ASSERT(it != finiteElements.end());

			const auto& fe = it->second;

			ASSERT_EQUAL(fe->SpaceDimension(), 1);

			ASSERT_EQUAL(fe->GetMaterial().GetName(), "Steel");

			if (fe->SpaceDimension() < 3)
			{
				ASSERT_DOUBLE_EQUAL(fe->GetParameter(), 1.);
			}
		}
	}
}
