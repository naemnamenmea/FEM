#include "stdafx.hpp"

#include <filesystem>
#include "auxiliary.hpp"
#include "fem.hpp"
#include "math_constants.hpp"
#include "test_runner.h"
#include "tests.hpp"

namespace
{
namespace fs = std::filesystem;

fs::path inputDir = "data/in";
fs::path outputDir = "data/out";
std::string fileName1 = "2dim_model_001.txt";
std::string fileName2 = "no_end_tag.txt";
}  // namespace

void TestFEMChekhovInputNodes()
{
	fs::path inputFile = inputDir / fileName1;

	{
		// test node:
		// 15 36.3 8.185 0. free free fixed

		FiniteElementModel<3> feModel;
		feModel.ReadData(inputFile);
		const auto& nodes = feModel.GetNodes();

		ASSERT_EQUAL(nodes.size(), 45);

		auto it = nodes.find(15);

		ASSERT(it != std::end(nodes));

		const auto& node = it->second;
		const auto& coord = node.GetCoord();
		const auto& freedomDegrees = node.GetFreedomDegrees();

		ASSERT(CompareContainersEqualityWithTolerance(coord, {36.3, 8.185, 0.}));

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

		FiniteElementModel<3> feModel;
		feModel.ReadData(inputFile);
		const auto& materials = feModel.GetMaterials();

		{
			auto it = materials.find({"Al_Alloy"});
			ASSERT(it != std::end(materials));

			const auto& material = *it;
			const auto& materialProp = material.GetProperties();

			ASSERT(isEqual(materialProp.GetDensity(), 0.0028));
		}

		{
			auto it = materials.find({"Steel"});
			ASSERT(it != std::end(materials));

			const auto& material = *it;
			const auto& materialProp = material.GetProperties();

			ASSERT(isEqual(materialProp.GetDensity(), 0.0078));
		}
	}
}

void TestFEMChekhovInputFiniteElements()
{
	fs::path inputFile = inputDir / fileName1;

	{
		// test finite element:

		FiniteElementModel<3> feModel;
		feModel.ReadData(inputFile);
		const auto& finiteElements = feModel.GetFiniteElements();

		ASSERT_EQUAL(finiteElements.size(), 56);

		// 17 Quad4 25 26 29 28 Al_Alloy .1
		{
			auto it = finiteElements.find(17);

			ASSERT(it != std::end(finiteElements));

			const auto& fe = it->second;

			ASSERT_EQUAL(fe->GetDim(), 2);

			const auto& nodes = fe->GetNodeNumbers();

			ASSERT_EQUAL(nodes.size(), 4);

			{
				decltype(nodes) expected = {25, 26, 29, 28};
				ASSERT(nodes == expected);
			}

			ASSERT_EQUAL(fe->GetMaterial().GetName(), "Al_Alloy");

			if (fe->GetDim() < 3)
			{
				ASSERT(isEqual(fe->GetParameter(), 0.1));
			}
		}

		// 37 Rod2 25 28 Steel 1.
		{
			auto it = finiteElements.find(37);

			ASSERT(it != std::end(finiteElements));

			const auto& fe = it->second;

			ASSERT_EQUAL(fe->GetDim(), 1);

			const auto& nodes = fe->GetNodeNumbers();

			ASSERT_EQUAL(nodes.size(), 2);

			{
				decltype(nodes) expected = {25, 28};
				ASSERT(nodes == expected);
			}

			ASSERT_EQUAL(fe->GetMaterial().GetName(), "Steel");

			if (fe->GetDim() < 3)
			{
				ASSERT(isEqual(fe->GetParameter(), 1.));
			}
		}
	}
}
