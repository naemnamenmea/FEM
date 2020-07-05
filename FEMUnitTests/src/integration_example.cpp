#include "FiniteElementModel.hpp"
#include "FiniteElementPool.hpp"
#include "GaussIntegrWrapper.hpp"
#include "point1d.hpp"
#include "point2d.hpp"
#include "point3d.hpp"
#include "math_constants.hpp"
#include "tests_helper.hpp"
#include <iostream>
#include <iomanip>
#include <vector>
#include <utility>
#include "tchar.h"

const std::string& GetTestDataDir();

static const std::string GetInData()
{
	return GetTestDataDir() + _T("/Folder/");
}

#define PATH GetInData()

template <int DIM, int ORD>
using GaussIntegrator = GaussIntegr::GaussIntegrWrapper<DIM, ORD>;

namespace
{
template <typename T>
T func1d(const point1d<T>& v)
{
	return 30. * pow(v[0], 6);
}

template <typename T>
T func2d(const point2d<T>& v)
{
	return exp(-v[0] - v[1]);
}

template <typename T>
T func3d(const point3d<T>& v)
{
	return exp(-v[0] - v[1] - v[2]);
}
}  // namespace

void TestGaussIntegr1dOrder()
{
	const size_t DIM = 1;
	const double expected = 60. / 7.;
	double prevDiff = std::numeric_limits<double>::max();
	std::vector<std::pair<int, double>> results;

	results.emplace_back(
		1, GaussIntegrator<DIM, 1>().Calculate<point1d<double>, double>(func1d<double>));
	results.emplace_back(
		2, GaussIntegrator<DIM, 2>().Calculate<point1d<double>, double>(func1d<double>));
	results.emplace_back(
		3, GaussIntegrator<DIM, 3>().Calculate<point1d<double>, double>(func1d<double>));
	results.emplace_back(
		4, GaussIntegrator<DIM, 4>().Calculate<point1d<double>, double>(func1d<double>));
	results.emplace_back(
		5, GaussIntegrator<DIM, 5>().Calculate<point1d<double>, double>(func1d<double>));
	// results.emplace_back(
	//	6, GaussIntegrator<DIM, 6>().Calculate<point1d<double>, double>(func1d<double>));

	std::cerr << std::setprecision(15);
	std::cerr << "expected = " << expected << std::endl;
	int enoughOrder = -1;
	for (auto& [order, result] : results)
	{
		const double diff = abs(result - expected);
		std::cerr << "order = " << order << "; result = " << result << "; abs_diff = " << diff
				  << std::endl;
		ASSERT(mathdef::is_lte(diff, prevDiff));
		if (enoughOrder == -1 && mathdef::is_zero(diff))
		{
			enoughOrder = order;
		}
		prevDiff = diff;
	}
	std::cerr << "enoughOrder: " << enoughOrder << std::endl;
}

void TestGaussIntegr2dOrder()
{
	const size_t DIM = 2;
	const double expected = 5.52439138216726;
	double prevDiff = std::numeric_limits<double>::max();
	std::vector<std::pair<int, double>> results;

	results.emplace_back(
		1, GaussIntegrator<DIM, 1>().Calculate<point2d<double>, double>(func2d<double>));
	results.emplace_back(
		2, GaussIntegrator<DIM, 2>().Calculate<point2d<double>, double>(func2d<double>));
	results.emplace_back(
		3, GaussIntegrator<DIM, 3>().Calculate<point2d<double>, double>(func2d<double>));
	results.emplace_back(
		4, GaussIntegrator<DIM, 4>().Calculate<point2d<double>, double>(func2d<double>));
	results.emplace_back(
		5, GaussIntegrator<DIM, 5>().Calculate<point2d<double>, double>(func2d<double>));
	// results.emplace_back(
	//	6, GaussIntegrator<DIM, 6>().Calculate<point2d<double>, double>(func2d<double>));

	std::cerr << std::setprecision(15);
	std::cerr << "expected = " << expected << std::endl;
	int enoughOrder = -1;
	for (auto& [order, result] : results)
	{
		const double diff = abs(result - expected);
		std::cerr << "order = " << order << "; result = " << result << "; abs_diff = " << diff
				  << std::endl;
		ASSERT(mathdef::is_lte(diff, prevDiff));
		if (enoughOrder == -1 && mathdef::is_zero(diff))
		{
			enoughOrder = order;
		}
		prevDiff = diff;
	}
	std::cerr << "enoughOrder: " << enoughOrder << std::endl;
}

void TestGaussIntegr3dOrder()
{
	const size_t DIM = 3;
	const double expected = 12.984542692957;
	double prevDiff = std::numeric_limits<double>::max();
	std::vector<std::pair<int, double>> results;

	results.emplace_back(
		1, GaussIntegrator<DIM, 1>().Calculate<point3d<double>, double>(func3d<double>));
	results.emplace_back(
		2, GaussIntegrator<DIM, 2>().Calculate<point3d<double>, double>(func3d<double>));
	results.emplace_back(
		3, GaussIntegrator<DIM, 3>().Calculate<point3d<double>, double>(func3d<double>));
	results.emplace_back(
		4, GaussIntegrator<DIM, 4>().Calculate<point3d<double>, double>(func3d<double>));
	results.emplace_back(
		5, GaussIntegrator<DIM, 5>().Calculate<point3d<double>, double>(func3d<double>));
	// results.emplace_back(
	//	6, GaussIntegrator<DIM, 6>().Calculate<point3d<double>, double>(func3d<double>));

	std::cerr << std::setprecision(15);
	std::cerr << "expected = " << expected << std::endl;
	int enoughOrder = -1;
	for (auto& [order, result] : results)
	{
		const double diff = abs(result - expected);
		std::cerr << "order = " << order << "; result = " << result << "; abs_diff = " << diff
				  << std::endl;
		ASSERT(mathdef::is_lte(diff, prevDiff));
		if (enoughOrder == -1 && mathdef::is_zero(diff))
		{
			enoughOrder = order;
		}
		prevDiff = diff;
	}
	std::cerr << "enoughOrder: " << enoughOrder << std::endl;
}

void Test3dFiniteElement()
{
	tFE_model fem;
	fem.ImportFromFile((PATH + "test3dFE.xml").c_str());

	ASSERT_DOUBLE_EQUAL(fem.GetFE(1).Volume(), fem.GetFE(2).Volume());
	ASSERT_DOUBLE_EQUAL(fem.GetFE(3).Volume(), 14.);  // why negative value occurs?
	ASSERT_DOUBLE_EQUAL(fem.GetFE(4).Volume(), 48.);  // why negative value occurs?
}
