#include "GaussIntegrWrapper.hpp"
#include "point1d.hpp"
#include "math_constants.hpp"
#include "tests_helper.hpp"
#include <iostream>

namespace
{
template <typename T>
T vfun1_1(const point1d<T>& v)
{
	return v[0] * v[0];
}

template <typename T>
T vfun1_2(const point1d<T>& v)
{
	return v[0] * (v[0] - 1);
}
}  // namespace

void TestGaussIntegr1d()
{
	const real_t expected_1 = 2. / 3;
	const real_t expected_2 = 2. / 3;

	{
		const GaussIntegr::GaussIntegrWrapper<1, 1> integr;

		const real_t result_1 = integr.Calculate<point1d<real_t>, real_t>(vfun1_1<real_t>);
		const real_t result_2 = integr.Calculate<point1d<real_t>, real_t>(vfun1_2<real_t>);

		ASSERT_DOUBLE_EQUAL(0., result_1);
		ASSERT_DOUBLE_EQUAL(0., result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<1, 2> integr;

		const real_t result_1 = integr.Calculate<point1d<real_t>, real_t>(vfun1_1<real_t>);
		const real_t result_2 = integr.Calculate<point1d<real_t>, real_t>(vfun1_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<1, 3> integr;

		const real_t result_1 = integr.Calculate<point1d<real_t>, real_t>(vfun1_1<real_t>);
		const real_t result_2 = integr.Calculate<point1d<real_t>, real_t>(vfun1_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}
}
