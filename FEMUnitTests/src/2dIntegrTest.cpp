#include "GaussIntegrWrapper.hpp"
#include "point2d.hpp"
#include "math_constants.hpp"
#include "tests_helper.hpp"
#include <iostream>

namespace
{
template <typename T>
T vfun2_1(const point2d<T>& v)
{
	return v[0] - v[1] * v[1];
}

template <typename T>
T vfun2_2(const point2d<T>& v)
{
	return v[1] * v[1] * (v[0] + 1);
}
}  // namespace

void TestGaussIntegr2d()
{
	const real_t expected_1 = -4. / 3;
	const real_t expected_2 = 4. / 3;

	{
		const GaussIntegr::GaussIntegrWrapper<2, 1> integr;

		const real_t result_1 = integr.Calculate<point2d<real_t>, real_t>(vfun2_1<real_t>);
		const real_t result_2 = integr.Calculate<point2d<real_t>, real_t>(vfun2_2<real_t>);

		ASSERT_DOUBLE_EQUAL(0., result_1);
		ASSERT_DOUBLE_EQUAL(0., result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<2, 2> integr;

		const real_t result_1 = integr.Calculate<point2d<real_t>, real_t>(vfun2_1<real_t>);
		const real_t result_2 = integr.Calculate<point2d<real_t>, real_t>(vfun2_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<2, 3> integr;

		const real_t result_1 = integr.Calculate<point2d<real_t>, real_t>(vfun2_1<real_t>);
		const real_t result_2 = integr.Calculate<point2d<real_t>, real_t>(vfun2_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}
}
