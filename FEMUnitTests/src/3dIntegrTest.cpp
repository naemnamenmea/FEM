#include "GaussIntegrWrapper.hpp"
#include "point3d.hpp"
#include "math_constants.hpp"
#include "tests_helper.hpp"
#include <iostream>

namespace
{
template <typename T>
T vfun3_1(const point3d<T>& v)
{
	return v[0] * v[1] + v[2] * v[2];
}

template <typename T>
T vfun3_2(const point3d<T>& v)
{
	return v[0] * v[0] * v[2] * v[2] + v[0] * v[0] * v[0] * v[1] * v[1] * v[2] + v[2];
}
}  // namespace

void TestGaussIntegr3d()
{
	const real_t expected_1 = 8. / 3;
	const real_t expected_2 = 8. / 9;

	{
		const GaussIntegr::GaussIntegrWrapper<3, 1> integr;

		const real_t result_1 = integr.Calculate<point3d<real_t>, real_t>(vfun3_1<real_t>);
		const real_t result_2 = integr.Calculate<point3d<real_t>, real_t>(vfun3_2<real_t>);

		ASSERT_DOUBLE_EQUAL(0., result_1);
		ASSERT_DOUBLE_EQUAL(0., result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<3, 2> integr;

		const real_t result_1 = integr.Calculate<point3d<real_t>, real_t>(vfun3_1<real_t>);
		const real_t result_2 = integr.Calculate<point3d<real_t>, real_t>(vfun3_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}

	{
		const GaussIntegr::GaussIntegrWrapper<3, 3> integr;

		const real_t result_1 = integr.Calculate<point3d<real_t>, real_t>(vfun3_1<real_t>);
		const real_t result_2 = integr.Calculate<point3d<real_t>, real_t>(vfun3_2<real_t>);

		ASSERT_DOUBLE_EQUAL(expected_1, result_1);
		ASSERT_DOUBLE_EQUAL(expected_2, result_2);
	}
}
