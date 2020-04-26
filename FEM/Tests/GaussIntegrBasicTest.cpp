#include "stdafx.hpp"

#include "GaussIntegrWrapper.hpp"
#include "point2d.hpp"
#include "test_runner.h"
#include <cmath>
#include <iostream>

namespace
{
template <typename T>
T fun2arr(const T* p)
{
	return p[0] * p[1];
}

template <typename T>
T fun2(const T& x, const T& y)
{
	return x * y;
}

template <typename T>
T vfun2(const point2d<T>& p)
{
	return p[0] * p[1];
}

template <typename T>
struct ffun2
{
	T operator()(T x, T y) const
	{
		return x * y;
	}
};

template <typename T>
point2d<T> vfunv(const point2d<T>& p)
{
	point2d<T> result(p);
	result[0] *= p[0];
	result[1] *= p[1];
	return result;
}
}  // namespace

void TestGaussIntegrBasic()
{
	typedef GaussIntegr::real_t real_t;

	const GaussIntegr::GaussIntegrWrapper<2, 1> integr;

	{
		real_t expected = 0;

		real_t result1;
		real_t result2 = 0;
		real_t result3;
		real_t result4 = 0;

		result1 = integr.ByPlus<real_t, real_t>(fun2<real_t>);
		integr.ByPlus<real_t>(fun2<real_t>, result2);
		result3 = integr.ByPlus_ArrArg<real_t>(fun2arr<real_t>);
		integr.ByPlus_ArrArg(fun2arr<real_t>, result4);

		ASSERT_EQUAL(expected, result1);
		ASSERT_EQUAL(expected, result2);
		ASSERT_EQUAL(expected, result3);
		ASSERT_EQUAL(expected, result4);
	}

	{
		real_t expected = 0;

		real_t result1 = integr.ByPlus_ArrArg<real_t>(fun2arr<real_t>);
		real_t result2 = integr.ByPlus<real_t, real_t>(fun2<real_t>);
		real_t result3 = integr.ByPlus<point2d<real_t>, real_t>(vfun2<real_t>);
		point2d<real_t> result4 = integr.ByPlus<point2d<real_t>, point2d<real_t>>(vfunv<real_t>);
		real_t result5 = integr.ByPlus<real_t, real_t>(ffun2<real_t>());

		ASSERT_EQUAL(expected, result1);
		ASSERT_EQUAL(expected, result2);
		ASSERT_EQUAL(expected, result3);
		// ASSERT_EQUAL(expected, result4);
		ASSERT_EQUAL(expected, result5);
	}
}
