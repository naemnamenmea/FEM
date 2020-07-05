#pragma once

#include "test_runner.hpp"
#include "math_constants.hpp"

struct DoubleComparator
{
	bool operator()(const double& lhs, const double& rhs)
	{
		return mathdef::is_eq(lhs, rhs);
	}
};

#define ASSERT_DOUBLE_EQUAL(x, y) ASSERT_EQUAL_CMP(x, y, DoubleComparator())
