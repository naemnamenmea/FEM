#pragma once

#include "math_constants.hpp"
#include <cmath>
#include <stdexcept>

template <class T>
T quick_pow(T x, unsigned int n)
{
	T res = 1;

	while (n)
	{
		if (n & 1)
			res *= x;
		x *= x;
		n >>= 1;
	}

	return res;
}

double round(double val, int precision);

template <typename Container>
bool CompareContainersEqualityWithTolerance(
	const Container& lhs, const Container& rhs, double tolerance = mathdef::MATH_TOL)
{
	if (lhs.size() != rhs.size())
	{
		throw std::runtime_error("Containers must be equal in size");
	}

	for (auto it = lhs.begin(); it != lhs.end(); ++it)
	{
		auto offset = std::distance(begin(lhs), it);

		if (!mathdef::is_eq(*it, *std::next(begin(rhs), offset), tolerance))
		{
			return false;
		}
	}

	return true;
}
