#pragma once

#include <cmath>
#include <limits>

namespace mathdef
{
const double MATH_TOL = 1e-12;
const float MATH_TOL_FLOAT = 1e-6f;
const double MATH_TOL_LD = 1e-16;

constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();

inline double math_tol(const double&)
{
	return MATH_TOL;
}

inline bool is_eq(const double& a, const double& b)
{
	return abs(a - b) < MATH_TOL;
}

inline const bool is_neq(const double& val1, const double& val2)
{
	return !is_eq(val1, val2);
}

inline bool is_eq(const long double& a, const long double& b)
{
	return abs(a - b) < MATH_TOL_LD;
}

inline bool is_eq(const double& a, const double& b, const double& tol)
{
	return abs(a - b) < tol;
}

inline const bool is_lt(const double& toCompare, const double& source)
{
	return toCompare < (source - MATH_TOL);
}

inline const bool is_gt(const double& toCompare, const double& source)
{
	return toCompare > (source + MATH_TOL);
}

inline const bool is_lte(const double& toCompare, const double& source)
{
	return toCompare <= (source + MATH_TOL);
}

inline const bool is_gte(const double& toCompare, const double& source)
{
	return toCompare >= (source - MATH_TOL);
}

inline const bool is_lte(const long double& toCompare, const long double& source)
{
	return toCompare <= (source + MATH_TOL);
}

inline const bool is_zero(const long double& value)
{
	return fabs(value) <= MATH_TOL_LD;
}

inline const bool is_not_zero(const long double& value)
{
	return !is_zero(value);
}
}  // namespace mathdef
