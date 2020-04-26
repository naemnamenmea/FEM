#pragma once
#include <cmath>
#include <limits>

namespace mathdef
{
const double MATH_TOL = 1e-12;

constexpr double __EPS = 1e-9;
constexpr double __MOD = 1e9 + 7;

constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();

template <typename T>
inline bool isZero(T val)
{
	return std::abs(val) < __EPS;
}

inline double math_tol(const double&)
{
	return MATH_TOL;
}

inline bool isEqual(const double& a, const double& b)
{
	return abs(a - b) < MATH_TOL;
}

inline bool isEqual(const double& a, const double& b, const double& tol)
{
	return abs(a - b) < tol;
}
}  // namespace mathdef
