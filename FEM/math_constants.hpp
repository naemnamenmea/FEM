#pragma once
#include <cmath>
#include <limits>

namespace mathdef
{
constexpr double __EPS = 1e-9;
constexpr double __MOD = 1e9 + 7;

constexpr double MAX_DOUBLE = std::numeric_limits<double>::max();
}  // namespace mathdef

template <typename T>
inline bool isZero(T val)
{
	return std::abs(val) < mathdef::__EPS;
}

template <typename T, typename R>
inline bool isEqual(T a, R b)
{
	return isZero(a - b);
}

template <typename T>
inline bool isEqual(T a, T b)
{
	return isZero(a - b);
}
