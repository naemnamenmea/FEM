#pragma once

#include <string_view>
#include <stdexcept>
#include "math_constants.hpp"

std::string_view& TrimFront(std::string_view& str, char delimiter = ' ');

std::string_view& TrimBack(std::string_view& str, char delimiter = ' ');

std::string_view& Trim(std::string_view& str, char delimiter = ' ');

template <typename Container>
bool CompareContainersEqualityWithTolerance(
	const Container& lhs, const Container& rhs, double tolerance = mathdef::__EPS)
{
	if (lhs.size() != rhs.size())
	{
		throw std::runtime_error("Containers must be equal in size");
	}

	for (auto it = lhs.begin(); it != lhs.end(); ++it)
	{
		auto offset = std::distance(begin(lhs), it);

		if (!mathdef::isEqual(*it, *std::next(begin(rhs), offset), tolerance))
		{
			return false;
		}
	}

	return true;
}
