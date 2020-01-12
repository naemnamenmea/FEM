#pragma once

#include <string_view>
#include "math_constants.hpp"

void TrimFront(std::string_view& sv, char symbol = ' ');

void TrimBack(std::string_view& sv, char symbol = ' ');

void Trim(std::string_view& sv, char symbol = ' ');

template <typename Container>
bool CompareContainersEqualityWithTolerance(
	const Container& lhs, const Container& rhs, double tolerance = mathdef::__EPS)
{
	if (lhs.size() != rhs.size())
	{
		throw std::runtime_error("Containers must be equal in size");
	}

	for (auto it = lhs.begin(); it != std::end(lhs); ++it)
	{
		auto offset = std::distance(begin(lhs), it);

		if (!isEqual(*it, *std::next(begin(rhs), offset)))
		{
			return false;
		}
	}

	return true;
}
