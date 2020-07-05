#pragma once

#include <iostream>
#include <vector>

template <typename T>
std::istream& operator>>(std::istream& is, std::vector<T>& v)
{
	for (size_t i = 0; i < v.size(); ++i)
	{
		is >> v[i];
	}
	return is;
}
