#include "stdafx.hpp"

#include "auxiliary.hpp"

std::string_view& TrimFront(std::string_view& str, char delimiter)
{
	while (!str.empty() && str.front() == delimiter)
	{
		str.remove_prefix(1);
	}

	return str;
}

std::string_view& TrimBack(std::string_view& str, char delimiter)
{
	while (!str.empty() && str.back() == delimiter)
	{
		str.remove_suffix(1);
	}

	return str;
}

std::string_view& Trim(std::string_view& str, char delimiter)
{
	TrimFront(str, delimiter);
	TrimBack(str, delimiter);

	return str;
}
