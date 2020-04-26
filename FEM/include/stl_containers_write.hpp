#pragma once

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

template <class Container>
std::string Join(const Container& container, char delimiter);

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
	return os << '[' << Join(v, ',') << ']';
}

template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m)
{
	return os << '{' << Join(m, ',') << '}';
}

template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<K, V>& m)
{
	return os << '{' << Join(m, ',') << '}';
}

template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::unordered_multimap<K, V>& m)
{
	return os << '{' << Join(m, ',') << '}';
}

template <typename L, typename R>
std::ostream& operator<<(std::ostream& os, const std::pair<L, R>& p)
{
	return os << '(' << p.first << ',' << p.second << ')';
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
	return os << '{' << Join(s, ',') << '}';
}

template <typename T, class Comp>
std::ostream& operator<<(std::ostream& os, const std::multiset<T, Comp>& s)
{
	return os << '{' << Join(s, ',') << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& s)
{
	return os << '(' << Join(s, ',') << ')';
}

template <typename Container>
std::string Join(const Container& container, char delimiter)
{
	std::stringstream ss;

	bool flag = false;
	for (const auto& el : container)
	{
		if (flag)
		{
			ss << delimiter;
		}
		flag = true;
		ss << el;
	}
	return ss.str();
}
