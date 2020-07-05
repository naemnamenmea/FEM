#pragma once

#include "NumTypes.hpp"

#include <string>
#include <vector>
#include <utility>
#include <algorithm>

/*template <typename T> struct fLessName
{
 bool operator()(std::pair<std::string,T> pr1, std::pair<std::string,T> pr2) const
		 {return pr1.first < pr2.first;}
};*/

template <typename T>
struct fTheSame2nd
{
	fTheSame2nd(T p) : m_second(p)
	{
	}
	bool operator()(std::pair<std::string, T> pr) const
	{
		return pr.second == m_second;
	}

	T m_second;
};

struct fTheSame1st
{
	fTheSame1st(const std::string& str) : m_first(str)
	{
	}
	template <typename T>
	bool operator()(std::pair<std::string, T> pr) const
	{
		return pr.first == m_first;
	}

	const std::string& m_first;
};

template <typename T>
class tNamesToFuns
{
public:
	typedef T* (*tPtrToFun)();

	tNamesToFuns(const char* names_[], /*const*/ tPtrToFun functions_[], size_t dim_) : m_array()
	{
		m_array.reserve(dim_);
		for (size_t i = 0; i < dim_; ++i)
			m_array.push_back(std::make_pair(std::string(names_[i]), functions_[i]));
		//      std::sort(Array.begin(),Array.end(),fLessName<tPtrToFun>());
	}
	const std::string& FindName(tPtrToFun p) const
	{
		return find_if(m_array.begin(), m_array.end(), fTheSame2nd<tPtrToFun>(p))->first;
	}
	T* CallFunction(const std::string& name_) const
	{
		return find_if(m_array.begin(), m_array.end(), fTheSame1st(name_))->second();
	}

private:
	std::vector<std::pair<std::string, tPtrToFun> > m_array;
};
