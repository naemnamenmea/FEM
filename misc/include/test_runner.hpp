#pragma once

#include "stl_containers_write.hpp"
#include "ThrowMessage.hpp"
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

template <class T, class U, class Comparator = std::equal_to<T>>
void AssertEqual(
	const T& t, const U& u, const std::string& hint = {}, Comparator cmp = Comparator())
{
#pragma warning(push)
#pragma warning(disable : 4389)
	if (!cmp(t, u))
#pragma warning(pop)
	{
		std::stringstream os;
		os << "Assertion failed: " << t << " != " << u;
		if (!hint.empty())
		{
			os << " hint: " << hint;
		}
		throw std::runtime_error(os.str());
	}
}

inline void Assert(bool b, const std::string& hint)
{
	AssertEqual(b, true, hint, std::equal_to<bool>());
}

class TestRunner
{
public:
	template <class TestFunc, class... Args>
	void RunTest(TestFunc func, std::string_view test_name, Args... args)
	{
		try
		{
			func(args...);
			std::cerr << test_name << " OK" << std::endl;
		}
		catch (const tMessage& e)
		{
			++m_failCount;
			std::cerr << test_name << " fail: " << e << std::endl;
		}
		catch (std::exception& e)
		{
			++m_failCount;
			std::cerr << test_name << " fail: " << e.what() << std::endl;
		}
		catch (...)
		{
			++m_failCount;
			std::cerr << test_name << " fail: Unknown exception caught" << std::endl;
		}
	}

	~TestRunner()
	{
		if (m_failCount > 0)
		{
			std::cerr << m_failCount << " unit tests failed. Terminate" << std::endl;
			exit(1);
		}
		else
		{
			std::cerr << "All tests passed." << std::endl;
		}
	}

private:
	int m_failCount = 0;
};

#define ASSERT_EQUAL_CMP(x, y, comparator)                                                      \
	{                                                                                           \
		std::stringstream __assert_equal_private_os;                                            \
		__assert_equal_private_os << #x << " != " << #y << ", " << __FILE__ << ":" << __LINE__; \
		AssertEqual(x, y, __assert_equal_private_os.str(), comparator);                         \
	}

#define ASSERT_EQUAL(x, y)                                                                      \
	{                                                                                           \
		std::stringstream __assert_equal_private_os;                                            \
		__assert_equal_private_os << #x << " != " << #y << ", " << __FILE__ << ":" << __LINE__; \
		AssertEqual(x, y, __assert_equal_private_os.str());                                     \
	}

#define ASSERT(x)                                                                        \
	{                                                                                    \
		std::stringstream __assert_equal_private_os;                                     \
		__assert_equal_private_os << #x << " is false, " << __FILE__ << ":" << __LINE__; \
		Assert(x, __assert_equal_private_os.str());                                      \
	}

#define RUN_TEST(tr, func, ...) tr.RunTest(func, #func, ##__VA_ARGS__)
