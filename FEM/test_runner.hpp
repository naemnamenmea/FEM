#pragma once

#include "stl_containers_write.hpp"
#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

template <class T, class U>
void AssertEqual(const T& t, const U& u, const std::string& hint = {})
{
	if (!(t == u))
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
	AssertEqual(b, true, hint);
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
		catch (std::exception& e)
		{
			++fail_count;
			std::cerr << test_name << " fail: " << e.what() << std::endl;
		}
		catch (...)
		{
			++fail_count;
			std::cerr << "Unknown exception caught" << std::endl;
		}
	}

	~TestRunner()
	{
		if (fail_count > 0)
		{
			std::cerr << fail_count << " unit tests failed. Terminate" << std::endl;
			exit(1);
		}
		else
		{
			std::cerr << "All tests passed." << std::endl;
		}
	}

private:
	int fail_count = 0;
};

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

#define RUN_TEST(tr, func, ...) tr.RunTest(func, #func, __VA_ARGS__)
