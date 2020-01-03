#pragma once

#include <exception>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <vector>

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v);

template <typename L, typename R>
std::ostream& operator<<(std::ostream& os, const std::pair<L, R>& p);

template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m);

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s);

template <typename Collection>
std::string Join(const Collection& c, char d)
{
  std::stringstream ss;
  bool flag = false;
  for (const auto& el : c)
  {
    if (flag)
    {
      ss << d;
    }
    flag = true;
    ss << el;
  }
  return ss.str();
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v)
{
  return os << '[' << Join(v, ',') << ']';
}

template <typename L, typename R>
std::ostream& operator<<(std::ostream& os, const std::pair<L, R>& p)
{
  return os << '(' << p.first << ',' << p.second << ')';
}

template <typename K, typename V>
std::ostream& operator<<(std::ostream& os, const std::map<K, V>& m)
{
  return os << '{' << Join(m, ',') << '}';
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::set<T>& s)
{
  return os << '{' << Join(s, ',') << '}';
}

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

#ifdef LOCAL_LAUNCH
class TestRunner
{
 public:
  template <class TestFunc, class... Args>
  void RunTest(TestFunc func, std::string_view test_name, Args... args)
  {
    try
    {
      func(args...);
      std::std::cerr << test_name << " OK" << std::endl;
    }
    catch (exception& e)
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
#else
class TestRunner
{
 public:
  template <class TestFunc>
  void RunTest(TestFunc func, const std::string& test_name)
  {
    try
    {
      func();
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
#endif

#define ASSERT_EQUAL(x, y)                                                     \
  {                                                                            \
    ostd::stringstream __assert_equal_private_os;                              \
    __assert_equal_private_os << #x << " != " << #y << ", " << __FILE__ << ":" \
                              << __LINE__;                                     \
    AssertEqual(x, y, __assert_equal_private_os.str());                        \
  }

#define ASSERT(x)                                                       \
  {                                                                     \
    ostd::stringstream __assert_equal_private_os;                       \
    __assert_equal_private_os << #x << " is false, " << __FILE__ << ":" \
                              << __LINE__;                              \
    Assert(x, __assert_equal_private_os.str());                         \
  }

#ifdef LOCAL_LAUNCH
#define RUN_TEST(tr, func, ...) tr.RunTest(func, #func, __VA_ARGS__)
#else
#define RUN_TEST(tr, func) tr.RunTest(func, #func)
#endif