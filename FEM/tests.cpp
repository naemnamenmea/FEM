#include "tests.hpp"
#include "test_runner.hpp"

void RunAllTests()
{
  TestRunner tr;
  RUN_TEST(tr, TestBasic2d);
  RUN_TEST(tr, TestBasic3d);
}