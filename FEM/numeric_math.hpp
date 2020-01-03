#include "cmath"

using namespace std;

template <class T>
T quick_pow(T x, unsigned int n)
{
  T res = 1;

  while (n)
  {
    if (n & 1) res *= x;
    x *= x;
    n >>= 1;
  }

  return res;
}

inline double round(double val, int precision)
{
  double mult = pow(10, precision);
  val *= mult;
  return (val < 0 ? ceil(val - 0.5) : floor(val + 0.5)) / mult;
}