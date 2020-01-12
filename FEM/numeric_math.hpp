#pragma once

#include "cmath"

using namespace std;

template <class T>
T quick_pow(T x, unsigned int n)
{
	T res = 1;

	while (n)
	{
		if (n & 1)
			res *= x;
		x *= x;
		n >>= 1;
	}

	return res;
}

double round(double val, int precision);
