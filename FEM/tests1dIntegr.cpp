// Author: Chekhov Vladimir Valerevich

#include "stdafx.hpp"

#include <cmath>
#include <iostream>
#include "GaussIntegrWrapper.hpp"
#include "point1d.hpp"
#include "tests.hpp"

// ��������� ������� ��������������� ��������� ��� ��������� ������:

template <typename T>
T fun1arr(const T* x)
{
	return (5 * x[0] * x[0] * x[0] - 1 * x[0]) / 0.5;
	// return std::abs(x[0]);
}

template <typename T>
T fun1(T x)
{
	return x;
}

template <typename T>
T vfun1(point1d<T> v)
{
	return v[0];
}

template <typename T>
point1d<T> vfunv(const point1d<T>& v)
{
	point1d<T> result(v);
	result += v;
	return result;
}

void TestGaussIntegr1d()
{
	typedef GaussIntegr::real_t real_t;

	static const GaussIntegr::GaussIntegrWrapper<1, 1> integr;

	// ���������� ���������� ��������������:

	{
		// � �������������� �������� +=

		real_t result(0.);

		std::cout << integr.ByPlus_ArrArg(fun1arr<real_t>, result) << '\n';

		/*arg*/  // - ��� ���� ���������� ��������� ��� ��������������� �������
		std::cout << integr.ByPlus<real_t>(vfun1<real_t>, result) << '\n';
		std::cout << integr.ByPlus<real_t>(fun1<real_t>, result) << '\n';
	}

	{
		// � �������������� �������� +

		/*ret*/
		std::cout << integr.ByPlus_ArrArg<real_t>(fun1arr<real_t>) << '\n';

		/*arg,ret*/
		std::cout << integr.ByPlus<real_t, real_t>(vfun1<real_t>) << '\n';
		std::cout << integr.ByPlus<real_t, real_t>(fun1<real_t>) << '\n';
	}

	{
		// �������������� ������-�������:
		point1d<real_t> v;

		std::cout << integr.ByPlus<point1d<real_t>>(vfunv<real_t>, v) << '\n';
		std::cout << integr.ByPlus<point1d<real_t>, point1d<real_t>>(vfunv<real_t>) << '\n';
	}
}
