#include "stdafx.hpp"

#include <cmath>
#include <iostream>
#include "GaussIntegrWrapper.hpp"
#include "point1d.hpp"
#include "tests.hpp"

// различные примеры подынтегральных выражений для трёхмерной задачи:

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

void TestGaussIntegrBasic()
{
	typedef GaussIntegr::real_t real_t;

	static const GaussIntegr::GaussIntegrWrapper<1, 1> integr;

	// выполнение численного интегрирования:

	/*
		x16 случаев рассмотреть: 4 вида функций х 4 метода
	*/

	{
		// с использованием операции +=

		{
			real_t numRes(0.);	
			point1d<real_t> vecRes(0.);	

			std::cout << integr.ByPlus_ArrArg(fun1arr<real_t>, numRes) << '\n';
			//std::cout << integr.ByPlus_ArrArg(fun1<real_t>, numRes) << '\n';
			//std::cout << integr.ByPlus_ArrArg(vfun1<real_t>, numRes) << '\n';
			//std::cout << integr.ByPlus_ArrArg(vfunv<real_t>, vecRes) << '\n';
		}

		{
			//real_t result(0.);	

			//std::cout << integr.ByPlus(fun1arr<real_t>, result) << '\n';
			//std::cout << integr.ByPlus(fun1<real_t>, result) << '\n';
			//std::cout << integr.ByPlus(vfun1<real_t>, result) << '\n';
			//std::cout << integr.ByPlus(vfunv<real_t>, result) << '\n';
		}
	}

	{
		// с использованием операции +

		/*ret*/
		std::cout << integr.ByPlus_ArrArg<real_t>(fun1arr<real_t>) << '\n';

		/*arg,ret*/
		std::cout << integr.ByPlus<point1d<real_t>, real_t>(vfun1<real_t>) << '\n';
		std::cout << integr.ByPlus<real_t, real_t>(fun1<real_t>) << '\n';
	}

	{
		// интегрирование вектор-функций:
		point1d<real_t> v;

		std::cout << integr.ByPlus<real_t>(vfunv<real_t>, v) << '\n';
		std::cout << integr.ByPlus<real_t, point1d<real_t>>(vfunv<real_t>) << '\n';
	}
}
