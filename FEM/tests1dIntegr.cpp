// Author: Chekhov Vladimir Valerevich

#include "stdafx.hpp"

#include <cmath>
#include <iostream>
#include "GaussIntegr.hpp"
#include "point1d.hpp"
#include "tests.hpp"

// различные примеры подынтегральных выражений дл€ трЄхмерной задачи:

template <typename T>
T fun1arr(const T* x)
{
	return (5 * x[0] * x[0] * x[0] - 3 * x[0]) / 0.5;
	// return std::abs(x[0]);
}

template <typename T>
T fun1(T x)
{
	return x;
}

template <typename T>
struct ffun1
{  // функциональный объект
	T operator()(T x) const
	{
		return x > 0 ? x : 0;
	}
};

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

//--------------------------------------------------------------------------

void TestBasic1d()
{
	typedef GaussIntegr::real_t real_t;

	// создание объекта, выполн€ющего интегрирование по √ауссу дл€ трЄхмерного
	// случа€ с пор€дком 4:
	static const GaussIntegr::fIntegrate<1, 4> integr;

	real_t result(0.);

	// выполнение численного интегрировани€:

	// с использованием операции +=

	// std::cout << integr.ByPlusAssgn_ArrArg(fun1arr<real_t>, result) << '\n';

	///*arg*/  // - тип чего приходитс€ указывать при инстанцировании объекта
	std::cout << integr.ByPlusAssgn<real_t>(ffun1<real_t>(), result) << '\n';
	// std::cout << integr.ByPlusAssgn<point1d<real_t>>(vfun1<real_t>, result)
	//          << '\n';

	//// с использованием операции +

	///*ret*/
	// std::cout << integr.ByPlus_ArrArg<real_t>(fun1arr<real_t>) << '\n';

	///*arg,ret*/
	// std::cout << integr.ByPlus<real_t, real_t>(ffun1<real_t>()) << '\n';
	// std::cout << integr.ByPlus<point1d<real_t>, real_t>(vfun1<real_t>) << '\n';

	//// интегрирование вектор-функций:
	// point1d<real_t> v;
	// integr.ByPlusAssgn<point1d<real_t>>(vfunv<real_t>, v);
	// std::cout << v << '\n';

	// v = integr.ByPlus<point1d<real_t>, point1d<real_t>>(vfunv<real_t>);
	// std::cout << v << '\n';
}
//--------------------------------------------------------------------------
