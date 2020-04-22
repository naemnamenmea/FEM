// Author: Chekhov Vladimir Valerevich

#include "stdafx.hpp"

#include <iostream>
#include "GaussIntegrWrapper.hpp"
#include "point3d.hpp"
#include "tests.hpp"

// различные примеры подынтегральных выражений для трёхмерной задачи:

template <typename T>
T fun3arr(const T* xyz)
{
	return (1 + xyz[0]) * xyz[1] * xyz[2];
}

template <typename T>
T fun3(T x, T y, T z)
{
	return x * y * z;
}

template <typename T>
struct ffun3
{  // функциональный объект
	T operator()(T x, T y, T z) const
	{
		return (x + 1) * y * z;
	}
};

template <typename T>
T vfun3(point3d<T> v)
{
	return v[0] / v[1] * v[2];
}

template <typename T>
point3d<T> vfunv(const point3d<T>& v)
{
	point3d<T> result(v);
	result += v;
	return result;
}

void TestGaussIntegr3d()
{
	typedef GaussIntegr::real_t real_t;

	static const GaussIntegr::GaussIntegrWrapper<3, 4> integr;

	real_t result(0.);

	// выполнение численного интегрирования:

	// с использованием операции +=

	std::cout << integr.ByPlus_ArrArg(fun3arr<real_t>, result) << '\n';

	/*arg*/  // - тип чего приходится указывать при инстанцировании объекта
	std::cout << integr.ByPlus<real_t>(ffun3<real_t>(), result) << '\n';
	std::cout << integr.ByPlus<point3d<real_t>>(vfun3<real_t>, result) << '\n';
	std::cout << integr.ByPlus<real_t>(fun3<real_t>, result) << '\n';

	// с использованием операции +

	/*ret*/
	std::cout << integr.ByPlus_ArrArg<real_t>(fun3arr<real_t>) << '\n';

	/*arg,ret*/
	std::cout << integr.ByPlus<real_t, real_t>(ffun3<real_t>()) << '\n';
	std::cout << integr.ByPlus<point3d<real_t>, real_t>(vfun3<real_t>) << '\n';
	std::cout << integr.ByPlus<real_t ,real_t>(fun3<real_t>) << '\n';

	// интегрирование вектор-функций:
	point3d<real_t> v;

	std::cout << integr.ByPlus<point3d<real_t>>(vfunv<real_t>, v) << '\n';
	std::cout << integr.ByPlus<point3d<real_t>, point3d<real_t>>(vfunv<real_t>) << '\n';
}
