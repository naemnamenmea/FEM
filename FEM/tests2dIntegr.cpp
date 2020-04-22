// Author: Chekhov Vladimir Valerevich

#include "stdafx.hpp"

#include <iostream>
#include "GaussIntegrWrapper.hpp"
#include "point2d.hpp"
#include "tests.hpp"

// ��������� ������� ��������������� ��������� ��� ��������� ������:

template <typename T>
T fun2arr(const T* xy)
{
	return xy[0] * xy[1];
}

template <typename T>
T fun2(T x, T y)
{
	return x * y;
}

template <typename T>
struct ffun2
{  // �������������� ������
	T operator()(T x, T y) const
	{
		return x * y;
	}
};

template <typename T>
T vfun2(point2d<T> v)
{
	return v[0] * v[1];
}

template <typename T>
point2d<T> vfunv(const point2d<T>& v)
{
	point2d<T> result(v);
	result += v;
	return result;
}

void TestGaussIntegr2d()
{
	typedef GaussIntegr::real_t real_t;

	static const GaussIntegr::GaussIntegrWrapper<2, 4> integr;

	real_t result(0.);

	// ���������� ���������� ��������������:

	// � �������������� �������� +=

	std::cout << integr.ByPlus_ArrArg(fun2arr<real_t>, result) << '\n';

	/*arg*/  // - ��� ���� ���������� ��������� ��� ��������������� �������
	std::cout << integr.ByPlus<real_t>(ffun2<real_t>(), result) << '\n';
	std::cout << integr.ByPlus<point2d<real_t>>(vfun2<real_t>, result) << '\n';
	std::cout << integr.ByPlus<real_t>(fun2<real_t>, result) << '\n';

	// � �������������� �������� +

	/*ret*/
	std::cout << integr.ByPlus_ArrArg<real_t>(fun2arr<real_t>) << '\n';

	/*arg,ret*/
	std::cout << integr.ByPlus<real_t, real_t>(ffun2<real_t>()) << '\n';
	std::cout << integr.ByPlus<point2d<real_t>, real_t>(vfun2<real_t>) << '\n';
	std::cout << integr.ByPlus<real_t, real_t>(fun2<real_t>) << '\n';

	// �������������� ������-�������:
	point2d<real_t> v;

	std::cout << integr.ByPlus<point2d<real_t>>(vfunv<real_t>, v) << '\n';
	std::cout << integr.ByPlus<point2d<real_t>, point2d<real_t>>(vfunv<real_t>) << '\n';
}
