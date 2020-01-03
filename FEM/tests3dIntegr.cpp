// Author: Chekhov Vladimir Valerevich

#include <iostream>
#include "GaussIntegr.hpp"
#include "point3d.hpp"
#include "tests.hpp"

// ��������� ������� ��������������� ��������� ��� ��������� ������:

template <typename T>
T fun3arr(const T* xyz)
{
  return xyz[0] * xyz[1] * xyz[2];
}

template <typename T>
T fun3(T x, T y, T z)
{
  return x * y * z;
}

template <typename T>
struct ffun3
{  // �������������� ������
  T operator()(T x, T y, T z) const { return x * y * z; }
};

template <typename T>
T vfun3(point3d<T> v)
{
  return v[0] * v[1] * v[2];
}

template <typename T>
point3d<T> vfunv(const point3d<T>& v)
{
  point3d<T> result(v);
  result += v;
  return result;
}

//--------------------------------------------------------------------------

void TestBasic3d()
{
  typedef GaussIntegr::real_t real_t;

  // �������� �������, ������������ �������������� �� ������ ��� ����������
  // ������ � �������� 4:
  static const GaussIntegr::fIntegrate<3, 4> integr;

  real_t result(0.);

  // ���������� ���������� ��������������:

  // � �������������� �������� +=

  std::cout << integr.ByPlusAssgn_ArrArg(fun3arr<real_t>, result) << '\n';

  /*arg*/  // - ��� ���� ���������� ��������� ��� ��������������� �������
  std::cout << integr.ByPlusAssgn<real_t>(ffun3<real_t>(), result) << '\n';
  std::cout << integr.ByPlusAssgn<point3d<real_t>>(vfun3<real_t>, result)
            << '\n';

  // � �������������� �������� +

  /*ret*/
  std::cout << integr.ByPlus_ArrArg<real_t>(fun3arr<real_t>) << '\n';

  /*arg,ret*/
  std::cout << integr.ByPlus<real_t, real_t>(ffun3<real_t>()) << '\n';
  std::cout << integr.ByPlus<point3d<real_t>, real_t>(vfun3<real_t>) << '\n';

  // �������������� ������-�������:
  point3d<real_t> v;
  integr.ByPlusAssgn<point3d<real_t>>(vfunv<real_t>, v);
  std::cout << v << '\n';

  v = integr.ByPlus<point3d<real_t>, point3d<real_t>>(vfunv<real_t>);
  std::cout << v << '\n';
}
//--------------------------------------------------------------------------