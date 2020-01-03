// Author: Chekhov Vladimir Valerevich

#include <iostream>
#include "GaussIntegr.hpp"
#include "point2d.hpp"
#include "tests.hpp"

// различные примеры подынтегральных выражений дл€ трЄхмерной задачи:

template <typename T>
T fun3arr(const T* xy)
{
  return xy[0] * xy[1];
}

template <typename T>
T fun3(T x, T y)
{
  return x * y;
}

template <typename T>
struct ffun3
{  // функциональный объект
  T operator()(T x, T y) const { return x * y; }
};

template <typename T>
T vfun3(point2d<T> v)
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

//--------------------------------------------------------------------------

void TestBasic2d()
{
  typedef GaussIntegr::real_t real_t;

  // создание объекта, выполн€ющего интегрирование по √ауссу дл€ трЄхмерного
  // случа€ с пор€дком 4:
  static const GaussIntegr::fIntegrate<2, 4> integr;

  real_t result(0.);

  // выполнение численного интегрировани€:

  // с использованием операции +=

  std::cout << integr.ByPlusAssgn_ArrArg(fun3arr<real_t>, result) << '\n';

  /*arg*/  // - тип чего приходитс€ указывать при инстанцировании объекта
  std::cout << integr.ByPlusAssgn<real_t>(ffun3<real_t>(), result) << '\n';
  std::cout << integr.ByPlusAssgn<point2d<real_t>>(vfun3<real_t>, result)
            << '\n';

  // с использованием операции +

  /*ret*/
  std::cout << integr.ByPlus_ArrArg<real_t>(fun3arr<real_t>) << '\n';

  /*arg,ret*/
  std::cout << integr.ByPlus<real_t, real_t>(ffun3<real_t>()) << '\n';
  std::cout << integr.ByPlus<point2d<real_t>, real_t>(vfun3<real_t>) << '\n';

  // интегрирование вектор-функций:
  point2d<real_t> v;
  integr.ByPlusAssgn<point2d<real_t>>(vfunv<real_t>, v);
  std::cout << v << '\n';

  v = integr.ByPlus<point2d<real_t>, point2d<real_t>>(vfunv<real_t>);
  std::cout << v << '\n';
}
//--------------------------------------------------------------------------
