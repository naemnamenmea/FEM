// Author: Chekhov Vladimir Valerevich

/* Numerical integration of an arbitrary function
   over the n-dimensional cube [-1..1] by the Gauss quadrature
   (n-dimensional array of the gauss coefficients is previously flatten to the
   one-dim array to be processed by one cycle/std.algo/unrolled cycle)
   --------------------------------------------------------------------------
    Template parameters:
    -------------------
    int ORD - order of the quadrature (must be >1)
    int DIM - space dimension (must be >0)
    typename FUNC -  integrand
    typedef real_t   -  real numbers type

    Control parameters:
    ------------------
   CYCLE_INSTEAD_OF_STD_ALGO   -   if defined then cycles are used else std.algo
   UNROLL_LOOPS_LIMIT      -  minimal array size to process by the
   cycle/std.algo (else by unrolled cycle) NESTED_CYCLE_FOR_DIM_LESS_THAN_4  -
   if defined then ordinary nested cycles are used for dimensions 1-3
*/
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
#pragma once
//---------------------------------------------------------------------------
#include <algorithm>
#include <cmath>
#include <numeric>
//---------------------------------------------------------------------------
namespace GaussIntegr
{
typedef double real_t;
//---------------------------------------------------------------------------
//#define NESTED_CYCLE_FOR_DIM_LESS_THAN_4
//#define CYCLE_INSTEAD_OF_STD_ALGO
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
//---------------Calculation of coefficients of Gauss quadrarure-------------
//---------------------------------------------------------------------------
real_t GaussX(int ord_, int num_, int newton_steps_);
real_t GaussW(int ord_, real_t x_);
//---------------------------------------------------------------------------
template <int ORD>
struct GaussCoef1D  // values of gauss coefficients for all orders - by
                    // newton-raphson's iters
{
 private:
  enum
  {
    NEWTON_ITERS_NUM = 4
  };

 public:
  template <int I_>
  static real_t x()
  {
    return GaussX(ORD, I_, NEWTON_ITERS_NUM);
  }
  template <int I_>
  static real_t w()
  {
    return GaussW(ORD, GaussX(ORD, I_, NEWTON_ITERS_NUM));
  }
  static real_t x(int i_)
  {
    return GaussX(ORD, i_, NEWTON_ITERS_NUM);
  }  // is necessary for rolled fInitGaussCoeffsArray
  static real_t w(int i_)
  {
    return GaussW(ORD, x(i_));
  }  // is necessary for rolled fInitGaussCoeffsArray
};   //---------------------
template <>
struct GaussCoef1D<1>  // specialization of explicit constant values for orders
                       // 1-5
{
  template <int>
  static real_t x()
  {
    return 0.l;
  }
  template <int>
  static real_t w()
  {
    return 2.l;
  }
  static real_t x(int) { return x<0>(); }
  static real_t w(int) { return w<0>(); }
};  //---------------------
template <>
struct GaussCoef1D<2>
{
 private:
  template <int I_>
  struct num;

 public:
  template <int I_>
  static real_t x()
  {
    return num<I_>::x();
  }
  template <int>
  static real_t w()
  {
    return 1.l;
  }
  static real_t x(int);
  static real_t w(int);
};
template <>
struct GaussCoef1D<2>::num<0>
{
  static real_t x() { return -1.l / std::sqrt(3.l); }
};
template <>
struct GaussCoef1D<2>::num<1>
{
  static real_t x() { return 1.l / std::sqrt(3.l); }
};
inline real_t GaussCoef1D<2>::x(int i_) { return i_ == 0 ? x<0>() : x<1>(); }
inline real_t GaussCoef1D<2>::w(int) { return w<0>(); }
//---------------------
template <>
struct GaussCoef1D<3>
{
 private:
  template <int I_>
  struct num;

 public:
  template <int I_>
  static real_t x()
  {
    return num<I_>::x();
  }
  template <int I_>
  static real_t w()
  {
    return num<I_>::w();
  }
  static real_t x(int);
  static real_t w(int);
};
template <>
struct GaussCoef1D<3>::num<0>
{
  static real_t x() { return -std::sqrt(0.6l); }
  static real_t w() { return 5.l / 9.l; }
};
template <>
struct GaussCoef1D<3>::num<1>
{
  static real_t x() { return 0.l; }
  static real_t w() { return 8.l / 9.l; }
};
template <>
struct GaussCoef1D<3>::num<2>
{
  static real_t x() { return std::sqrt(0.6l); }
  static real_t w() { return 5.l / 9.l; }
};
inline real_t GaussCoef1D<3>::x(int i_)
{
  return i_ == 0 ? x<0>() : i_ == 1 ? x<1>() : x<2>();
}
inline real_t GaussCoef1D<3>::w(int i_)
{
  return i_ == 0 ? w<0>() : i_ == 1 ? w<1>() : w<2>();
}
//---------------------
template <>
struct GaussCoef1D<4>
{
 private:
  template <int I_>
  struct num;

 public:
  template <int I_>
  static real_t x()
  {
    return num<I_>::x();
  }
  template <int I_>
  static real_t w()
  {
    return num<I_>::w();
  }
  static real_t x(int);
  static real_t w(int);
};
template <>
struct GaussCoef1D<4>::num<0>
{
  static real_t x() { return -std::sqrt((3.l + std::sqrt(4.8l)) / 7.l); }
  static real_t w() { return 0.5l - std::sqrt(30.l) / 36.l; }
};
template <>
struct GaussCoef1D<4>::num<1>
{
  static real_t x() { return -std::sqrt((3.l - std::sqrt(4.8l)) / 7.l); }
  static real_t w() { return 0.5l + std::sqrt(30.l) / 36.l; }
};
template <>
struct GaussCoef1D<4>::num<2>
{
  static real_t x() { return std::sqrt((3.l - std::sqrt(4.8l)) / 7.l); }
  static real_t w() { return 0.5l + std::sqrt(30.l) / 36.l; }
};
template <>
struct GaussCoef1D<4>::num<3>
{
  static real_t x() { return std::sqrt((3.l + std::sqrt(4.8l)) / 7.l); }
  static real_t w() { return 0.5l - std::sqrt(30.l) / 36.l; }
};
inline real_t GaussCoef1D<4>::x(int i_)
{
  return i_ == 0 ? x<0>() : i_ == 1 ? x<1>() : i_ == 2 ? x<2>() : x<3>();
}
inline real_t GaussCoef1D<4>::w(int i_)
{
  return i_ == 0 ? w<0>() : i_ == 1 ? w<1>() : i_ == 2 ? w<2>() : w<3>();
}
//---------------------
template <>
struct GaussCoef1D<5>
{
 private:
  template <int I_>
  struct num;

 public:
  template <int I_>
  static real_t x()
  {
    return num<I_>::x();
  }
  template <int I_>
  static real_t w()
  {
    return num<I_>::w();
  }
  static real_t x(int);
  static real_t w(int);
};
template <>
struct GaussCoef1D<5>::num<0>
{
  static real_t x() { return -std::sqrt(5.l + std::sqrt(40.l / 7.l)) / 3.l; }
  static real_t w() { return (16.1l - std::sqrt(29.575l)) / 45.l; }
};
template <>
struct GaussCoef1D<5>::num<1>
{
  static real_t x() { return -std::sqrt(5.l - std::sqrt(40.l / 7.l)) / 3.l; }
  static real_t w() { return (16.1l + std::sqrt(29.575l)) / 45.l; }
};
template <>
struct GaussCoef1D<5>::num<2>
{
  static real_t x() { return 0.l; }
  static real_t w() { return 128.l / 225.l; }
};
template <>
struct GaussCoef1D<5>::num<3>
{
  static real_t x() { return std::sqrt(5.l - std::sqrt(40.l / 7.l)) / 3.l; }
  static real_t w() { return (16.1l + std::sqrt(29.575l)) / 45.l; }
};
template <>
struct GaussCoef1D<5>::num<4>
{
  static real_t x() { return std::sqrt(5.l + std::sqrt(40.l / 7.l)) / 3.l; }
  static real_t w() { return (16.1l - std::sqrt(29.575l)) / 45.l; }
};
inline real_t GaussCoef1D<5>::x(int i_)
{
  return i_ == 0
             ? x<0>()
             : i_ == 1 ? x<1>() : i_ == 2 ? x<2>() : i_ == 3 ? x<3>() : x<4>();
}
inline real_t GaussCoef1D<5>::w(int i_)
{
  return i_ == 0
             ? w<0>()
             : i_ == 1 ? w<1>() : i_ == 2 ? w<2>() : i_ == 3 ? w<3>() : w<4>();
}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------
template <int DIM>
struct IntegrCoef  // dataStruct for one point for numerical integration of
                   // arbitrary dimension and order
{
  real_t wxyz[1 + DIM];  // contains w, x, y, z, ...
  real_t operator[](int i_) const { return wxyz[i_]; }
  real_t& operator[](int i_) { return wxyz[i_]; }
  //  template <int I> real_t  Coord()const {return wxyz[I];}
  //  template <int I> real_t& Coord()      {return wxyz[I];}
  const real_t* coord() const { return wxyz + 1; }
  real_t coord(int i_) const
  {  // 1-x, 2-y, 3-z, 4-...
    return wxyz[i_];
  }
  real_t weight() const { return wxyz[0]; }
};
//---------------------------------------------------------------------------
template <int M, int N>
struct Pow
{
  enum
  {
    result = M * Pow<M, N - 1>::result
  };
};
template <int M>
struct Pow<M, 1>
{
  enum
  {
    result = M
  };
};
//---------------------------------------------------------------------------
template <int DIM, int ORD>
struct fInitGaussCoeffsArray  // initialization of array of Gauss coefficients
                              // for any order and dimension
{
  static void Do(IntegrCoef<DIM>* coeffs_)
  {
    for (int i = 0; i < ORD; ++i, coeffs_ += Pow<ORD, DIM - 1>::result)
      fInitGaussCoeffsArray<DIM - 1, ORD>::Do(coeffs_, GaussCoef1D<ORD>::w(i),
                                              GaussCoef1D<ORD>::x(i));
  }
  template <int DIM1>  // nontemplate DIM1=DIM+1 is sufficient if DIM <= 3
  static void Do(IntegrCoef<DIM1>* coeffs_, real_t w_on_y_,
                 real_t next_coord_component_)
  {
    for (int i = 0; i < ORD; ++i)
    {
      fInitGaussCoeffsArray<DIM - 1, ORD>::Do(
          coeffs_, GaussCoef1D<ORD>::w(i) * w_on_y_, GaussCoef1D<ORD>::x(i));
      for (int j = 0; j < Pow<ORD, DIM - 1>::result; ++j, ++coeffs_)
        (*coeffs_)[DIM + 1] = next_coord_component_;
    }
  }
};
//---------------------------------------------------------------------------
template <int ORD, int I_>
struct fInitGaussCoeffsArray1D
{
  static void Do(IntegrCoef<1>* coeffs_)
  {
    fInitGaussCoeffsArray1D<ORD, I_ - 1>::Do(coeffs_);
    coeffs_[I_ - 1][0] = GaussCoef1D<ORD>::w(I_ - 1);  //<I_-1>();
    coeffs_[I_ - 1][1] = GaussCoef1D<ORD>::x(I_ - 1);  //<I_-1>();
  }
};
template <int ORD>
struct fInitGaussCoeffsArray1D<ORD, 0>
{
  static void Do(IntegrCoef<1>*) {}
};
//---------------------------------------------------------------------------
template <int ORD, int DIM, int I_ = ORD>
struct fInitGaussCoeffsSlice
{
  void operator()(IntegrCoef<DIM>* coeffs_, real_t w_on_y_, real_t y_) const
  {
    fInitGaussCoeffsSlice<ORD, DIM, I_ - 1>()(coeffs_, w_on_y_, y_);
    coeffs_[I_ - 1][0] = GaussCoef1D<ORD>::w(I_ - 1) /*<I_-1>()*/ * w_on_y_;
    coeffs_[I_ - 1][1] = GaussCoef1D<ORD>::x(I_ - 1);  //<I_-1>();
    coeffs_[I_ - 1][2] = y_;
  }
};
template <int ORD, int DIM>
struct fInitGaussCoeffsSlice<ORD, DIM, 0>
{
  void operator()(IntegrCoef<DIM>*, real_t, real_t) const {}
};
//---------------------------------------------------------------------------
template <int ORD>
struct fInitGaussCoeffsArray<1, ORD>
{
  static void Do(IntegrCoef<1>* coeffs_)
  {
    fInitGaussCoeffsArray1D<ORD, ORD>::Do(coeffs_);
  }
  template <int DIM>
  static void Do(IntegrCoef<DIM>* coeffs_, real_t w_on_y_, real_t y_)
  {
    fInitGaussCoeffsSlice<ORD, DIM, ORD>()(coeffs_, w_on_y_, y_);
  }
};

//---------------------------------------------------------------------------
//--------------------Integrating:-------------------------------------------
//---------------------------------------------------------------------------

//------- Adapter for subintegrand to convert any function one arg of type const
// real_t*
//------- due to array of gauss coefs has type real_t[]
template <int DIM, typename Targ, typename Tresult, typename FUNC>
struct fFunArgAdapter;  //    Tresult FUNC(Targ) --> Tresult FUNC(const real_t*)
//-----
template <typename Targ, typename Tresult, typename FUNC>
struct fFunArgAdapter<1, Targ, Tresult, FUNC>
{
  FUNC f;
  mutable Targ coord;
  fFunArgAdapter(FUNC f_) : f(f_), coord(0.l) {}
  Tresult operator()(const real_t* args_) const
  {
    coord[0] = args_[0];
    return f(coord);
  }
};
//-----
template <typename Tresult, typename FUNC>
struct fFunArgAdapter<1, real_t, Tresult, FUNC>
{
  FUNC f;
  fFunArgAdapter(FUNC f_) : f(f_) {}
  Tresult operator()(const real_t* args_) const { return f(args_[0]); }
};
//-----
template <typename Targ, typename Tresult, typename FUNC>
struct fFunArgAdapter<2, Targ, Tresult, FUNC>
{
  FUNC f;
  mutable Targ coord;
  fFunArgAdapter(FUNC f_) : f(f_) /*, coord(0.l)*/ {}
  Tresult operator()(const real_t* args_) const
  {
    coord[0] = args_[0];
    coord[1] = args_[1];
    return f(coord);
  }
};
//-----
template <typename Tresult, typename FUNC>
struct fFunArgAdapter<2, real_t, Tresult, FUNC>
{
  FUNC f;
  fFunArgAdapter(FUNC f_) : f(f_) {}
  Tresult operator()(const real_t* args_) const
  {
    return f(args_[0], args_[1]);
  }
};
//-----
template <typename Targ, typename Tresult, typename FUNC>
struct fFunArgAdapter<3, Targ, Tresult, FUNC>
{
  FUNC f;
  mutable Targ coord;
  fFunArgAdapter(FUNC f_) : f(f_) {}
  Tresult operator()(const real_t* args_) const
  {
    coord[0] = args_[0];
    coord[1] = args_[1];
    coord[2] = args_[2];
    return f(coord);
  }
};
//-----
template <typename Tresult, typename FUNC>
struct fFunArgAdapter<3, real_t, Tresult, FUNC>
{
  FUNC f;
  fFunArgAdapter(FUNC f_) : f(f_) {}
  Tresult operator()(const real_t* args_) const
  {
    return f(args_[0], args_[1], args_[2]);
  }
};
//------------ adapter-function to implicitly define types of template arguments
template <int DIM, typename Targ, typename Tresult, typename FUNC>
fFunArgAdapter<DIM, Targ, Tresult, FUNC> AdaptArgsToArray(FUNC f_)
{
  return fFunArgAdapter<DIM, Targ, Tresult, FUNC>(f_);
}
//---------------------------------------------------------------------
//-----------functional objects to be used in algorithms---------------
//---------------------------------------------------------------------
template <typename FUNC>
struct fweighted_plus
{
  FUNC f;
  fweighted_plus(FUNC f_) : f(f_) {}

  template <typename T, int DIM>
  T operator()(const T& prev_sum_, const IntegrCoef<DIM>& coef_) const
  {
    return prev_sum_ + f(coef_.coord()) * coef_.weight();
  }
};

template <typename FUNC>
inline fweighted_plus<FUNC> weighted_plus(FUNC f_)
{
  return fweighted_plus<FUNC>(f_);
}
//---------------------------------------------------------------------
template <typename T, typename FUNC>
struct fweighted_plus_with_assgn
{
  fweighted_plus_with_assgn(FUNC f_, T& result0_) : f(f_), result(result0_) {}

  template <int DIM>
  void operator()(const IntegrCoef<DIM>& coef_) const
  //     {  result  +=  f(coef_.coord()) * coef_.weight(); } // чомусь не робе
  //     пiд g++
  {
    T funret = f(coef_.coord());
    funret *= coef_.weight();
    result += funret;
  }

  FUNC f;
  T& result;
};

template <typename T, typename FUNC>
inline fweighted_plus_with_assgn<T, FUNC> weighted_plus_with_assgn(FUNC f_,
                                                                   T& result0_)
{
  return fweighted_plus_with_assgn<T, FUNC>(f_, result0_);
}
//---------------------------------------------------------------------
//-----------Unrolled algorithms---------------------------------------
//---------------------------------------------------------------------
#ifndef UNROLL_LOOPS_LIMIT
#define UNROLL_LOOPS_LIMIT 10
#endif
//---------------------------------------------------------------------
template <typename T, int I_, typename Tret, typename FUNC>
struct fAccumulate_unrolled
{
  static Tret Do(const T* BEGIN_, const Tret& prev_sum_, FUNC fun_)
  {
    return fun_(fAccumulate_unrolled<T, I_ - 1, Tret, FUNC>::Do(
                    BEGIN_, prev_sum_, fun_),
                *(BEGIN_ + I_ - 1));
  }
};
//--------------
template <typename T, typename Tret, typename FUNC>
struct fAccumulate_unrolled<T, 0, Tret, FUNC>
{
  static Tret Do(const T* BEGIN_, const Tret& prev_sum_, FUNC fun_)
  {
    return prev_sum_;
  }
};
//--------------
template <int I_, typename Tret, typename T, typename FUNC>
inline Tret Accumulate_unrolled(const T* BEGIN_, const Tret& prev_sum_,
                                FUNC fun_)
{
  return fAccumulate_unrolled<T, I_, Tret, FUNC>::Do(BEGIN_, prev_sum_, fun_);
}
//---------------------------------------------------------------------
template <bool ARRAY_IS_BIG>
struct fRunAccumulate;
//-------
template <>
struct fRunAccumulate<true>
{
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
  template <typename TArray, typename Tret, typename FUNC>
  static Tret Do(const TArray& arr_, const Tret& sum0_, FUNC f_)
  {
    return std::accumulate(arr_.begin(), arr_.end(), sum0_, f_);
  }
#else
  template <typename TArray, typename Tret, typename FUNC, typename TPtr>
  static Tret Do(const TArray& arr_, const Tret& sum0_, FUNC f_)
  {
    for (TPtr p = arr_.begin(); p < arr_.end(); ++p) sum0_ = sum0_ + f_(*p);
    return sum0_;
  }
#endif
};
//-------
template <>
struct fRunAccumulate<false>
{
  template <typename TArray, typename Tret, typename FUNC>
  static Tret Do(const TArray& arr_, const Tret& sum0_, FUNC f_)
  {
    return Accumulate_unrolled<TArray::size, Tret>(arr_.begin(), sum0_, f_);
  }
};
//---------------------------------------------------------------------
template <typename T, int I_, typename FUNC>
struct fForEach_unrolled
{
  static void Do(const T* BEGIN_, FUNC fun_)
  {
    fForEach_unrolled<T, I_ - 1, FUNC>::Do(BEGIN_, fun_);
    fun_(*(BEGIN_ + I_ - 1));
  }
};
//--------------
template <typename T, typename FUNC>
struct fForEach_unrolled<T, 0, FUNC>
{
  static void Do(const T*, FUNC) {}
};
//--------------
template <int I_, typename T, typename FUNC>
inline void ForEach_unrolled(const T* BEGIN_, FUNC fun_)
{
  fForEach_unrolled<T, I_, FUNC>::Do(BEGIN_, fun_);
}
//---------------------------------------------------------------------
template <bool ARRAY_IS_BIG>
struct fRunForEach;
//-------
template <>
struct fRunForEach<true>
{
  template <typename TArray, typename FUNC>
  static void Do(const TArray& arr_, FUNC f_)
  {
    std::for_each(arr_.begin(), arr_.end(), f_);
  }
};  //-------
template <>
struct fRunForEach<false>
{
  template <typename TArray, typename FUNC>
  static void Do(const TArray& arr_, FUNC f_)
  {
    ForEach_unrolled<TArray::size>(arr_.begin(), f_);
  }
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------
template <int DIM, int ORD>
struct fIntegrate  // 1D-array of Gauss coefs for any dimension/order with
                   // integration
{
  enum
  {
    size = Pow<ORD, DIM>::result
  };

  /*static*/ IntegrCoef<DIM> Coeffs[size];

  fIntegrate() { fInitGaussCoeffsArray<DIM, ORD>::Do(Coeffs); }
  const IntegrCoef<DIM>* begin() const { return Coeffs; }
  const IntegrCoef<DIM>* end() const { return Coeffs + size; }
  //-------
  template <typename Tret, typename FUNC>
  Tret ByPlus_ArrArg(FUNC f_) const
  {
    return fRunAccumulate<(size >= UNROLL_LOOPS_LIMIT)>::Do(*this, Tret(0.l),
                                                            weighted_plus(f_));
  }  //-------
  template <typename Tret, typename FUNC>
  Tret& ByPlusAssgn_ArrArg(FUNC f_, Tret& result_) const
  {
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
    fRunForEach<(size >= UNROLL_LOOPS_LIMIT)>::Do(
        *this, weighted_plus_with_assgn(f_, result_));
#else
    Tret tmp;
    for (const IntegrCoef<DIM>* p = begin(); p < end(); ++p)
      result_ += (tmp = f_(p->coord())) *= p->weight();
#endif
    return result_;
  }  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret ByPlus(FUNC f_) const
  {
    return fRunAccumulate<(size >= UNROLL_LOOPS_LIMIT)>::Do(
        *this, Tret(),
        /*Tret(0.l), */ weighted_plus(AdaptArgsToArray<DIM, Targ, Tret>(f_)));
  }  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret& ByPlusAssgn(FUNC f_, Tret& result_) const
  {
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
    fRunForEach<(size >= UNROLL_LOOPS_LIMIT)>::Do(
        *this, weighted_plus_with_assgn(AdaptArgsToArray<DIM, Targ, Tret>(f_),
                                        result_));
#else
    fFunArgAdapter<DIM, Targ, Tret, FUNC> f(f_);
    Tret f_from_coord;
    for (const IntegrCoef<DIM>* p = begin(); p < end(); ++p)
      result_ += ((f_from_coord = f(p->coord())) *= p->weight());
#endif
    return result_;
  }  //-------
};
//---------------------------------------------------------------------
//---------------------------------------------------------------------
#ifdef NESTED_CYCLE_FOR_DIM_LESS_THAN_4
template <int ORD>
struct fIntegrate<1, ORD>
{
  /*static*/ IntegrCoef<1> Coeffs[ORD];

  fIntegrate() { fInitGaussCoeffsArray<1, ORD>::Do(Coeffs); }
  //-------
  /* template <typename Targ, typename Tret, typename FUNC>
   Tret ByPlus_ArrArg(FUNC f_) const
     {
     }//-------
  */
  template <typename Targ, typename Tret, typename FUNC>
  Tret ByPlus(FUNC f_) const
  {
    Tret result(0.l);
    for (int i = 0; i < ORD; ++i)
      result = result + f_(Coeffs[i].coord(1)) * Coeffs[i].weight();
    return result;
  }  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret& ByPlusAssgn(FUNC f_, Tret& result_) const
  {
    static Targ arg(0.);
    static Tret tmp;
    for (int i = 0; i < ORD; ++i)
    //       result_ += f_(Coeffs[i].coord(1))*Coeffs[i].weight();
    {
      arg[0] = Coeffs[i].coord(1);
      (tmp = f_(arg)) *= Coeffs[i].weight();
      result_ += tmp;
    }
    return result_;
  }
};
//---------------------------------------------------------------------
template <int ORD>
struct fIntegrate<2, ORD>
{
  /*static*/ IntegrCoef<1> Coeffs[ORD];

  fIntegrate() { fInitGaussCoeffsArray<1, ORD>::Do(Coeffs); }
  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret ByPlus(FUNC f_) const
  {
    Tret result(0.l);
    //    AdaptArgsToArray<2,Targ,Tret> f(f_);  // - error in g++: expected ';'
    //    before 'f'
    for (int i = 0, j; i < ORD; ++i)
      for (j = 0; j < ORD; ++j)
        //       result = result +
        //       f(Coeffs[i].coord(1),Coeffs[j].coord(1))*Coeffs[i].weight()*Coeffs[j].weight();
        result = result + AdaptArgsToArray<2, Targ, Tret>(f_)(
                              Coeffs[i].coord(1), Coeffs[j].coord(1)) *
                              Coeffs[i].weight() * Coeffs[j].weight();
    return result;
  }  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret& ByPlusAssgn(FUNC f_, Tret& result_) const
  {
    static Targ arg(0.);
    //    static Tret tmp;
    for (int i = 0, j; i < ORD; ++i)
    {
      arg[0] = Coeffs[i].coord(1);
      for (j = 0; j < ORD; ++j)
      {
        arg[1] = Coeffs[j].coord(1);
        /*       tmp = f_(arg);
               tmp *= (Coeffs[i].weight()*Coeffs[j].weight());   // = error
               result_ += tmp;*/
        result_ += (f_(arg) *= (Coeffs[i].weight() * Coeffs[j].weight()));
      }
    }
    return result_;
  }  //-------
};
//---------------------------------------------------------------------
template <int ORD>
struct fIntegrate<3, ORD>
{
  /*static*/ IntegrCoef<1> Coeffs[ORD];

  fIntegrate() { fInitGaussCoeffsArray<1, ORD>::Do(Coeffs); }
  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret ByPlus(FUNC f_) const
  {
    Tret result(0.l);
    //    AdaptArgsToArray<3,Targ,Tret> f(f_);
    for (int i = 0, j, k; i < ORD; ++i)
      for (j = 0; j < ORD; ++j)
        for (k = 0; k < ORD; ++k)
          //       result = result +
          //       f(Coeffs[i].coord(1),Coeffs[j].coord(2),Coeffs[k].coord(3))*Coeffs[i].weight()*Coeffs[j].weight()*Coeffs[k].weight();
          result =
              result +
              AdaptArgsToArray<3, Targ, Tret>(f_)(
                  Coeffs[i].coord(1), Coeffs[j].coord(1), Coeffs[k].coord(1)) *
                  Coeffs[i].weight() * Coeffs[j].weight() * Coeffs[k].weight();
    return result;
  }  //-------
  template <typename Targ, typename Tret, typename FUNC>
  Tret& ByPlusAssgn(FUNC f_, Tret& result_) const
  {
    //    AdaptArgsToArray<3,Targ,Tret> f(f_);  // - error in g++: expected ';'
    //    before 'f'
    static Targ arg;
    for (int i = 0, j, k; i < ORD; ++i)
    {
      arg[0] = Coeffs[i].coord(1);
      for (j = 0; j < ORD; ++j)
      {
        arg[1] = Coeffs[j].coord(1);
        for (k = 0; k < ORD; ++k)
        {
          //       result_ +=
          //       f(Coeffs[i].coord(1),Coeffs[j].coord(2),Coeffs[k].coord(3))*Coeffs[i].weight()*Coeffs[j].weight()*Coeffs[k].weight();
          //       result_ +=
          //       AdaptArgsToArray<3,Targ,Tret>(f_)(Coeffs[i].coord(1),Coeffs[j].coord(1),Coeffs[k].coord(1))*(Coeffs[i].weight()*Coeffs[j].weight()*Coeffs[k].weight());
          arg[2] = Coeffs[k].coord(1);
          result_ += (f_(arg) *= (Coeffs[i].weight() * Coeffs[j].weight() *
                                  Coeffs[k].weight()));
        }
      }
    }
    return result_;
  }  //-------
};
//---------------------------------------------------------------------------
#endif  // NESTED_CYCLE_FOR_DIM_LESS_THAN_4
}  // namespace GaussIntegr
