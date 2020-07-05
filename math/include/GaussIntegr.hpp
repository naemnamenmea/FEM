/* Numerical integration of an arbitrary function
over the n-dimensional cube [-1..1] by the Gauss quadrature
(n-dimensional array of the gauss coefficients is previously flatten to the one-dim array to be
processed by one cycle/std.algo/unrolled cycle)

Author: Vladimir V. Chekhov (http://vch.crimea.ua)
This code is described in: the paper submitted to ZNU Herald
(http://web.znu.edu.ua/cms/index.php?action=category/browse&site_id=5&lang=eng&category_id=1146)
Version of Sept 14, 2013

		Template parameters:
		-
		int ORD - order of the quadrature (must be >1)
		int DIM - space dimension (must be >0)
		typename FUN - integrand

		Control parameters:

CYCLE_INSTEAD_OF_STD_ALGO - if defined then cycles are used else std.algo
UNROLL_LOOPS_LIMIT - minimal array size to process by the cycle/std.algo (else by unrolled cycle)
NESTED_CYCLE_FOR_DIM_LESS_THAN_4 - if defined then ordinary nested cycles are used for dimensions
1-3
*/

#pragma once

#include "NumTypes.hpp"
#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif
#include <algorithm>
#include <numeric>
#include <cmath>

//#define NESTED_CYCLE_FOR_DIM_LESS_THAN_4
//#define CYCLE_INSTEAD_OF_STD_ALGO

// Calculation of coefficients of Gauss quadrarure

real_t GaussX(int ord_, int num_, int newton_steps_);
real_t GaussW(int ord_, real_t x_);

template <int ORD>
struct GaussCoef1D	// values of gauss coefficients for all orders - by newton-raphson's iters
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
};
template <>
struct GaussCoef1D<1>  // specialization of explicit constant values for orders 1-5
{
	template <int>
	static real_t x()
	{
		return 0.;
	}
	template <int>
	static real_t w()
	{
		return 2.;
	}
	static real_t x(int)
	{
		return x<0>();
	}
	static real_t w(int)
	{
		return w<0>();
	}
};
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
		return 1.;
	}
	static real_t x(int);
	static real_t w(int);
};
template <>
struct GaussCoef1D<2>::num<0>
{
	static real_t x()
	{
		return -1. / std::sqrt(3.);
	}
};
template <>
struct GaussCoef1D<2>::num<1>
{
	static real_t x()
	{
		return 1. / std::sqrt(3.);
	}
};
inline real_t GaussCoef1D<2>::x(int i_)
{
	return i_ == 0 ? x<0>() : x<1>();
}
inline real_t GaussCoef1D<2>::w(int)
{
	return w<0>();
}

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
	static real_t x()
	{
		return -std::sqrt(0.6);
	}
	static real_t w()
	{
		return 5. / 9.;
	}
};
template <>
struct GaussCoef1D<3>::num<1>
{
	static real_t x()
	{
		return 0.;
	}
	static real_t w()
	{
		return 8. / 9.;
	}
};
template <>
struct GaussCoef1D<3>::num<2>
{
	static real_t x()
	{
		return std::sqrt(0.6);
	}
	static real_t w()
	{
		return 5. / 9.;
	}
};
inline real_t GaussCoef1D<3>::x(int i_)
{
	return i_ == 0 ? x<0>() : i_ == 1 ? x<1>() : x<2>();
}
inline real_t GaussCoef1D<3>::w(int i_)
{
	return i_ == 0 ? w<0>() : i_ == 1 ? w<1>() : w<2>();
}

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
	static real_t x()
	{
		return -std::sqrt((3. + std::sqrt(4.8)) / 7.);
	}
	static real_t w()
	{
		return 0.5 - std::sqrt(30.) / 36.;
	}
};
template <>
struct GaussCoef1D<4>::num<1>
{
	static real_t x()
	{
		return -std::sqrt((3. - std::sqrt(4.8)) / 7.);
	}
	static real_t w()
	{
		return 0.5 + std::sqrt(30.) / 36.;
	}
};
template <>
struct GaussCoef1D<4>::num<2>
{
	static real_t x()
	{
		return std::sqrt((3. - std::sqrt(4.8)) / 7.);
	}
	static real_t w()
	{
		return 0.5 + std::sqrt(30.) / 36.;
	}
};
template <>
struct GaussCoef1D<4>::num<3>
{
	static real_t x()
	{
		return std::sqrt((3. + std::sqrt(4.8)) / 7.);
	}
	static real_t w()
	{
		return 0.5 - std::sqrt(30.) / 36.;
	}
};
inline real_t GaussCoef1D<4>::x(int i_)
{
	return i_ == 0 ? x<0>() : i_ == 1 ? x<1>() : i_ == 2 ? x<2>() : x<3>();
}
inline real_t GaussCoef1D<4>::w(int i_)
{
	return i_ == 0 ? w<0>() : i_ == 1 ? w<1>() : i_ == 2 ? w<2>() : w<3>();
}

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
	static real_t x()
	{
		return -std::sqrt(5. + std::sqrt(40. / 7.)) / 3.;
	}
	static real_t w()
	{
		return (16.1 - std::sqrt(29.575)) / 45.;
	}
};
template <>
struct GaussCoef1D<5>::num<1>
{
	static real_t x()
	{
		return -std::sqrt(5. - std::sqrt(40. / 7.)) / 3.;
	}
	static real_t w()
	{
		return (16.1 + std::sqrt(29.575)) / 45.;
	}
};
template <>
struct GaussCoef1D<5>::num<2>
{
	static real_t x()
	{
		return 0.;
	}
	static real_t w()
	{
		return 128. / 225.;
	}
};
template <>
struct GaussCoef1D<5>::num<3>
{
	static real_t x()
	{
		return std::sqrt(5. - std::sqrt(40. / 7.)) / 3.;
	}
	static real_t w()
	{
		return (16.1 + std::sqrt(29.575)) / 45.;
	}
};
template <>
struct GaussCoef1D<5>::num<4>
{
	static real_t x()
	{
		return std::sqrt(5. + std::sqrt(40. / 7.)) / 3.;
	}
	static real_t w()
	{
		return (16.1 - std::sqrt(29.575)) / 45.;
	}
};
inline real_t GaussCoef1D<5>::x(int i_)
{
	return i_ == 0 ? x<0>() : i_ == 1 ? x<1>() : i_ == 2 ? x<2>() : i_ == 3 ? x<3>() : x<4>();
}
inline real_t GaussCoef1D<5>::w(int i_)
{
	return i_ == 0 ? w<0>() : i_ == 1 ? w<1>() : i_ == 2 ? w<2>() : i_ == 3 ? w<3>() : w<4>();
}

template <int DIM>
struct IntegrCoef  // dataStruct for one point for numerical integration of arbitrary dimension and
				   // order
{
	real_t wxyz[1 + DIM];  // contains w, x, y, z, ...
	real_t operator[](int i_) const
	{
		return wxyz[i_];
	}
	real_t& operator[](int i_)
	{
		return wxyz[i_];
	}
	//  template <int I> real_t  Coord()const {return wxyz[I];}
	//  template <int I> real_t& Coord()      {return wxyz[I];}
	const real_t* coord() const
	{
		return wxyz + 1;
	}
	real_t coord(int i_) const
	{  // 1-x, 2-y, 3-z, 4-...
#ifdef STRONGCHECK
		Assert(
			i_ > 0 && i_ <= DIM,
			tMessage("invalid coordinate index=")
				<< i_ << "must be>0 and <=" << DIM << " in IntegrCoef<DIM>::coord(int): ");
#endif
		return wxyz[i_];
	}
	real_t weight() const
	{
		return wxyz[0];
	}
};

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

/*template <int DIM, int ORD, int I_=ORD>
struct fSetAllCoeff // развёртка цикла
{
// void operator()(IntegrCoef<DIM>* coeffs_) const
 static void Do(IntegrCoef<DIM>* coeffs_)
	{//fSetAllCoeff<DIM,ORD,I_-1>()(coeffs_);
	 fSetAllCoeff<DIM,ORD,I_-1>::Do(coeffs_);
	 fInitGaussCoeffsArray<DIM-1,ORD>::Do(coeffs_+(I_-1)*Pow<ORD,DIM-1>::result,GaussCoef1D<ORD>::w<I_-1>(),GaussCoef1D<ORD>::x<I_-1>());
//
fInitGaussCoeffsArray<DIM-1,ORD>::Do(coeffs_+=(I_-1)*Pow<ORD,DIM-1>::result,GaussCoef1D<ORD>::w<I_-1>(),GaussCoef1D<ORD>::x<I_-1>());
	}
};
template <int DIM, int ORD>  struct fSetAllCoeff<DIM,ORD,0>
//{ void operator()(IntegrCoef<DIM>*) const {} };
{ static void Do(IntegrCoef<DIM>*) {} };

template <int DIM, int DIM1, int ORD, int I_=ORD>
struct fSetAllCoeff2 // развёртка цикла
{
 void operator()(IntegrCoef<DIM1>* coeffs_, real_t w_on_y_, real_t next_coord_component_) const
	{fSetAllCoeff2<DIM,DIM1,ORD,I_-1>()(coeffs_,w_on_y_,next_coord_component_);
	 for (IntegrCoef<DIM1>* p=coeffs_; p<coeffs_+Pow<ORD,DIM-1>::result; ++p)
			(*p)[DIM+1] = next_coord_component_;
	 fInitGaussCoeffsArray<DIM-1,ORD>::Do(coeffs_+=(I_-1)*Pow<ORD,DIM-1>::result,GaussCoef1D<ORD>::w<I_-1>()*w_on_y_,GaussCoef1D<ORD>::x<I_-1>());
	}
};
template <int DIM, int DIM1, int ORD>  struct fSetAllCoeff2<DIM,DIM1,ORD,0>
{ void operator()(IntegrCoef<DIM1>*,real_t,real_t) const {} };*/

template <int DIM, int ORD>
struct fInitGaussCoeffsArray  // initialization of array of Gauss coefficients for any order and
							  // dimension
{
	static void Do(IntegrCoef<DIM>* coeffs_)
	{
		//     fSetAllCoeff<DIM,ORD,ORD>::Do(coeffs_);
		for (small_t i = 0; i < ORD; ++i, coeffs_ += Pow<ORD, DIM - 1>::result)
			fInitGaussCoeffsArray<DIM - 1, ORD>::Do(
				coeffs_, GaussCoef1D<ORD>::w(i), GaussCoef1D<ORD>::x(i));
	}
	template <int DIM1>	 // nontemplate DIM1=DIM+1 is sufficient if DIM <= 3
	static void Do(IntegrCoef<DIM1>* coeffs_, real_t w_on_y_, real_t next_coord_component_)
	{
		//     fSetAllCoeff2<DIM,DIM1,ORD,ORD>()(coeffs_,w_on_y_,next_coord_component_);
		for (small_t i = 0; i < ORD; ++i)
		{
			fInitGaussCoeffsArray<DIM - 1, ORD>::Do(
				coeffs_, GaussCoef1D<ORD>::w(i) * w_on_y_, GaussCoef1D<ORD>::x(i));
			for (small_t j = 0; j < Pow<ORD, DIM - 1>::result; ++j, ++coeffs_)
				(*coeffs_)[DIM + 1] = next_coord_component_;
		}
	}
};

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
	static void Do(IntegrCoef<1>*)
	{
	}
};

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
	void operator()(IntegrCoef<DIM>*, real_t, real_t) const
	{
	}
};

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

// Integrating:
// Adapter for subintegrand to convert any function one arg of type const real_t*
// due to array of gauss coefs has type real_t[]

template <int DIM, typename Targ, typename Tresult, typename FUN>
struct fFunArgAdapter;	//    Tresult FUN(Targ) --> Tresult FUN(const real_t*)

template <typename Targ, typename Tresult, typename FUN>
struct fFunArgAdapter<1, Targ, Tresult, FUN>
{
	FUN f;
	mutable Targ coord;
	fFunArgAdapter(FUN f_) : f(f_), coord(0.)
	{
	}
	Tresult operator()(const real_t* args_) const
	{
		coord[0] = args_[0];
		return f(coord);
	}
};

template <typename Tresult, typename FUN>
struct fFunArgAdapter<1, real_t, Tresult, FUN>
{
	FUN f;
	fFunArgAdapter(FUN f_) : f(f_)
	{
	}
	Tresult operator()(const real_t* args_) const
	{
		return f(args_[0]);
	}
};

template <typename Targ, typename Tresult, typename FUN>
struct fFunArgAdapter<2, Targ, Tresult, FUN>
{
	FUN f;
	mutable Targ coord;
	fFunArgAdapter(FUN f_) : f(f_), coord(0.)
	{
	}
	Tresult operator()(const real_t* args_) const
	{
		coord[0] = args_[0];
		coord[1] = args_[1];
		return f(coord);
	}
};

template <typename Targ, typename Tresult, typename FUN>
struct fFunArgAdapter<3, Targ, Tresult, FUN>
{
	FUN f;
	mutable Targ coord;
	fFunArgAdapter(FUN f_) : f(f_)
	{
	}
	Tresult operator()(const real_t* args_) const
	{
		coord[0] = args_[0];
		coord[1] = args_[1];
		coord[2] = args_[2];
		return f(coord);
	}
};

template <typename Tresult, typename FUN>
struct fFunArgAdapter<2, real_t, Tresult, FUN>
{
	FUN f;
	fFunArgAdapter(FUN f_) : f(f_)
	{
	}
	Tresult operator()(const real_t* args_) const
	{
		return f(args_[0], args_[1]);
	}
};
// adapter-function to implicitly define types of template arguments
template <int DIM, typename Targ, typename Tresult, typename FUN>
fFunArgAdapter<DIM, Targ, Tresult, FUN> AdaptArgsToArray(FUN f_)
{
	return fFunArgAdapter<DIM, Targ, Tresult, FUN>(f_);
}

// functional objects to be used in algorithms

template <typename FUN>
struct fweighted_plus
{
	FUN f;
	fweighted_plus(FUN f_) : f(f_)
	{
	}

	template <typename T, int DIM>
	T operator()(const T& prev_sum_, const IntegrCoef<DIM>& coef_) const
	{
		return prev_sum_ + f(coef_.coord()) * coef_.weight();
	}
};
template <typename FUN>
inline fweighted_plus<FUN> weighted_plus(FUN f_)
{
	return fweighted_plus<FUN>(f_);
}

template <typename T, typename FUN>
struct fweighted_plus_with_assgn
{
	FUN f;
	T& result;
	fweighted_plus_with_assgn(FUN f_, T& result0_) : f(f_), result(result0_)
	{
	}
	template <int DIM>
	void operator()(const IntegrCoef<DIM>& coef_) const
	//     {  result  +=  f(coef_.coord()) * coef_.weight(); } // чомусь не робе пiд g++
	{
		T funret = f(coef_.coord());
		funret *= coef_.weight();
		result += funret;
	}
};
template <typename T, typename FUN>
inline fweighted_plus_with_assgn<T, FUN> weighted_plus_with_assgn(FUN f_, T& result0_)
{
	return fweighted_plus_with_assgn<T, FUN>(f_, result0_);
}

// Unrolled algorithms

#ifndef UNROLL_LOOPS_LIMIT
#define UNROLL_LOOPS_LIMIT 10
#endif

template <typename T, int I_, typename Tret, typename FUN>
struct fAccumulate_unrolled
{
	static Tret Do(const T* BEGIN_, const Tret& prev_sum_, FUN fun_)
	{
		return fun_(
			fAccumulate_unrolled<T, I_ - 1, Tret, FUN>::Do(BEGIN_, prev_sum_, fun_),
			*(BEGIN_ + I_ - 1));
	}
};
template <typename T, typename Tret, typename FUN>
struct fAccumulate_unrolled<T, 0, Tret, FUN>
{
	static Tret Do(const T* BEGIN_, const Tret& prev_sum_, FUN fun_)
	{
		return prev_sum_;
	}
};
template <int I_, typename Tret, typename T, typename FUN>
inline Tret Accumulate_unrolled(const T* BEGIN_, const Tret& prev_sum_, FUN fun_)
{
	return fAccumulate_unrolled<T, I_, Tret, FUN>::Do(BEGIN_, prev_sum_, fun_);
}

template <bool ARRAY_IS_BIG>
struct fRunAccumulate;

template <>
struct fRunAccumulate<true>
{
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
	template <typename TArray, typename Tret, typename FUN>
	static Tret Do(const TArray& arr_, const Tret& sum0_, FUN f_)
	{
		return std::accumulate(arr_.begin(), arr_.end(), sum0_, f_);
	}
#else
	template <typename TArray, typename Tret, typename FUN, typename TPtr>
	static Tret Do(const TArray& arr_, const Tret& sum0_, FUN f_)
	{
		for (TPtr p = arr_.begin(); p < arr_.end(); ++p) sum0_ = sum0_ + f_(*p);
		return sum0_;
	}
#endif
};
template <>
struct fRunAccumulate<false>
{
	template <typename TArray, typename Tret, typename FUN>
	static Tret Do(const TArray& arr_, const Tret& sum0_, FUN f_)
	{
		return Accumulate_unrolled<TArray::size, Tret>(arr_.begin(), sum0_, f_);
	}
};

template <typename T, int I_, typename FUN>
struct fForEach_unrolled
{
	static void Do(const T* BEGIN_, FUN fun_)
	{
		fForEach_unrolled<T, I_ - 1, FUN>::Do(BEGIN_, fun_);
		fun_(*(BEGIN_ + I_ - 1));
	}
};
template <typename T, typename FUN>
struct fForEach_unrolled<T, 0, FUN>
{
	static void Do(const T*, FUN)
	{
	}
};

template <int I_, typename T, typename FUN>
inline void ForEach_unrolled(const T* BEGIN_, FUN fun_)
{
	fForEach_unrolled<T, I_, FUN>::Do(BEGIN_, fun_);
}

template <bool ARRAY_IS_BIG>
struct fRunForEach;

template <>
struct fRunForEach<true>
{
	template <typename TArray, typename FUN>
	static void Do(const TArray& arr_, FUN f_)
	{
		std::for_each(arr_.begin(), arr_.end(), f_);
	}
};
template <>
struct fRunForEach<false>
{
	template <typename TArray, typename FUN>
	static void Do(const TArray& arr_, FUN f_)
	{
		ForEach_unrolled<TArray::size>(arr_.begin(), f_);
	}
};

template <int DIM, int ORD>
struct fIntegrate  // 1D-array of Gauss coefs for any dimension/order with integration
{
	enum
	{
		size = Pow<ORD, DIM>::result
	};

	/*static*/ IntegrCoef<DIM> Coeffs[size];

	fIntegrate()
	{
		fInitGaussCoeffsArray<DIM, ORD>::Do(Coeffs);
	}
	const IntegrCoef<DIM>* begin() const
	{
		return Coeffs;
	}
	const IntegrCoef<DIM>* end() const
	{
		return Coeffs + size;
	}
	template <typename Tret, typename FUN>
	Tret ByPlus_ArrArg(FUN f_) const
	{
		return fRunAccumulate<(size >= UNROLL_LOOPS_LIMIT)>::Do(*this, Tret(0.), weighted_plus(f_));
	}
	template <typename Tret, typename FUN>
	Tret& ByPlusAssgn_ArrArg(FUN f_, Tret& result_) const
	{
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
		fRunForEach<(size >= UNROLL_LOOPS_LIMIT)>::Do(*this, weighted_plus_with_assgn(f_, result_));
#else
		for (const IntegrCoef<DIM>* p = begin(); p < end(); ++p)
			result_ += f_(p->coord()) *= p->weight();
#endif
		//    for (const IntegrCoef<DIM>* p=cfbegin_; p<cfend_; ++p)
		//       result_ +=  f_(p->coord()) *= p->weight();
		//   std::for_each(begin(),end(),weighted_plus_with_assgn(f_,result_));
		//   ForEach_unrolled<size>(begin(),weighted_plus_with_assgn(f_,result_));
		return result_;
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret ByPlus(FUN f_) const
	{
		return fRunAccumulate<(size >= UNROLL_LOOPS_LIMIT)>::Do(
			*this, Tret(0.), weighted_plus(AdaptArgsToArray<DIM, Targ, Tret>(f_)));
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret& ByPlusAssgn(FUN f_, Tret& result_) const
	{
#ifndef CYCLE_INSTEAD_OF_STD_ALGO
		fRunForEach<(size >= UNROLL_LOOPS_LIMIT)>::Do(
			*this, weighted_plus_with_assgn(AdaptArgsToArray<DIM, Targ, Tret>(f_), result_));
#else
		fFunArgAdapter<DIM, Targ, Tret, FUN> f(f_);
		Tret f_from_coord;
		for (const IntegrCoef<DIM>* p = begin(); p < end(); ++p)
			result_ += ((f_from_coord = f(p->coord())) *= p->weight());
#endif
		return result_;
	}
};

#ifdef NESTED_CYCLE_FOR_DIM_LESS_THAN_4
template <int ORD>
struct fIntegrate<1, ORD>
{
	/*static*/ IntegrCoef<1> Coeffs[ORD];

	fIntegrate()
	{
		fInitGaussCoeffsArray<1, ORD>::Do(Coeffs);
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret ByPlus(FUN f_) const
	{
		Tret result(0.);
		for (int i = 0; i < ORD; ++i) result = result + f_(Coeffs[i].coord(1)) * Coeffs[i].weight();
		return result;
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret& ByPlusAssgn(FUN f_, Tret& result_) const
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

template <int ORD>
struct fIntegrate<2, ORD>
{
	/*static*/ IntegrCoef<1> Coeffs[ORD];

	fIntegrate()
	{
		fInitGaussCoeffsArray<1, ORD>::Do(Coeffs);
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret ByPlus(FUN f_) const
	{
		Tret result(0.);
		//    AdaptArgsToArray<2,Targ,Tret> f(f_);  // - error in g++: expected ';' before 'f'
		for (int i = 0, j; i < ORD; ++i)
			for (j = 0; j < ORD; ++j)
				//       result = result +
				//       f(Coeffs[i].coord(1),Coeffs[j].coord(1))*Coeffs[i].weight()*Coeffs[j].weight();
				result = result + AdaptArgsToArray<2, Targ, Tret>(f_)(
									  Coeffs[i].coord(1), Coeffs[j].coord(1)) *
									  Coeffs[i].weight() * Coeffs[j].weight();
		return result;
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret& ByPlusAssgn(FUN f_, Tret& result_) const
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
	}
};

/*template <int ORD>// sample from mooney_fem
struct fIntegr2D
{
 IntegrCoef<1> coeffs[ORD];
 fIntegr2D()  { fInitGaussCoeffsArray<1,ORD>::Do(coeffs); }

 template <typename Tres, typename FUN>
 Tres& ByPlusAssgn(FUN f_, Tres& result_)
  {
   for (int i=0, j; i<ORD; ++i)
	for (j=0; j<ORD; ++j)
	 {
	  result_ += (f_(coeffs[i].coord(1),coeffs[j].coord(1)) * (coeffs[i].weight() *
coeffs[j].weight()));
	 }
   return result_;
 }
};*/

template <int ORD>
struct fIntegrate<3, ORD>
{
	/*static*/ IntegrCoef<1> Coeffs[ORD];

	fIntegrate()
	{
		fInitGaussCoeffsArray<1, ORD>::Do(Coeffs);
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret ByPlus(FUN f_) const
	{
		Tret result(0.);
		//    AdaptArgsToArray<3,Targ,Tret> f(f_);
		for (int i = 0, j, k; i < ORD; ++i)
			for (j = 0; j < ORD; ++j)
				for (k = 0; k < ORD; ++k)
					//       result = result +
					//       f(Coeffs[i].coord(1),Coeffs[j].coord(2),Coeffs[k].coord(3))*Coeffs[i].weight()*Coeffs[j].weight()*Coeffs[k].weight();
					result =
						result + AdaptArgsToArray<3, Targ, Tret>(f_)(
									 Coeffs[i].coord(1), Coeffs[j].coord(1), Coeffs[k].coord(1)) *
									 Coeffs[i].weight() * Coeffs[j].weight() * Coeffs[k].weight();
		return result;
	}
	template <typename Targ, typename Tret, typename FUN>
	Tret& ByPlusAssgn(FUN f_, Tret& result_) const
	{
		//    AdaptArgsToArray<3,Targ,Tret> f(f_);  // - error in g++: expected ';' before 'f'
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
					result_ +=
						(f_(arg) *= (Coeffs[i].weight() * Coeffs[j].weight() * Coeffs[k].weight()));
				}
			}
		}
		return result_;
	}
};

#endif	// NESTED_CYCLE_FOR_DIM_LESS_THAN_4
//
// Set of 1D-arrays of Gauss coefs for any dimension and all possible
// orders (dynamic switching) with integration
//
// template <int DIM, int MAXORD>
// struct fIntegrateAllOrds_core: public fIntegrateAllOrds_core<DIM,MAXORD-1>, public
// fIntegrate<DIM,MAXORD>
//{
// int Size(int ord_) const
//   {return pow(float(ord_),DIM);}
// const IntegrCoef<DIM>* begin(int ord_) const
//   {return ord_ == MAXORD ? fIntegrate<DIM,MAXORD>::begin() :
//   fIntegrateAllOrds_core<DIM,MAXORD-1>::begin(ord_);}
// const IntegrCoef<DIM>* end(int ord_) const
//   {return begin(ord_) + Size(ord_);}
//
// template <typename Targ, typename Tret, typename FUN>
// Tret ByPlus(int ord_, FUN f_) const
//   {
//    #ifdef STRONGCHECK
//     Assert(ord_>0 && ord_<=MAXORD, "integration order exceeds max allowable value in
//     fIntegrateAllOrds_core<DIM,MAXORD>::IntByPlus(int,FUN)");
//    #endif
//    return ord_ == MAXORD  ?  fIntegrate<DIM,MAXORD>::ByPlus<Targ/*gcc error: expected
//    primary-expression before ","*/,Tret/*gcc error: expected primary-expression before ">"*/>(f_)
//                           :  fIntegrateAllOrds_core<DIM,MAXORD-1>::ByPlus<Targ/*gcc error:
//                           expected primary-expression before ","*/,Tret/*gcc error: expected
//                           primary-expression before ">"*/>(ord_,f_);
//   }
//
// template <typename Targ, typename Tret, typename FUN>
// Tret& ByPlusAssgn(int ord_, FUN f_, Tret& result_) const
//   {
//    #ifdef STRONGCHECK
//     Assert(ord_>0 && ord_<=MAXORD, "integration order exceeds max allowable value in
//     fIntegrateAllOrds_core<DIM,MAXORD>::IntByPlusAssgn(int,FUN,Tret&)");
//    #endif
//    return ord_ == MAXORD  ?  fIntegrate<DIM,MAXORD>::ByPlusAssgn<Targ/*gcc error: expected
//    primary-expression before ">"*/>(f_,result_)
//                           :  fIntegrateAllOrds_core<DIM,MAXORD-1>::ByPlusAssgn<Targ/*gcc error:
//                           expected primary-expression before ">"*/>(ord_,f_,result_);
//   }
//};
//
// template <int DIM>
// struct fIntegrateAllOrds_core<DIM,1>: public fIntegrate<DIM,1>
//{
// int Size(int) const {return 1;}
// const IntegrCoef<DIM>* begin(int) const
//   {return fIntegrate<DIM,1>::begin();}
// const IntegrCoef<DIM>* end(int) const
//   {return fIntegrate<DIM,1>::begin() + 1;}
//
// template <typename Targ, typename Tret, typename FUN>
// Tret ByPlus(int, FUN f_) const
//   {
//    return fIntegrate<DIM,1>::ByPlus<Targ,Tret>(f_); // gcc errors: expected primary-expression
//    before ",",  expected primary-expression before ">"
//   }
// template <typename Targ, typename Tret, typename FUN>
// Tret& ByPlusAssgn(int, FUN f_, Tret& result_) const
//   {
//    return fIntegrate<DIM,1>::ByPlusAssgn<Targ>(f_,result_); // gcc error: expected
//    primary-expression before ">"
//   }
//};
//
// template <int DIM, int MAXORD>
// struct fIntegrateAllOrds: public fIntegrateAllOrds_core<DIM,MAXORD>
//{
//// small_t CurrentOrder;
// const IntegrCoef<DIM> *CurrentBegin, *CurrentEnd;
// void SetOrder(small_t ord_)
//    {
//     #ifdef STRONGCHECK
//      Assert(ord_>0 && ord_<=MAXORD,"Improper gauss order in fIntegrateAllOrds::SetOrder");
//     #endif
////     CurrentOrder = ord_;
//     CurrentBegin = fIntegrateAllOrds_core<DIM,MAXORD>::begin(ord_);
//     CurrentEnd = fIntegrateAllOrds_core<DIM,MAXORD>::end(ord_);
//    }
// fIntegrateAllOrds(int ord_=MAXORD) {SetOrder(ord_);}
//
// template <typename Targ, typename Tret, typename FUN>
// Tret ByPlus_ArrArg_algo(FUN f_) const
//   {
//    return std::accumulate(CurrentBegin,CurrentEnd,Tret(0.),weighted_plus(f_));
//   }
// template <typename Targ, typename Tret, typename FUN>
// Tret ByPlus_algo(FUN f_) const
//   {
//    return
//    std::accumulate(CurrentBegin,CurrentEnd,Tret(0.),weighted_plus(AdaptArgsToArray<DIM,Targ,Tret>(f_)));
//   }
// template <typename Targ, typename Tret, typename FUN>
// Tret& ByPlusAssgn_ArrArg_algo(FUN f_, Tret& result_) const
//   {
//    std::for_each(CurrentBegin,CurrentEnd,weighted_plus_with_assgn(f_,result_));
//    return result_;
//   }
// template <typename Targ, typename Tret, typename FUN>
// Tret& ByPlusAssgn_algo(FUN f_, Tret& result_) const
//   {
//    std::for_each(CurrentBegin,CurrentEnd,weighted_plus_with_assgn(AdaptArgsToArray<DIM,Targ,Tret>(f_),result_));
//    return result_;
//   }
//};
