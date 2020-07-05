#pragma once

#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif

#include <cmath>
#include <cstddef>	//nullptr

#if defined __BORLANDC__ && __BORLANDC__ >= 0x550
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::ceil;
using std::cos;
using std::cosh;
using std::exp;
using std::fabs;
using std::floor;
using std::fmod;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sin;
using std::sinh;
using std::sqrt;
using std::tan;
using std::tanh;
#endif	//__BORLANDC__

// Settings for compilation

#ifndef GAUSS_ORDER_1D
#define GAUSS_ORDER_1D 3
#endif
#ifndef GAUSS_ORDER_2D
#define GAUSS_ORDER_2D 3
#endif
#ifndef GAUSS_ORDER_3D
#define GAUSS_ORDER_3D 3
#endif
#ifndef UNROLL_LOOPS_LIMIT
#define UNROLL_LOOPS_LIMIT 0
#endif

// Aliases for numeric types

typedef double real_t;
typedef unsigned short int small_t;

//       real_t abs  (real_t);
inline real_t sign(real_t x_)
{
	return x_ < 0. ? -1. : 1.;
}
inline real_t sign(real_t y_, real_t x_)
{
	return x_ < 0. ? -y_ : y_;
}

// functionals for for_each

template <typename T>
struct address_is
{
	address_is(const T& ref_) : m_address(&ref_)
	{
	}
	address_is(const T* ptr_) : m_address(ptr_)
	{
	}
	bool operator()(const T& ref_) const
	{
		return m_address == &ref_;
	}
	// bool operator()(const T* ptr_) const {return address == ptr_;}

	const T* const m_address;
};

template <typename T>
struct fDelete_object
{
	void operator()(T*& p_) const
	{
		delete p_;
		p_ = nullptr;
	}
};

template <typename T1, typename T2>
struct fMilt
{
	fMilt(T2 v_) : m_value(v_)
	{
	}
	void operator()(T1& o_) const
	{
		o_ *= m_value;
	}

	const T2 m_value;
};
