#pragma once

#include "NumTypes.hpp"
#include "math_constants.hpp"
#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif

#include <cstddef>	//nullptr
#include <algorithm>
#include <functional>
#include <fstream>
#include <cmath>

enum component_t
{
	X = 1,
	Y,
	Z,
	R,
	Fi,
	Ro,
	Psi,
	XX,
	YY,
	ZZ,
	XY,
	YX,
	YZ,
	ZY,
	ZX,
	XZ
};
enum tensor_kind_t
{
	eZero,
	eUnit
};
class Tensor1a;
class Tensor2s;
class Tensor2a;
class Tensor3s;
class Tensor4s;

class Tensor1s
{
public:
	const real_t& operator[](small_t i_) const
	{
		return m_data[i_];
	}  // is public because need for fFunArgAdapter<1,Targ,Tresult,FUN> for integration
	real_t& operator[](small_t i_)
	{
		return m_data[i_];
	}
	// public:
	Tensor1s() : m_data()
	{
	}
	Tensor1s(tensor_kind_t, small_t, real_t = 1.);	// Basis vector multiplied by factor_
	Tensor1s(float x_, float y_, float z_)
	{
		m_data[0] = x_;
		m_data[1] = y_;
		m_data[2] = z_;
	}
	Tensor1s(double x_, double y_, double z_)
	{
		m_data[0] = x_;
		m_data[1] = y_;
		m_data[2] = z_;
	}
	Tensor1s(long double x_, long double y_, long double z_)
	{
		m_data[0] = static_cast<real_t>(x_);
		m_data[1] = static_cast<real_t>(y_);
		m_data[2] = static_cast<real_t>(z_);
	}
	explicit Tensor1s(real_t v_)
	{
		__Assign(v_);
	}
	explicit Tensor1s(const real_t* const);
	Tensor1s(const Tensor1s& o_)
	{
		__Assign(o_);
	}
	explicit Tensor1s(const Tensor1a&);
	const real_t& operator()(small_t) const;
	real_t& operator()(small_t);
	Tensor1s& operator=(const Tensor1a& o_);
	Tensor1s& operator=(const Tensor1s& o_)
	{
		__Assign(o_);
		return *this;
	}
	Tensor1s& operator-=(const Tensor1s& o_)
	{
		m_data[0] -= o_[0];
		m_data[1] -= o_[1];
		m_data[2] -= o_[2];
		return *this;
	}
	Tensor1s& operator+=(const Tensor1s& o_)
	{
		m_data[0] += o_[0];
		m_data[1] += o_[1];
		m_data[2] += o_[2];
		return *this;
	}
	Tensor1s& operator*=(real_t v_)
	{
		m_data[0] *= v_;
		m_data[1] *= v_;
		m_data[2] *= v_;
		return *this;
	}
	Tensor1s& operator/=(real_t);
	Tensor1s operator-() const
	{
		return Tensor1s(-m_data[0], -m_data[1], -m_data[2]);
	}
	Tensor1s operator+(const Tensor1s& o_) const
	{
		return Tensor1s(m_data[0] + o_[0], m_data[1] + o_[1], m_data[2] + o_[2]);
	}
	Tensor1s operator-(const Tensor1s& o_) const
	{
		return Tensor1s(m_data[0] - o_[0], m_data[1] - o_[1], m_data[2] - o_[2]);
	}
	Tensor1s operator*(real_t v_) const
	{
		return Tensor1s(m_data[0] * v_, m_data[1] * v_, m_data[2] * v_);
	}
	bool operator==(const Tensor1s& o_) const
	{
		return mathdef::is_eq(m_data[0], o_[0]) && mathdef::is_eq(m_data[1], o_[1]) &&
			   mathdef::is_eq(m_data[2], o_[2]);
	}  // is used in "DxfStream.cpp"
	Tensor1s& Negate()
	{
		m_data[0] = -m_data[0];
		m_data[1] = -m_data[1];
		m_data[2] = -m_data[2];
		return *this;
	}
	real_t Length() const
	{
		return sqrt(m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2]);
	}
	Tensor1s& Normalize()
	{
		return operator/=(Length());
	}
	Tensor1s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor1s& Assign1(small_t i_, real_t v_ = 1.)
	{
		__Assign(0.);
		m_data[--i_] = v_;
		return *this;
	}
	Tensor1s& Assign(real_t x_, real_t y_, real_t z_)
	{
		m_data[0] = x_;
		m_data[1] = y_;
		m_data[2] = z_;
		return *this;
	}
	Tensor1s& Assign(small_t i_, real_t v_)
	{
		m_data[--i_] = v_;
		return *this;
	}  // temporarily - to be removed
	   //   Tensor1& Assign(component_t, real_t);
	   //   Tensor1& Add(component_t, real_t);
	Tensor1s& CrossMultiplyByBasisVector(small_t);
	Tensor1s& CrossMultiplyByBasisVector_left(small_t);
	Tensor1s& CrossProductWithBasisVector(small_t, Tensor1s&) const;
	Tensor1s& CrossProductWithBasisVector_left(small_t, Tensor1s&) const;
	real_t DotProduct(const Tensor1s&) const;
	Tensor1s& CrossMultiply(const Tensor1s&);
	Tensor1s& CrossMultiply_left(const Tensor1s&);
	Tensor1s& CrossProduct(const Tensor1s&, Tensor1s&) const;
	Tensor2s& DirectProduct(const Tensor1s&, Tensor2s&) const;
	Tensor2a& DirectProduct(const Tensor1s&, Tensor2a&) const;
	Tensor1s& DotMultiply(const Tensor2s&);
	Tensor1s& DotMultiply_left(const Tensor2s&);
	Tensor1s& DotProduct(const Tensor2s&, Tensor1s&) const;
	Tensor2s& CrossProduct(const Tensor2s&, Tensor2s&) const;
	Tensor3s& DirectProduct(const Tensor2s&, Tensor3s&) const;
	Tensor2s& DotProduct(const Tensor3s&, Tensor2s&) const;
	Tensor3s& CrossProduct(const Tensor3s&, Tensor3s&) const;
	Tensor4s& DirectProduct(const Tensor3s&, Tensor4s&) const;
	Tensor3s& DotProduct(const Tensor4s&, Tensor3s&) const;
	Tensor4s& CrossProduct(const Tensor4s&, Tensor4s&) const;
	// friend std::ostream& operator<< (std::ostream& out_, const Tensor1s& t_);

private:
	friend class Tensor1a;
	friend class Tensor2s;
	friend class Tensor3s;
	friend class Tensor4s;

	void __Assign(real_t v_)
	{
		m_data[0] = m_data[1] = m_data[2] = v_;
	}
	void __Assign(const real_t* const d_)
	{
		m_data[0] = d_[0];
		m_data[1] = d_[1];
		m_data[2] = d_[2];
	}
	void __Assign(const Tensor1s& o_)
	{
		m_data[0] = o_.m_data[0];
		m_data[1] = o_.m_data[1];
		m_data[2] = o_.m_data[2];
	}

	real_t m_data[3];
};

class RetTensor1;
class RetTensor2;

class Tensor2a;

class Tensor1a
{
public:
	//////////////////////// temporarily!!! - only for tNode::LinkDisplacement
	const Tensor1s* Data() const
	{
		return m_pData;
	}
	////////////////////////
	Tensor1a() : m_pData(nullptr)
	{
	}
	Tensor1a(real_t x_, real_t y_, real_t z_) : m_pData(new Tensor1s(x_, y_, z_))
	{
	}
	explicit Tensor1a(const Tensor1s& o_) : m_pData(new Tensor1s(o_))
	{
	}
	Tensor1a(const Tensor1a& o_)
		: m_pData(o_.is0() ? (Tensor1s*)nullptr : new Tensor1s(*o_.m_pData))
	{
	}
	Tensor1a(RetTensor1&);
	virtual ~Tensor1a()
	{
		delete m_pData;
	}
	bool is0() const
	{
		return m_pData == nullptr;
	}
	bool is0(small_t) const
	{
		return is0();
	}
	const real_t& operator()(small_t) const;
	real_t& operator()(small_t);
	Tensor1a& operator=(RetTensor1&);
	Tensor1a& operator=(const Tensor1a& o_)
	{
		if (__AdjustAdd(o_))
			m_pData->operator=(*o_.m_pData);
		return *this;
	}
	Tensor1a& operator-=(const Tensor1a&);
	Tensor1a& operator+=(const Tensor1a& o_)
	{
		if (__AdjustAdd(o_))
			m_pData->operator+=(*o_.m_pData);
		return *this;
	}
	Tensor1a& operator*=(real_t v_)
	{
		return is0() ? *this : (m_pData->operator*=(v_), *this);
	}
	Tensor1a& operator/=(real_t);
	RetTensor1 operator-() const;
	RetTensor1 operator-(const Tensor1a&) const;
	RetTensor1 operator+(const Tensor1a&) const;
	RetTensor1 operator*(real_t) const;
	RetTensor1 operator/(real_t) const;
	Tensor1a& Negate()
	{
		return is0() ? *this : (m_pData->Negate(), *this);
	}
	real_t Length() const
	{
		return is0() ? 0. : m_pData->Length();
	}
	Tensor1a& Normalize();
	Tensor1a& Assign0()
	{
		delete m_pData;
		m_pData = nullptr;
		return *this;
	}
	Tensor1a& Assign1(small_t i_, real_t v_ = 1.)
	{
		if (is0())
			m_pData = new Tensor1s(true, i_, v_);
		else
			m_pData->Assign1(i_, v_);
		return *this;
	}
	Tensor1a& Assign(small_t i_, real_t v_)
	{
		if (is0())
			m_pData = new Tensor1s(0.);
		m_pData->m_data[--i_] = v_;
		return *this;
	}
	Tensor1a& Assign(real_t x_, real_t y_, real_t z_)
	{
		is0() ? m_pData = new Tensor1s(x_, y_, z_) : &m_pData->Assign(x_, y_, z_);
		return *this;
	}
	Tensor1a& Add(small_t i_, real_t v_)
	{
		if (is0())
			m_pData = new Tensor1s(0.);
		m_pData->m_data[--i_] += v_;
		return *this;
	}
	Tensor1a& CrossMultiplyByBasisVector(small_t i_)
	{
		return is0() ? *this : (m_pData->CrossMultiplyByBasisVector(i_), *this);
	}
	Tensor1a& CrossMultiplyByBasisVector_left(small_t i_)
	{
		return is0() ? *this : (m_pData->CrossMultiplyByBasisVector_left(i_), *this);
	}
	Tensor1a& CrossProductWithBasisVector(small_t i_, Tensor1a& r_) const
	{
		return (r_ = *this).CrossMultiplyByBasisVector(i_);
	}
	Tensor1a& CrossProductWithBasisVector_left(small_t i_, Tensor1a& r_) const
	{
		return (r_ = *this).CrossMultiplyByBasisVector_left(i_);
	}
	RetTensor1 CrossProductWithBasisVector(small_t) const;
	RetTensor1 CrossProductWithBasisVector_left(small_t) const;
	real_t DotProduct(const Tensor1a& o_) const
	{
		return is0() || o_.is0() ? 0. : m_pData->DotProduct(*o_.m_pData);
	}
	Tensor1a& CrossMultiply(const Tensor1a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->CrossMultiply(*o_.m_pData);
		return *this;
	}
	Tensor1a& CrossMultiply_left(const Tensor1a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->CrossMultiply_left(*o_.m_pData);
		return *this;
	}
	Tensor1a& CrossProduct(const Tensor1a& o_, Tensor1a& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).CrossMultiply(o_);
	}
	RetTensor1 CrossProduct(const Tensor1a&) const;
	Tensor2a& DirectProduct(const Tensor1a&, Tensor2a&) const;
	RetTensor2 DirectProduct(const Tensor1a&) const;
	Tensor1a& DotMultiply(const Tensor2a&);
	Tensor1a& DotMultiply_left(const Tensor2a&);
	Tensor1a& DotProduct(const Tensor2a&, Tensor1a&) const;
	RetTensor1 DotProduct(const Tensor2a&) const;
	Tensor2a& CrossProduct(const Tensor2a&, Tensor2a&) const;
	RetTensor2 CrossProduct(const Tensor2a&) const;
	Tensor3s& DirectProduct(const Tensor2a&, Tensor3s&) const;
	Tensor2a& DotProduct(const Tensor3s&, Tensor2a&) const;
	RetTensor2 DotProduct(const Tensor3s&) const;
	Tensor3s& CrossProduct(const Tensor3s&, Tensor3s&) const;
	Tensor4s& DirectProduct(const Tensor3s&, Tensor4s&) const;
	Tensor3s& DotProduct(const Tensor4s&, Tensor3s&) const;
	Tensor4s& CrossProduct(const Tensor4s&, Tensor4s&) const;

protected:
	Tensor1s* m_pData;

private:
	friend class Tensor1s;
	friend class RetTensor1;
	friend class Tensor2a;
	// friend std::ostream& operator<< (std::ostream& out_, const Tensor1a& t_);

	Tensor1s* __pData()
	{
		return is0() ? (m_pData = new Tensor1s) : m_pData;
	}
	bool __AdjustAdd(const Tensor1a& o_)
	{
		return is0() ? (o_.is0() ? false : (m_pData = new Tensor1s(*o_.m_pData), false))
					 : (o_.is0() ? false : true);
	}
	template <typename T>
	bool __AdjustMult(const T& o_)
	{
		return is0() ? false : (o_.is0() ? Assign0(), false : true);
	}
	template <typename Tfactor, typename Tresult>
	bool __AdjustMult(const Tfactor& o_, Tresult& r_) const
	{
		return is0() || o_.is0() ? r_.Assign0(), false : true;
	}
	const real_t operator[](small_t i_) const
	{
		return is0() ? 0. : m_pData->operator[](i_);
	}
};
// inline std::ostream& operator<< (std::ostream& out_, const Tensor1s& t_) {return out_<<"|
// "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" |\n";} inline std::ostream& operator<< (std::ostream& out_,
// const Tensor1a& t_) {return out_<<"| "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" |\n";}

class RetTensor1
{
public:
	~RetTensor1()
	{
		delete m_pData;
	}

private:
	friend class Tensor1a;
	friend class Tensor2a;

	RetTensor1() : m_pData(nullptr)
	{
	}
	RetTensor1(Tensor1a& o_) : m_pData(o_.m_pData)
	{
		o_.m_pData = nullptr;
	}
	//   RetTensor1(const RetTensor1& o_): pData(o_.pData) {o_.pData = nullptr;}
	//   explicit RetTensor1(const real_t* const p_): pData(p_==nullptr ? nullptr : new
	//   Tensor1s(p_)) {}
	RetTensor1(real_t x_, real_t y_, real_t z_) : m_pData(new Tensor1s(x_, y_, z_))
	{
	}

	mutable Tensor1s* m_pData;
};

// std::ostream& operator<< (std::ostream&, const BasicTensor2&);

// class Tensor4;
class SymmetricTensor4s;
class SymmetricTensor2s;

class Tensor2s
{
private:
#define __SET_ROWS          \
	*m_rows = m_data;       \
	m_rows[1] = m_data + 3; \
	m_rows[2] = m_data + 6;

public:
	Tensor2s() : m_data()
	{
		__SET_ROWS;
	}
	Tensor2s(::tensor_kind_t, real_t v_ = 0.)
	{
		__SET_ROWS;
		Assign1(v_);
	}
	explicit Tensor2s(real_t v_)
	{
		__SET_ROWS;
		__Assign(v_);
	}
	explicit Tensor2s(const real_t (*v_)[3])
	{
		__SET_ROWS;
		__Assign2D(v_);
	}
	explicit Tensor2s(const real_t** v_)
	{
		__SET_ROWS;
		__Assign2D(v_);
	}
	Tensor2s(const real_t (*v_)[3], bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(v_);
	}
	Tensor2s(const real_t** v_, bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(v_);
	}
	Tensor2s(const Tensor1s&, const Tensor1s&, const Tensor1s&);
	Tensor2s(const Tensor1s&, const Tensor1s&, const Tensor1s&, bool);
	Tensor2s(const Tensor2s& o_)
	{
		__SET_ROWS;
		__Assign(o_.m_data);
	}
	Tensor2s(const Tensor2s& o_, bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(o_.m_rows);
	}
	const real_t& operator()(small_t, small_t) const;
	real_t& operator()(small_t, small_t);
	Tensor2s& operator=(const Tensor2s&);
	Tensor2s& operator-=(const Tensor2s&);
	Tensor2s& operator+=(const Tensor2s&);
	Tensor2s& operator*=(real_t v_)
	{
		std::for_each(m_data, m_data + 9, fMilt<real_t, real_t>(v_));
		return *this;
	}
	Tensor2s operator*(real_t v_) const
	{
		return Tensor2s(*this) *= v_;
	}  // tmp to remove
	Tensor2s& operator/=(real_t);
	Tensor2s& Negate()
	{
		std::transform(m_data, m_data + 9, m_data, std::negate<real_t>());
		return *this;
	}
	Tensor2s& Transpose();
	Tensor2s& AddTransposed(const Tensor2s&);
	SymmetricTensor2s& Symmetrization_twofold(SymmetricTensor2s&) const;
	Tensor2s& AssignRows(const Tensor1s&, const Tensor1s&, const Tensor1s&);
	Tensor2s& AssignCols(const Tensor1s&, const Tensor1s&, const Tensor1s&);
	Tensor2s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor2s& Assign1(real_t v_ = 1.)
	{
		__Assign(0.);
		m_data[0] = m_data[4] = m_data[8] = v_;
		return *this;
	}
	real_t Contraction() const
	{
		return m_data[0] + m_data[4] + m_data[8];
	}
	real_t Det() const;
	real_t I1() /*Linear    invariant*/ const
	{
		return Contraction();
	}
	real_t I2() /*Quadratic invariant*/ const;
	real_t I3() /*Cubic     invariant*/ const
	{
		return Det();
	}
	real_t Minor(small_t, small_t) const;
	Tensor2s& Invert();
	Tensor1s& AttachedVector(Tensor1s&) const;
	Tensor1s AttachedVector() const;
	//   Tensor1s& RotationVector(Tensor1s&) const;//for orthogonal tensors only
	Tensor1s& ScalarProductWithBasisVector(small_t, Tensor1s&) const;
	Tensor1s& ScalarProductWithBasisVector_left(small_t, Tensor1s&) const;
	Tensor1s ScalarProductWithBasisVector(small_t) const;
	Tensor1s ScalarProductWithBasisVector_left(small_t) const;
	void EigenValues(real_t[3]) const;
	small_t EigenVectors(
		real_t /*current e.value*/, Tensor2s& /*e.vects*/, small_t = 1 /*begin position*/) const;
	void Eigen(real_t[3], Tensor2s&) const;
	Tensor2s& EigenVectors(Tensor2s& r_) const
	{
		real_t evs[3];
		Eigen(evs, r_);
		return r_;
	}
	Tensor2s& RotateToEigenSystem(Tensor2s&);
	Tensor2s& AssignDev();
	Tensor2s& Dev(Tensor2s& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	Tensor1s& DotProduct(const Tensor1s&, Tensor1s&) const;
	Tensor2s& CrossMultiply(const Tensor1s&);
	Tensor2s& CrossMultiply_left(const Tensor1s&);
	Tensor2s& CrossProduct(const Tensor1s&, Tensor2s&) const;
	Tensor3s& DirectProduct(const Tensor1s&, Tensor3s&) const;
	real_t Dot2Product(const Tensor2s&) const;
	Tensor2s& DotMultiply(const Tensor2s&);
	Tensor2s& DotMultiply_left(const Tensor2s&);
	Tensor2s& DotProduct(const Tensor2s&, Tensor2s&) const;
	Tensor3s& CrossProduct(const Tensor2s&, Tensor3s&) const;
	Tensor4s& DirectProduct(const Tensor2s&, Tensor4s&) const;
	Tensor1s& Dot2Product(const Tensor3s&, Tensor1s&) const;
	Tensor3s& DotProduct(const Tensor3s&, Tensor3s&) const;
	Tensor4s& CrossProduct(const Tensor3s&, Tensor4s&) const;
	//   Tensor2s& Dot2Multiply      (const Tensor4s&,small_t,small_t);
	//   Tensor2s& Dot2Multiply      (const Tensor4&,small_t,small_t);
	Tensor2s& Dot2Multiply(const SymmetricTensor4s&, small_t, small_t);
	Tensor2s& Dot2Multiply(const Tensor4s&);
	Tensor2s& Dot2Multiply_left(const Tensor4s&);
	Tensor2s& Dot2Product(const Tensor4s&, Tensor2s&) const;
	Tensor4s& DotProduct(const Tensor4s&, Tensor4s&) const;
	//   Tensor2s& ApplyFunction(const fClassicFunction&);

private:
	friend class Tensor1s;
	friend class Tensor2a;
	friend class Tensor3s;
	friend class Tensor4s;

	struct tIndex
	{
	public:
		tIndex() : m_first(), m_second()
		{
		}
		tIndex(small_t i_, small_t j_) : m_first(i_), m_second(j_)
		{
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);
		//        Index&   operator()(small_t i_, small_t j_){first=i_; second=j_; return *this;}

	private:
		small_t m_first, m_second;
	};

	//   void __set_rows(){const_cast<real_t*>(*Rows) = Data; const_cast<real_t*>(Rows[1]) = Data+3;
	//   const_cast<real_t*>(Rows[2]) = Data+6;}
	const real_t* operator[](small_t i_) const
	{
		return m_rows[i_];
	}
	real_t* operator[](small_t i_)
	{
		return m_rows[i_];
	}
	void __Assign(real_t v_)
	{
		std::fill(m_data, m_data + 9, v_);
	}
	void __Assign(const real_t* v_)
	{
		std::copy(v_, v_ + 9, m_data);
	}
	template <typename T>
	void __Assign2D(const T v_)
	{
		for (small_t i = 0, j; i < 3; ++i)
			for (j = 0; j < 3; ++j) m_rows[i][j] = v_[i][j];
	}
	template <typename T>
	void __Assign2DTranspose(const T v_)
	{
		for (small_t i = 0, j; i < 3; ++i)
			for (j = 0; j < 3; ++j) m_rows[i][j] = v_[j][i];
	}
	const real_t& operator[](tIndex) const;
	real_t& operator[](tIndex);
	const real_t& operator()(tIndex) const;
	real_t& operator()(tIndex);
	void MaxAbsIndex(small_t&, small_t&) const;

	real_t m_data[9];
	real_t* /*const*/ m_rows[3];
};

class RetTensor2;
class Tensor2CaptureDataAdapter;

class Tensor2a
{
public:
	Tensor2a() : m_pData(nullptr)
	{
	}
	explicit Tensor2a(const Tensor2s& o_) : m_pData(new Tensor2s(o_))
	{
	}
	Tensor2a(const Tensor2a& o_)
		: m_pData(o_.is0() ? (Tensor2s*)nullptr : new Tensor2s(*o_.m_pData))
	{
	}
	Tensor2a(RetTensor2&);
	virtual ~Tensor2a()
	{
		delete m_pData;
	}
	bool is0() const
	{
		return m_pData == nullptr;
	}
	bool is0(small_t, small_t) const
	{
		return m_pData == nullptr;
	}
	const real_t& operator()(small_t, small_t) const;
	real_t& operator()(small_t, small_t);
	Tensor2a& operator=(RetTensor2&);
	Tensor2a& operator=(const Tensor2s& o_)
	{
		if (is0())
			m_pData = new Tensor2s(o_);
		else
			m_pData->operator=(o_);
		return *this;
	}
	Tensor2a& operator=(const Tensor2a& o_)
	{
		if (__AdjustAdd(o_))
			m_pData->operator=(*o_.m_pData);
		return *this;
	}
	Tensor2a& operator-=(const Tensor2a&);
	Tensor2a& operator+=(const Tensor2a& o_)
	{
		if (__AdjustAdd(o_))
			m_pData->operator+=(*o_.m_pData);
		return *this;
	}
	Tensor2a& operator+=(Tensor2CaptureDataAdapter);
	Tensor2a& operator*=(real_t v_)
	{
		return is0() ? *this : (m_pData->operator*=(v_), *this);
	}
	Tensor2a& operator/=(real_t);
	RetTensor2 operator-() const;
	RetTensor2 operator-(const Tensor2a&) const;
	RetTensor2 operator+(const Tensor2a&) const;
	RetTensor2 operator*(real_t) const;
	RetTensor2 operator/(real_t) const;
	real_t Contraction() const
	{
		return is0() ? 0. : m_pData->Contraction();
	}
	real_t Det() const
	{
		return is0() ? 0. : m_pData->Det();
	}
	real_t I1() const
	{
		return Contraction();
	}
	real_t I2() const
	{
		return is0() ? 0. : m_pData->I2();
	}
	real_t I3() const
	{
		return Det();
	}
	real_t Minor(small_t i_, small_t j_) const
	{
		return is0() ? 0. : m_pData->Minor(i_, j_);
	}
	Tensor2a& Negate()
	{
		return is0() ? *this : (m_pData->Negate(), *this);
	}
	Tensor2a& Transpose()
	{
		return is0() ? *this : (m_pData->Transpose(), *this);
	}
	Tensor2a& AddTransposed(const Tensor2a& o_)
	{
		return is0() ? (o_.is0() ? *this : (m_pData = new Tensor2s(*o_.m_pData, false), *this))
					 : (o_.is0() ? *this : (m_pData->AddTransposed(*o_.m_pData), *this));
	}
	Tensor2a& AddTransposed(Tensor2CaptureDataAdapter);
	Tensor2a& Invert();
	Tensor2a& Assign0()
	{
		delete m_pData;
		m_pData = nullptr;
		return *this;
	}
	Tensor2a& Assign0toData()
	{
		m_pData->Assign0();
		return *this;
	}
	Tensor2a& Assign1(real_t v_ = 1.)
	{
		if (is0())
			m_pData = new Tensor2s(eUnit, v_);
		else
			m_pData->Assign1(v_);
		return *this;
	}
	//   Tensor2& AssignRows(const Tensor1& r1_, const Tensor1& r2_, const Tensor1& r3_) {return
	//   r1_.is0()&&r2_.is0()&&r3_.is0()? (Assign0(),*this) :
	//   (AllocData(),pData->AssignRows(r1_,r2_,r3_),*this);} Tensor2& AssignCols(const Tensor1&
	//   r1_, const Tensor1& r2_, const Tensor1& r3_) {return r1_.is0()&&r2_.is0()&&r3_.is0()?
	//   (Assign0(),*this) : (AllocData(),pData->AssignCols(r1_,r2_,r3_),*this);} Tensor1s&
	//   EigenValues(Tensor1s& r_) const {return is0()? r_.Assign0() : pData->EigenValues(r_);}
	//   small_t EigenVectors(real_t ev_, Tensor2& evc_, small_t p_=1) const {return is0()?
	//   (evc_.Assign0(),3) : pData->EigenVectors(ev_,evc_,p_);} Tensor2&
	//   RotateToEigenSystem(Tensor2& r_) {return is0()? (r_.Assign1(),*this) :
	//   (pData->RotateToEigenSystem(r_),*this);} Tensor1& AttachedVector(Tensor1& r_) const {return
	//   is0()? r_.Assign0() : pData->AttachedVector(r_);}
	Tensor2a& AssignDev()
	{
		return is0() ? *this : (m_pData->AssignDev(), *this);
	}
	Tensor2a& Dev(Tensor2a& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	//   Tensor1& ScalarProductWithBasisVector     (small_t i_, Tensor1& r_) const {return is0()?
	//   r_.Assign0() : pData->ScalarProductWithBasisVector(i_,r_);} Tensor1&
	//   ScalarProductWithBasisVector_left(small_t i_, Tensor1& r_) const {return is0()?
	//   r_.Assign0() : pData->ScalarProductWithBasisVector_left(i_,r_);}
	Tensor1a& DotProduct(const Tensor1a& o_, Tensor1a& r_) const
	{
		return is0() ? r_.Assign0() : (r_ = o_).DotMultiply_left(*this);
	}
	RetTensor1 DotProduct(const Tensor1a& o_) const
	{
		Tensor1a result;
		return DotProduct(o_, result);
	}
	Tensor2a& CrossMultiply(const Tensor1a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->CrossMultiply(*o_.m_pData);
		return *this;
	}
	Tensor2a& CrossMultiply_left(const Tensor1a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->CrossMultiply_left(*o_.m_pData);
		return *this;
	}
	Tensor2a& CrossProduct(const Tensor1a& o_, Tensor2a& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).CrossMultiply(o_);
	}
	RetTensor2 CrossProduct(const Tensor1a&) const;
	Tensor3s& DirectProduct(const Tensor1a&, Tensor3s&) const;
	real_t Dot2Product(const Tensor2a& o_) const
	{
		return is0() || o_.is0() ? 0. : m_pData->Dot2Product(*o_.m_pData);
	}
	Tensor2a& DotMultiply(const Tensor2a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->DotMultiply(*o_.m_pData);
		return *this;
	}
	Tensor2a& DotMultiply_left(const Tensor2a& o_)
	{
		if (__AdjustMult(o_))
			m_pData->DotMultiply_left(*o_.m_pData);
		return *this;
	}
	Tensor2a& DotProduct(const Tensor2a& o_, Tensor2a& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).DotMultiply(o_);
	}
	RetTensor2 DotProduct(const Tensor2a&) const;
	Tensor3s& CrossProduct(const Tensor2a&, Tensor3s&) const;
	Tensor4s& DirectProduct(const Tensor2a&, Tensor4s&) const;
	Tensor1a& Dot2Product(const Tensor3s&, Tensor1a&) const;
	RetTensor1 Dot2Product(const Tensor3s& o_) const
	{
		Tensor1a result;
		return Dot2Product(o_, result);
	}
	Tensor3s& DotProduct(const Tensor3s&, Tensor3s&) const;
	Tensor4s& CrossProduct(const Tensor3s&, Tensor4s&) const;
	//   Tensor2a& Dot2Multiply      (const Tensor4&,small_t,small_t);
	Tensor2a& Dot2Multiply(const SymmetricTensor4s&, small_t, small_t);
	Tensor2a& Dot2Multiply(const Tensor4s&);
	Tensor2a& Dot2Multiply_left(const Tensor4s&);
	Tensor2a& Dot2Product(const Tensor4s&, Tensor2a&) const;
	RetTensor2 Dot2Product(const Tensor4s&) const;
	Tensor4s& DotProduct(const Tensor4s&, Tensor4s&) const;
	//   Tensor2& ApplyFunction(const fClassicFunction& fun_) {AllocData(0.);
	//   pData->ApplyFunction(fun_); return *this;}

protected:
	Tensor2s* m_pData;
	//   Tensor2a(Tensor2s* addr_): pData(addr_) {}

private:
	friend class Tensor1s;
	friend class Tensor1a;
	friend class RetTensor2;

	Tensor2s* __pData()
	{
		return is0() ? (m_pData = new Tensor2s) : m_pData;
	}
	bool __AdjustAdd(const Tensor2a& o_)
	{
		return is0() ? (o_.is0() ? false : (m_pData = new Tensor2s(*o_.m_pData), false))
					 : (o_.is0() ? false : true);
	}
	template <typename T>
	bool __AdjustMult(const T& o_)
	{
		return is0() ? false : (o_.is0() ? Assign0(), false : true);
	}
	template <typename Tfactor, typename Tresult>
	bool __AdjustMult(const Tfactor& o_, Tresult& r_) const
	{
		return is0() || o_.is0() ? r_.Assign0(), false : true;
	}
	const real_t* operator[](small_t i_) const
	{
		return m_pData->operator[](i_);
	}
};
// std::ostream& operator<< (std::ostream&, const Tensor2a&);

class RetTensor2
{
public:
	bool is0() const
	{
		return m_pData == nullptr;
	}
	~RetTensor2()
	{
		delete m_pData;
	}
	void SwapData(Tensor2s*);

private:
	friend class Tensor1a;
	friend class Tensor2a;

	RetTensor2() : m_pData(nullptr)
	{
	}
	RetTensor2(Tensor2a& o_) : m_pData(o_.m_pData)
	{
		o_.m_pData = nullptr;
	}
	//   RetTensor2(const RetTensor2& o_): pData(o_.pData) {o_.pData = nullptr;}
	//   RetTensor2(const real_t* const p_): pData(p_==nullptr ? nullptr : new Tensor1s(p_)) {}
	//   RetTensor2(real_t x_, real_t y_, real_t z_): pData(new Tensor1s(x_,y_,z_)) {}

	mutable Tensor2s* m_pData;
};

class Tensor2CaptureDataAdapter
{
public:
	explicit Tensor2CaptureDataAdapter(Tensor2a& o_) : m_ref(o_)
	{
	}

private:
	friend class Tensor2a;

	Tensor2a& m_ref;
};

class SymmetricTensor2s
{
public:
	SymmetricTensor2s() : m_data()
	{
	}
	////   Tensor2s(tensor_kind_t, real_t v_=0.) {__SET_ROWS; Assign1(v_);}
	explicit SymmetricTensor2s(real_t v_)
	{
		__Assign(v_);
	}
	SymmetricTensor2s(const SymmetricTensor2s& o_)
	{
		__Assign(o_.m_data);
	}
	real_t operator()(small_t, small_t) const;
	//         real_t& operator()(small_t, small_t);
	//   SymmetricTensor2s& operator=(const SymmetricTensor2s&);
	//   SymmetricTensor2s& operator-=(const SymmetricTensor2s&);
	SymmetricTensor2s& operator+=(const SymmetricTensor2s& o_)
	{
		m_data[0] += o_.m_data[0];
		m_data[1] += o_.m_data[1];
		m_data[2] += o_.m_data[2];
		m_data[3] += o_.m_data[3];
		m_data[4] += o_.m_data[4];
		m_data[5] += o_.m_data[5];
		return *this;
	}
	SymmetricTensor2s& operator*=(real_t v_)
	{
		std::for_each(m_data, m_data + 6, fMilt<real_t, real_t>(v_));
		return *this;
	}
	//   SymmetricTensor2s& operator/=(real_t);
	SymmetricTensor2s operator*(real_t v_)
	{
		SymmetricTensor2s res(*this);
		return res *= v_;
	}
	//   SymmetricTensor2s& Negate() {std::transform(Data, Data+6, Data, std::negate<real_t>());
	//   return *this;}
	////   SymmetricTensor2s& Transpose() {return *this;}
	SymmetricTensor2s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	SymmetricTensor2s& Assign1(real_t v_ = 1.)
	{
		m_data[0] = m_data[1] = m_data[2] = v_;
		m_data[3] = m_data[4] = m_data[5] = 0.;
		return *this;
	}
	real_t Contraction() const
	{
		return m_data[0] + m_data[1] + m_data[2];
	}
	//   real_t Det() const;
	real_t I1() /*Linear    invariant*/ const
	{
		return Contraction();
	}
	//   real_t I2()/*Quadratic invariant*/ const;
	//   real_t I3()/*Cubic     invariant*/ const {return Det();}
	//   real_t Minor(small_t,small_t) const;
	//   SymmetricTensor2s& Invert();
	////   Tensor1s& AttachedVector(Tensor1s&) const;
	////   Tensor1s  AttachedVector() const;
	////   Tensor1s& ScalarProductWithBasisVector     (small_t, Tensor1s&) const;
	////   Tensor1s& ScalarProductWithBasisVector_left(small_t, Tensor1s&) const;
	////   Tensor1s  ScalarProductWithBasisVector     (small_t) const;
	////   Tensor1s  ScalarProductWithBasisVector_left(small_t) const;
	//   void EigenValues(real_t[3]) const;
	////   small_t EigenVectors(real_t/*current e.value*/, SymmetricTensor2s&/*e.vects*/,
	/// small_t=1/*begin position*/) const; /   void Eigen(real_t[3],SymmetricTensor2s&) const; /
	/// SymmetricTensor2s& EigenVectors(SymmetricTensor2s& r_) const {real_t evs[3]; Eigen(evs,r_);
	/// return r_;} /   SymmetricTensor2s& RotateToEigenSystem(SymmetricTensor2s&);
	SymmetricTensor2s& AssignDev();
	SymmetricTensor2s& Dev(SymmetricTensor2s& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	//   Tensor1s& DotProduct      (const Tensor1s&, Tensor1s&) const;
	////   SymmetricTensor2s& CrossMultiply     (const Tensor1s&);
	////   SymmetricTensor2s& CrossMultiply_left(const Tensor1s&);
	////   SymmetricTensor2s& CrossProduct      (const Tensor1s&, SymmetricTensor2s&) const;

private:
	friend class Tensor2s;

	void __Assign(real_t v_)
	{
		std::fill(m_data, m_data + 6, v_);
	}
	void __Assign(const real_t* v_)
	{
		std::copy(v_, v_ + 6, m_data);
	}

	real_t m_data[6];  // 11 22 33 12 23 31
};

inline real_t SymmetricTensor2s::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3,
		"invalid index in SymmetricTensor2s::operator() const");
#endif
	return (i_ == j_) ? m_data[--i_] : m_data[(25 - 3 * (i_ + j_)) * (i_ + j_) / 2 - 21];
}

inline SymmetricTensor2s& SymmetricTensor2s::AssignDev()
{
	const real_t p = I1() / 3.;
	m_data[0] -= p;
	m_data[1] -= p;
	m_data[2] -= p;
	return *this;
}

class Tensor3s
{
private:
#define __SET_ROWS3                                                    \
	*m_rows = static_cast<pArr33_t>(static_cast<void*>(m_data));       \
	m_rows[1] = static_cast<pArr33_t>(static_cast<void*>(m_data + 9)); \
	m_rows[2] = static_cast<pArr33_t>(static_cast<void*>(m_data + 18));

	struct tIndex
	{
	public:
		tIndex() : m_first(0), m_second(0), m_third(0)
		{
		}
		tIndex(small_t i_, small_t j_, small_t k_) : m_first(i_), m_second(j_), m_third(k_)
		{
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);
		tIndex& operator()(small_t i_, small_t j_, small_t k_)
		{
			m_first = i_;
			m_second = j_;
			m_third = k_;
			return *this;
		}

	private:
		small_t m_first, m_second, m_third;
	};

public:
	Tensor3s()
		: m_data(),
		  m_rows(){__SET_ROWS3} Tensor3s(tensor_kind_t, real_t v_ = 1.)
		: m_rows(), m_data()
	{
		__SET_ROWS3;
		Assign1(v_);
	}  // Levy-Chivita
	explicit Tensor3s(real_t v_)
	{
		__SET_ROWS3;
		__Assign(v_);
	}
	Tensor3s(const Tensor3s& o_)
	{
		__SET_ROWS3;
		__Assign(o_.m_data);
	}
	real_t operator[](const tIndex&) const;
	// real_t& operator[](const tIndex&);
	real_t operator()(small_t, small_t, small_t) const;
	real_t& operator()(small_t, small_t, small_t);
	Tensor3s& operator=(const Tensor3s&);
	Tensor3s& operator-=(const Tensor3s&);
	Tensor3s& operator+=(const Tensor3s&);
	Tensor3s& operator*=(real_t v_)
	{
		std::for_each(m_data, m_data + 27, fMilt<real_t, real_t>(v_));
		return *this;
	}
	Tensor3s& operator/=(real_t);
	Tensor3s& Negate()
	{
		std::transform(m_data, m_data + 27, m_data, std::negate<real_t>());
		return *this;
	}
	bool is0() const
	{
		return false;
	}
	Tensor3s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor3s& Assign1(real_t = 1.);	 // Levy-Chivita
	Tensor1s& Contraction(small_t, small_t, Tensor1s&) const;
	Tensor2s& ScalarProductWithBasisVector(small_t, Tensor2s&) const;
	Tensor2s& ScalarProductWithBasisVector_left(small_t, Tensor2s&) const;
	Tensor2s& DotProduct(const Tensor1s&, Tensor2s&) const;
	Tensor3s& CrossMultiply(const Tensor1s&);
	Tensor3s& CrossMultiply_left(const Tensor1s&);
	Tensor3s& CrossProduct(const Tensor1s&, Tensor3s&) const;
	Tensor4s& DirectProduct(const Tensor1s&, Tensor4s&) const;
	Tensor1s& Dot2Product(const Tensor2s&, Tensor1s&) const;
	Tensor3s& DotMultiply(const Tensor2s&);
	Tensor3s& DotMultiply_left(const Tensor2s&);
	Tensor3s& DotProduct(const Tensor2s&, Tensor3s&) const;
	real_t Dot3Product(const Tensor3s&) const;
	Tensor2s& Dot2Product(const Tensor3s&, Tensor2s&) const;
	Tensor4s& DotProduct(const Tensor3s&, Tensor4s&) const;
	Tensor1s& Dot3Product(const Tensor4s&, Tensor1s&) const;
	Tensor3s& Dot2Multiply(const Tensor4s&);
	Tensor3s& Dot2Multiply_left(const Tensor4s&);
	Tensor3s& Dot2Product(const Tensor4s&, Tensor3s&) const;

private:
	friend class Tensor1s;
	friend class Tensor2s;
	friend class Tensor4s;

	typedef real_t (*pArr33_t)[3][3];
	typedef real_t Arr33_t[3][3];

	void __Assign(real_t v_)
	{
		std::fill(m_data, m_data + 27, v_);
	}
	void __Assign(const real_t* v_)
	{
		std::copy(v_, v_ + 27, m_data);
	}
	const Arr33_t& operator[](small_t i_) const
	{
		return *m_rows[i_];
	}
	Arr33_t& operator[](small_t i_)
	{
		return *m_rows[i_];
	}

	real_t m_data[27];
	//   real_t (*Rows[3])[3][3];
	pArr33_t m_rows[3];
};

/*ostream& operator<< (ostream& out_, Tensor3s& t_)
{
 Tensor3::Index index;
 for (small_t i=0,j; i<3; ++i)
   {
	out_<<"| ";
	for (j=0; j<3; ++j)
	   out_<<t_[index(i,j,0)]<<' '<<t_[index(i,j,1)]<<' '<<t_[index(i,j,2)]<<"\t\t";
	out_<<"|\n";
   }
 return out_;
}*/

class Tensor4
{
public:
	bool is0() const
	{
		return false;
	}

protected:
	struct tIndex
	{
	public:
		tIndex() : m_indices()
		{
		}
		tIndex(const small_t i_, const small_t j_, const small_t k_, const small_t l_)
		{
			m_indices[0] = i_;
			m_indices[1] = j_;
			m_indices[2] = k_;
			m_indices[3] = l_;
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);

	private:
		small_t m_indices[4];
	};
	//  virtual real_t operator[](tIndex) const =0;

private:
	friend class Tensor2s;
};

class SymmetricTensor4s;

class Tensor4s : public Tensor4
{
private:
#define __SET_ROWS4                                                      \
	*m_rows = static_cast<pArr333_t>(static_cast<void*>(m_data));        \
	m_rows[1] = static_cast<pArr333_t>(static_cast<void*>(m_data + 27)); \
	m_rows[2] = static_cast<pArr333_t>(static_cast<void*>(m_data + 54));

public:
	Tensor4s()
		: m_data(),
		  m_rows(){__SET_ROWS4} Tensor4s(tensor_kind_t, small_t kind_, real_t v_ = 1.)
		: m_data(), m_rows()
	{
		__SET_ROWS4;
		Assign1(kind_, v_);
	}  // four isotropic tensors
	explicit Tensor4s(real_t v_)
	{
		__SET_ROWS4;
		__Assign(v_);
	}
	explicit Tensor4s(const SymmetricTensor4s&);
	Tensor4s(const Tensor4s& o_)
	{
		__SET_ROWS4;
		__Assign(o_.m_data);
	}
	const real_t& operator()(small_t, small_t, small_t, small_t) const;
	real_t& operator()(small_t, small_t, small_t, small_t);
	//   Tensor4s& operator=  (const Tensor4&);
	Tensor4s& operator=(const Tensor4s&);
	Tensor4s& operator=(const SymmetricTensor4s&);
	Tensor4s& operator-=(const Tensor4s&);
	Tensor4s& operator+=(const Tensor4s&);
	Tensor4s& operator*=(real_t v_)
	{
		std::for_each(m_data, m_data + 81, fMilt<real_t, real_t>(v_));
		return *this;
	}
	Tensor4s& operator/=(real_t);
	Tensor4s& Negate()
	{
		std::transform(m_data, m_data + 81, m_data, std::negate<real_t>());
		return *this;
	}
	Tensor4s& Assign0(small_t, small_t, small_t, small_t);
	Tensor4& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor4s& Assign1(small_t, real_t = 1.);  // three isotropic tensors
	Tensor2s& Contraction(small_t, small_t, Tensor2s&) const;
	Tensor3s& ScalarProductWithBasisVector(small_t, Tensor3s&) const;
	Tensor3s& ScalarProductWithBasisVector_left(small_t, Tensor3s&) const;
	Tensor3s& DotProduct(const Tensor1s&, Tensor3s&) const;
	Tensor4s& CrossMultiply(const Tensor1s&);
	Tensor4s& CrossMultiply_left(const Tensor1s&);
	Tensor4s& CrossProduct(const Tensor1s&, Tensor4s&) const;
	Tensor2s& Dot2Product(const Tensor2s&, Tensor2s&) const;
	Tensor4s& DotMultiply(const Tensor2s&, small_t, small_t);
	Tensor4s& DotMultiply(const Tensor2s&);
	Tensor4s& DotMultiply_left(const Tensor2s&);
	Tensor4s& DotProduct(const Tensor2s&, Tensor4s&) const;
	Tensor1s& Dot3Product(const Tensor3s&, Tensor1s&) const;
	Tensor3s& Dot2Product(const Tensor3s&, Tensor3s&) const;
	real_t Dot4Product(const Tensor4s&) const;
	Tensor2s& Dot3Product(const Tensor4s&, Tensor2s&) const;
	Tensor4s& Dot2Multiply(const Tensor4s&);
	Tensor4s& Dot2Multiply_left(const Tensor4s&);
	Tensor4s& Dot2Product(const Tensor4s&, Tensor4s&) const;
	Tensor4s& TransformAll(const Tensor2s&);

private:
	friend class Tensor1s;
	friend class Tensor2s;
	friend class Tensor3s;

	typedef real_t (*pArr333_t)[3][3][3];
	typedef real_t Arr333_t[3][3][3];

	void __Assign(real_t v_)
	{
		std::fill(m_data, m_data + 81, v_);
	}
	void __Assign(const real_t* v_)
	{
		std::copy(v_, v_ + 81, m_data);
	}
	const Arr333_t& operator[](small_t i_) const
	{
		return *m_rows[i_];
	}
	Arr333_t& operator[](small_t i_)
	{
		return *m_rows[i_];
	}
	/*virtual*/ real_t operator[](tIndex) const;
	/*virtual*/ real_t& operator[](tIndex);

	real_t m_data[81];
	//   real_t (*Rows[3])[3][3][3];
	pArr333_t m_rows[3];
};
// std::ostream& operator<< (std::ostream&, Tensor4&);

class SymmetricTensor4s : public Tensor4
{
public:
	real_t operator()(small_t, small_t) const;
	real_t& operator()(small_t, small_t);
	real_t operator()(small_t, small_t, small_t, small_t) const;
	SymmetricTensor4s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	SymmetricTensor4s& TransformAll(const Tensor2s&);

private:
	friend class Tensor4s;
	friend class Tensor2s;

	void __Assign(real_t v_)
	{
		std::fill(m_data, m_data + 21, v_);
	}
	void __Assign(const real_t* v_)
	{
		std::copy(v_, v_ + 21, m_data);
	}
	/*virtual*/ real_t operator[](tIndex) const;

	real_t m_data[21];	// 11..16,22..26,33..36,44..46,55,56,66
};

//==========================INLINES:============================================

//==========================Tensor1s============================================

inline Tensor1s::Tensor1s(
	tensor_kind_t, small_t number_, real_t v_)	// Basis vector multiplied by v_
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= 3,
		"invalid number of coordinate in Tensor1s::Tensor1s(bool,small_t,real_t)");
#endif
	m_data[--number_] = v_;
	m_data[number_ == 0 ? 1 : 0] = m_data[number_ == 2 ? 1 : 2] = 0.;
}

inline Tensor1s::Tensor1s(const real_t* const p_)
{
#ifdef STRONGCHECK
	Assert(p_ != nullptr, "non-allocated data in Tensor1s::Tensor1s(const real_t* const)");
#endif
	__Assign(p_);
}

inline Tensor1s::Tensor1s(const Tensor1a& o_)
{
	operator=(o_);
}

inline Tensor1s& Tensor1s::operator=(const Tensor1a& o_)
{
	if (o_.is0())
		Assign0();
	else
		__Assign(*o_.m_pData);
	return *this;
}

inline real_t& Tensor1s::operator()(small_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "invalid index in Tensor1s::()");
#endif
	return m_data[--i_];
}

inline const real_t& Tensor1s::operator()(small_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "index is beyond the range in Tensor1s::()const");
#endif
	return m_data[--i_];
}

inline Tensor1s& Tensor1s::operator/=(real_t v_)
{
	// std::cout << "\n\tRun Tensor1s::operator/=\n";

#ifdef STRONGCHECK
	Assert(mathdef::is_not_zero(v_), "zero divisor in Tensor1s::operator/=(real_t)");
#endif
	return operator*=(1. / v_);
}

inline Tensor1s& Tensor1s::CrossMultiplyByBasisVector(small_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "invalid index in Tensor1s::CrossMultiplyByBasisVector");
#endif
	const small_t j = (--i_) == 0 ? 1 : (i_ == 1 ? 2 : 0), k = 3 - i_ - j;
	m_data[i_] = m_data[j];
	m_data[j] = m_data[k];
	m_data[k] = -m_data[i_];
	m_data[i_] = 0.;
	return *this;
}

inline Tensor1s& Tensor1s::CrossMultiply(const Tensor1s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "vector multiplication by itself in Tensor1s::CrossMultiply");
#endif
	// if (&o_ == this) {Assign0(); return *this;}
	real_t data_0 = m_data[0], data_1 = m_data[1];
	m_data[0] = data_1 * o_[2] - m_data[2] * o_[1];
	m_data[1] = m_data[2] * o_[0] - data_0 * o_[2];
	m_data[2] = data_0 * o_[1] - data_1 * o_[0];
	return *this;
}

inline Tensor1s& Tensor1s::CrossMultiply_left(const Tensor1s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "vector multiplication by itself in Tensor1s::CrossMultiply_left");
#endif
	// if (&o_ == this) {Assign0(); return *this;}
	real_t data_0 = m_data[0], data_1 = m_data[1];
	m_data[0] = o_[1] * m_data[2] - o_[2] * data_1;
	m_data[1] = o_[2] * data_0 - o_[0] * m_data[2];
	m_data[2] = o_[0] * data_1 - o_[1] * data_0;
	return *this;
}

inline Tensor1s& Tensor1s::CrossMultiplyByBasisVector_left(small_t i_)
{
	CrossMultiplyByBasisVector(i_);
	m_data[0] = -m_data[0];
	m_data[1] = -m_data[1];
	m_data[2] = -m_data[2];
	return *this;
}

inline Tensor1s& Tensor1s::CrossProductWithBasisVector(small_t i_, Tensor1s& result_) const
{
	return (result_ = *this).CrossMultiplyByBasisVector(i_);
}

inline Tensor1s& Tensor1s::CrossProductWithBasisVector_left(small_t i_, Tensor1s& result_) const
{
	return (result_ = *this).CrossMultiplyByBasisVector_left(i_);
}

inline real_t Tensor1s::DotProduct(const Tensor1s& o_) const
{
	return m_data[0] * o_.m_data[0] + m_data[1] * o_.m_data[1] + m_data[2] * o_.m_data[2];
}

inline Tensor1s& Tensor1s::DotProduct(const Tensor2s& o_, Tensor1s& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}

inline Tensor1s& Tensor1s::CrossProduct(const Tensor1s& o_, Tensor1s& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}

inline Tensor2s& Tensor1s::CrossProduct(const Tensor2s& o_, Tensor2s& result_) const
{  //!!!
	return (result_ = o_).CrossMultiply_left(*this);
}

inline Tensor3s& Tensor1s::CrossProduct(const Tensor3s& o_, Tensor3s& result_) const
{
	return (result_ = o_).CrossMultiply_left(*this);
}

inline Tensor4s& Tensor1s::CrossProduct(const Tensor4s& o_, Tensor4s& result_) const
{
	return (result_ = o_).CrossMultiply_left(*this);
}

inline Tensor2a& Tensor1s::DirectProduct(const Tensor1s& o_, Tensor2a& result_) const
{
	DirectProduct(o_, *result_.__pData());
	return result_;
}

//==========================Tensor1a=============================================

inline Tensor1a::Tensor1a(RetTensor1& o_) : m_pData(o_.m_pData)
{
	o_.m_pData = nullptr;
}

inline const real_t& Tensor1a::operator()(small_t i_) const
{
#ifdef STRONGCHECK
	Assert(!is0(), "nullptr data in Tensor1a::operator()const");
#endif
	return m_pData->operator()(i_);
}

inline real_t& Tensor1a::operator()(small_t i_)
{
#ifdef STRONGCHECK
	Assert(!is0(), "nullptr data in Tensor1a::operator()");
#endif
	return m_pData->operator()(i_);
}

inline Tensor1a& Tensor1a::operator=(RetTensor1& o_)
{
	std::swap(m_pData, o_.m_pData);
	return *this;
}

inline Tensor1a& Tensor1a::Normalize()
{
#ifdef STRONGCHECK
	Assert(!is0(), "attempt to normalize zero vector in Tensor1a::Normalize()");
#endif
	m_pData->Normalize();
	return *this;
}

inline Tensor1a& Tensor1a::operator/=(real_t v_)
{
#ifdef STRONGCHECK
	Assert(mathdef::is_not_zero(v_), "zero divisor in Tensor1a::operator/=(real_t)");
#endif
	return is0() ? *this : (m_pData->operator/=(v_), *this);
}

inline RetTensor1 Tensor1a::operator-() const
{
	return is0() ? RetTensor1()
				 : RetTensor1(-m_pData->m_data[0], -m_pData->m_data[1], -m_pData->m_data[2]);
}

inline RetTensor1 Tensor1a::operator-(const Tensor1a& o_) const
{
	return Tensor1a(*this) -= o_;
}

inline RetTensor1 Tensor1a::operator+(const Tensor1a& o_) const
{
	return Tensor1a(*this) += o_;
}

inline RetTensor1 Tensor1a::operator*(real_t v_) const
{
	return Tensor1a(*this) *= v_;
}

inline RetTensor1 Tensor1a::operator/(real_t v_) const
{
	return Tensor1a(*this) /= v_;
}

inline RetTensor1 Tensor1a::CrossProductWithBasisVector(small_t i_) const
{
	Tensor1a result;
	return CrossProductWithBasisVector(i_, result);
}

inline RetTensor1 Tensor1a::CrossProductWithBasisVector_left(small_t i_) const
{
	Tensor1a result;
	return CrossProductWithBasisVector_left(i_, result);
}

inline RetTensor1 Tensor1a::CrossProduct(const Tensor1a& o_) const
{
	Tensor1a result;
	return CrossProduct(o_, result);
}

inline RetTensor1 Tensor1a::DotProduct(const Tensor2a& o_) const
{
	Tensor1a result;
	return DotProduct(o_, result);
}

inline Tensor2a& Tensor1a::DirectProduct(const Tensor1a& o_, Tensor2a& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DirectProduct(*o_.m_pData, *result_.__pData());
	return result_;
}

inline RetTensor2 Tensor1a::DirectProduct(const Tensor1a& o_) const
{
	Tensor2a result;
	return DirectProduct(o_, result);
}

inline Tensor1a& Tensor1a::DotMultiply(const Tensor2a& o_)
{
	if (__AdjustMult(o_))
		m_pData->DotMultiply(*o_.m_pData);
	return *this;
}

inline Tensor1a& Tensor1a::DotMultiply_left(const Tensor2a& o_)
{
	if (__AdjustMult(o_))
		m_pData->DotMultiply_left(*o_.m_pData);
	return *this;
}

inline Tensor1a& Tensor1a::DotProduct(const Tensor2a& o_, Tensor1a& result_) const
{
	return o_.is0() ? result_.Assign0() : (result_ = *this).DotMultiply(o_);
}

inline Tensor2a& Tensor1a::CrossProduct(const Tensor2a& o_, Tensor2a& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->CrossProduct(*o_.m_pData, *result_.__pData());
	return result_;
}

inline Tensor3s& Tensor1a::DirectProduct(const Tensor2a& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DirectProduct(*o_.m_pData, result_);
	return result_;
}

inline Tensor2a& Tensor1a::DotProduct(const Tensor3s& o_, Tensor2a& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DotProduct(o_, *result_.__pData());
	return result_;
}

inline Tensor3s& Tensor1a::CrossProduct(const Tensor3s& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->CrossProduct(o_, result_);
	return result_;
}

inline Tensor4s& Tensor1a::DirectProduct(const Tensor3s& o_, Tensor4s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DirectProduct(o_, result_);
	return result_;
}

inline Tensor3s& Tensor1a::DotProduct(const Tensor4s& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DotProduct(o_, result_);
	return result_;
}

inline Tensor4s& Tensor1a::CrossProduct(const Tensor4s& o_, Tensor4s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->CrossProduct(o_, result_);
	return result_;
}

inline RetTensor2 Tensor1a::CrossProduct(const Tensor2a& o_) const
{
	Tensor2a result;
	return CrossProduct(o_, result);
}

inline RetTensor2 Tensor1a::DotProduct(const Tensor3s& o_) const
{
	Tensor2a result;
	return DotProduct(o_, result);
}

//==========================Tensor2s=======================================

inline small_t Tensor2s::tIndex::operator[](small_t no_) const
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 2, "invalid Index2");
#endif
	return (no_ == 1) ? m_first : m_second;
}

inline small_t& Tensor2s::tIndex::operator[](small_t no_)
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 2, "invalid Index2 const");
#endif
	return (no_ == 1) ? m_first : m_second;
}

inline const real_t& Tensor2s::operator[](Tensor2s::tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[1] < 3 && i_[2] >= 0 && i_[2] < 3,
		"invalid index in Tensor2s::operator[]const");
#endif
	return m_rows[i_[1]][i_[2]];
}

inline real_t& Tensor2s::operator[](Tensor2s::tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[1] < 3 && i_[2] >= 0 && i_[2] < 3,
		"invalid index in Tensor2s::operator[]");
#endif
	return m_rows[i_[1]][i_[2]];
}

inline const real_t& Tensor2s::operator()(Tensor2s::tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] > 0 && i_[1] <= 3 && i_[2] > 0 && i_[2] <= 3,
		"invalid index in Tensor2s::operator()const");
#endif
	return m_rows[--i_[1]][--i_[2]];
}

inline real_t& Tensor2s::operator()(tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] > 0 && i_[1] <= 3 && i_[2] > 0 && i_[2] <= 3,
		"invalid index in Tensor2s::operator()");
#endif
	return m_rows[--i_[1]][--i_[2]];
}

inline const real_t& Tensor2s::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3,
		"invalid index in real_t Tensor2s::operator() const");
#endif
	return m_rows[--i_][--j_];
}

inline real_t& Tensor2s::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3, "invalid index in real_t Tensor2s::operator()");
#endif
	return m_rows[--i_][--j_];
}

inline Tensor2s& Tensor2s::operator=(const Tensor2s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor2s::operator= (const Tensor2s&)");
#endif
	// if (&o_!=this)
	__Assign(o_.m_data);
	return *this;
}

inline Tensor2s& Tensor2s::operator/=(real_t v_)
{
#ifdef STRONGCHECK
	Assert(mathdef::is_not_zero(v_), "zero divisor in Tensor2s::operator/=(real_t)");
#endif
	return operator*=(1. / v_);
}

inline Tensor2s& Tensor2s::AddTransposed(const Tensor2s& o_)
{
	if (&o_ == this)
	{
		m_data[0] += m_data[0];
		m_data[1] += m_data[3];
		m_data[2] += m_data[6];
		m_data[3] = m_data[1];
		m_data[4] += m_data[4];
		m_data[5] += m_data[7];
		m_data[6] = m_data[2];
		m_data[7] = m_data[5];
		m_data[8] += m_data[8];
	}
	else
	{
		m_data[0] += o_.m_data[0];
		m_data[1] += o_.m_data[3];
		m_data[2] += o_.m_data[6];
		m_data[3] += o_.m_data[1];
		m_data[4] += o_.m_data[4];
		m_data[5] += o_.m_data[7];
		m_data[6] += o_.m_data[2];
		m_data[7] += o_.m_data[5];
		m_data[8] += o_.m_data[8];
	}
	return *this;
}

inline SymmetricTensor2s& Tensor2s::Symmetrization_twofold(SymmetricTensor2s& result_) const
{
	// result_.Data[6] ~ 11 22 33 12 23 31
	result_.m_data[0] = m_data[0] + m_data[0];
	result_.m_data[3] = m_data[1] + m_data[3];
	result_.m_data[5] = m_data[2] + m_data[6];
	result_.m_data[1] = m_data[4] + m_data[4];
	result_.m_data[4] = m_data[5] + m_data[7];
	result_.m_data[2] = m_data[8] + m_data[8];
	return result_;
}

inline Tensor2s& Tensor2s::Transpose()
{
	std::iter_swap(m_data + 1, m_data + 3);
	std::iter_swap(m_data + 2, m_data + 6);
	std::iter_swap(m_data + 5, m_data + 7);
	return *this;
}

inline real_t Tensor2s::Det() const
{
	return m_data[0] * m_data[4] * m_data[8] + m_data[1] * m_data[5] * m_data[6] +
		   m_data[2] * m_data[3] * m_data[7] - m_data[2] * m_data[4] * m_data[6] -
		   m_data[0] * m_data[5] * m_data[7] - m_data[1] * m_data[3] * m_data[8];
}

inline real_t Tensor2s::I2() const
{
	return m_data[0] * m_data[4] + m_data[4] * m_data[8] + m_data[8] * m_data[0] -
		   m_data[1] * m_data[3] - m_data[5] * m_data[7] - m_data[6] * m_data[2];
}

inline Tensor1s& Tensor2s::AttachedVector(Tensor1s& result_) const
{
	return result_.Assign(
		0.5 * (m_data[7] - m_data[5]),
		0.5 * (m_data[2] - m_data[6]),
		0.5 * (m_data[3] - m_data[1]));
}

inline Tensor1s Tensor2s::AttachedVector() const
{
	return Tensor1s(
		0.5 * (m_data[7] - m_data[5]),
		0.5 * (m_data[2] - m_data[6]),
		0.5 * (m_data[3] - m_data[1]));
}

inline Tensor1s& Tensor2s::ScalarProductWithBasisVector(small_t i_, Tensor1s& result_) const
{
	--i_;
	return result_.Assign(m_rows[0][i_], m_rows[1][i_], m_rows[2][i_]);
}

inline Tensor1s& Tensor2s::ScalarProductWithBasisVector_left(small_t i_, Tensor1s& result_) const
{
	--i_;
	return result_.Assign(m_rows[i_][0], m_rows[i_][1], m_rows[i_][2]);
}

inline Tensor1s Tensor2s::ScalarProductWithBasisVector(small_t i_) const
{
	--i_;
	return Tensor1s(m_rows[0][i_], m_rows[1][i_], m_rows[2][i_]);
}

inline Tensor1s Tensor2s::ScalarProductWithBasisVector_left(small_t i_) const
{
	--i_;
	return Tensor1s(m_rows[i_][0], m_rows[i_][1], m_rows[i_][2]);
}

inline Tensor2s& Tensor2s::AssignDev()
{
	const real_t p = I1() / 3.;
	m_data[0] -= p;
	m_data[4] -= p;
	m_data[8] -= p;
	return *this;
}

inline Tensor1s& Tensor2s::DotProduct(const Tensor1s& o_, Tensor1s& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}

inline Tensor2s& Tensor2s::CrossProduct(const Tensor1s& o_, Tensor2s& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}

inline Tensor2s& Tensor2s::DotProduct(const Tensor2s& o_, Tensor2s& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}

inline Tensor3s& Tensor2s::DotProduct(const Tensor3s& o_, Tensor3s& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}

inline Tensor2s& Tensor2s::Dot2Product(const Tensor4s& o_, Tensor2s& result_) const
{
	return (result_ = *this).Dot2Multiply(o_);
}

inline Tensor4s& Tensor2s::DotProduct(const Tensor4s& o_, Tensor4s& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}

//==========================Tensor2a=============================================

inline Tensor2a::Tensor2a(RetTensor2& o_) : m_pData(o_.m_pData)
{
	o_.m_pData = nullptr;
}

inline const real_t& Tensor2a::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(!is0(), "nullptr data in Tensor2a::operator()const");
#endif
	return m_pData->operator()(i_, j_);
}

inline real_t& Tensor2a::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(!is0(), "nullptr data in Tensor2a::operator()");
#endif
	return m_pData->operator()(i_, j_);
}

inline Tensor2a& Tensor2a::operator=(RetTensor2& o_)
{
	std::swap(m_pData, o_.m_pData);
	return *this;
}

inline Tensor2a& Tensor2a::operator+=(Tensor2CaptureDataAdapter o_)
{
	if (is0())
	{
		m_pData = o_.m_ref.m_pData;
		o_.m_ref.m_pData = nullptr;
	}
	else if (!o_.m_ref.is0())
		m_pData->operator+=(*o_.m_ref.m_pData);
	return *this;
}

inline Tensor2a& Tensor2a::AddTransposed(Tensor2CaptureDataAdapter o_)
{
	if (!o_.m_ref.is0())
		if (is0())
		{
			m_pData = o_.m_ref.m_pData;
			o_.m_ref.m_pData = nullptr;
			m_pData->Transpose();
		}
		else
			m_pData->AddTransposed(*o_.m_ref.m_pData);
	return *this;
}

inline RetTensor2 Tensor2a::operator-() const
{
	return Tensor2a(*this).Negate();
}

inline RetTensor2 Tensor2a::operator-(const Tensor2a& o_) const
{
	return Tensor2a(*this) -= o_;
}

inline RetTensor2 Tensor2a::operator+(const Tensor2a& o_) const
{
	return Tensor2a(*this) += o_;
}

inline RetTensor2 Tensor2a::operator*(real_t v_) const
{
	return Tensor2a(*this) *= v_;
}

inline RetTensor2 Tensor2a::operator/(real_t v_) const
{
	return Tensor2a(*this) /= v_;
}

inline Tensor2a& Tensor2a::Invert()
{
#ifdef STRONGCHECK
	Assert(!is0(), "nullptr data in Tensor2a::Invert()");
#endif
	m_pData->Invert();
	return *this;
}

inline RetTensor2 Tensor2a::DotProduct(const Tensor2a& o_) const
{
	return Tensor2a(*this).DotMultiply(o_);
}

inline Tensor2a& Tensor2a::operator/=(real_t v_)
{
#ifdef STRONGCHECK
	Assert(mathdef::is_not_zero(v_), "zero divisor in Tensor2a::operator/=(real_t)");
#endif
	return is0() ? *this : (m_pData->operator/=(v_), *this);
}

inline Tensor3s& Tensor2a::DirectProduct(const Tensor1a& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DirectProduct(*o_.m_pData, result_);
	return result_;
}

inline Tensor3s& Tensor2a::CrossProduct(const Tensor2a& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->CrossProduct(*o_.m_pData, result_);
	return result_;
}

inline Tensor4s& Tensor2a::DirectProduct(const Tensor2a& o_, Tensor4s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DirectProduct(*o_.m_pData, result_);
	return result_;
}

inline Tensor1a& Tensor2a::Dot2Product(const Tensor3s& o_, Tensor1a& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->Dot2Product(o_, *result_.__pData());
	return result_;
}

inline Tensor3s& Tensor2a::DotProduct(const Tensor3s& o_, Tensor3s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->DotProduct(o_, result_);
	return result_;
}

inline Tensor4s& Tensor2a::CrossProduct(const Tensor3s& o_, Tensor4s& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->CrossProduct(o_, result_);
	return result_;
}

// inline Tensor2a& Tensor2a::Dot2Multiply (const Tensor4& o_, small_t i_, small_t j_)
inline Tensor2a& Tensor2a::Dot2Multiply(const SymmetricTensor4s& o_, small_t i_, small_t j_)
{
	if (__AdjustMult(o_))
		m_pData->Dot2Multiply(o_, i_, j_);
	return *this;
}

inline Tensor2a& Tensor2a::Dot2Multiply(const Tensor4s& o_)
{
	if (__AdjustMult(o_))
		m_pData->Dot2Multiply(o_);
	return *this;
}

inline Tensor2a& Tensor2a::Dot2Multiply_left(const Tensor4s& o_)
{
	if (__AdjustMult(o_))
		m_pData->Dot2Multiply_left(o_);
	return *this;
}

inline Tensor2a& Tensor2a::Dot2Product(const Tensor4s& o_, Tensor2a& result_) const
{
	if (__AdjustMult(o_, result_))
		m_pData->Dot2Product(o_, *result_.__pData());
	return result_;
}

inline Tensor4s& Tensor2a::DotProduct(const Tensor4s& o_, Tensor4s& result_) const
{  //!!!!!!!!!!! может, лучше использовать DotMultiply?
	if (__AdjustMult(o_, result_))
		m_pData->DotProduct(o_, result_);
	return result_;
}

inline RetTensor2 Tensor2a::CrossProduct(const Tensor1a& o_) const
{
	Tensor2a result;
	return CrossProduct(o_, result);
}

inline RetTensor2 Tensor2a::Dot2Product(const Tensor4s& o_) const
{
	Tensor2a result;
	return Dot2Product(o_, result);
}

//==========================Tensor3s==========================================

inline small_t Tensor3s::tIndex::operator[](small_t indexNo_) const
{
#ifdef STRONGCHECK
	Assert(indexNo_ > 0 && indexNo_ <= 3, "invalid Index3 const");
#endif
	return (indexNo_ == 1) ? m_first : (indexNo_ == 2 ? m_second : m_third);
}

inline small_t& Tensor3s::tIndex::operator[](small_t indexNo_)
{
#ifdef STRONGCHECK
	Assert(indexNo_ > 0 && indexNo_ <= 3, "invalid Index3");
#endif
	return (indexNo_ == 1) ? m_first : (indexNo_ == 2 ? m_second : m_third);
}

inline real_t Tensor3s::operator[](const Tensor3s::tIndex& i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] < 3 && i_[2] < 3 && i_[3] < 3, "invalid index in Tensor3s::operator[](tIndex3)const");
#endif
	return (*m_rows[i_[1]])[i_[2]][i_[3]];
}

inline real_t Tensor3s::operator()(small_t i_, small_t j_, small_t k_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3,
		"invalid index in Tensor3s::operator()const");
#endif
	return (*m_rows[--i_])[--j_][--k_];
}

inline real_t& Tensor3s::operator()(small_t i_, small_t j_, small_t k_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3,
		"invalid index in Tensor3s::operator()");
#endif
	return (*m_rows[--i_])[--j_][--k_];
}

inline Tensor3s& Tensor3s::operator=(const Tensor3s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor2s::operator= (const Tensor2s&)");
#endif
	// if (&o_!=this)
	__Assign(o_.m_data);
	return *this;
}

inline Tensor3s& Tensor3s::operator/=(real_t v_)
{
#ifdef STRONGCHECK
	Assert(
		mathdef::is_not_zero(v_),
		"suspicious absence of exception in Tensor3s::operator/=(real_t)");
#endif
	return operator*=(1. / v_);
}

inline Tensor3s& Tensor3s::CrossProduct(const Tensor1s& o_, Tensor3s& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}

inline Tensor3s& Tensor3s::DotProduct(const Tensor2s& o_, Tensor3s& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}

inline Tensor3s& Tensor3s::Dot2Product(const Tensor4s& o_, Tensor3s& result_) const
{
	return (result_ = *this).Dot2Multiply(o_);
}

//==========================Tensor4s==========================================

inline Tensor4s::Tensor4s(const SymmetricTensor4s& o_)
{
	operator=(o_);
}

inline small_t Tensor4::tIndex::operator[](small_t no_) const
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 4, "invalid Tensor4::tIndex[]const");
#endif
	return m_indices[--no_];
}

inline small_t& Tensor4::tIndex::operator[](small_t no_)
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 4, "invalid Tensor4::tIndex[]");
#endif
	return m_indices[--no_];
}

inline const real_t& Tensor4s::operator()(small_t i_, small_t j_, small_t k_, small_t l_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s::operator[](small_t,small_t,small_t,small_t)const");
#endif
	return (*m_rows[--i_])[--j_][--k_][--l_];
}

inline real_t& Tensor4s::operator()(small_t i_, small_t j_, small_t k_, small_t l_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s::operator[](small_t,small_t,small_t,small_t)");
#endif
	return (*m_rows[--i_])[--j_][--k_][--l_];
}

inline Tensor4s& Tensor4s::operator=(const Tensor4s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor4s::operator= (const Tensor4s&)");
#endif
	// if (&o_ != this)
	__Assign(o_.m_data);
	return *this;
}

inline Tensor4s& Tensor4s::operator/=(real_t v_)
{
#ifdef STRONGCHECK
	Assert(
		mathdef::is_not_zero(v_),
		"suspicious absence of exception in Tensor4s::operator/=(real_t)");
#endif
	return operator*=(1. / v_);
}

inline Tensor4s& Tensor4s::Assign0(small_t i_, small_t j_, small_t k_, small_t l_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s::Assign0(small_t,small_t,small_t,small_t)");
#endif
	(*m_rows[--i_])[--j_][--k_][--l_] = 0.;
	return *this;
}

inline Tensor4s& Tensor4s::CrossProduct(const Tensor1s& o_, Tensor4s& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}

inline Tensor2s& Tensor4s::Dot2Product(const Tensor2s& o_, Tensor2s& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}

inline Tensor4s& Tensor4s::DotProduct(const Tensor2s& o_, Tensor4s& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}

inline Tensor3s& Tensor4s::Dot2Product(const Tensor3s& o_, Tensor3s& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}

inline Tensor4s& Tensor4s::Dot2Product(const Tensor4s& o_, Tensor4s& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}

inline Tensor4s& Tensor4s::TransformAll(const Tensor2s& tr_)
{
	DotMultiply(tr_, 1, 1);
	DotMultiply(tr_, 2, 1);
	DotMultiply(tr_, 3, 1);
	DotMultiply(tr_);  //<=>4,1
	/* обратно (если tr_ - тензор поворота, т.е. tr_^T==tr_^-1):
//
	 * DotMultiply_left(tr_);//<=>right(1,2)
 DotMultiply(tr_,1,2);
 DotMultiply(tr_,2,2);

	 * DotMultiply(tr_,3,2);
 DotMultiply(tr_,4,2);
*/
	return *this;
}

//==========================SymmetricTensor4s===================================
//отображение индексов: 11->1, 22->2, 33->3, 12->4, 23->5, 13->6
inline unsigned short int map1(unsigned short int i_)
{
	switch (i_)
	{
		case 1:
			return 1;
		case 2:
			return 2;
		case 3:
			return 3;
		case 4:
			return 1;
		case 5:
			return 2;
		default:
			return 1;
	}
}
inline unsigned short int map2(unsigned short int i_)
{
	switch (i_)
	{
		case 1:
			return 1;
		case 2:
			return 2;
		case 3:
			return 3;
		case 4:
			return 2;
		case 5:
			return 3;
		default:
			return 3;
	}
}

inline SymmetricTensor4s& SymmetricTensor4s::TransformAll(const Tensor2s& /*ll*/)
{
	Tensor2s l;
	l(1, 1) = 1. / sqrt(2.);
	l(1, 2) = -1. / sqrt(2.);
	l(1, 3) = 0.;
	l(2, 1) = 1. / sqrt(2.);
	l(2, 2) = 1. / sqrt(2.);
	l(2, 3) = 0.;
	l(3, 1) = 0.;
	l(3, 2) = 0.;
	l(3, 3) = 1.;
	const Tensor2s& q(l);
	operator()(1, 1) = 165.64;
	operator()(1, 2) = 63.94;
	operator()(1, 3) = 63.94;
	operator()(1, 4) = 0.;
	operator()(1, 5) = 0.;
	operator()(1, 6) = 0.;

	operator()(2, 2) = 165.64;
	operator()(2, 3) = 63.94;
	operator()(2, 4) = 0.;
	operator()(2, 5) = 0.;
	operator()(2, 6) = 0.;

	operator()(3, 3) = 165.64;
	operator()(3, 4) = 0.;
	operator()(3, 5) = 0.;
	operator()(3, 6) = 0.;

	operator()(4, 4) = 79.51;
	operator()(4, 5) = 0.;
	operator()(4, 6) = 0.;

	operator()(5, 5) = 79.51;
	operator()(5, 6) = 0.;

	operator()(6, 6) = 79.51;

	std::ofstream out("transform_all.txt");
	out << "Rotation:\n";
	for (small_t i = 1, j; i <= 3; ++i)
	{
		for (j = 1; j <= 3; ++j) out << q(i, j) << '\t';
		out << '\n';
	}
	out << "Before rotation:\n";
	// out.close();
	for (small_t i = 1, j; i <= 6; ++i)
	{
		// out.open("transform_all.txt",std::ios::ate|std::ios::app);
		for (j = 1; j < i; ++j) out << '\t';
		for (j = i; j <= 6; ++j)
		{
			out << operator()(i, j);
			//    out.close();
			//    out.open("transform_all.txt",std::ios::ate|std::ios::app);
			out << '\t';
		}
		out << '\n';
	}
	// out.close();
	const SymmetricTensor4s A(*this);
	/* real_t
      t1 = q(3,1)*q(3,1),
      t2 = t1*t1,
		  t4 = A(3,6)*q(1,1),
		  t5 = A(3,5)*q(2,1),
		  t6 = t4+t5,
		  t7 = t1*q(3,1),
		  t11 = 4.0*A(5,5)+2.0*A(2,3),
		  t12 = q(2,1)*q(2,1),
		  t15 = A(3,4)+2.0*A(5,6),
		  t16 = t15*q(1,1),
		  t19 = q(1,1)*q(1,1),
		  t21 = A(1,3)+2.0*A(6,6),
		  t26 = t12*q(2,1),
		  t29 = 2.0*A(4,5)+A(2,6),
		  t30 = t29*q(1,1),
		  t33 = A(1,5)+2.0*A(4,6),
		  t36 = t19*q(1,1),
		  t40 = t12*t12,
		  t46 = A(1,2)+2.0*A(4,4),
		  t53 = t19*t19,
		  t56 = q(3,2)*q(3,2),
		  t57 = t56*A(3,3),
		  t58 = A(3,6)*q(1,2),
		  t59 = A(3,5)*q(2,2),
		  t60 = t58+t59,
		  t61 = 2.0*q(3,2)*t60,
		  t62 = q(1,2)*q(1,2),
		  t64 = q(2,2)*q(2,2),
		  t66 = A(3,4)*q(1,2),
		  t69 = t57+t61+t62*A(1,3)+A(2,3)*t64+2.0*q(2,2)*t66,
		  t72 = 2.0*t56*A(3,5),
		  t73 = A(5,5)*q(2,2),
		  t74 = A(5,6)*q(1,2),
		  t77 = A(4,5)*q(1,2),
		  t81 = 2.0*t64*A(2,5),
		  t84 = t72+4.0*q(3,2)*(t73+t74)+4.0*q(2,2)*t77+t81+2.0*t62*A(1,5),
		  t86 = t56*A(3,6),
		  t87 = A(6,6)*q(1,2),
		  t88 = A(5,6)*q(2,2),
		  t92 = A(4,6)*q(1,2),
		  t95 = t62*A(1,6),
		  t96 = t86+2.0*q(3,2)*(t87+t88)+t64*A(2,6)+2.0*q(2,2)*t92+t95,
		  t102 = A(2,5)*q(2,2),
		  t103 = A(2,6)*q(1,2),
		  t107 = t64*A(2,2),
		  t108 = A(2,4)*q(1,2),
		  t110 = 2.0*q(2,2)*t108,
		  t111 = t56*A(2,3)+2.0*q(3,2)*(t102+t103)+t62*A(1,2)+t107+t110,
		  t114 = A(4,5)*q(2,2),
		  t117 = A(4,4)*q(1,2),
		  t120 = t62*A(1,4),
		  t121 = t64*A(2,4),
		  t122 = t56*A(3,4)+2.0*q(3,2)*(t114+t92)+2.0*q(2,2)*t117+t120+t121,
		  t127 = q(2,2)*A(1,5),
		  t128 = A(1,6)*q(1,2),
		  t132 = A(1,4)*q(1,2),
		  t134 = 2.0*q(2,2)*t132,
		  t135 = t62*A(1,1),
		  t136 = A(1,3)*t56+2.0*q(3,2)*(t127+t128)+A(1,2)*t64+t134+t135,
		  t139 = t1*A(3,3),
		  t140 = 2.0*q(3,1)*t6,
		  t143 = A(3,4)*q(1,1),
		  t147 = q(3,3)*q(3,3),
		  t150 = 2.0*t1*A(3,5),
		  t151 = A(5,5)*q(2,1),
		  t152 = A(5,6)*q(1,1),
		  t155 = A(4,5)*q(1,1),
		  t159 = 2.0*t12*A(2,5),
		  t164 = t1*A(3,6),
		  t165 = A(6,6)*q(1,1),
		  t166 = A(5,6)*q(2,1),
		  t170 = A(4,6)*q(1,1),
		  t173 = t19*A(1,6),
		  t180 = A(2,5)*q(2,1),
		  t181 = A(2,6)*q(1,1),
		  t185 = t12*A(2,2),
		  t186 = A(2,4)*q(1,1),
		  t188 = 2.0*q(2,1)*t186,
		  t190 = q(2,3)*q(2,3),
		  t193 = q(2,1)*A(4,5),
		  t196 = A(4,4)*q(1,1),
		  t199 = t19*A(1,4),
		  t200 = t12*A(2,4),
		  t205 = q(1,3)*q(1,3),
		  t207 = q(2,1)*A(1,5),
		  t208 = A(1,6)*q(1,1),
		  t212 = A(1,4)*q(1,1),
		  t214 = 2.0*q(2,1)*t212,
		  t215 = t19*A(1,1),
		  t220 = t58+t59+A(3,3)*q(3,2),
		  t222 = A(3,5)*q(3,2),
		  t225 = A(2,3)+2.0*A(5,5),
		  t227 = t15*q(1,2),
		  t228 = 3.0*t222+q(2,2)*t225+t227,
		  t230 = A(3,6)*q(3,2),
		  t234 = 3.0*t230+t15*q(2,2)+t21*q(1,2),
		  t240 = t29*q(1,2),
		  t241 = q(3,2)*t225+3.0*t102+t240,
		  t246 = q(3,2)*t15+t29*q(2,2)+t33*q(1,2),
		  t253 = q(3,2)*t21+t33*q(2,2)+3.0*t128,
		  t259 = t108+q(2,2)*A(2,2)+A(2,5)*q(3,2),
		  t262 = q(2,2)*A(2,4),
		  t265 = q(3,2)*t29+3.0*t262+t46*q(1,2),
		  t271 = q(3,2)*t33+t46*q(2,2)+3.0*t132,
		  t277 = A(1,6)*q(3,2)+A(1,4)*q(2,2)+A(1,1)*q(1,2),
		  t280 = q(3,3)*t220,
		  t282 = q(2,2)*A(2,3)+t66+t222,
		  t286 = A(1,3)*q(1,2)+q(2,2)*A(3,4)+t230,
		  t290 = t74+t73+t222,
		  t293 = t77+t102+A(5,5)*q(3,2),
		  t296 = q(3,2)*A(5,6),
		  t297 = A(1,5)*q(1,2)+t114+t296,
		  t302 = t230+t88+t87,
		  t305 = q(2,2)*A(2,6)+t296+t92,
		  t308 = q(2,2)*A(4,6),
		  t309 = A(6,6)*q(3,2)+t308+t128,
		  t317 = t102+A(2,3)*q(3,2)+t103,
		  t319 = q(2,3)*t259,
		  t322 = q(3,2)*A(2,6)+t262+A(1,2)*q(1,2),
		  t327 = t114+A(3,4)*q(3,2)+t92,
		  t330 = A(4,5)*q(3,2)+t262+t117,
		  t334 = A(4,4)*q(2,2)+t132+A(4,6)*q(3,2),
		  t341 = t127+q(3,2)*A(1,3)+t128,
		  t345 = t132+A(1,2)*q(2,2)+A(1,5)*q(3,2),
		  t347 = t277*q(1,3),
		  t351 = A(3,6)*q(1,3),
		  t352 = A(3,5)*q(2,3),
		  t354 = t351+t352+A(3,3)*q(3,3),
		  t359 = t15*q(1,3),
		  t360 = 3.0*q(3,3)*A(3,5)+q(2,3)*t225+t359,
		  t361 = q(2,1)*t360,
		  t362 = A(3,6)*q(3,3),
		  t367 = t21*q(1,3),
		  t377 = t29*q(1,3),
		  t378 = q(3,3)*t225+3.0*A(2,5)*q(2,3)+t377,
		  t383 = q(3,3)*t15+t29*q(2,3)+t33*q(1,3),
		  t391 = q(3,3)*t21+t33*q(2,3)+3.0*q(1,3)*A(1,6),
		  t394 = q(3,1)*(t12*t378+2.0*q(2,1)*q(1,1)*t383+t391*t19),
		  t398 = A(2,5)*q(3,3)+q(1,3)*A(2,4)+A(2,2)*q(2,3),
		  t399 = t26*t398,
		  t401 = A(2,4)*q(2,3),
		  t403 = t46*q(1,3),
		  t411 = q(3,3)*t33+t46*q(2,3)+3.0*q(1,3)*A(1,4),
		  t413 = q(2,1)*t19*t411,
		  t417 = q(1,3)*A(1,1)+A(1,6)*q(3,3)+q(2,3)*A(1,4),
		  t418 = t417*t36,
		  t422 = t56*q(3,2),
		  t431 = t64*q(2,2),
		  t436 = t62*q(1,2),
		  t440 = t64*t64,
		  t450 = t62*t62,
		  t466 = t5+t4+A(3,3)*q(3,1),
		  t468 = A(3,5)*q(3,1),
		  t471 = 3.0*t468+q(2,1)*t225+t16,
		  t473 = q(3,1)*A(3,6),
		  t477 = 3.0*t473+q(2,1)*t15+t21*q(1,1),
		  t483 = q(3,1)*t225+3.0*t180+t30,
		  t488 = q(3,1)*t15+q(2,1)*t29+t33*q(1,1),
		  t495 = q(3,1)*t21+q(2,1)*t33+3.0*t208,
		  t501 = A(2,5)*q(3,1)+A(2,2)*q(2,1)+t186,
		  t505 = (A(2,6)+2.0*A(4,5))/3.0,
		  t507 = A(2,4)*q(2,1),
		  t508 = t46*q(1,1),
		  t517 = q(3,1)*t33+q(2,1)*t46+3.0*t212,
		  t523 = q(3,1)*A(1,6)+q(1,1)*A(1,1)+A(1,4)*q(2,1),
		  t530 = 3.0*t362+t15*q(2,3)+t367,
		  t544 = q(3,3)*t505+t401+t403/3.0,
		  t552 = q(3,3)*t466,
		  t562 = t152+t468+t151,
		  t565 = t155+A(5,5)*q(3,1)+t180,
		  t567 = q(3,1)*A(5,6),
		  t574 = t473+t165+t166,
		  t579 = A(4,6)*q(2,1),
		  t581 = t579+q(3,1)*A(6,6)+t208,
		  t591 = q(2,3)*t501,
		  t602 = t507+A(4,5)*q(3,1)+t196,
		  t606 = A(4,4)*q(2,1)+q(3,1)*A(4,6)+t212,
		  t619 = t523*q(1,3),
		  t626 = t147*q(3,3),
		  t635 = t190*q(2,3),
		  t640 = t205*q(1,3),
		  t644 = t190*t190,
		  t655 = t205*t205,
		  t658 = q(3,1)*t220,
		  t678 = q(2,1)*t259,
		  t691 = t277*q(1,1),
		  t740 = t57+t61+2.0*q(2,2)*t74+A(5,5)*t64+t62*A(6,6),
		  t742 = A(5,5)+A(2,3),
		  t744 = A(3,4)+A(5,6),
		  t745 = t744*q(1,2),
		  t749 = A(2,6)+A(4,5),
		  t750 = t749*q(1,2),
		  t755 = t72+q(3,2)*(2.0*q(2,2)*t742+2.0*t745)+t81+2.0*q(2,2)*t750+2.0*t62*A(4,6),
		  t757 = t744*q(2,2),
		  t758 = A(1,3)+A(6,6),
		  t759 = t758*q(1,2),
		  t763 = A(4,6)+A(1,5),
		  t764 = t763*q(1,2),
		  t766 = t86+q(3,2)*(t757+t759)+t64*A(4,5)+q(2,2)*t764+t95,
		  t775 = t56*A(5,5)+2.0*q(3,2)*(t102+t77)+t110+t107+t62*A(4,4),
		  t778 = t749*q(2,2),
		  t781 = A(4,4)+A(1,2),
		  t782 = t781*q(1,2),
		  t784 = t56*A(5,6)+q(3,2)*(t778+t764)+t121+q(2,2)*t782+t120,
		  t792 = t56*A(6,6)+2.0*q(3,2)*(t128+t308)+t134+A(4,4)*t64+t135,
		  t798 = t56*(t552+q(2,3)*t562+t574*q(1,3)),
		  t801 = t744*q(1,1),
		  t806 = t749*q(1,1),
		  t809 = q(3,1)*t744,
		  t810 = q(2,1)*t749,
		  t815 =
	 q(2,2)*(q(3,3)*(2.0*t468+q(2,1)*t742+t801)+q(2,3)*(q(3,1)*t742+2.0*t180+t806)+q(1,3)*(t809+t810+2.0*t170)),
		  t817 = q(2,1)*t744,
		  t818 = t758*q(1,1),
		  t822 = t763*q(1,1),
		  t829 = (q(3,1)*t758+q(2,1)*t763+2.0*t208)*q(1,3),
		  t837 = t64*(q(3,3)*t565+t591+t602*q(1,3)),
		  t843 = t781*q(1,1),
		  t850 = q(1,3)*(q(3,1)*t763+q(2,1)*t781+2.0*t212),
		  t857 = t62*(q(3,3)*t581+q(2,3)*t606+t619),
		  t865 = 2.0*t222+q(2,2)*t742+t745,
		  t869 = q(3,2)*t742+2.0*t102+t750,
		  t871 = q(3,2)*t744,
		  t873 = t871+t778+2.0*t92,
		  t878 = 2.0*t230+t757+t759,
		  t881 = t871+2.0*t114+t764,
		  t886 = q(3,2)*t758+t763*q(2,2)+2.0*t128,
		  t897 = 2.0*t296+t778+t764,
		  t901 = q(3,2)*t749+2.0*t262+t782,
		  t906 = q(3,2)*t763+t781*q(2,2)+2.0*t132;

	 operator()(1,1) =
		  t2*A(3,3)+4.0*t7*t6+t1*(t12*t11+4.0*q(2,1)*t16+2.0*t21*t19)
		  +4.0*q(3,1)*(t26*A(2,5)+t12*t30+q(2,1)*t33*t19+t36*A(1,6))+t40*A(2,2)
		  +4.0*q(1,1)*t26*A(2,4)+2.0*t12*t46*t19+4.0*q(2,1)*t36*A(1,4)+t53*A(1,1);
	 operator()(1,2) =
		  t1*t69+q(3,1)*(q(2,1)*t84+2.0*t96*q(1,1))+t12*t111+2.0*q(2,1)*t122*q(1,1)+t136*t19;
	 operator()(1,3) =
		  t147*(t139+t140+t19*A(1,3)+A(2,3)*t12+2.0*q(2,1)*t143)+q(3,3)
		  *(q(2,3)*(t150+4.0*q(3,1)*(t151+t152)+4.0*q(2,1)*t155+t159+2.0*t19*A(1,5))
		  +2.0*(t164+2.0*q(3,1)*(t165+t166)+t12*A(2,6)+2.0*q(2,1)*t170+t173)*q(1,3))
		  +t190*(t1*A(2,3)+2.0*q(3,1)*(t180+t181)+t19*A(1,2)+t185+t188)+2.0*q(2,3)
		  *(t1*A(3,4)+2.0*q(3,1)*(t193+t170)+2.0*q(2,1)*t196+t199+t200)*q(1,3)+(A(1,3)
		  *t1+2.0*q(3,1)*(t207+t208)+A(1,2)*t12+t214+t215)*t205;
	 operator()(1,4) =
		  t7*t220+t1*(q(2,1)*t228+q(1,1)*t234)
		  +q(3,1)*(t12*t241+2.0*q(2,1)*q(1,1)*t246+t19*t253)+t26*t259+t12*t265*q(1,1)+q(2,1)*t19*t271+t277*t36;
	 operator()(1,5) =
		  t1*(t280+q(2,3)*t282+t286*q(1,3))+q(3,1)*(q(2,1)*(2.0*q(3,3)
		  *t290+2.0*q(2,3)*t293+2.0*t297*q(1,3))+2.0*q(1,1)*(q(3,3)*t302+q(2,3)*t305
		  +t309*q(1,3)))+t12*(q(3,3)*t317+t319+t322*q(1,3))+2.0*q(2,1)*(q(3,3)*t327
		  +q(2,3)*t330+t334*q(1,3))*q(1,1)+t19*(q(3,3)*t341+q(2,3)*t345+t347);
	 operator()(1,6) =
		  t7*t354+t1*(t361+t530*q(1,1))+t394+t399+3.0*t12*q(1,1)*t544+t413+t418;
	 operator()(2,2) =
		  t56*t56*A(3,3)+4.0*t422*t60+t56*(t64*t11+4.0*q(2,2)*t227+2.0*t21*t62)
		  +4.0*q(3,2)*(t431*A(2,5)+t64*t240+q(2,2)*t33*t62+t436*A(1,6))+t440*A(2,2)
		  +4.0*t431*t108+2.0*t64*t46*t62+4.0*q(2,2)*t436*A(1,4)+t450*A(1,1);
	 operator()(2,3) =
		  t147*t69+q(3,3)*(q(2,3)*t84+2.0*q(1,3)*t96)+t190*t111+2.0*q(2,3)*q(1,3)*t122+t205*t136;
	 operator()(2,4) =
		  t422*t466+t56*(q(2,2)*t471+t477*q(1,2))+q(3,2)*(t64*t483+2.0*q(2,2)*q(1,2)*t488+t62*t495)
		  +t431*t501+3.0*t64*(q(3,1)*t505+t507+t508/3.0)*q(1,2)+q(2,2)*t517*t62+t523*t436;
	 operator()(2,5) =
		  t422*t354+t56*(q(2,2)*t360+t530*q(1,2))+q(3,2)*(t64*t378+2.0*q(2,2)*q(1,2)*t383+t62*t391)
		  +t431*t398+3.0*t64*t544*q(1,2)+q(2,2)*t411*t62+t417*t436;
	 operator()(2,6) =
		  t56*(t552+q(2,3)*(t468+A(2,3)*q(2,1)+t143)+(A(3,4)*q(2,1)+q(1,1)
		  *A(1,3)+t473)*q(1,3))+q(3,2)*(q(2,2)*(2.0*q(3,3)*t562+2.0*q(2,3)*t565
		  +2.0*(t567+q(1,1)*A(1,5)+t193)*q(1,3))+2.0*(q(3,3)*t574+q(2,3)
		  *(A(2,6)*q(2,1)+t170+t567)+t581*q(1,3))*q(1,2))+t64*(q(3,3)*(t181+t180+A(2,3)*q(3,1))
		  +t591+(t507+q(3,1)*A(2,6)+q(1,1)*A(1,2))*q(1,3))+2.0*q(2,2)
		  *(q(3,3)*(t170+t193+q(3,1)*A(3,4))+q(2,3)*t602+t606*q(1,3))*q(1,2)+t62*(q(3,3)
		  *(t208+t207+q(3,1)*A(1,3))+q(2,3)*(A(1,2)*q(2,1)+A(1,5)*q(3,1)+t212)+t619);
	 operator()(3,3) =
		  t147*t147*A(3,3)+4.0*t626*(t351+t352)+t147*(t190*t11+4.0*q(2,3)*t359+2.0*t21*t205)
		  +4.0*q(3,3)*(t635*A(2,5)+t190*t377+q(2,3)*t33*t205+t640*A(1,6))+t644*A(2,2)
		  +4.0*q(1,3)*t635*A(2,4)+2.0*t190*t46*t205+4.0*q(2,3)*t640*A(1,4)+t655*A(1,1);
	 operator()(3,4) =
		  t147*(t658+q(2,1)*t282+t286*q(1,1))+q(3,3)*(q(2,3)*(2.0*q(3,1)*t290
		  +2.0*q(2,1)*t293+2.0*t297*q(1,1))+2.0*q(1,3)*(q(3,1)*t302+q(2,1)
		  *t305+t309*q(1,1)))+t190*(q(3,1)*t317+t678+t322*q(1,1))+2.0*q(2,3)*q(1,3)
		  *(q(3,1)*t327+q(2,1)*t330+t334*q(1,1))+t205*(q(3,1)*t341+q(2,1)*t345+t691);
	 operator()(3,5) =
		  t626*t220+t147*(q(2,3)*t228+q(1,3)*t234)+q(3,3)*(t190*t241+2.0
		  *q(2,3)*q(1,3)*t246+t253*t205)+t635*t259+t190*q(1,3)*t265+q(2,3)*t271*t205+t277*t640;
	 operator()(3,6) =
		  t626*t466+t147*(q(2,3)*t471+q(1,3)*t477)+q(3,3)*(t190*t483+2.0*q(2,3)*q(1,3)*t488+t495*t205)
		  +t635*t501+t190*q(1,3)*(q(3,1)*t29+3.0*t507+t508)+q(2,3)*t517*t205+t523*t640;
	 operator()(4,4) =
		  t1*t740+q(3,1)*(q(2,1)*t755+2.0*t766*q(1,1))+t12*t775+2.0*q(2,1)*q(1,1)*t784+t19*t792;
	 operator()(4,5) =
		  t798+q(3,2)*(t815+q(1,2)*(q(3,3)*(2.0*t473+t817+t818)+q(2,3)*(t809+2.0*t193+t822)+t829))
		  +t837+q(2,2)*q(1,2)*(q(3,3)*(2.0*t567+t810+t822)+q(2,3)*(q(3,1)*t749+2.0*t507+t843)+t850)+t857;
	 operator()(4,6) =
		  t1*(t280+q(2,3)*t290+t302*q(1,3))+q(3,1)*(q(2,1)*(q(3,3)*t865
		  +q(2,3)*t869+q(1,3)*t873)+q(1,1)*(q(3,3)*t878+q(2,3)*t881+q(1,3)*t886))
		  +t12*(q(3,3)*t293+t319+t330*q(1,3))+q(2,1)*q(1,1)*(q(3,3)*t897+q(2,3)*t901
		  +q(1,3)*t906)+t19*(q(3,3)*t309+q(2,3)*t334+t347);
	 operator()(5,5) =
		  t147*t740+q(3,3)*(q(2,3)*t755+2.0*q(1,3)*t766)+t190*t775+2.0*q(2,3)*q(1,3)*t784+t205*t792;
	 operator()(5,6) =
		  t147*(t658+q(2,1)*t290+t302*q(1,1))+q(3,3)*(q(2,3)*(q(3,1)*t865+q(2,1)*t869
		  +t873*q(1,1))+(q(3,1)*t878+q(2,1)*t881+q(1,1)*t886)*q(1,3))+t190*(q(3,1)*t293
		  +t678+t330*q(1,1))+q(2,3)*(q(3,1)*t897+q(2,1)*t901+q(1,1)*t906)*q(1,3)+t205*(q(3,1)*t309+q(2,1)*t334+t691);
	 operator()(6,6) =
		  t147*(t139+t140+2.0*q(1,1)*t166+t19*A(6,6)+A(5,5)*t12)
		  +q(3,3)*(q(2,3)*(t150+q(3,1)*(2.0*q(2,1)*t742+2.0*t801)+t159+2.0*q(2,1)*t806
		  +2.0*t19*A(4,6))+2.0*q(1,3)*(t164+q(3,1)*(t817+t818)+t12*A(4,5)+q(2,1)*t822+t173))
		  +t190*(t1*A(5,5)+2.0*q(3,1)*(t155+t180)+t185+t19*A(4,4)+t188)
		  +2.0*q(2,3)*q(1,3)*(t1*A(5,6)+q(3,1)*(t810+t822)+t200+q(2,1)*t843+t199)
		  +t205*(t1*A(6,6)+2.0*q(3,1)*(t579+t208)+t214+A(4,4)*t12+t215);*/
	static real_t g[6][6];
	g[0][0] = l(1, 1) * l(1, 1);
	g[0][1] = l(1, 2) * l(1, 2);
	g[0][2] = l(1, 3) * l(1, 3);
	g[0][3] = 2. * l(1, 2) * l(1, 3);
	g[0][4] = 2. * l(1, 3) * l(1, 1);
	g[0][5] = 2. * l(1, 2) * l(1, 1);
	g[1][0] = l(2, 1) * l(2, 1);
	g[1][1] = l(2, 2) * l(2, 2);
	g[1][2] = l(2, 3) * l(2, 3);
	g[1][3] = 2. * l(2, 3) * l(2, 2);
	g[1][4] = 2. * l(2, 3) * l(2, 1);
	g[1][5] = 2. * l(2, 2) * l(2, 1);
	g[2][0] = l(3, 1) * l(3, 1);
	g[2][1] = l(3, 2) * l(3, 2);
	g[2][2] = l(3, 3) * l(3, 3);
	g[2][3] = 2. * l(3, 3) * l(3, 2);
	g[2][4] = 2. * l(3, 3) * l(3, 1);
	g[2][5] = 2. * l(3, 2) * l(3, 1);
	g[3][0] = l(3, 1) * l(2, 1);
	g[3][1] = l(3, 2) * l(2, 2);
	g[3][2] = l(3, 3) * l(2, 3);
	g[3][3] = l(3, 3) * l(2, 2) + l(3, 2) * l(2, 3);
	g[3][4] = l(3, 3) * l(2, 1) + l(3, 1) * l(2, 3);
	g[3][5] = l(3, 1) * l(2, 2) + l(3, 2) * l(2, 1);
	g[4][0] = l(3, 1) * l(1, 1);
	g[4][1] = l(3, 2) * l(1, 2);
	g[4][2] = l(3, 3) * l(1, 3);
	g[4][3] = l(3, 3) * l(1, 2) + l(3, 2) * l(1, 3);
	g[4][4] = l(3, 3) * l(1, 1) + l(3, 1) * l(1, 3);
	g[4][5] = l(3, 1) * l(1, 2) + l(3, 2) * l(1, 1);
	g[5][0] = l(2, 1) * l(1, 1);
	g[5][1] = l(1, 2) * l(2, 2);
	g[5][2] = l(1, 3) * l(2, 3);
	g[5][3] = l(1, 3) * l(2, 2) + l(1, 2) * l(2, 3);
	g[5][4] = l(1, 1) * l(2, 3) + l(1, 3) * l(2, 1);
	g[5][5] = l(1, 1) * l(2, 2) + l(1, 2) * l(2, 1);

	/* g[0][0]=l(1,1)*l(1,1); g[0][1]=l(1,2)*l(1,2); g[0][2]=l(1,3)*l(1,3);
	 * g[0][4]=          2.*l(1,2)*l(1,3); g[0][5]=          2.*l(1,3)*l(1,1);
	 * g[0][3]=          2.*l(1,2)*l(1,1);
 g[1][0]=l(2,1)*l(2,1); g[1][1]=l(2,2)*l(2,2);
	 * g[1][2]=l(2,3)*l(2,3); g[1][4]=          2.*l(2,3)*l(2,2);
	 * g[1][5]=          2.*l(2,3)*l(2,1); g[1][3]=          2.*l(2,2)*l(2,1);

	 * g[2][0]=l(3,1)*l(3,1); g[2][1]=l(3,2)*l(3,2); g[2][2]=l(3,3)*l(3,3);
	 * g[2][4]=          2.*l(3,3)*l(3,2); g[2][5]=          2.*l(3,3)*l(3,1);
	 * g[2][3]=          2.*l(3,2)*l(3,1);
 g[4][0]=l(3,1)*l(2,1); g[4][1]=l(3,2)*l(2,2);
	 * g[4][2]=l(3,3)*l(2,3); g[4][4]=l(3,3)*l(2,2)+l(3,2)*l(2,3);
	 * g[4][5]=l(3,3)*l(2,1)+l(3,1)*l(2,3); g[4][3]=l(3,1)*l(2,2)+l(3,2)*l(2,1);

	 * g[5][0]=l(3,1)*l(1,1); g[5][1]=l(3,2)*l(1,2); g[5][2]=l(3,3)*l(1,3);
	 * g[5][4]=l(3,3)*l(1,2)+l(3,2)*l(1,3); g[5][5]=l(3,3)*l(1,1)+l(3,1)*l(1,3);
	 * g[5][3]=l(3,1)*l(1,2)+l(3,2)*l(1,1);
 g[3][0]=l(2,1)*l(1,1); g[3][1]=l(1,2)*l(2,2);
	 * g[3][2]=l(1,3)*l(2,3); g[3][4]=l(1,3)*l(2,2)+l(1,2)*l(2,3);
	 * g[3][5]=l(1,1)*l(2,3)+l(1,3)*l(2,1); g[3][3]=l(1,1)*l(2,2)+l(1,2)*l(2,1);
*/
	for (small_t i_ = 0, i1_ = 1; i_ < 6; ++i_, ++i1_)
		for (small_t l_ = i_, l1_ = i1_; l_ < 6; ++l_, ++l1_)
			for (small_t j_ = 0, j1_ = 1; j_ < 6; ++j_, ++j1_)
				for (small_t k_ = 0, k1_ = 1; k_ < 6; ++k_, ++k1_)
					operator()(i1_, l1_) = A(j1_, k1_) * g[i_][j_] * g[l_][k_];

	/*Tensor4s T4;//(*this);
for (small_t i=1,j,k,l; i<=3; ++i)
for (j=1; j<=3; ++j)
for (k=1; k<=3;
	 * ++k)
for (l=1; l<=3; ++l)
  T4(i,j,k,l) = operator()(i,j,k,l);
T4.TransformAll(q);
for
	 * (small_t i=1,j; i<=6; ++i)
for (j=i; j<=6; ++j)
  operator()(i,j) =
	 * T4(map1(i),map2(i),map1(j),map2(j));
*/
	// out.open("transform_all.txt",std::ios::ate|std::ios::app);
	out << "After rotation:\n";
	// out.close();
	for (small_t i = 1, j; i <= 6; ++i)
	{
		// out.open("transform_all.txt",std::ios::ate|std::ios::app);
		for (j = 1; j < i; ++j) out << '\t';
		for (j = i; j <= 6; ++j)
		{
			out << operator()(i, j);
			//    out.close();
			//    out.open("transform_all.txt",std::ios::ate|std::ios::app);
			out << '\t';
		}
		out << '\n';
	}
	out.close();
	return *this;
}

inline real_t SymmetricTensor4s::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 6 && j_ <= 6,
		"invalid index in SymmetricTensor4s::operator(small_t,small_t)const");
#endif
	if (--i_ > --j_)
		std::swap(i_, j_);
	return m_data[(13 - i_) * i_ / 2 + j_ - i_];
}

inline real_t& SymmetricTensor4s::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ <= 6 && i_ <= j_,
		tMessage("invalid index in SymmetricTensor4s::operator(small_t,small_t), must be 0<")
			<< i_ << "<=" << j_ << "<=6");
#endif
	return m_data[(13 - (i_ - 1)) * (i_ - 1) / 2 + j_ - i_];
}

inline real_t SymmetricTensor4s::operator()(small_t i_, small_t j_, small_t k_, small_t l_) const
{  //отображение индексов: 11->1, 22->2, 33->3, 12->4, 23->5, 13->6
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in SymmetricTensor4s::operator(small_t i_,small_t,small_t i_,small_t)const");
#endif
	return operator()(
		i_ == j_ ? i_ : (25 - 3 * (i_ + j_)) * (i_ + j_) / 2 - 20,
		k_ == l_ ? k_ : (25 - 3 * (k_ + l_)) * (k_ + l_) / 2 - 20);
}

inline real_t SymmetricTensor4s::operator[](Tensor4::tIndex i_) const
{  //отображение индексов: 00->1, 11->2, 22->3, 01->4, 12->5, 02->6
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[2] >= 0 && i_[3] >= 0 && i_[4] >= 0 && i_[1] < 3 && i_[2] < 3 &&
			i_[3] < 3 && i_[4] < 3,
		"invalid index in SymmetricTensor4s::operator[]const");
#endif
	return operator()(
		i_[1] == i_[2] ? (i_[1] + 1) : (13 - 3 * (i_[1] + i_[2])) * (i_[1] + i_[2]) / 2 - 1,
		i_[3] == i_[4] ? (i_[3] + 1) : (13 - 3 * (i_[3] + i_[4])) * (i_[3] + i_[4]) / 2 - 1);
}

inline real_t Tensor4s::operator[](Tensor4::tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[2] >= 0 && i_[3] >= 0 && i_[4] >= 0 && i_[1] < 3 && i_[2] < 3 &&
			i_[3] < 3 && i_[4] < 3,
		"invalid index in Tensor4s::operator[](Index4)const");
#endif
	return (*m_rows[i_[1]])[i_[2]][i_[3]][i_[4]];
}

inline real_t& Tensor4s::operator[](Tensor4::tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[2] >= 0 && i_[3] >= 0 && i_[4] >= 0 && i_[1] < 3 && i_[2] < 3 &&
			i_[3] < 3 && i_[4] < 3,
		"invalid index in Tensor4s::operator[](Index4)");
#endif
	return (*m_rows[i_[1]])[i_[2]][i_[3]][i_[4]];
}

//#undef __SET_ROWS
#undef __SET_ROWS3
#undef __SET_ROWS4
