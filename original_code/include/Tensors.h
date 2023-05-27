// Author: Chekhov Vladimir Valerevich

#ifndef TensorsH
#define TensorsH
//------------------------------------------------------------------------------
#include <cstddef>	//NULL
#include <algorithm>
#include <functional>
#include <fstream>
#include <cmath>
//------------------------------------------------------------------------------
#include "NumTypes.h"
#include "math_constants.hpp"

#ifdef STRONGCHECK
#include "ThrowMessage.h"
#endif
//------------------------------------------------------------------------------
enum COMPONENT
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
enum class TENSOR_KIND
{
	ZERO,
	UNIT
};
template <typename T>
class Tensor1a;
template <typename T>
class Tensor2s;
template <typename T>
class Tensor2a;
template <typename T>
class Tensor3s;
template <typename T>
class Tensor4s;
//------------------------------------------------------------------------------
template <typename T>
class Tensor1s
{
	friend class Tensor1a<T>;
	friend class Tensor2s<T>;
	friend class Tensor3s<T>;
	friend class Tensor4s<T>;

private:
	T Data[3];
	void __Assign(T v_)
	{
		Data[0] = Data[1] = Data[2] = v_;
	}
	void __Assign(const T* const d_)
	{
		Data[0] = d_[0];
		Data[1] = d_[1];
		Data[2] = d_[2];
	}
	void __Assign(const Tensor1s<T>& o_)
	{
		Data[0] = o_.Data[0];
		Data[1] = o_.Data[1];
		Data[2] = o_.Data[2];
	}

public:
	const T& operator[](small_t i_) const
	{
		return Data[i_];
	}  // is public because need for fFunArgAdapter<1,Targ,Tresult,FUN> for integration
	T& operator[](small_t i_)
	{
		return Data[i_];
	}
	// public:
	Tensor1s<T>() : Data{}
	{
	}
	Tensor1s<T>(TENSOR_KIND, small_t, T = 1.);	 // Basis vector multiplied by factor_
	Tensor1s<T>(float x_, float y_, float z_)
	{
		Data[0] = x_;
		Data[1] = y_;
		Data[2] = z_;
	}
	Tensor1s<T>(double x_, double y_, double z_)
	{
		Data[0] = x_;
		Data[1] = y_;
		Data[2] = z_;
	}
	Tensor1s<T>(long double x_, long double y_, long double z_)
	{
		Data[0] = x_;
		Data[1] = y_;
		Data[2] = z_;
	}
	explicit Tensor1s<T>(T v_)
	{
		__Assign(v_);
	}
	explicit Tensor1s<T>(const T* const);
	Tensor1s<T>(const Tensor1s<T>& o_)
	{
		__Assign(o_);
	}
	explicit Tensor1s<T>(const Tensor1a<T>&);
	const T& operator()(small_t) const;
	T& operator()(small_t);
	Tensor1s<T>& operator=(const Tensor1a<T>& o_);
	Tensor1s<T>& operator=(const Tensor1s<T>& o_)
	{
		__Assign(o_);
		return *this;
	}
	Tensor1s<T>& operator-=(const Tensor1s<T>& o_)
	{
		Data[0] -= o_[0];
		Data[1] -= o_[1];
		Data[2] -= o_[2];
		return *this;
	}
	Tensor1s<T>& operator+=(const Tensor1s<T>& o_)
	{
		Data[0] += o_[0];
		Data[1] += o_[1];
		Data[2] += o_[2];
		return *this;
	}
	Tensor1s<T>& operator*=(T v_)
	{
		Data[0] *= v_;
		Data[1] *= v_;
		Data[2] *= v_;
		return *this;
	}
	Tensor1s<T>& operator/=(T);
	Tensor1s<T> operator-() const
	{
		return Tensor1s<T>(-Data[0], -Data[1], -Data[2]);
	}
	Tensor1s<T> operator+(const Tensor1s<T>& o_) const
	{
		return Tensor1s<T>(Data[0] + o_[0], Data[1] + o_[1], Data[2] + o_[2]);
	}
	Tensor1s<T> operator-(const Tensor1s<T>& o_) const
	{
		return Tensor1s<T>(Data[0] - o_[0], Data[1] - o_[1], Data[2] - o_[2]);
	}
	Tensor1s<T> operator*(T v_) const
	{
		return Tensor1s<T>(Data[0] * v_, Data[1] * v_, Data[2] * v_);
	}
	bool operator==(const Tensor1s<T>& o_) const
	{
		return mathdef::is_eq(Data[0], o_[0]) && mathdef::is_eq(Data[1], o_[1]) &&
			   mathdef::is_eq(Data[2], o_[2]);
	}  // is used in "DxfStream.cpp"
	Tensor1s<T>& Negate()
	{
		Data[0] = -Data[0];
		Data[1] = -Data[1];
		Data[2] = -Data[2];
		return *this;
	}
	T Length() const
	{
		return sqrt(Data[0] * Data[0] + Data[1] * Data[1] + Data[2] * Data[2]);
	}
	Tensor1s<T>& Normalize()
	{
		return operator/=(Length());
	}
	Tensor1s<T>& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor1s<T>& Assign1(small_t i_, T v_ = 1.)
	{
		__Assign(0.);
		Data[--i_] = v_;
		return *this;
	}
	Tensor1s<T>& Assign(T x_, T y_, T z_)
	{
		Data[0] = x_;
		Data[1] = y_;
		Data[2] = z_;
		return *this;
	}
	Tensor1s<T>& Assign(small_t i_, T v_)
	{
		Data[--i_] = v_;
		return *this;
	}  // temporarily - to be removed
	   //   Tensor1& Assign(component_t, T);
	   //   Tensor1& Add(component_t, T);
	Tensor1s<T>& CrossMultiplyByBasisVector(small_t);
	Tensor1s<T>& CrossMultiplyByBasisVector_left(small_t);
	Tensor1s<T>& CrossProductWithBasisVector(small_t, Tensor1s<T>&) const;
	Tensor1s<T>& CrossProductWithBasisVector_left(small_t, Tensor1s<T>&) const;
	T DotProduct(const Tensor1s<T>&) const;
	Tensor1s<T>& CrossMultiply(const Tensor1s<T>&);
	Tensor1s<T>& CrossMultiply_left(const Tensor1s<T>&);
	Tensor1s<T>& CrossProduct(const Tensor1s<T>&, Tensor1s<T>&) const;
	Tensor2s<T>& DirectProduct(const Tensor1s<T>&, Tensor2s<T>&) const;
	Tensor2a<T>& DirectProduct(const Tensor1s<T>&, Tensor2a<T>&) const;
	Tensor1s<T>& DotMultiply(const Tensor2s<T>&);
	Tensor1s<T>& DotMultiply_left(const Tensor2s<T>&);
	Tensor1s<T>& DotProduct(const Tensor2s<T>&, Tensor1s<T>&) const;
	Tensor2s<T>& CrossProduct(const Tensor2s<T>&, Tensor2s<T>&) const;
	Tensor3s<T>& DirectProduct(const Tensor2s<T>&, Tensor3s<T>&) const;
	Tensor2s<T>& DotProduct(const Tensor3s<T>&, Tensor2s<T>&) const;
	Tensor3s<T>& CrossProduct(const Tensor3s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& DirectProduct(const Tensor3s<T>&, Tensor4s<T>&) const;
	Tensor3s<T>& DotProduct(const Tensor4s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& CrossProduct(const Tensor4s<T>&, Tensor4s<T>&) const;
	friend std::istream& operator>>(std::istream& is_, Tensor1s<T>& t_);
	friend std::ostream& operator<<(std::ostream& os_, const Tensor1s<T>& t_);
};

template <class T>
std::istream& operator>>(std::istream& is, Tensor1s<T>& t)
{
	is >> t.Data[0] >> t.Data[1] >> t.Data[2];
	return is;
}

template <class T>
std::ostream& operator<<(std::ostream& os, const Tensor1s<T>& t)
{
	os << t.Data[0] << ' ' << t.Data[1] << ' ' << t.Data[2];
	return os;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template <typename T>
class RetTensor1;
template <typename T>
class RetTensor2;
//------------------------------------------------------------------------------
template <typename T>
class Tensor2a;

template <typename T>
class Tensor1a
{
	friend class Tensor1s<T>;
	friend class RetTensor1<T>;
	friend class Tensor2a<T>;
	// friend std::ostream& operator<< (std::ostream& out_, const Tensor1a<T>& t_);

protected:
	Tensor1s<T>* pData;

private:
	Tensor1s<T>* __pData()
	{
		return is0() ? (pData = new Tensor1s<T>) : pData;
	}
	bool __AdjustAdd(const Tensor1a<T>& o_)
	{
		return is0() ? (o_.is0() ? false : (pData = new Tensor1s<T>(*o_.pData), false))
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
	const T operator[](small_t i_) const
	{
		return is0() ? 0. : pData->operator[](i_);
	}

public:
	//////////////////////// temporarily!!! - only for Node::LinkDisplacement
	const Tensor1s<T>* Data() const
	{
		return pData;
	}
	////////////////////////
	Tensor1a<T>() : pData(NULL)
	{
	}
	Tensor1a<T>(T x_, T y_, T z_) : pData(new Tensor1s<T>(x_, y_, z_))
	{
	}
	explicit Tensor1a<T>(const Tensor1s<T>& o_) : pData(new Tensor1s<T>(o_))
	{
	}
	Tensor1a<T>(const Tensor1a<T>& o_)
		: pData(o_.is0() ? (Tensor1s<T>*)NULL : new Tensor1s<T>(*o_.pData))
	{
	}
	Tensor1a<T>(RetTensor1<T>&);
	virtual ~Tensor1a<T>()
	{
		delete pData;
	}
	bool is0() const
	{
		return pData == NULL;
	}
	bool is0(small_t) const
	{
		return is0();
	}
	const T& operator()(small_t) const;
	T& operator()(small_t);
	Tensor1a<T>& operator=(RetTensor1<T>&);
	Tensor1a<T>& operator=(const Tensor1a<T>& o_)
	{
		if (__AdjustAdd(o_))
			pData->operator=(*o_.pData);
		return *this;
	}
	Tensor1a<T>& operator-=(const Tensor1a<T>&);
	Tensor1a<T>& operator+=(const Tensor1a<T>& o_)
	{
		if (__AdjustAdd(o_))
			pData->operator+=(*o_.pData);
		return *this;
	}
	Tensor1a<T>& operator*=(T v_)
	{
		return is0() ? *this : (pData->operator*=(v_), *this);
	}
	Tensor1a<T>& operator/=(T);
	RetTensor1<T> operator-() const;
	RetTensor1<T> operator-(const Tensor1a<T>&) const;
	RetTensor1<T> operator+(const Tensor1a<T>&) const;
	RetTensor1<T> operator*(T) const;
	RetTensor1<T> operator/(T) const;
	Tensor1a<T>& Negate()
	{
		return is0() ? *this : (pData->Negate(), *this);
	}
	T Length() const
	{
		return is0() ? 0. : pData->Length();
	}
	Tensor1a<T>& Normalize();
	Tensor1a<T>& Assign0()
	{
		delete pData;
		pData = NULL;
		return *this;
	}
	Tensor1a<T>& Assign1(small_t i_, T v_ = 1.)
	{
		if (is0())
			pData = new Tensor1s<T>(true, i_, v_);
		else
			pData->Assign1(i_, v_);
		return *this;
	}
	Tensor1a<T>& Assign(small_t i_, T v_)
	{
		if (is0())
			pData = new Tensor1s<T>(0.);
		pData->Data[--i_] = v_;
		return *this;
	}
	Tensor1a<T>& Assign(T x_, T y_, T z_)
	{
		is0() ? pData = new Tensor1s<T>(x_, y_, z_) : &pData->Assign(x_, y_, z_);
		return *this;
	}
	Tensor1a<T>& Add(small_t i_, T v_)
	{
		if (is0())
			pData = new Tensor1s<T>(0.);
		pData->Data[--i_] += v_;
		return *this;
	}
	Tensor1a<T>& CrossMultiplyByBasisVector(small_t i_)
	{
		return is0() ? *this : (pData->CrossMultiplyByBasisVector(i_), *this);
	}
	Tensor1a<T>& CrossMultiplyByBasisVector_left(small_t i_)
	{
		return is0() ? *this : (pData->CrossMultiplyByBasisVector_left(i_), *this);
	}
	Tensor1a<T>& CrossProductWithBasisVector(small_t i_, Tensor1a<T>& r_) const
	{
		return (r_ = *this).CrossMultiplyByBasisVector(i_);
	}
	Tensor1a<T>& CrossProductWithBasisVector_left(small_t i_, Tensor1a<T>& r_) const
	{
		return (r_ = *this).CrossMultiplyByBasisVector_left(i_);
	}
	RetTensor1<T> CrossProductWithBasisVector(small_t) const;
	RetTensor1<T> CrossProductWithBasisVector_left(small_t) const;
	T DotProduct(const Tensor1a<T>& o_) const
	{
		return is0() || o_.is0() ? 0. : pData->DotProduct(*o_.pData);
	}
	Tensor1a<T>& CrossMultiply(const Tensor1a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->CrossMultiply(*o_.pData);
		return *this;
	}
	Tensor1a<T>& CrossMultiply_left(const Tensor1a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->CrossMultiply_left(*o_.pData);
		return *this;
	}
	Tensor1a<T>& CrossProduct(const Tensor1a<T>& o_, Tensor1a<T>& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).CrossMultiply(o_);
	}
	RetTensor1<T> CrossProduct(const Tensor1a<T>&) const;
	Tensor2a<T>& DirectProduct(const Tensor1a<T>&, Tensor2a<T>&) const;
	RetTensor2<T> DirectProduct(const Tensor1a<T>&) const;
	Tensor1a<T>& DotMultiply(const Tensor2a<T>&);
	Tensor1a<T>& DotMultiply_left(const Tensor2a<T>&);
	Tensor1a<T>& DotProduct(const Tensor2a<T>&, Tensor1a<T>&) const;
	RetTensor1<T> DotProduct(const Tensor2a<T>&) const;
	Tensor2a<T>& CrossProduct(const Tensor2a<T>&, Tensor2a<T>&) const;
	RetTensor2<T> CrossProduct(const Tensor2a<T>&) const;
	Tensor3s<T>& DirectProduct(const Tensor2a<T>&, Tensor3s<T>&) const;
	Tensor2a<T>& DotProduct(const Tensor3s<T>&, Tensor2a<T>&) const;
	RetTensor2<T> DotProduct(const Tensor3s<T>&) const;
	Tensor3s<T>& CrossProduct(const Tensor3s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& DirectProduct(const Tensor3s<T>&, Tensor4s<T>&) const;
	Tensor3s<T>& DotProduct(const Tensor4s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& CrossProduct(const Tensor4s<T>&, Tensor4s<T>&) const;
};
// inline std::ostream& operator<< (std::ostream& out_, const Tensor1s<T>& t_) {return out_<<"|
// "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" |\n";} inline std::ostream& operator<< (std::ostream& out_,
// const Tensor1a<T>& t_) {return out_<<"| "<<t_[0]<<" "<<t_[1]<<" "<<t_[2]<<" |\n";}
//------------------------------------------------------------------------------
template <typename T>
class RetTensor1
{
	friend class Tensor1a<T>;
	friend class Tensor2a<T>;
	mutable Tensor1s<T>* pData;
	RetTensor1<T>() : pData(NULL)
	{
	}
	RetTensor1<T>(Tensor1a<T>& o_) : pData(o_.pData)
	{
		o_.pData = NULL;
	}
	//   RetTensor1<T>(const RetTensor1<T>& o_): pData(o_.pData) {o_.pData = NULL;}
	//   explicit RetTensor1<T>(const T* const p_): pData(p_==NULL ? NULL : new Tensor1s<T>(p_)) {}
	RetTensor1<T>(T x_, T y_, T z_) : pData(new Tensor1s<T>(x_, y_, z_))
	{
	}

public:
	~RetTensor1<T>()
	{
		delete pData;
	}
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// std::ostream& operator<< (std::ostream&, const BasicTensor2&);
//------------------------------------------------------------------------------
// class Tensor4;
template <typename T>
class SymmetricTensor4s;
template <typename T>
class SymmetricTensor2s;
//------------------------------------------------------------------------------
template <typename T>
class Tensor2s
{
	friend class Tensor1s<T>;
	friend class Tensor2a<T>;
	friend class Tensor3s<T>;
	friend class Tensor4s<T>;

private:
	struct tIndex
	{
	private:
		small_t first, second;

	public:
		tIndex() : first(0), second(0)
		{
		}
		tIndex(small_t i_, small_t j_) : first(i_), second(j_)
		{
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);
		//        Index&   operator()(small_t i_, small_t j_){first=i_; second=j_; return *this;}
	};

	T Data[9];
	T* /*const*/ Rows[3];

#define __SET_ROWS      \
	*Rows = Data;       \
	Rows[1] = Data + 3; \
	Rows[2] = Data + 6;
	//   void __set_rows(){const_cast<T*>(*Rows) = Data; const_cast<T*>(Rows[1]) = Data+3;
	//   const_cast<T*>(Rows[2]) = Data+6;}
	const T* operator[](small_t i_) const
	{
		return Rows[i_];
	}
	T* operator[](small_t i_)
	{
		return Rows[i_];
	}
	void __Assign(T v_)
	{
		std::fill(Data, Data + 9, v_);
	}
	void __Assign(const T* v_)
	{
		std::copy(v_, v_ + 9, Data);
	}
	template <typename T>
	void __Assign2D(const T v_)
	{
		for (small_t i = 0, j; i < 3; ++i)
			for (j = 0; j < 3; ++j) Rows[i][j] = v_[i][j];
	}
	template <typename T>
	void __Assign2DTranspose(const T v_)
	{
		for (small_t i = 0, j; i < 3; ++i)
			for (j = 0; j < 3; ++j) Rows[i][j] = v_[j][i];
	}
	const T& operator[](tIndex) const;
	T& operator[](tIndex);
	const T& operator()(tIndex) const;
	T& operator()(tIndex);
	void MaxAbsIndex(small_t&, small_t&) const;

public:
	Tensor2s<T>() : Data{}, Rows{}
	{
		__SET_ROWS;
	}
	Tensor2s<T>(TENSOR_KIND, T v_ = 0.) : Data{}, Rows{}
	{
		__SET_ROWS;
		Assign1(v_);
	}
	explicit Tensor2s<T>(T v_)
	{
		__SET_ROWS;
		__Assign(v_);
	}
	explicit Tensor2s<T>(const T (*v_)[3])
	{
		__SET_ROWS;
		__Assign2D(v_);
	}
	explicit Tensor2s<T>(const T** v_)
	{
		__SET_ROWS;
		__Assign2D(v_);
	}
	Tensor2s<T>(const T (*v_)[3], bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(v_);
	}
	Tensor2s<T>(const T** v_, bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(v_);
	}
	Tensor2s<T>(const Tensor1s<T>&, const Tensor1s<T>&, const Tensor1s<T>&);
	Tensor2s<T>(const Tensor1s<T>&, const Tensor1s<T>&, const Tensor1s<T>&, bool);
	Tensor2s<T>(const Tensor2s<T>& o_)
	{
		__SET_ROWS;
		__Assign(o_.Data);
	}
	Tensor2s<T>(const Tensor2s<T>& o_, bool)
	{
		__SET_ROWS;
		__Assign2DTranspose(o_.Rows);
	}
	const T& operator()(small_t, small_t) const;
	T& operator()(small_t, small_t);
	Tensor2s<T>& operator=(const Tensor2s<T>&);
	Tensor2s<T>& operator-=(const Tensor2s<T>&);
	Tensor2s<T>& operator+=(const Tensor2s<T>&);
	Tensor2s<T>& operator*=(T v_);
	// Tensor2s<T>& operator*=(T v_)
	//{
	//	std::for_each(Data, Data + 9, fMilt<T, T>(v_));
	//	return *this;
	//}
	Tensor2s<T> operator*(T v_) const
	{
		return Tensor2s<T>(*this) *= v_;
	}  // tmp to remove
	Tensor2s<T>& operator/=(T);
	Tensor2s<T>& Negate()
	{
		std::transform(Data, Data + 9, Data, std::negate<T>());
		return *this;
	}
	Tensor2s<T>& Transpose();
	Tensor2s<T>& AddTransposed(const Tensor2s<T>&);
	SymmetricTensor2s<T>& Symmetrization_twofold(SymmetricTensor2s<T>&) const;
	Tensor2s<T>& AssignRows(const Tensor1s<T>&, const Tensor1s<T>&, const Tensor1s<T>&);
	Tensor2s<T>& AssignCols(const Tensor1s<T>&, const Tensor1s<T>&, const Tensor1s<T>&);
	Tensor2s<T>& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor2s<T>& Assign1(T v_ = 1.)
	{
		__Assign(0.);
		Data[0] = Data[4] = Data[8] = v_;
		return *this;
	}
	T Contraction() const
	{
		return Data[0] + Data[4] + Data[8];
	}
	T Det() const;
	T I1() /*Linear    invariant*/ const
	{
		return Contraction();
	}
	T I2() /*Quadratic invariant*/ const;
	T I3() /*Cubic     invariant*/ const
	{
		return Det();
	}
	T Minor(small_t, small_t) const;
	Tensor2s<T>& Invert();
	Tensor1s<T>& AttachedVector(Tensor1s<T>&) const;
	Tensor1s<T> AttachedVector() const;
	//   Tensor1s<T>& RotationVector(Tensor1s<T>&) const;//for orthogonal tensors only
	Tensor1s<T>& ScalarProductWithBasisVector(small_t, Tensor1s<T>&) const;
	Tensor1s<T>& ScalarProductWithBasisVector_left(small_t, Tensor1s<T>&) const;
	Tensor1s<T> ScalarProductWithBasisVector(small_t) const;
	Tensor1s<T> ScalarProductWithBasisVector_left(small_t) const;
	void EigenValues(T[3]) const;
	small_t EigenVectors(
		T /*current e.value*/, Tensor2s<T>& /*e.vects*/, small_t = 1 /*begin position*/) const;
	void Eigen(T[3], Tensor2s<T>&) const;
	Tensor2s<T>& EigenVectors(Tensor2s<T>& r_) const
	{
		T evs[3];
		Eigen(evs, r_);
		return r_;
	}
	Tensor2s<T>& RotateToEigenSystem(Tensor2s<T>&);
	Tensor2s<T>& AssignDev();
	Tensor2s<T>& Dev(Tensor2s<T>& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	Tensor1s<T>& DotProduct(const Tensor1s<T>&, Tensor1s<T>&) const;
	Tensor2s<T>& CrossMultiply(const Tensor1s<T>&);
	Tensor2s<T>& CrossMultiply_left(const Tensor1s<T>&);
	Tensor2s<T>& CrossProduct(const Tensor1s<T>&, Tensor2s<T>&) const;
	Tensor3s<T>& DirectProduct(const Tensor1s<T>&, Tensor3s<T>&) const;
	T Dot2Product(const Tensor2s<T>&) const;
	Tensor2s<T>& DotMultiply(const Tensor2s<T>&);
	Tensor2s<T>& DotMultiply_left(const Tensor2s<T>&);
	Tensor2s<T>& DotProduct(const Tensor2s<T>&, Tensor2s<T>&) const;
	Tensor3s<T>& CrossProduct(const Tensor2s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& DirectProduct(const Tensor2s<T>&, Tensor4s<T>&) const;
	Tensor1s<T>& Dot2Product(const Tensor3s<T>&, Tensor1s<T>&) const;
	Tensor3s<T>& DotProduct(const Tensor3s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& CrossProduct(const Tensor3s<T>&, Tensor4s<T>&) const;
	//   Tensor2s<T>& Dot2Multiply      (const Tensor4s<T>&,small_t,small_t);
	//   Tensor2s<T>& Dot2Multiply      (const Tensor4&,small_t,small_t);
	Tensor2s<T>& Dot2Multiply(const SymmetricTensor4s<T>&, small_t, small_t);
	Tensor2s<T>& Dot2Multiply(const Tensor4s<T>&);
	Tensor2s<T>& Dot2Multiply_left(const Tensor4s<T>&);
	Tensor2s<T>& Dot2Product(const Tensor4s<T>&, Tensor2s<T>&) const;
	Tensor4s<T>& DotProduct(const Tensor4s<T>&, Tensor4s<T>&) const;
	//   Tensor2s<T>& ApplyFunction(const fClassicFunction&);
};
//------------------------------------------------------------------------------
template <typename T>
class RetTensor2;
template <typename T>
class Tensor2CaptureDataAdapter;
//------------------------------------------------------------------------------
template <typename T>
class Tensor2a
{
	friend class Tensor1s<T>;
	friend class Tensor1a<T>;
	friend class RetTensor2<T>;

protected:
	Tensor2s<T>* pData;

private:
	Tensor2s<T>* __pData()
	{
		return is0() ? (pData = new Tensor2s<T>) : pData;
	}
	bool __AdjustAdd(const Tensor2a<T>& o_)
	{
		return is0() ? (o_.is0() ? false : (pData = new Tensor2s<T>(*o_.pData), false))
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
	const T* operator[](small_t i_) const
	{
		return pData->operator[](i_);
	}

protected:
	//   Tensor2a<T>(Tensor2s<T>* addr_): pData(addr_) {}
public:
	Tensor2a<T>() : pData(NULL)
	{
	}
	explicit Tensor2a<T>(const Tensor2s<T>& o_) : pData(new Tensor2s<T>(o_))
	{
	}
	Tensor2a<T>(const Tensor2a<T>& o_)
		: pData(o_.is0() ? (Tensor2s<T>*)NULL : new Tensor2s<T>(*o_.pData))
	{
	}
	Tensor2a<T>(RetTensor2<T>&);
	virtual ~Tensor2a<T>()
	{
		delete pData;
	}
	bool is0() const
	{
		return pData == NULL;
	}
	bool is0(small_t, small_t) const
	{
		return pData == NULL;
	}
	const T& operator()(small_t, small_t) const;
	T& operator()(small_t, small_t);
	Tensor2a<T>& operator=(RetTensor2<T>&);
	Tensor2a<T>& operator=(const Tensor2s<T>& o_)
	{
		if (is0())
			pData = new Tensor2s<T>(o_);
		else
			pData->operator=(o_);
		return *this;
	}
	Tensor2a<T>& operator=(const Tensor2a<T>& o_)
	{
		if (__AdjustAdd(o_))
			pData->operator=(*o_.pData);
		return *this;
	}
	Tensor2a<T>& operator-=(const Tensor2a<T>&);
	Tensor2a<T>& operator+=(const Tensor2a<T>& o_)
	{
		if (__AdjustAdd(o_))
			pData->operator+=(*o_.pData);
		return *this;
	}
	Tensor2a<T>& operator+=(Tensor2CaptureDataAdapter<T>);
	Tensor2a<T>& operator*=(T v_)
	{
		return is0() ? *this : (pData->operator*=(v_), *this);
	}
	Tensor2a<T>& operator/=(T);
	RetTensor2<T> operator-() const;
	RetTensor2<T> operator-(const Tensor2a<T>&) const;
	RetTensor2<T> operator+(const Tensor2a<T>&) const;
	RetTensor2<T> operator*(T) const;
	RetTensor2<T> operator/(T) const;
	T Contraction() const
	{
		return is0() ? 0. : pData->Contraction();
	}
	T Det() const
	{
		return is0() ? 0. : pData->Det();
	}
	T I1() const
	{
		return Contraction();
	}
	T I2() const
	{
		return is0() ? 0. : pData->I2();
	}
	T I3() const
	{
		return Det();
	}
	T Minor(small_t i_, small_t j_) const
	{
		return is0() ? 0. : pData->Minor(i_, j_);
	}
	Tensor2a<T>& Negate()
	{
		return is0() ? *this : (pData->Negate(), *this);
	}
	Tensor2a<T>& Transpose()
	{
		return is0() ? *this : (pData->Transpose(), *this);
	}
	Tensor2a<T>& AddTransposed(const Tensor2a<T>& o_)
	{
		return is0() ? (o_.is0() ? *this : (pData = new Tensor2s<T>(*o_.pData, false), *this))
					 : (o_.is0() ? *this : (pData->AddTransposed(*o_.pData), *this));
	}
	Tensor2a<T>& AddTransposed(Tensor2CaptureDataAdapter<T>);
	Tensor2a<T>& Invert();
	Tensor2a<T>& Assign0()
	{
		delete pData;
		pData = NULL;
		return *this;
	}
	Tensor2a<T>& Assign0toData()
	{
		pData->Assign0();
		return *this;
	}
	Tensor2a<T>& Assign1(T v_ = 1.)
	{
		if (is0())
			pData = new Tensor2s<T>(TENSOR_KIND::UNIT, v_);
		else
			pData->Assign1(v_);
		return *this;
	}
	//   Tensor2& AssignRows(const Tensor1& r1_, const Tensor1& r2_, const Tensor1& r3_) {return
	//   r1_.is0()&&r2_.is0()&&r3_.is0()? (Assign0(),*this) :
	//   (AllocData(),pData->AssignRows(r1_,r2_,r3_),*this);} Tensor2& AssignCols(const Tensor1&
	//   r1_, const Tensor1& r2_, const Tensor1& r3_) {return r1_.is0()&&r2_.is0()&&r3_.is0()?
	//   (Assign0(),*this) : (AllocData(),pData->AssignCols(r1_,r2_,r3_),*this);} Tensor1s<T>&
	//   EigenValues(Tensor1s<T>& r_) const {return is0()? r_.Assign0() : pData->EigenValues(r_);}
	//   small_t EigenVectors(T ev_, Tensor2& evc_, small_t p_=1) const {return is0()?
	//   (evc_.Assign0(),3) : pData->EigenVectors(ev_,evc_,p_);} Tensor2&
	//   RotateToEigenSystem(Tensor2& r_) {return is0()? (r_.Assign1(),*this) :
	//   (pData->RotateToEigenSystem(r_),*this);} Tensor1& AttachedVector(Tensor1& r_) const {return
	//   is0()? r_.Assign0() : pData->AttachedVector(r_);}
	Tensor2a<T>& AssignDev()
	{
		return is0() ? *this : (pData->AssignDev(), *this);
	}
	Tensor2a<T>& Dev(Tensor2a<T>& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	//   Tensor1& ScalarProductWithBasisVector     (small_t i_, Tensor1& r_) const {return is0()?
	//   r_.Assign0() : pData->ScalarProductWithBasisVector(i_,r_);} Tensor1&
	//   ScalarProductWithBasisVector_left(small_t i_, Tensor1& r_) const {return is0()?
	//   r_.Assign0() : pData->ScalarProductWithBasisVector_left(i_,r_);}
	Tensor1a<T>& DotProduct(const Tensor1a<T>& o_, Tensor1a<T>& r_) const
	{
		return is0() ? r_.Assign0() : (r_ = o_).DotMultiply_left(*this);
	}
	RetTensor1<T> DotProduct(const Tensor1a<T>& o_) const
	{
		Tensor1a<T> result;
		return DotProduct(o_, result);
	}
	Tensor2a<T>& CrossMultiply(const Tensor1a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->CrossMultiply(*o_.pData);
		return *this;
	}
	Tensor2a<T>& CrossMultiply_left(const Tensor1a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->CrossMultiply_left(*o_.pData);
		return *this;
	}
	Tensor2a<T>& CrossProduct(const Tensor1a<T>& o_, Tensor2a<T>& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).CrossMultiply(o_);
	}
	RetTensor2<T> CrossProduct(const Tensor1a<T>&) const;
	Tensor3s<T>& DirectProduct(const Tensor1a<T>&, Tensor3s<T>&) const;
	T Dot2Product(const Tensor2a<T>& o_) const
	{
		return is0() || o_.is0() ? 0. : pData->Dot2Product(*o_.pData);
	}
	Tensor2a<T>& DotMultiply(const Tensor2a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->DotMultiply(*o_.pData);
		return *this;
	}
	Tensor2a<T>& DotMultiply_left(const Tensor2a<T>& o_)
	{
		if (__AdjustMult(o_))
			pData->DotMultiply_left(*o_.pData);
		return *this;
	}
	Tensor2a<T>& DotProduct(const Tensor2a<T>& o_, Tensor2a<T>& r_) const
	{
		return o_.is0() ? r_.Assign0() : (r_ = *this).DotMultiply(o_);
	}
	RetTensor2<T> DotProduct(const Tensor2a<T>&) const;
	Tensor3s<T>& CrossProduct(const Tensor2a<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& DirectProduct(const Tensor2a<T>&, Tensor4s<T>&) const;
	Tensor1a<T>& Dot2Product(const Tensor3s<T>&, Tensor1a<T>&) const;
	RetTensor1<T> Dot2Product(const Tensor3s<T>& o_) const
	{
		Tensor1a<T> result;
		return Dot2Product(o_, result);
	}
	Tensor3s<T>& DotProduct(const Tensor3s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& CrossProduct(const Tensor3s<T>&, Tensor4s<T>&) const;
	//   Tensor2a<T>& Dot2Multiply      (const Tensor4&,small_t,small_t);
	Tensor2a<T>& Dot2Multiply(const SymmetricTensor4s<T>&, small_t, small_t);
	Tensor2a<T>& Dot2Multiply(const Tensor4s<T>&);
	Tensor2a<T>& Dot2Multiply_left(const Tensor4s<T>&);
	Tensor2a<T>& Dot2Product(const Tensor4s<T>&, Tensor2a<T>&) const;
	RetTensor2<T> Dot2Product(const Tensor4s<T>&) const;
	Tensor4s<T>& DotProduct(const Tensor4s<T>&, Tensor4s<T>&) const;
	//   Tensor2& ApplyFunction(const fClassicFunction& fun_) {AllocData(0.);
	//   pData->ApplyFunction(fun_); return *this;}
};
// std::ostream& operator<< (std::ostream&, const Tensor2a<T>&);
//------------------------------------------------------------------------------
template <typename T>
class RetTensor2
{
	friend class Tensor1a<T>;
	friend class Tensor2a<T>;
	mutable Tensor2s<T>* pData;
	RetTensor2<T>() : pData(NULL)
	{
	}
	RetTensor2<T>(Tensor2a<T>& o_) : pData(o_.pData)
	{
		o_.pData = NULL;
	}
	//   RetTensor2<T>(const RetTensor2<T>& o_): pData(o_.pData) {o_.pData = NULL;}
	//   RetTensor2<T>(const T* const p_): pData(p_==NULL ? NULL : new Tensor1s<T>(p_)) {}
	//   RetTensor2<T>(T x_, T y_, T z_): pData(new Tensor1s<T>(x_,y_,z_)) {}
public:
	bool is0() const
	{
		return pData == NULL;
	}
	~RetTensor2<T>()
	{
		delete pData;
	}
	void SwapData(Tensor2s<T>*);
};
//------------------------------------------------------------------------------
template <typename T>
class Tensor2CaptureDataAdapter
{
	friend class Tensor2a<T>;
	Tensor2a<T>& Ref;

public:
	explicit Tensor2CaptureDataAdapter(Tensor2a<T>& o_) : Ref(o_)
	{
	}
};
//------------------------------------------------------------------------------
template <typename T>
class SymmetricTensor2s
{
	friend class Tensor2s<T>;
	T Data[6];	// 11 22 33 12 23 31

	void __Assign(T v_)
	{
		std::fill(Data, Data + 6, v_);
	}
	void __Assign(const T* v_)
	{
		std::copy(v_, v_ + 6, Data);
	}

public:
	SymmetricTensor2s<T>() : Data{}
	{
	}
	////   Tensor2s<T>(tensor_kind_t, T v_=0.) {__SET_ROWS; Assign1(v_);}
	explicit SymmetricTensor2s<T>(T v_)
	{
		__Assign(v_);
	}
	SymmetricTensor2s<T>(const SymmetricTensor2s<T>& o_)
	{
		__Assign(o_.Data);
	}
	T operator()(small_t, small_t) const;
	//         T& operator()(small_t, small_t);
	//   SymmetricTensor2s<T>& operator=(const SymmetricTensor2s<T>&);
	//   SymmetricTensor2s<T>& operator-=(const SymmetricTensor2s<T>&);
	SymmetricTensor2s<T>& operator+=(const SymmetricTensor2s<T>& o_)
	{
		Data[0] += o_.Data[0];
		Data[1] += o_.Data[1];
		Data[2] += o_.Data[2];
		Data[3] += o_.Data[3];
		Data[4] += o_.Data[4];
		Data[5] += o_.Data[5];
		return *this;
	}
	SymmetricTensor2s<T>& operator*=(T v_);
	// SymmetricTensor2s<T>& operator*=(T v_)
	//{
	//	std::for_each(Data, Data + 6, fMilt<T, T>(v_));
	//	return *this;
	//}
	//   SymmetricTensor2s<T>& operator/=(T);
	SymmetricTensor2s<T> operator*(T v_)
	{
		SymmetricTensor2s<T> res(*this);
		return res *= v_;
	}
	//   SymmetricTensor2s<T>& Negate() {std::transform(Data, Data+6, Data, std::negate<T>());
	//   return *this;}
	////   SymmetricTensor2s<T>& Transpose() {return *this;}
	SymmetricTensor2s<T>& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	SymmetricTensor2s<T>& Assign1(T v_ = 1.)
	{
		Data[0] = Data[1] = Data[2] = v_;
		Data[3] = Data[4] = Data[5] = 0.;
		return *this;
	}
	T Contraction() const
	{
		return Data[0] + Data[1] + Data[2];
	}
	//   T Det() const;
	T I1() /*Linear    invariant*/ const
	{
		return Contraction();
	}
	//   T I2()/*Quadratic invariant*/ const;
	//   T I3()/*Cubic     invariant*/ const {return Det();}
	//   T Minor(small_t,small_t) const;
	//   SymmetricTensor2s<T>& Invert();
	////   Tensor1s<T>& AttachedVector(Tensor1s<T>&) const;
	////   Tensor1s<T>  AttachedVector() const;
	////   Tensor1s<T>& ScalarProductWithBasisVector     (small_t, Tensor1s<T>&) const;
	////   Tensor1s<T>& ScalarProductWithBasisVector_left(small_t, Tensor1s<T>&) const;
	////   Tensor1s<T>  ScalarProductWithBasisVector     (small_t) const;
	////   Tensor1s<T>  ScalarProductWithBasisVector_left(small_t) const;
	//   void EigenValues(T[3]) const;
	////   small_t EigenVectors(T/*current e.value*/, SymmetricTensor2s<T>&/*e.vects*/,
	/// small_t=1/*begin position*/) const; /   void Eigen(T[3],SymmetricTensor2s<T>&) const; /
	/// SymmetricTensor2s<T>& EigenVectors(SymmetricTensor2s<T>& r_) const {T evs[3]; Eigen(evs,r_);
	/// return r_;} /   SymmetricTensor2s<T>& RotateToEigenSystem(SymmetricTensor2s<T>&);
	SymmetricTensor2s<T>& AssignDev();
	SymmetricTensor2s<T>& Dev(SymmetricTensor2s<T>& dev_) const
	{
		return (dev_ = *this).AssignDev();
	}
	//   Tensor1s<T>& DotProduct      (const Tensor1s<T>&, Tensor1s<T>&) const;
	////   SymmetricTensor2s<T>& CrossMultiply     (const Tensor1s<T>&);
	////   SymmetricTensor2s<T>& CrossMultiply_left(const Tensor1s<T>&);
	////   SymmetricTensor2s<T>& CrossProduct      (const Tensor1s<T>&, SymmetricTensor2s<T>&) const;
};
//------------------------------------------------------------------------------
template <typename T>
inline T SymmetricTensor2s<T>::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3,
		"invalid index in SymmetricTensor2s<T>::operator() const");
#endif
	return (i_ == j_) ? Data[--i_] : Data[(25 - 3 * (i_ + j_)) * (i_ + j_) / 2 - 21];
}
//------------------------------------------------------------------------------
template <typename T>
inline SymmetricTensor2s<T>& SymmetricTensor2s<T>::AssignDev()
{
	const T p = I1() / 3.l;
	Data[0] -= p;
	Data[1] -= p;
	Data[2] -= p;
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
class Tensor3s
{
	friend class Tensor1s<T>;
	friend class Tensor2s<T>;
	friend class Tensor4s<T>;

private:
	T Data[27];
	//   T (*Rows[3])[3][3];
	typedef T (*pArr33_t)[3][3];
	typedef T Arr33_t[3][3];
	pArr33_t Rows[3];
#define __SET_ROWS3                                                \
	*Rows = static_cast<pArr33_t>(static_cast<void*>(Data));       \
	Rows[1] = static_cast<pArr33_t>(static_cast<void*>(Data + 9)); \
	Rows[2] = static_cast<pArr33_t>(static_cast<void*>(Data + 18));
	void __Assign(T v_)
	{
		std::fill(Data, Data + 27, v_);
	}
	void __Assign(const T* v_)
	{
		std::copy(v_, v_ + 27, Data);
	}
	const Arr33_t& operator[](small_t i_) const
	{
		return *Rows[i_];
	}
	Arr33_t& operator[](small_t i_)
	{
		return *Rows[i_];
	}
	struct tIndex
	{
	private:
		small_t first, second, third;

	public:
		tIndex() : first(0), second(0), third(0)
		{
		}
		tIndex(small_t i_, small_t j_, small_t k_) : first(i_), second(j_), third(k_)
		{
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);
		tIndex& operator()(small_t i_, small_t j_, small_t k_)
		{
			first = i_;
			second = j_;
			third = k_;
			return *this;
		}
	};

public:
	Tensor3s<T>() : Data{}, Rows{} {__SET_ROWS3};
	Tensor3s<T>(TENSOR_KIND, T v_ = 1.) : Data{}, Rows{}
	{
		__SET_ROWS3;
		Assign1(v_);
	}  // Levy-Chivita
	explicit Tensor3s<T>(T v_)
	{
		__SET_ROWS3;
		__Assign(v_);
	}
	Tensor3s<T>(const Tensor3s<T>& o_)
	{
		__SET_ROWS3;
		__Assign(o_.Data);
	}
	T operator[](const tIndex&) const;
	// T& operator[](const tIndex&);
	T operator()(small_t, small_t, small_t) const;
	T& operator()(small_t, small_t, small_t);
	Tensor3s<T>& operator=(const Tensor3s<T>&);
	Tensor3s<T>& operator-=(const Tensor3s<T>&);
	Tensor3s<T>& operator+=(const Tensor3s<T>&);
	Tensor3s<T>& operator*=(T v_);
	// Tensor3s<T>& operator*=(T v_)
	//{
	//	std::for_each(Data, Data + 27, fMilt<T, T>(v_));
	//	return *this;
	//}
	Tensor3s<T>& operator/=(T);
	Tensor3s<T>& Negate()
	{
		std::transform(Data, Data + 27, Data, std::negate<T>());
		return *this;
	}
	bool is0() const
	{
		return false;
	}
	Tensor3s<T>& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor3s<T>& Assign1(T = 1.);	// Levy-Chivita
	Tensor1s<T>& Contraction(small_t, small_t, Tensor1s<T>&) const;
	Tensor2s<T>& ScalarProductWithBasisVector(small_t, Tensor2s<T>&) const;
	Tensor2s<T>& ScalarProductWithBasisVector_left(small_t, Tensor2s<T>&) const;
	Tensor2s<T>& DotProduct(const Tensor1s<T>&, Tensor2s<T>&) const;
	Tensor3s<T>& CrossMultiply(const Tensor1s<T>&);
	Tensor3s<T>& CrossMultiply_left(const Tensor1s<T>&);
	Tensor3s<T>& CrossProduct(const Tensor1s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& DirectProduct(const Tensor1s<T>&, Tensor4s<T>&) const;
	Tensor1s<T>& Dot2Product(const Tensor2s<T>&, Tensor1s<T>&) const;
	Tensor3s<T>& DotMultiply(const Tensor2s<T>&);
	Tensor3s<T>& DotMultiply_left(const Tensor2s<T>&);
	Tensor3s<T>& DotProduct(const Tensor2s<T>&, Tensor3s<T>&) const;
	T Dot3Product(const Tensor3s<T>&) const;
	Tensor2s<T>& Dot2Product(const Tensor3s<T>&, Tensor2s<T>&) const;
	Tensor4s<T>& DotProduct(const Tensor3s<T>&, Tensor4s<T>&) const;
	Tensor1s<T>& Dot3Product(const Tensor4s<T>&, Tensor1s<T>&) const;
	Tensor3s<T>& Dot2Multiply(const Tensor4s<T>&);
	Tensor3s<T>& Dot2Multiply_left(const Tensor4s<T>&);
	Tensor3s<T>& Dot2Product(const Tensor4s<T>&, Tensor3s<T>&) const;
};
//------------------------------------------------------------------------------
/*ostream& operator<< (ostream& out_, Tensor3s<T>& t_)
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
//------------------------------------------------------------------------------
template <typename T>
class Tensor4
{
	friend class Tensor2s<T>;

protected:
	struct tIndex
	{
	private:
		small_t indices[4];

	public:
		tIndex() : indices{}
		{
		}
		tIndex(const small_t i_, const small_t j_, const small_t k_, const small_t l_)
		{
			indices[0] = i_;
			indices[1] = j_;
			indices[2] = k_;
			indices[3] = l_;
		}
		small_t operator[](small_t) const;
		small_t& operator[](small_t);
	};
	//  virtual T operator[](tIndex) const =0;
public:
	bool is0() const
	{
		return false;
	}
};
//------------------------------------------------------------------------------
template <typename T>
class SymmetricTensor4s;
//------------------------------------------------------------------------------
template <typename T>
class Tensor4s : public Tensor4<T>
{
	typedef typename Tensor4<T>::tIndex tIndex;

	friend class Tensor1s<T>;
	friend class Tensor2s<T>;
	friend class Tensor3s<T>;

private:
	T Data[81];
	//   T (*Rows[3])[3][3][3];
	typedef T (*pArr333_t)[3][3][3];
	typedef T Arr333_t[3][3][3];
	pArr333_t Rows[3];
#define __SET_ROWS4                                                  \
	*Rows = static_cast<pArr333_t>(static_cast<void*>(Data));        \
	Rows[1] = static_cast<pArr333_t>(static_cast<void*>(Data + 27)); \
	Rows[2] = static_cast<pArr333_t>(static_cast<void*>(Data + 54));
	void __Assign(T v_)
	{
		std::fill(Data, Data + 81, v_);
	}
	void __Assign(const T* v_)
	{
		std::copy(v_, v_ + 81, Data);
	}
	const Arr333_t& operator[](small_t i_) const
	{
		return *Rows[i_];
	}
	Arr333_t& operator[](small_t i_)
	{
		return *Rows[i_];
	}
	/*virtual*/ T operator[](tIndex) const;
	/*virtual*/ T& operator[](tIndex);

public:
	Tensor4s<T>() : Data{}, Rows{} {__SET_ROWS4} Tensor4s<T>(TENSOR_KIND, small_t kind_, T v_ = 1.)
	{
		__SET_ROWS4;
		Assign1(kind_, v_);
	}  // four isotropic tensors
	explicit Tensor4s<T>(T v_)
	{
		__SET_ROWS4;
		__Assign(v_);
	}
	explicit Tensor4s<T>(const SymmetricTensor4s<T>&);
	Tensor4s<T>(const Tensor4s<T>& o_)
	{
		__SET_ROWS4;
		__Assign(o_.Data);
	}
	const T& operator()(small_t, small_t, small_t, small_t) const;
	T& operator()(small_t, small_t, small_t, small_t);
	//   Tensor4s<T>& operator=  (const Tensor4&);
	Tensor4s<T>& operator=(const Tensor4s<T>&);
	Tensor4s<T>& operator=(const SymmetricTensor4s<T>&);
	Tensor4s<T>& operator-=(const Tensor4s<T>&);
	Tensor4s<T>& operator+=(const Tensor4s<T>&);
	Tensor4s<T>& operator*=(T v_);
	// Tensor4s<T>& operator*=(T v_)
	//{
	//	std::for_each(Data, Data + 81, fMilt<T, T>(v_));
	//	return *this;
	//}
	Tensor4s<T>& operator/=(T);
	Tensor4s<T>& Negate()
	{
		std::transform(Data, Data + 81, Data, std::negate<T>());
		return *this;
	}
	Tensor4s<T>& Assign0(small_t, small_t, small_t, small_t);
	Tensor4<T>& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	Tensor4s<T>& Assign1(small_t, T = 1.);	 // three isotropic tensors
	Tensor2s<T>& Contraction(small_t, small_t, Tensor2s<T>&) const;
	Tensor3s<T>& ScalarProductWithBasisVector(small_t, Tensor3s<T>&) const;
	Tensor3s<T>& ScalarProductWithBasisVector_left(small_t, Tensor3s<T>&) const;
	Tensor3s<T>& DotProduct(const Tensor1s<T>&, Tensor3s<T>&) const;
	Tensor4s<T>& CrossMultiply(const Tensor1s<T>&);
	Tensor4s<T>& CrossMultiply_left(const Tensor1s<T>&);
	Tensor4s<T>& CrossProduct(const Tensor1s<T>&, Tensor4s<T>&) const;
	Tensor2s<T>& Dot2Product(const Tensor2s<T>&, Tensor2s<T>&) const;
	Tensor4s<T>& DotMultiply(const Tensor2s<T>&, small_t, small_t);
	Tensor4s<T>& DotMultiply(const Tensor2s<T>&);
	Tensor4s<T>& DotMultiply_left(const Tensor2s<T>&);
	Tensor4s<T>& DotProduct(const Tensor2s<T>&, Tensor4s<T>&) const;
	Tensor1s<T>& Dot3Product(const Tensor3s<T>&, Tensor1s<T>&) const;
	Tensor3s<T>& Dot2Product(const Tensor3s<T>&, Tensor3s<T>&) const;
	T Dot4Product(const Tensor4s<T>&) const;
	Tensor2s<T>& Dot3Product(const Tensor4s<T>&, Tensor2s<T>&) const;
	Tensor4s<T>& Dot2Multiply(const Tensor4s<T>&);
	Tensor4s<T>& Dot2Multiply_left(const Tensor4s<T>&);
	Tensor4s<T>& Dot2Product(const Tensor4s<T>&, Tensor4s<T>&) const;
	Tensor4s<T>& TransformAll(const Tensor2s<T>&);
};
// std::ostream& operator<< (std::ostream&, Tensor4&);
//------------------------------------------------------------------------------
template <typename T>
class SymmetricTensor4s : public Tensor4<T>
{
	typedef typename Tensor4<T>::tIndex tIndex;

	friend class Tensor4s<T>;
	friend class Tensor2s<T>;

private:
	T Data[21];	 // 11..16,22..26,33..36,44..46,55,56,66

	void __Assign(T v_)
	{
		std::fill(Data, Data + 21, v_);
	}
	void __Assign(const T* v_)
	{
		std::copy(v_, v_ + 21, Data);
	}
	/*virtual*/ T operator[](tIndex) const;

public:
	T operator()(small_t, small_t) const;
	T& operator()(small_t, small_t);
	T operator()(small_t, small_t, small_t, small_t) const;
	SymmetricTensor4s& Assign0()
	{
		__Assign(0.);
		return *this;
	}
	SymmetricTensor4s& TransformAll(const Tensor2s<T>&);
};
//------------------------------------------------------------------------------
//==========================INLINES:============================================
//------------------------------------------------------------------------------
//==========================Tensor1s<T>============================================
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>::Tensor1s(TENSOR_KIND, small_t number_, T v_)  // Basis vector multiplied by v_
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= 3,
		"invalid number of coordinate in Tensor1s<T>::Tensor1s<T>(bool,small_t,T)");
#endif
	Data[--number_] = v_;
	Data[number_ == 0 ? 1 : 0] = Data[number_ == 2 ? 1 : 2] = 0.;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>::Tensor1s(const T* const p_)
{
#ifdef STRONGCHECK
	Assert(p_ != NULL, "non-allocated data in Tensor1s<T>::Tensor1s<T>(const T* const)");
#endif
	__Assign(p_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>::Tensor1s(const Tensor1a<T>& o_)
{
	operator=(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::operator=(const Tensor1a<T>& o_)
{
	if (o_.is0())
		Assign0();
	else
		__Assign(*o_.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor1s<T>::operator()(small_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "invalid index in Tensor1s<T>::()");
#endif
	return Data[--i_];
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor1s<T>::operator()(small_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "index is beyond the range in Tensor1s<T>::()const");
#endif
	return Data[--i_];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::operator/=(T v_)
{
	// std::cout << "\n\tRun Tensor1s<T>::operator/=\n";

#ifdef STRONGCHECK
	Assert(NE0(v_), "zero divisor in Tensor1s<T>::operator/=(T)");
#endif
	return operator*=(1. / v_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossMultiplyByBasisVector(small_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= 3, "invalid index in Tensor1s<T>::CrossMultiplyByBasisVector");
#endif
	const small_t j = (--i_) == 0 ? 1 : (i_ == 1 ? 2 : 0), k = 3 - i_ - j;
	Data[i_] = Data[j];
	Data[j] = Data[k];
	Data[k] = -Data[i_];
	Data[i_] = 0.;
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossMultiply(const Tensor1s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "vector multiplication by itself in Tensor1s<T>::CrossMultiply");
#endif
	// if (&o_ == this) {Assign0(); return *this;}
	T data_0 = Data[0], data_1 = Data[1];
	Data[0] = data_1 * o_[2] - Data[2] * o_[1];
	Data[1] = Data[2] * o_[0] - data_0 * o_[2];
	Data[2] = data_0 * o_[1] - data_1 * o_[0];
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossMultiply_left(const Tensor1s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "vector multiplication by itself in Tensor1s<T>::CrossMultiply_left");
#endif
	// if (&o_ == this) {Assign0(); return *this;}
	T data_0 = Data[0], data_1 = Data[1];
	Data[0] = o_[1] * Data[2] - o_[2] * data_1;
	Data[1] = o_[2] * data_0 - o_[0] * Data[2];
	Data[2] = o_[0] * data_1 - o_[1] * data_0;
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossMultiplyByBasisVector_left(small_t i_)
{
	CrossMultiplyByBasisVector(i_);
	Data[0] = -Data[0];
	Data[1] = -Data[1];
	Data[2] = -Data[2];
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossProductWithBasisVector(small_t i_, Tensor1s<T>& result_) const
{
	return (result_ = *this).CrossMultiplyByBasisVector(i_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossProductWithBasisVector_left(
	small_t i_, Tensor1s<T>& result_) const
{
	return (result_ = *this).CrossMultiplyByBasisVector_left(i_);
}
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor1s<T>::DotProduct(const Tensor1s<T>& o_) const
{
	return Data[0] * o_.Data[0] + Data[1] * o_.Data[1] + Data[2] * o_.Data[2];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::DotProduct(const Tensor2s<T>& o_, Tensor1s<T>& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor1s<T>::CrossProduct(const Tensor1s<T>& o_, Tensor1s<T>& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor1s<T>::CrossProduct(const Tensor2s<T>& o_, Tensor2s<T>& result_) const
{  //!!!
	return (result_ = o_).CrossMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor1s<T>::CrossProduct(const Tensor3s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = o_).CrossMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor1s<T>::CrossProduct(const Tensor4s<T>& o_, Tensor4s<T>& result_) const
{
	return (result_ = o_).CrossMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor1s<T>::DirectProduct(const Tensor1s<T>& o_, Tensor2a<T>& result_) const
{
	DirectProduct(o_, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
//==========================Tensor1a<T>=============================================
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>::Tensor1a(RetTensor1<T>& o_) : pData(o_.pData)
{
	o_.pData = NULL;
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor1a<T>::operator()(small_t i_) const
{
#ifdef STRONGCHECK
	Assert(!is0(), "NULL data in Tensor1a<T>::operator()const");
#endif
	return pData->operator()(i_);
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor1a<T>::operator()(small_t i_)
{
#ifdef STRONGCHECK
	Assert(!is0(), "NULL data in Tensor1a<T>::operator()");
#endif
	return pData->operator()(i_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::operator=(RetTensor1<T>& o_)
{
	std::swap(pData, o_.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::Normalize()
{
#ifdef STRONGCHECK
	Assert(!is0(), "attempt to normalize zero vector in Tensor1a<T>::Normalize()");
#endif
	pData->Normalize();
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::operator/=(T v_)
{
#ifdef STRONGCHECK
	Assert(NE0(v_), "zero divisor in Tensor1a<T>::operator/=(T)");
#endif
	return is0() ? *this : (pData->operator/=(v_), *this);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::operator-() const
{
	return is0() ? RetTensor1<T>()
				 : RetTensor1<T>(-pData->Data[0], -pData->Data[1], -pData->Data[2]);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::operator-(const Tensor1a<T>& o_) const
{
	return Tensor1a<T>(*this) -= o_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::operator+(const Tensor1a<T>& o_) const
{
	return Tensor1a<T>(*this) += o_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::operator*(T v_) const
{
	return Tensor1a<T>(*this) *= v_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::operator/(T v_) const
{
	return Tensor1a<T>(*this) /= v_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::CrossProductWithBasisVector(small_t i_) const
{
	Tensor1a<T> result;
	return CrossProductWithBasisVector(i_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::CrossProductWithBasisVector_left(small_t i_) const
{
	Tensor1a<T> result;
	return CrossProductWithBasisVector_left(i_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::CrossProduct(const Tensor1a<T>& o_) const
{
	Tensor1a<T> result;
	return CrossProduct(o_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor1<T> Tensor1a<T>::DotProduct(const Tensor2a<T>& o_) const
{
	Tensor1a<T> result;
	return DotProduct(o_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor1a<T>::DirectProduct(const Tensor1a<T>& o_, Tensor2a<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DirectProduct(*o_.pData, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor1a<T>::DirectProduct(const Tensor1a<T>& o_) const
{
	Tensor2a<T> result;
	return DirectProduct(o_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::DotMultiply(const Tensor2a<T>& o_)
{
	if (__AdjustMult(o_))
		pData->DotMultiply(*o_.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::DotMultiply_left(const Tensor2a<T>& o_)
{
	if (__AdjustMult(o_))
		pData->DotMultiply_left(*o_.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor1a<T>::DotProduct(const Tensor2a<T>& o_, Tensor1a<T>& result_) const
{
	return o_.is0() ? result_.Assign0() : (result_ = *this).DotMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor1a<T>::CrossProduct(const Tensor2a<T>& o_, Tensor2a<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->CrossProduct(*o_.pData, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor1a<T>::DirectProduct(const Tensor2a<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DirectProduct(*o_.pData, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor1a<T>::DotProduct(const Tensor3s<T>& o_, Tensor2a<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DotProduct(o_, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor1a<T>::CrossProduct(const Tensor3s<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->CrossProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor1a<T>::DirectProduct(const Tensor3s<T>& o_, Tensor4s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DirectProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor1a<T>::DotProduct(const Tensor4s<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DotProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor1a<T>::CrossProduct(const Tensor4s<T>& o_, Tensor4s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->CrossProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor1a<T>::CrossProduct(const Tensor2a<T>& o_) const
{
	Tensor2a<T> result;
	return CrossProduct(o_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor1a<T>::DotProduct(const Tensor3s<T>& o_) const
{
	Tensor2a<T> result;
	return DotProduct(o_, result);
}
//------------------------------------------------------------------------------
//==========================Tensor2s<T>=======================================
//------------------------------------------------------------------------------
template <typename T>
inline small_t Tensor2s<T>::tIndex::operator[](small_t no_) const
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 2, "invalid Index2");
#endif
	return (no_ == 1) ? first : second;
}
//------------------------------------------------------------------------------
template <typename T>
inline small_t& Tensor2s<T>::tIndex::operator[](small_t no_)
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 2, "invalid Index2 const");
#endif
	return (no_ == 1) ? first : second;
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor2s<T>::operator[](Tensor2s<T>::tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[1] < 3 && i_[2] >= 0 && i_[2] < 3,
		"invalid index in Tensor2s<T>::operator[]const");
#endif
	return Rows[i_[1]][i_[2]];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor2s<T>::operator[](Tensor2s<T>::tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[1] < 3 && i_[2] >= 0 && i_[2] < 3,
		"invalid index in Tensor2s<T>::operator[]");
#endif
	return Rows[i_[1]][i_[2]];
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor2s<T>::operator()(Tensor2s<T>::tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] > 0 && i_[1] <= 3 && i_[2] > 0 && i_[2] <= 3,
		"invalid index in Tensor2s<T>::operator()const");
#endif
	return Rows[--i_[1]][--i_[2]];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor2s<T>::operator()(tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] > 0 && i_[1] <= 3 && i_[2] > 0 && i_[2] <= 3,
		"invalid index in Tensor2s<T>::operator()");
#endif
	return Rows[--i_[1]][--i_[2]];
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor2s<T>::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3, "invalid index in T Tensor2s<T>::operator() const");
#endif
	return Rows[--i_][--j_];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor2s<T>::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3, "invalid index in T Tensor2s<T>::operator()");
#endif
	return Rows[--i_][--j_];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::operator=(const Tensor2s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor2s<T>::operator= (const Tensor2s<T>&)");
#endif
	// if (&o_!=this)
	__Assign(o_.Data);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::operator/=(T v_)
{
#ifdef STRONGCHECK
	Assert(NE0(v_), "zero divisor in Tensor2s<T>::operator/=(T)");
#endif
	return operator*=(1. / v_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::AddTransposed(const Tensor2s<T>& o_)
{
	if (&o_ == this)
	{
		Data[0] += Data[0];
		Data[1] += Data[3];
		Data[2] += Data[6];
		Data[3] = Data[1];
		Data[4] += Data[4];
		Data[5] += Data[7];
		Data[6] = Data[2];
		Data[7] = Data[5];
		Data[8] += Data[8];
	}
	else
	{
		Data[0] += o_.Data[0];
		Data[1] += o_.Data[3];
		Data[2] += o_.Data[6];
		Data[3] += o_.Data[1];
		Data[4] += o_.Data[4];
		Data[5] += o_.Data[7];
		Data[6] += o_.Data[2];
		Data[7] += o_.Data[5];
		Data[8] += o_.Data[8];
	}
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline SymmetricTensor2s<T>& Tensor2s<T>::Symmetrization_twofold(
	SymmetricTensor2s<T>& result_) const
{
	// result_.Data[6] ~ 11 22 33 12 23 31
	result_.Data[0] = Data[0] + Data[0];
	result_.Data[3] = Data[1] + Data[3];
	result_.Data[5] = Data[2] + Data[6];
	result_.Data[1] = Data[4] + Data[4];
	result_.Data[4] = Data[5] + Data[7];
	result_.Data[2] = Data[8] + Data[8];
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::Transpose()
{
	std::iter_swap(Data + 1, Data + 3);
	std::iter_swap(Data + 2, Data + 6);
	std::iter_swap(Data + 5, Data + 7);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor2s<T>::Det() const
{
	return Data[0] * Data[4] * Data[8] + Data[1] * Data[5] * Data[6] + Data[2] * Data[3] * Data[7] -
		   Data[2] * Data[4] * Data[6] - Data[0] * Data[5] * Data[7] - Data[1] * Data[3] * Data[8];
}
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor2s<T>::I2() const
{
	return Data[0] * Data[4] + Data[4] * Data[8] + Data[8] * Data[0] - Data[1] * Data[3] -
		   Data[5] * Data[7] - Data[6] * Data[2];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor2s<T>::AttachedVector(Tensor1s<T>& result_) const
{
	return result_.Assign(
		0.5 * (Data[7] - Data[5]), 0.5 * (Data[2] - Data[6]), 0.5 * (Data[3] - Data[1]));
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T> Tensor2s<T>::AttachedVector() const
{
	return Tensor1s<T>(
		0.5 * (Data[7] - Data[5]), 0.5 * (Data[2] - Data[6]), 0.5 * (Data[3] - Data[1]));
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor2s<T>::ScalarProductWithBasisVector(
	small_t i_, Tensor1s<T>& result_) const
{
	--i_;
	return result_.Assign(Rows[0][i_], Rows[1][i_], Rows[2][i_]);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor2s<T>::ScalarProductWithBasisVector_left(
	small_t i_, Tensor1s<T>& result_) const
{
	--i_;
	return result_.Assign(Rows[i_][0], Rows[i_][1], Rows[i_][2]);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T> Tensor2s<T>::ScalarProductWithBasisVector(small_t i_) const
{
	--i_;
	return Tensor1s<T>(Rows[0][i_], Rows[1][i_], Rows[2][i_]);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T> Tensor2s<T>::ScalarProductWithBasisVector_left(small_t i_) const
{
	--i_;
	return Tensor1s<T>(Rows[i_][0], Rows[i_][1], Rows[i_][2]);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::AssignDev()
{
	const T p = I1() / 3.l;
	Data[0] -= p;
	Data[4] -= p;
	Data[8] -= p;
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1s<T>& Tensor2s<T>::DotProduct(const Tensor1s<T>& o_, Tensor1s<T>& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::CrossProduct(const Tensor1s<T>& o_, Tensor2s<T>& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::DotProduct(const Tensor2s<T>& o_, Tensor2s<T>& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor2s<T>::DotProduct(const Tensor3s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor2s<T>::Dot2Product(const Tensor4s<T>& o_, Tensor2s<T>& result_) const
{
	return (result_ = *this).Dot2Multiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor2s<T>::DotProduct(const Tensor4s<T>& o_, Tensor4s<T>& result_) const
{
	return (result_ = o_).DotMultiply_left(*this);
}
//------------------------------------------------------------------------------
//==========================Tensor2a<T>=============================================
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>::Tensor2a(RetTensor2<T>& o_) : pData(o_.pData)
{
	o_.pData = NULL;
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor2a<T>::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(!is0(), "NULL data in Tensor2a<T>::operator()const");
#endif
	return pData->operator()(i_, j_);
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor2a<T>::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(!is0(), "NULL data in Tensor2a<T>::operator()");
#endif
	return pData->operator()(i_, j_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::operator=(RetTensor2<T>& o_)
{
	std::swap(pData, o_.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::operator+=(Tensor2CaptureDataAdapter<T> o_)
{
	if (is0())
	{
		pData = o_.Ref.pData;
		o_.Ref.pData = NULL;
	}
	else if (!o_.Ref.is0())
		pData->operator+=(*o_.Ref.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::AddTransposed(Tensor2CaptureDataAdapter<T> o_)
{
	if (!o_.Ref.is0())
		if (is0())
		{
			pData = o_.Ref.pData;
			o_.Ref.pData = NULL;
			pData->Transpose();
		}
		else
			pData->AddTransposed(*o_.Ref.pData);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::operator-() const
{
	return Tensor2a<T>(*this).Negate();
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::operator-(const Tensor2a<T>& o_) const
{
	return Tensor2a<T>(*this) -= o_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::operator+(const Tensor2a<T>& o_) const
{
	return Tensor2a<T>(*this) += o_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::operator*(T v_) const
{
	return Tensor2a<T>(*this) *= v_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::operator/(T v_) const
{
	return Tensor2a<T>(*this) /= v_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::Invert()
{
#ifdef STRONGCHECK
	Assert(!is0(), "NULL data in Tensor2a<T>::Invert()");
#endif
	pData->Invert();
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::DotProduct(const Tensor2a<T>& o_) const
{
	return Tensor2a<T>(*this).DotMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::operator/=(T v_)
{
#ifdef STRONGCHECK
	Assert(NE0(v_), "zero divisor in Tensor2a<T>::operator/=(T)");
#endif
	return is0() ? *this : (pData->operator/=(v_), *this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor2a<T>::DirectProduct(const Tensor1a<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DirectProduct(*o_.pData, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor2a<T>::CrossProduct(const Tensor2a<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->CrossProduct(*o_.pData, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor2a<T>::DirectProduct(const Tensor2a<T>& o_, Tensor4s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DirectProduct(*o_.pData, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor1a<T>& Tensor2a<T>::Dot2Product(const Tensor3s<T>& o_, Tensor1a<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->Dot2Product(o_, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor2a<T>::DotProduct(const Tensor3s<T>& o_, Tensor3s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->DotProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor2a<T>::CrossProduct(const Tensor3s<T>& o_, Tensor4s<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->CrossProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
// inline Tensor2a<T>& Tensor2a<T>::Dot2Multiply (const Tensor4& o_, small_t i_, small_t j_)
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::Dot2Multiply(
	const SymmetricTensor4s<T>& o_, small_t i_, small_t j_)
{
	if (__AdjustMult(o_))
		pData->Dot2Multiply(o_, i_, j_);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::Dot2Multiply(const Tensor4s<T>& o_)
{
	if (__AdjustMult(o_))
		pData->Dot2Multiply(o_);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::Dot2Multiply_left(const Tensor4s<T>& o_)
{
	if (__AdjustMult(o_))
		pData->Dot2Multiply_left(o_);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2a<T>& Tensor2a<T>::Dot2Product(const Tensor4s<T>& o_, Tensor2a<T>& result_) const
{
	if (__AdjustMult(o_, result_))
		pData->Dot2Product(o_, *result_.__pData());
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor2a<T>::DotProduct(const Tensor4s<T>& o_, Tensor4s<T>& result_) const
{  //!!!!!!!!!!! может, лучше использовать DotMultiply?
	if (__AdjustMult(o_, result_))
		pData->DotProduct(o_, result_);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::CrossProduct(const Tensor1a<T>& o_) const
{
	Tensor2a<T> result;
	return CrossProduct(o_, result);
}
//------------------------------------------------------------------------------
template <typename T>
inline RetTensor2<T> Tensor2a<T>::Dot2Product(const Tensor4s<T>& o_) const
{
	Tensor2a<T> result;
	return Dot2Product(o_, result);
}
//------------------------------------------------------------------------------
//==========================Tensor3s<T>==========================================
//------------------------------------------------------------------------------
template <typename T>
inline small_t Tensor3s<T>::tIndex::operator[](small_t indexNo_) const
{
#ifdef STRONGCHECK
	Assert(indexNo_ > 0 && indexNo_ <= 3, "invalid Index3 const");
#endif
	return (indexNo_ == 1) ? first : (indexNo_ == 2 ? second : third);
}
//------------------------------------------------------------------------------
template <typename T>
inline small_t& Tensor3s<T>::tIndex::operator[](small_t indexNo_)
{
#ifdef STRONGCHECK
	Assert(indexNo_ > 0 && indexNo_ <= 3, "invalid Index3");
#endif
	return (indexNo_ == 1) ? first : (indexNo_ == 2 ? second : third);
}
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor3s<T>::operator[](const Tensor3s<T>::tIndex& i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] < 3 && i_[2] < 3 && i_[3] < 3,
		"invalid index in Tensor3s<T>::operator[](tIndex3)const");
#endif
	return (*Rows[i_[1]])[i_[2]][i_[3]];
}
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor3s<T>::operator()(small_t i_, small_t j_, small_t k_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3,
		"invalid index in Tensor3s<T>::operator()const");
#endif
	return (*Rows[--i_])[--j_][--k_];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor3s<T>::operator()(small_t i_, small_t j_, small_t k_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3,
		"invalid index in Tensor3s<T>::operator()");
#endif
	return (*Rows[--i_])[--j_][--k_];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor3s<T>::operator=(const Tensor3s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor2s<T>::operator= (const Tensor2s<T>&)");
#endif
	// if (&o_!=this)
	__Assign(o_.Data);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor3s<T>::operator/=(T v_)
{
#ifdef STRONGCHECK
	Assert(NE0(v_), "suspicious absence of exception in Tensor3s<T>::operator/=(T)");
#endif
	return operator*=(1. / v_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor3s<T>::CrossProduct(const Tensor1s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor3s<T>::DotProduct(const Tensor2s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor3s<T>::Dot2Product(const Tensor4s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = *this).Dot2Multiply(o_);
}
//------------------------------------------------------------------------------
//==========================Tensor4s<T>==========================================
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>::Tensor4s(const SymmetricTensor4s<T>& o_)
{
	operator=(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline small_t Tensor4<T>::tIndex::operator[](small_t no_) const
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 4, "invalid Tensor4::tIndex[]const");
#endif
	return indices[--no_];
}
//------------------------------------------------------------------------------
template <typename T>
inline small_t& Tensor4<T>::tIndex::operator[](small_t no_)
{
#ifdef STRONGCHECK
	Assert(no_ > 0 && no_ <= 4, "invalid Tensor4::tIndex[]");
#endif
	return indices[--no_];
}
//------------------------------------------------------------------------------
template <typename T>
inline const T& Tensor4s<T>::operator()(small_t i_, small_t j_, small_t k_, small_t l_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s<T>::operator[](small_t,small_t,small_t,small_t)const");
#endif
	return (*Rows[--i_])[--j_][--k_][--l_];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor4s<T>::operator()(small_t i_, small_t j_, small_t k_, small_t l_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s<T>::operator[](small_t,small_t,small_t,small_t)");
#endif
	return (*Rows[--i_])[--j_][--k_][--l_];
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::operator=(const Tensor4s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "assignment to itself in Tensor4s<T>::operator= (const Tensor4s<T>&)");
#endif
	// if (&o_ != this)
	__Assign(o_.Data);
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::operator/=(T v_)
{
#ifdef STRONGCHECK
	Assert(NE0(v_), "suspicious absence of exception in Tensor4s<T>::operator/=(T)");
#endif
	return operator*=(1. / v_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::Assign0(small_t i_, small_t j_, small_t k_, small_t l_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && k_ > 0 && l_ > 0 && i_ <= 3 && j_ <= 3 && k_ <= 3 && l_ <= 3,
		"invalid index in Tensor4s<T>::Assign0(small_t,small_t,small_t,small_t)");
#endif
	(*Rows[--i_])[--j_][--k_][--l_] = 0.;
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::CrossProduct(const Tensor1s<T>& o_, Tensor4s<T>& result_) const
{
	return (result_ = *this).CrossMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor2s<T>& Tensor4s<T>::Dot2Product(const Tensor2s<T>& o_, Tensor2s<T>& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::DotProduct(const Tensor2s<T>& o_, Tensor4s<T>& result_) const
{
	return (result_ = *this).DotMultiply(o_);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor3s<T>& Tensor4s<T>::Dot2Product(const Tensor3s<T>& o_, Tensor3s<T>& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::Dot2Product(const Tensor4s<T>& o_, Tensor4s<T>& result_) const
{
	return (result_ = o_).Dot2Multiply_left(*this);
}
//------------------------------------------------------------------------------
template <typename T>
inline Tensor4s<T>& Tensor4s<T>::TransformAll(const Tensor2s<T>& tr_)
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
//------------------------------------------------------------------------------
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
	}
	return 1;
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
	}
	return 3;
}
//------------------------------------------------------------------------------
template <typename T>
inline SymmetricTensor4s<T>& SymmetricTensor4s<T>::TransformAll(const Tensor2s<T>& /*ll*/)
{
	using std::sqrt;
	Tensor2s<T> tensor2s;
	tensor2s(1, 1) = 1. / sqrt(2.);
	tensor2s(1, 2) = -1. / sqrt(2.);
	tensor2s(1, 3) = 0.;
	tensor2s(2, 1) = 1. / sqrt(2.);
	tensor2s(2, 2) = 1. / sqrt(2.);
	tensor2s(2, 3) = 0.;
	tensor2s(3, 1) = 0.;
	tensor2s(3, 2) = 0.;
	tensor2s(3, 3) = 1.;
	const Tensor2s<T>& q(tensor2s);
	operator()(1, 1) = 165.64l;
	operator()(1, 2) = 63.94l;
	operator()(1, 3) = 63.94l;
	operator()(1, 4) = 0.;
	operator()(1, 5) = 0.;
	operator()(1, 6) = 0.;

	operator()(2, 2) = 165.64l;
	operator()(2, 3) = 63.94l;
	operator()(2, 4) = 0.;
	operator()(2, 5) = 0.;
	operator()(2, 6) = 0.;

	operator()(3, 3) = 165.64l;
	operator()(3, 4) = 0.;
	operator()(3, 5) = 0.;
	operator()(3, 6) = 0.;

	operator()(4, 4) = 79.51l;
	operator()(4, 5) = 0.;
	operator()(4, 6) = 0.;

	operator()(5, 5) = 79.51l;
	operator()(5, 6) = 0.;

	operator()(6, 6) = 79.51l;

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
	/* T
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
	static T g[6][6];
	g[0][0] = tensor2s(1, 1) * tensor2s(1, 1);
	g[0][1] = tensor2s(1, 2) * tensor2s(1, 2);
	g[0][2] = tensor2s(1, 3) * tensor2s(1, 3);
	g[0][3] = 2. * tensor2s(1, 2) * tensor2s(1, 3);
	g[0][4] = 2. * tensor2s(1, 3) * tensor2s(1, 1);
	g[0][5] = 2. * tensor2s(1, 2) * tensor2s(1, 1);
	g[1][0] = tensor2s(2, 1) * tensor2s(2, 1);
	g[1][1] = tensor2s(2, 2) * tensor2s(2, 2);
	g[1][2] = tensor2s(2, 3) * tensor2s(2, 3);
	g[1][3] = 2. * tensor2s(2, 3) * tensor2s(2, 2);
	g[1][4] = 2. * tensor2s(2, 3) * tensor2s(2, 1);
	g[1][5] = 2. * tensor2s(2, 2) * tensor2s(2, 1);
	g[2][0] = tensor2s(3, 1) * tensor2s(3, 1);
	g[2][1] = tensor2s(3, 2) * tensor2s(3, 2);
	g[2][2] = tensor2s(3, 3) * tensor2s(3, 3);
	g[2][3] = 2. * tensor2s(3, 3) * tensor2s(3, 2);
	g[2][4] = 2. * tensor2s(3, 3) * tensor2s(3, 1);
	g[2][5] = 2. * tensor2s(3, 2) * tensor2s(3, 1);
	g[3][0] = tensor2s(3, 1) * tensor2s(2, 1);
	g[3][1] = tensor2s(3, 2) * tensor2s(2, 2);
	g[3][2] = tensor2s(3, 3) * tensor2s(2, 3);
	g[3][3] = tensor2s(3, 3) * tensor2s(2, 2) + tensor2s(3, 2) * tensor2s(2, 3);
	g[3][4] = tensor2s(3, 3) * tensor2s(2, 1) + tensor2s(3, 1) * tensor2s(2, 3);
	g[3][5] = tensor2s(3, 1) * tensor2s(2, 2) + tensor2s(3, 2) * tensor2s(2, 1);
	g[4][0] = tensor2s(3, 1) * tensor2s(1, 1);
	g[4][1] = tensor2s(3, 2) * tensor2s(1, 2);
	g[4][2] = tensor2s(3, 3) * tensor2s(1, 3);
	g[4][3] = tensor2s(3, 3) * tensor2s(1, 2) + tensor2s(3, 2) * tensor2s(1, 3);
	g[4][4] = tensor2s(3, 3) * tensor2s(1, 1) + tensor2s(3, 1) * tensor2s(1, 3);
	g[4][5] = tensor2s(3, 1) * tensor2s(1, 2) + tensor2s(3, 2) * tensor2s(1, 1);
	g[5][0] = tensor2s(2, 1) * tensor2s(1, 1);
	g[5][1] = tensor2s(1, 2) * tensor2s(2, 2);
	g[5][2] = tensor2s(1, 3) * tensor2s(2, 3);
	g[5][3] = tensor2s(1, 3) * tensor2s(2, 2) + tensor2s(1, 2) * tensor2s(2, 3);
	g[5][4] = tensor2s(1, 1) * tensor2s(2, 3) + tensor2s(1, 3) * tensor2s(2, 1);
	g[5][5] = tensor2s(1, 1) * tensor2s(2, 2) + tensor2s(1, 2) * tensor2s(2, 1);

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
	for (small_t i = 0, j, k, l, i1 = 1, j1, k1, l1; i < 6; ++i, ++i1)
		for (l = i, l1 = i1; l < 6; ++l, ++l1)
			for (j = 0, j1 = 1; j < 6; ++j, ++j1)
				for (k = 0, k1 = 1; k < 6; ++k, ++k1)
					operator()(i1, l1) = A(j1, k1) * g[i][j] * g[l][k];

	/*Tensor4s<T> T4;//(*this);
  for (small_t i=1,j,k,l; i<=3; ++i)
  for (j=1; j<=3; ++j)
  for
	 * (k=1; k<=3; ++k)
  for (l=1; l<=3; ++l)
    T4(i,j,k,l) = operator()(i,j,k,l);

	 * T4.TransformAll(q);
  for (small_t i=1,j; i<=6; ++i)
  for (j=i; j<=6; ++j)

	 * operator()(i,j) = T4(map1(i),map2(i),map1(j),map2(j));
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
//------------------------------------------------------------------------------
template <typename T>
inline T SymmetricTensor4s<T>::operator()(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= 6 && j_ <= 6,
		"invalid index in SymmetricTensor4s::operator(small_t,small_t)const");
#endif
	if (--i_ > --j_)
		std::swap(i_, j_);
	return Data[(13 - i_) * i_ / 2 + j_ - i_];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& SymmetricTensor4s<T>::operator()(small_t i_, small_t j_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ <= 6 && i_ <= j_,
		tMessage("invalid index in SymmetricTensor4s::operator(small_t,small_t), must be 0<")
			<< i_ << "<=" << j_ << "<=6");
#endif
	return Data[(13 - (i_ - 1)) * (i_ - 1) / 2 + j_ - i_];
}
//------------------------------------------------------------------------------
template <typename T>
inline T SymmetricTensor4s<T>::operator()(small_t i_, small_t j_, small_t k_, small_t l_) const
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
//------------------------------------------------------------------------------
template <typename T>
inline T SymmetricTensor4s<T>::operator[](tIndex i_) const
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
//------------------------------------------------------------------------------
template <typename T>
inline T Tensor4s<T>::operator[](tIndex i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[2] >= 0 && i_[3] >= 0 && i_[4] >= 0 && i_[1] < 3 && i_[2] < 3 &&
			i_[3] < 3 && i_[4] < 3,
		"invalid index in Tensor4s<T>::operator[](Index4)const");
#endif
	return (*Rows[i_[1]])[i_[2]][i_[3]][i_[4]];
}
//------------------------------------------------------------------------------
template <typename T>
inline T& Tensor4s<T>::operator[](tIndex i_)
{
#ifdef STRONGCHECK
	Assert(
		i_[1] >= 0 && i_[2] >= 0 && i_[3] >= 0 && i_[4] >= 0 && i_[1] < 3 && i_[2] < 3 &&
			i_[3] < 3 && i_[4] < 3,
		"invalid index in Tensor4s<T>::operator[](Index4)");
#endif
	return (*Rows[i_[1]])[i_[2]][i_[3]][i_[4]];
}
//------------------------------------------------------------------------------
//#undef __SET_ROWS
#undef __SET_ROWS3
#undef __SET_ROWS4
//------------------------------------------------------------------------------

#endif	// ndef TensorsH
