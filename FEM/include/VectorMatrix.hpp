#pragma once

#include "HybridTensors.hpp"
#include "NumTypes.hpp"

#include <vector>
#include <string>
#include <functional>
#include <algorithm>

class Tensor2a;
// class Tensor2s;
class tNode;
class tFinitElement;
class tFE_model;
class tNodalTensor1;
class tNodalTensor2;
class tFEsTensor2;
// class tFEsSetOfTensor2;

// Sparse Vectors/Matrices of nodal tensors

class tNodalTensor1Column
{
public:
	tNodalTensor1Column(const tFE_model&, const std::string& name_ = std::string());
	tNodalTensor1Column()
	{
	}
	//   ~tNodalTensor1Column() {}
	//   tNodalTensor1Column(const tNodalTensor1Column& o_): Caption(o_.Caption),
	//   Components(o_.Components) {}
	tNodalTensor1Column& Link(const tFE_model&);
	tNodalTensor1Column& SetName(const char* name_)
	{
		m_caption = name_;
		return *this;
	}
	const std::string& Name() const
	{
		return m_caption;
	}
	size_t Dimension() const
	{
		return m_components.size();
	}
	size_t NumberOfDOFs() const;
	tNodalTensor1Column& Assign0();
	const tNode& Node(size_t) const;
	tNodalTensor1Column& Assign0(size_t);
	bool is0(size_t) const;
	const tNodalTensor1& operator[](size_t) const;
	tNodalTensor1& operator[](size_t);
	const tNodalTensor1& operator()(size_t) const;
	tNodalTensor1& operator()(size_t);
	const tNodalTensor1& operator()(const tNode&) const;
	tNodalTensor1& operator()(const tNode&);
	tNodalTensor1Column& AssignComponent(const tNodalTensor1& o_)
	{
		operator()(o_.Node()) = o_;
		return *this;
	}
	tNodalTensor1Column& operator+=(const tNodalTensor1& o_)
	{
		operator()(o_.Node()) += o_;
		return *this;
	}
	tNodalTensor1Column& operator+=(const tNodalTensor1Column&);
	tNodalTensor1Column& operator-=(const tNodalTensor1Column&);
	const tNodalTensor1Column& LinkAsDisplacements() const;

private:
	std::string m_caption;
	std::vector<tNodalTensor1> m_components;
	mutable std::vector<tNodalTensor1>::const_iterator
		m_constNodeCursor;								// only for operator()(const tNode&)const
	std::vector<tNodalTensor1>::iterator m_nodeCursor;	// only for operator()(const tNode&)
};

/*class tFEsTensor2Column
{
 private:
   std::string Caption;
   std::vector<tFEsTensor2> Components;
 public:
   tFEsTensor2Column(const tFE_model&, const std::string& name_=std::string());
   tFEsTensor2Column(): Caption(), Components() {}
//   ~tFEsTensor2Column() {}
   tFEsTensor2Column& Link(const tFE_model&);
   tFEsTensor2Column& SetName(const std::string& name_) {Caption = name_; return *this;}
   const std::string& Name() const {return Caption;}
   size_t Dimension() const {return Components.size();}
   tFEsTensor2Column& Assign0();
   struct fTransformer {
						virtual void operator()(tFEsTensor2&) const =0;
					   };
   tFEsTensor2Column& Transform(const fTransformer&);
//   tFEsTensor2Column& RotateToLocalAxes();
   const tFinitElement& FE(size_t i_) const;
   tFEsTensor2Column& Assign0(size_t i_);
   bool is0 (size_t i_) const;
   const tFEsTensor2& operator[](size_t i_) const;
		 tFEsTensor2& operator[](size_t i_);
   const tFEsTensor2& operator()(size_t i_) const;
		 tFEsTensor2& operator()(size_t i_);
   const tFEsTensor2& operator()(const tFinitElement&) const;
		 tFEsTensor2& operator()(const tFinitElement&);
   tFEsTensor2Column& operator+=(const tFEsTensor2&);
};*/

// class tFEsSetOfTensor2Column
//{
// public:
//   enum Tensor2Kind_t {linear_strain, Cauchy_stress};
///*   class tTensor2Kind {private:
//                         Tensor2Kind_t Kind;
//                         bool GlobalDir;
//                         tTensor2Kind();//: Kind(linear_strain), GlobalDir(true) {}//is not used
//                       public:
//                         tTensor2Kind(Tensor2Kind_t k_,bool g_=true): Kind(k_), GlobalDir(g_) {}
//                         tTensor2Kind& operator=(Tensor2Kind_t k_) {Kind=k_; return *this;}
//                         bool operator==(const tTensor2Kind& c_) const {return Kind==c_.Kind &&
//                         GlobalDir==c_.GlobalDir;} bool operator!=(const tTensor2Kind& c_) const
//                         {return !operator==(c_);} bool operator<(const tTensor2Kind& c_) const
//                         {return (Kind<c_.Kind || (Kind==c_.Kind &&
//                         (GlobalDir&&!c_.BelongsToGlobalCoordSystem())))?true:false;} bool
//                         BelongsToGlobalCoordSystem() const {return GlobalDir;} tTensor2Kind&
//                         SetBelongingToLocal() {GlobalDir=false; return *this;} tTensor2Kind&
//                         SetBelongingToGlobal(){GlobalDir=true; return *this;} tTensor2Kind
//                         ConjByMatLow() const {return
//                         tTensor2Kind(Kind==linear_strain?Cauchy_stress:linear_strain,GlobalDir);}
//                         tTensor2Kind ConjByRotation() const {return
//                         tTensor2Kind(Kind,!GlobalDir);} bool isStrain() const {return
//                         Kind==linear_strain;} bool isStress() const {return Kind==Cauchy_stress;}
//                         tTensor2Kind& Assign(Tensor2Kind_t k_) {Kind=k_;return *this;};
//                      };*/
// private:
//   std::vector<tTensor2Kind> Kinds;
//   std::vector<tFEsSetOfTensor2> Components;
// public:
//   tFEsSetOfTensor2Column(const tFE_model&, size_t=1);
//   tFEsSetOfTensor2Column& Include(tTensor2Kind);
//   size_t ColNo(tTensor2Kind) const;
//   struct fTransformer {
//                        virtual void operator()(tFEsSetOfTensor2&) const =0;
//                       };
//   tFEsSetOfTensor2Column& Transform(const fTransformer&);
//   tFEsSetOfTensor2Column& Swap(tTensor2Kind,tFEsTensor2Column&);
//};

class tNodalTensor2SymMatrix
{
public:
	//   tNodalTensor2SymMatrix(const tFE_model&, const std::string&);
	//   tNodalTensor2SymMatrix(): Caption(), pNodes(), pComponents() {}
	//   ~tNodalTensor2SymMatrix();
	tNodalTensor2SymMatrix& SetName(const char* name_)
	{
		m_caption = name_;
		return *this;
	}
	tNodalTensor2SymMatrix& Link(const tFE_model&);
	size_t Dimension() const
	{
		return m_components.size();
	}
	const std::string& Name() const
	{
		return m_caption;
	}
	size_t NumberOfDOFs() const;
	tNodalTensor2SymMatrix& Assign0();
	size_t BandWidth(size_t) const;
	size_t MaxBandWidth() const;
	const tNode& Node(size_t) const;
	bool is0(size_t, size_t) const;
	const Tensor2a& operator()(size_t, size_t) const;
	tNodalTensor2SymMatrix& operator+=(tNodalTensor2&);

private:
	std::string m_caption;
	std::vector<const tNode*> m_pNodes;
	//   std::vector<Tensor2s**> pComponents;
	std::vector<std::vector<Tensor2a> > m_components;
};

inline const tNode& tNodalTensor2SymMatrix::Node(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_pNodes.size(), "Invalid index in tNodalTensor2SymMatrix::Node");
#endif
	return *m_pNodes[--i_];
}

inline bool tNodalTensor2SymMatrix::is0(size_t i_, size_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= m_components.size() && j_ <= m_components.size(),
		"Invalid index in tNodalTensor2SymMatrix::is0");
	Assert(
		i_ >= j_,
		"Indices in non-existing upper half of symmetric matrix in tNodalTensor2SymMatrix::is0");
#endif
	return m_components[--i_][--j_].is0();
	// return ((--j_)>(--i_)? pComponents[j_][i_] : pComponents[i_][j_]) == nullptr;
}

inline const Tensor2a& tNodalTensor2SymMatrix::operator()(size_t i_, size_t j_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && j_ > 0 && i_ <= m_components.size() && j_ <= m_components.size(),
		"Invalid index in tNodalTensor2SymMatrix::operator()");
	Assert(
		i_ >= j_,
		"Indices in non-existing upper half of symmetric matrix in "
		"tNodalTensor2SymMatrix::operator()");
	Assert(!is0(i_, j_), "Absent (nullptr) component in tNodalTensor2SymMatrix::operator()");
#endif
	return m_components[--i_][--j_];
	// return *((--j_)>(--i_)? pComponents[j_][i_] : pComponents[i_][j_]);
}

// Sparse Vectors/Matrices of reals

class tConstColumnVector
{
public:
	virtual ~tConstColumnVector()
	{
	}
	virtual size_t Dimension() const = 0;
	virtual bool is0(size_t) const = 0;
	virtual const real_t& operator[](size_t) const = 0;
	virtual const real_t& operator()(size_t) const = 0;
	virtual tConstColumnVector& Swap(size_t, size_t) = 0;
};

class tColumnVector : public tConstColumnVector
{
public:
	virtual ~tColumnVector()
	{
	}

	virtual tColumnVector& Assign0() = 0;
	virtual tColumnVector& Assign(size_t, real_t) = 0;
	virtual tColumnVector& Add(size_t, real_t) = 0;
	virtual real_t& operator[](size_t) = 0;
	virtual real_t& operator()(size_t) = 0;
};

class tSparseColumn : public tColumnVector
{
public:
	tSparseColumn(const tConstColumnVector&, bool);
	~tSparseColumn();
	size_t Dimension() const override
	{
		return m_pComponents.size();
	}
	tColumnVector& Assign0() override;
	tColumnVector& Assign(size_t, real_t) override;
	tColumnVector& Add(size_t, real_t) override;
#ifdef STRONGCHECK
	bool is0(size_t i_) const override
	{
		Assert(i_ > 0 && i_ <= m_pComponents.size(), "invalid index in is0");
		return m_pComponents[--i_] == nullptr;
	}

	const real_t& operator[](size_t i_) const override
	{
		Assert(i_ >= 0 && i_ < m_pComponents.size(), "invalid index in operator[] const");
		Assert(
			m_pComponents[i_] != nullptr,
			"result is non-allocated (==nullptr) in operator[] const");
		return *m_pComponents[i_];
	}

	real_t& operator[](size_t i_) override
	{
		Assert(i_ >= 0 && i_ < m_pComponents.size(), "invalid index in operator[]");
		Assert(m_pComponents[i_] != nullptr, "result is non-allocated (==nullptr) in operator[]");
		return *m_pComponents[i_];
	}

	const real_t& operator()(size_t i_) const override
	{
		Assert(i_ > 0 && i_ <= m_pComponents.size(), "invalid index in operator() const");
		Assert(
			m_pComponents[--i_] != nullptr,
			"result is non-allocated (==nullptr) in operator() const");
		return *m_pComponents[i_];
	}

	real_t& operator()(size_t i_) override
	{
		Assert(i_ > 0 && i_ <= m_pComponents.size(), "invalid index in operator()");
		Assert(m_pComponents[--i_] != nullptr, "result is non-allocated (==nullptr) in operator()");
		return *m_pComponents[i_];
	}

	tConstColumnVector& Swap(size_t i_, size_t j_) override
	{
		Assert(
			i_ > 0 && i_ <= m_pComponents.size() && j_ > 0 && j_ <= m_pComponents.size(),
			"invalid index in Swap");
		// real_t* tmp = pComponents[--i_];  pComponents[i_] = pComponents[--j_];  pComponents[j_] =
		// tmp;
		std::swap(m_pComponents[--i_], m_pComponents[--j_]);
		return *this;
	}
#else	// ifndef STRONGCHECK
	bool is0(size_t i_) const override
	{
		return m_pComponents[--i_] == nullptr;
	}
	const real_t& operator[](size_t i_) const override
	{
		return *m_pComponents[i_];
	}
	real_t& operator[](size_t i_) override
	{
		return *m_pComponents[i_];
	}
	const real_t& operator()(size_t i_) const override
	{
		return *m_pComponents[--i_];
	}
	real_t& operator()(size_t i_) override
	{
		return *m_pComponents[--i_];
	}
	tConstColumnVector& Swap(size_t i_, size_t j_) override
	{
		std::swap(m_pComponents[--i_], m_pComponents[--j_]);
		return *this;
	}
#endif	// def STRONGCHECK

private:
	std::vector<real_t*> m_pComponents;
};

class tMarkedGraph;

class tConstRefColumn : public tConstColumnVector
{
public:
	tConstRefColumn(const tNodalTensor1Column&, const tMarkedGraph&);
	virtual ~tConstRefColumn()
	{
		delete[] m_pComponents;
	}
	size_t Dimension() const override
	{
		return m_size;
	}
#ifdef STRONGCHECK
	bool is0(size_t i_) const override
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tConstRefColumn::is0");
		return m_pComponents[--i_] == nullptr;
	}

	const real_t& operator[](size_t i_) const override
	{
		Assert(i_ >= 0 && i_ < m_size, "invalid index in tConstRefColumn::operator[] const");
		Assert(
			m_pComponents[i_] != nullptr,
			"result is non-allocated (==nullptr) in tConstRefColumn::operator[] const");
		return *m_pComponents[i_];
	}

	const real_t& operator()(size_t i_) const override
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tConstRefColumn::operator() const");
		Assert(
			m_pComponents[--i_] != nullptr,
			"result is non-allocated (==nullptr) in tConstRefColumn::operator() const");
		return *m_pComponents[i_];
	}

	tConstColumnVector& Swap(size_t i_, size_t j_) override
	{
		Assert(
			i_ > 0 && i_ <= m_size && j_ > 0 && j_ <= m_size,
			"invalid index in tConstRefColumn::Swap");
		// const real_t* tmp = pComponents[--i_];  pComponents[i_] = pComponents[--j_];
		// pComponents[j_] = tmp;
		std::swap(m_pComponents[--i_], m_pComponents[--j_]);
		return *this;
	}
#else	// ifndef STRONGCHECK
	virtual bool is0(size_t i_) const
	{
		return m_pComponents[--i_] == nullptr;
	}
	virtual const real_t& operator[](size_t i_) const
	{
		return *m_pComponents[i_];
	}
	virtual const real_t& operator()(size_t i_) const
	{
		return *m_pComponents[--i_];
	}
	virtual tConstColumnVector& Swap(size_t i_, size_t j_)
	{
		std::swap(m_pComponents[--i_], m_pComponents[--j_]);
		return *this;
	}
#endif	// def STRONGCHECK

private:
	typedef const real_t* prtToConstReal;

	const size_t m_size;
	prtToConstReal* m_pComponents;
};

class tRefTensor1ComponentsColumn : public tColumnVector
{
public:
	tRefTensor1ComponentsColumn(tNodalTensor1Column&, const tMarkedGraph&);
	~tRefTensor1ComponentsColumn()
	{
		UpdateSource();
		Assign0();	// delete *pValue for all elements
		delete[] m_components;
	};
	size_t Dimension() const override
	{
		return m_size;
	}
	tColumnVector& Assign0() override;
	tColumnVector& Assign(size_t, real_t) override;
	tColumnVector& Add(size_t, real_t) override;
#ifdef STRONGCHECK
	bool is0(size_t i_) const override
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tRefTensor1ComponentsColumn::is0");
		return m_components[--i_].is0();
	}

	const real_t& operator[](size_t i_) const override
	{
		Assert(
			i_ >= 0 && i_ < m_size,
			"invalid index in tRefTensor1ComponentsColumn::operator[] const");
		Assert(
			!m_components[i_].is0(),
			"result is non-allocated (==nullptr) in tRefTensor1ComponentsColumn::operator[] const");
		return m_components[i_]();
	}

	real_t& operator[](size_t i_) override
	{
		Assert(i_ >= 0 && i_ < m_size, "invalid index in tRefTensor1ComponentsColumn::operator[]");
		Assert(
			!m_components[i_].is0(),
			"result is non-allocated (==nullptr) in tRefTensor1ComponentsColumn::operator[]");
		return m_components[i_]();
	}

	const real_t& operator()(size_t i_) const override
	{
		Assert(
			i_ > 0 && i_ <= m_size,
			"invalid index in tRefTensor1ComponentsColumn::operator() const");
		Assert(
			!m_components[--i_].is0(),
			"result is non-allocated (==nullptr) in tRefTensor1ComponentsColumn::operator() const");
		return m_components[i_]();
	}

	real_t& operator()(size_t i_) override
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tRefTensor1ComponentsColumn::operator()");
		Assert(
			!m_components[--i_].is0(),
			"result is non-allocated (==nullptr) in tRefTensor1ComponentsColumn::operator()");
		return m_components[i_]();
	}

	tConstColumnVector& Swap(size_t i_, size_t j_) override
	{
		Assert(
			i_ > 0 && i_ <= m_size && j_ > 0 && j_ <= m_size,
			"invalid index in tRefTensor1ComponentsColumn::Swap");
		// tElement tmp = Components[--i_];  Components[i_] = Components[--j_];  Components[j_] =
		// tmp;
		std::swap(m_components[--i_], m_components[--j_]);
		return *this;
	}
#else	// ifndef STRONGCHECK
	bool is0(size_t i_) const override
	{
		return m_components[--i_].is0();
	}
	const real_t& operator[](size_t i_) const override
	{
		return m_components[i_]();
	}
	real_t& operator[](size_t i_) override
	{
		return m_components[i_]();
	}
	const real_t& operator()(size_t i_) const override
	{
		return m_components[--i_]();
	}
	real_t& operator()(size_t i_) override
	{
		return m_components[--i_]();
	}
	tConstColumnVector& Swap(size_t i_, size_t j_) override
	{
		std::swap(m_components[--i_], m_components[--j_]);
		return *this;
	}
#endif	// def STRONGCHECK
	tRefTensor1ComponentsColumn& UpdateSource();

private:
	struct tElement
	{
		tElement() : m_pValue(nullptr), m_index(), m_pTensor()
		{
		}
		bool is0() const
		{
			return m_pValue == nullptr;
		}
		real_t& operator()()
		{
			return *m_pValue;
		}
		const real_t& operator()() const
		{
			return *m_pValue;
		}

		tNodalTensor1* m_pTensor;
		small_t m_index;
		real_t* m_pValue;
	};

	const size_t m_size;
	tElement* m_components;
};

class tSymConstRefMatrix
{
public:
	tSymConstRefMatrix(const tNodalTensor2SymMatrix&, const tMarkedGraph&);
	~tSymConstRefMatrix()
	{
		delete[] m_pComponents;
	}
	size_t Dimension() const
	{
		return m_size;
	}
	size_t BandWidth(size_t) const;
	const tSymConstRefMatrix& ProfileAndMaxBandWidth(size_t&, size_t&) const;
#ifdef STRONGCHECK
	bool is0(size_t i_, size_t j_) const
	{
		if (i_ < j_)  //{size_t k=i_; i_=j_; j_=k;}
			std::swap(i_, j_);
		Assert(j_ > 0 && i_ <= m_size, "invalid index in tSymConstRefMatrix::is0(size_t,size_t)");
		--i_;
		--j_;
		Assert(
			i_ - j_ < m_pComponents[i_].size(),
			"index beyond band in tSymConstRefMatrix::is0(size_t,size_t)");
		return m_pComponents[i_][i_ - j_] == nullptr;
	}

	const real_t& DiagonalComponent(size_t i_) const
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tSymConstRefMatrix::operator()(size_t)");
		return *(m_pComponents[--i_].front());
	}

	const real_t& operator()(size_t i_, size_t j_) const
	{
		if (i_ < j_)  //{size_t k=i_; i_=j_; j_=k;}
			std::swap(i_, j_);
		Assert(
			j_ > 0 && i_ <= m_size,
			"invalid index in tSymConstRefMatrix::operator()(size_t,size_t)");
		--i_;
		--j_;
		Assert(
			i_ - j_ < m_pComponents[i_].size(),
			"index beyond band in tSymConstRefMatrix::operator()(size_t,size_t)");
		return *(m_pComponents[i_][i_ - j_]);
	}
#else										   // ifndef STRONGCHECK
	bool is0(size_t i_, size_t j_) const
	{
		return ((--i_) < (--j_) ? m_pComponents[j_][j_ - i_] : m_pComponents[i_][i_ - j_]) ==
			   nullptr;
	}
	const real_t& DiagonalComponent(size_t i_) const
	{
		return *(m_pComponents[--i_].front());
	}
	const real_t& operator()(size_t i_, size_t j_) const
	{
		return *((--i_) < (--j_) ? m_pComponents[j_][j_ - i_] : m_pComponents[i_][i_ - j_]);
	}
#endif										   // def STRONGCHECK
	tSymConstRefMatrix& Pack0();			   // redundant
	tSymConstRefMatrix& Swap(size_t, size_t);  // redundant
	tColumnVector& SolveLinAlgSystem(const tConstColumnVector&, tColumnVector&) const;

private:
	const size_t m_size;
	std::vector<const real_t*>* m_pComponents;
};

class tTriangularMatrix
{
public:
	tTriangularMatrix(const tSymConstRefMatrix&);
	~tTriangularMatrix()
	{
		delete[] m_components;
	};
	size_t Dimension() const
	{
		return m_size;
	}
	size_t BandWidth(size_t rowNo_) const
	{
		return m_components[--rowNo_].size() - 1;
	}
#ifdef STRONGCHECK
	const real_t& DiagonalComponent(size_t i_) const
	{
		Assert(i_ > 0 && i_ <= m_size, "invalid index in tTriangularMatrix::operator()(size_t)");
		return m_components[--i_].back();
	}
	const real_t& operator()(size_t i_, size_t j_) const
	{
		Assert(
			i_ > 0 && j_ > 0 && i_ <= m_size && j_ <= m_size,
			"invalid index in tTriangularMatrix::operator()(size_t,size_t)");
		size_t i = m_isLower ? i_ : j_, j = m_isLower ? j_ : i_;
		Assert(
			j <= i && j >= (i - BandWidth(i)),
			"index beyond the band in tTriangularMatrix::operator()(size_t,size_t)");
		--i;
		--j;
		return m_components[i][j + m_components[i].size() - i - 1];
	}
	bool is0(size_t i_, size_t j_) const
	{
		Assert(
			i_ > 0 && j_ > 0 && i_ <= m_size && j_ <= m_size,
			"invalid index in tTriangularMatrix::is0(size_t,size_t)");
		size_t i = m_isLower ? i_ : j_, j = m_isLower ? j_ : i_;
		return (j > i || j < (i - BandWidth(i)));
	}
#else	// ifndef STRONGCHECK
	const real_t& DiagonalComponent(size_t i_) const
	{
		return m_components[--i_].back();
	}
	const real_t& operator()(size_t i_, size_t j_) const
	{
		return (
			--i_,
			--j_,
			(m_isLower ? m_components[i_][j_ + m_components[i_].size() - i_ - 1]
					   : m_components[j_][i_ + m_components[j_].size() - j_ - 1]));
	}
	bool is0(size_t i_, size_t j_) const
	{
		return m_isLower ? (j_ > i_ || j_ < (i_ - BandWidth(i_)))
						 : (i_ > j_ || i_ < (j_ - BandWidth(j_)));
	}
#endif	// def STRONGCHECK
	tTriangularMatrix& Transpose()
	{
		m_isLower = !m_isLower;
		return *this;
	}

	tColumnVector& SolveLinAlgSystem(const tConstColumnVector& right_, tColumnVector& sol_) const
	{
		return m_isLower ? SolveLowerSystem(right_, sol_) : SolveUpperSystem(right_, sol_);
	}

private:
	tColumnVector& SolveLowerSystem(const tConstColumnVector&, tColumnVector&) const;
	tColumnVector& SolveUpperSystem(const tConstColumnVector&, tColumnVector&) const;

	const size_t m_size;
	std::vector<real_t>* m_components;
	bool m_isLower;
};

inline const tNode& tNodalTensor1Column::Node(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_components.size(), "invalid index in tNodalTensor1Column::Node");
#endif
	return m_components[--i_].Node();
}

inline tNodalTensor1Column& tNodalTensor1Column::Assign0(size_t i_)
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && i_ <= m_components.size(),
		"invalid index in tNodalTensor1Column::Assign0(size_t)");
#endif
	m_components[--i_].Assign0();
	return *this;
}

inline bool tNodalTensor1Column::is0(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_components.size(), "invalid index in tNodalTensor1Column::is0");
#endif
	return m_components[--i_].is0();
}

inline tNodalTensor1& tNodalTensor1Column::operator[](size_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ >= 0 && i_ < m_components.size(), "invalid index in tNodalTensor1Column::operator[]");
#endif
	return m_components[i_];
}

inline const tNodalTensor1& tNodalTensor1Column::operator[](size_t i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ >= 0 && i_ < m_components.size(),
		"invalid index in tNodalTensor1Column::operator[] const");
#endif
	return m_components[i_];
}

inline tNodalTensor1& tNodalTensor1Column::operator()(size_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_components.size(), "invalid index in tNodalTensor1Column::operator()");
#endif
	return m_components[--i_];
}

inline const tNodalTensor1& tNodalTensor1Column::operator()(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && i_ <= m_components.size(),
		"invalid index in tNodalTensor1Column::operator() const");
#endif
	return m_components[--i_];
}

inline tNodalTensor1Column& tNodalTensor1Column::Assign0()
{
	std::for_each(m_components.begin(), m_components.end(), std::mem_fn(&tNodalTensor1::Assign0));
	return *this;
}

inline const tNodalTensor1Column& tNodalTensor1Column::LinkAsDisplacements() const
{
	std::for_each(
		m_components.begin(), m_components.end(), std::mem_fn(&tNodalTensor1::LinkAsDispl));
	return *this;
}

/*inline tFEsTensor2Column& tFEsTensor2Column::Assign0()
{
 std::for_each(Components.begin(),Components.end(),std::mem_fn(&tFEsTensor2::Assign0));
// for (vector<tFEsTensor2>::iterator p=Components.begin(); p<Components.end(); ++p)
//           p->Assign0();
 return *this;
}

inline const tFinitElement& tFEsTensor2Column::FE (size_t i_) const
{
#ifdef STRONGCHECK
 Assert(i_> 0  &&  i_<= Components.size(), "invalid index in tFEsTensor2Column::Node");
#endif
 return Components[--i_].FE();
}

inline tFEsTensor2Column& tFEsTensor2Column::Assign0 (size_t i_)
{
#ifdef STRONGCHECK
 Assert(i_> 0  &&  i_<= Components.size(), "invalid index in
tFEsTensor2Column::Assign0(size_t)"); #endif Components[--i_].Assign0(); return *this;
}

inline bool tFEsTensor2Column::is0 (size_t i_) const
{
#ifdef STRONGCHECK
 Assert(i_> 0  &&  i_<= Components.size(), "invalid index in tFEsTensor2Column::is0");
#endif
 return Components[--i_].is0();
}

inline tFEsTensor2& tFEsTensor2Column::operator[](size_t i_)
{
#ifdef STRONGCHECK
 Assert(i_>= 0  &&  i_< Components.size(), "invalid index in tFEsTensor2Column::operator[]");
#endif
 return Components[i_];
}

inline const tFEsTensor2& tFEsTensor2Column::operator[](size_t i_) const
{
#ifdef STRONGCHECK
 Assert(i_>= 0  &&  i_< Components.size(), "invalid index in tFEsTensor2Column::operator[] const");
#endif
 return Components[i_];
}

inline tFEsTensor2& tFEsTensor2Column::operator()(size_t i_)
{
#ifdef STRONGCHECK
 Assert(i_> 0  &&  i_<= Components.size(), "invalid index in tFEsTensor2Column::operator()");
#endif
 return Components[--i_];
}

inline const tFEsTensor2& tFEsTensor2Column::operator()(size_t i_) const
{
#ifdef STRONGCHECK
 Assert(i_> 0  &&  i_<= Components.size(), "invalid index in tFEsTensor2Column::operator() const");
#endif
 return Components[--i_];
}

inline tFEsTensor2Column& tFEsTensor2Column::operator+= (const tFEsTensor2& o_)
{
 operator()(o_.FE()) += o_;
 return *this;
}*/
