#pragma once

#include "NumTypes.hpp"
#include "Tensors.hpp"
#include "FiniteElement.hpp"
#include "HybridTensors.hpp"
#include "NameList.hpp"
#include "Node.hpp"
#include "GaussIntegr.hpp"

#include <vector>
#include <string>
#include <map>
#include <istream>

class tNode;
class tMaterial;
class tAnalysis;

class t1D_FE : virtual public tFinitElement
{
public:
	virtual ~t1D_FE()
	{
	}

	class fCrossSecArea
	{
	public:
		fCrossSecArea(const t1D_FE& fe_) : m_FE(fe_), m_result(0)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return m_result = m_FE.CrosSecArea(x_);
		}

	private:
		mutable real_t m_result;
		const t1D_FE& m_FE;
	};

	Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}
	Tensor2a& Integrate(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const override
	{
		real_t result(0.);
		return Integrate<Tensor1s>(fun_, result);
	}
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const override	// for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}

	size_t SpaceDimension() const override
	{
		return 1;
	}
	Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const override;

	const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const override;
	real_t Jacobian_l2g(const Tensor1s&) const override;

	Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const override;
	Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const override;
	virtual Tensor1s& Tangent(const Tensor1s&, Tensor1s&) const;
	t1D_FE& Set1dStress(bool setToStress_)
	{
		m_stateIs1dStress = setToStress_;
		return *this;
	}
	bool Has1dStressState() const
	{
		return m_stateIs1dStress;
	}
	virtual t1D_FE& SetCrosSecArea(real_t) = 0;
	virtual real_t CrosSecArea(real_t = 0.) const = 0;
	real_t CrosSecArea(const Tensor1s& lc_) const
	{
		return CrosSecArea(lc_(1));
	}
	virtual real_t Length() const
	{
		real_t result(0.);
		const fOne<Tensor1s> one;
		return Integrate1D<Tensor1s>(one, result);
	}

	t1D_FE() : m_stateIs1dStress()
	{
	}

	const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const override;

protected:
	bool m_stateIs1dStress;

private:
	static const fIntegrate<1, GAUSS_ORDER_1D> Integrator;

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D_local(TFun fun_, TFunRet& result_) const  // result_ isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D(TFun fun_, TFunRet& result_) const	 // result_ isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D_local<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D<TFunArg>(MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}

	mutable std::map<real_t, Tensor2s> m_jacobyMatrixDB;
	mutable std::map<real_t, real_t> m_jacobianDB;
};

class t2D_FE : virtual public tFinitElement
{
public:
	virtual ~t2D_FE()
	{
	}
	class fThickness
	{
	public:
		fThickness(const t2D_FE& fe_) : m_FE(fe_), m_result(0)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return m_result = m_FE.Thickness(x_(1), x_(2));
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_, const TArg& y_) const
		{
			return m_result = m_FE.Thickness(x_, y_);
		}

	private:
		mutable real_t m_result;
		const t2D_FE& m_FE;
	};

	Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const override;
	Tensor1s& UnitNormal(const Tensor1s& locoord_, Tensor1s& result_) const
	{
		Tensor1s tmp1, tmp2;
		LocalBasis(locoord_, tmp1, tmp2, result_);
		return result_;
	}
	Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const override;
	Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const override;

	const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const override;
	real_t Jacobian_l2g(const Tensor1s& lc_) const override;

	virtual real_t Integrate2D_local(const fIntegrand<Tensor1s, real_t>& fun_) const = 0;
	virtual real_t Integrate2D(const fIntegrand<Tensor1s, real_t>& fun_) const = 0;

	t2D_FE& Set2dStress(bool setToStress_)
	{
		m_stateIs2dStress = setToStress_;
		return *this;
	}
	bool Has2dStressState() const
	{
		return m_stateIs2dStress;
	}
	virtual t2D_FE& SetThickness(real_t) = 0;
	virtual real_t Thickness(real_t = 0., real_t = 0.) const = 0;
	virtual real_t ComputeArea() const
	{
		const fOne<Tensor1s> one;
		return Integrate2D(one);
	}

	const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const override;

protected:
	t2D_FE() : m_stateIs2dStress(false)
	{
	}

	bool m_stateIs2dStress;

private:
	mutable std::map<std::pair<real_t, real_t>, Tensor2s> m_jacobyMatrixDB;
	mutable std::map<std::pair<real_t, real_t>, real_t> m_jacobianDB;
};

class t3D_FE : virtual public tFinitElement
{
public:
	const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const override;
	real_t Jacobian_l2g(const Tensor1s&) const override;
	Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const override;
	/*     virtual real_t&  Integrate3D(const fIntegrand<real_t>& fun_, real_t&  result_)const
							{return Integrate3D<real_t>(fun_,(result_=0.));}*/
	//     virtual real_t&  Integrate (const fIntegrand<real_t>&  fun_, real_t&  result_) const
	//        {return Integrate3D<real_t>(fun_,(result_=0.));}
	//     virtual Tensor1a& Integrate (const fIntegrand<Tensor1a>& fun_, Tensor1a& result_)const
	//        {return Integrate3D<Tensor1a>(fun_,result_);}
	//     virtual Tensor2a& Integrate_local (const fIntegrand<Tensor1s,Tensor2a>& fun_, Tensor2a&
	//     result_)const
	//        {return result_;/*Integrate3D<Tensor2a>(fun_,result_);*/}

	size_t SpaceDimension() const override
	{
		return 3;
	}

private:
	mutable std::map<std::pair<std::pair<real_t, real_t>, real_t>, Tensor2s> m_jacobyMatrixDB;
	mutable std::map<std::pair<std::pair<real_t, real_t>, real_t>, real_t> m_jacobianDB;
	//     template<typename T> T& Integrate3D (const tFinitElement::fIntegrand<T>&, T&) const;//2nd
	//     parameter must be zero!!! template<typename T> T& Integrate (const
	//     tFinitElement::fIntegrand<T>&, T&) const;//2nd parameter must be zero!!!
};

class tRectangle : public t2D_FE
{
public:
	Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}
	Tensor2a& Integrate(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const override
	{
		real_t result(0.);
		return Integrate<Tensor1s>(fun_, result);
	}
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const override	// for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	real_t Integrate2D_local(const fIntegrand<Tensor1s, real_t>& /*fun_*/) const override
	{
		real_t result(0.);
		const fOne<Tensor1s> one;
		return Integrate2D_local<Tensor1s>(one, result);
	}
	real_t Integrate2D(const fIntegrand<Tensor1s, real_t>& /*fun_*/) const override
	{
		real_t result(0.);
		const fOne<Tensor1s> one;
		return Integrate2D<Tensor1s>(one, result);
	}

	size_t SpaceDimension() const override
	{
		return 2;
	}

private:
	static const fIntegrate<2, GAUSS_ORDER_2D> Integrator;

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D_local(
		TFun fun_, TFunRet& result_) const	// 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D(TFun fun_, TFunRet& result_) const	 // 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate2D_local<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fThickness(*this)), result_);
		//                    return Integrate2D_local<TFunArg>(fun_,result_); //isn't right, for
		//                    experiment only
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate2D<TFunArg>(MultBy<TFunArg, TFunRet>(fun_, fThickness(*this)), result_);
	}
};

#pragma warning(push)
#pragma warning(disable : 4250)
class tIsoRod2ConstSec : public tIsoparametricFE, public t1D_FE
{
public:
	static tFinitElement* NewFE()
	{
		return new tIsoRod2ConstSec();
	}

	const std::string& Kind() const override
	{
		return KindName;
	}
	real_t CoordShapeFun(size_t, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(size_t, const Tensor1s&, Tensor1s&) const override;
	//   virtual real_t   DisplShapeFun             (size_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(size_t, const Tensor1s&, Tensor1s&) const;
	real_t CoordShapeFun(const tNode&, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const override;
	//   virtual real_t   DisplShapeFun             (const tNode&, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;

	//   virtual Tensor1s& Tangent(const Tensor1s&, Tensor1s& result_) const {return
	//   (result_=Node(2).Coord();result_-=Node(1).Coord());} virtual real_t Length() const {return
	//   (Tensor1s d=Node(2).Coord();d-=Node(1).Coord()).Length();}
	t1D_FE& SetCrosSecArea(real_t newval_) override
	{
		m_area = newval_;
		return *this;
	}
	real_t CrosSecArea(real_t = 0.) const override
	{
		return m_area;
	}
	//   virtual real_t Volume() const {return Length()*CrosSecArea();}

private:
	friend class tPtrShapeFunArray;
	friend class tPtrShapeFunGradArray;

	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray() : tIsoparametricFE::tArrayOfPtrsToShapeFunctions(2)
		{
			operator[](0) = ShapeFun1;
			operator[](1) = ShapeFun2;
		}
	};
	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray() : tIsoparametricFE::tArrayOfPtrsToShapeFunGradients(2)
		{
			operator[](0) = ShapeGrad1;
			operator[](1) = ShapeGrad2;
		}
	};

	static const std::string& KindName;
	static const tPtrShapeFunArray pShapeFunctions;
	static const tPtrShapeFunGradArray pShapeFunGrads;

	static real_t ShapeFun1(const Tensor1s&);
	static real_t ShapeFun2(const Tensor1s&);
	static Tensor1s& ShapeGrad1(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad2(const Tensor1s&, Tensor1s&);

	explicit tIsoRod2ConstSec() : tFinitElement(2), m_area(0)
	{
	}

	real_t m_area;
};

class tIsoQuad4ConsThick : public tIsoparametricFE, public tRectangle
{
public:
	static tFinitElement* NewFE()
	{
		return new tIsoQuad4ConsThick();
	}

	real_t CoordShapeFun(size_t, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(size_t, const Tensor1s&, Tensor1s&) const override;
	//   virtual real_t   DisplShapeFun             (size_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(size_t, const Tensor1s&, Tensor1s&) const;
	real_t CoordShapeFun(const tNode&, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const override;
	//   virtual real_t   DisplShapeFun             (const tNode&, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;

	const std::string& Kind() const override
	{
		return KindName;
	}
	//   virtual real_t Volume() const {return Area() * ThicknessValue;}
	t2D_FE& SetThickness(real_t newval_) override
	{
		m_thicknessValue = newval_;
		return *this;
	}
	real_t Thickness(real_t = 0., real_t = 0.) const override
	{
		return m_thicknessValue;
	}

private:
	friend class tPtrShapeFunArray;
	friend class tPtrShapeFunGradArray;

	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};

	static const std::string& KindName;
	static const tPtrShapeFunArray pShapeFunctions;
	static const tPtrShapeFunGradArray pShapeFunGrads;

	static real_t ShapeFun1(const Tensor1s&);
	static real_t ShapeFun2(const Tensor1s&);
	static real_t ShapeFun3(const Tensor1s&);
	static real_t ShapeFun4(const Tensor1s&);
	static Tensor1s& ShapeGrad1(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad2(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad3(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad4(const Tensor1s&, Tensor1s&);

	explicit tIsoQuad4ConsThick() : tFinitElement(4), m_thicknessValue(0)
	{
	}

	real_t m_thicknessValue;
};

class tIsoParallelepiped8 : public tIsoparametricFE, public t3D_FE
{
public:
	static tFinitElement* NewFE()
	{
		return new tIsoParallelepiped8();
	}

	Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}
	Tensor2a& Integrate(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const override
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const override
	{
		real_t result(0.);
		return Integrate<Tensor1s>(fun_, result);
	}
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const override	// for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}

	real_t CoordShapeFun(size_t, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(size_t, const Tensor1s&, Tensor1s&) const override;
	real_t CoordShapeFun(const tNode&, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const override;

	const std::string& Kind() const override
	{
		return KindName;
	}

private:
	friend class tPtrShapeFunArray;
	friend class tPtrShapeFunGradArray;

	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};

	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};

	static const std::string& KindName;
	static const fIntegrate<3, GAUSS_ORDER_3D> Integrator;
	static const tPtrShapeFunArray pShapeFunctions;
	static const tPtrShapeFunGradArray pShapeFunGrads;

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D_local(
		TFun fun_, TFunRet& result_) const	// 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D(TFun fun_, TFunRet& result_) const	 // 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}

	/*   class fOne
		 {
			private:
				mutable real_t Result;
			public:
//        fOne() {}

		 template <typename TArg>
					real_t& operator()(const TArg& x_) const {return Result
		 = 1.;}
				template <typename TArg>
					real_t& operator()(const TArg& x_, const
		 TArg& y_) const {return Result = 1.;}
		 }; */
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D_local<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fOne<Tensor1s>()), result_);
		//                    return Integrate3D_local<TFunArg>(fun_,result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D<TFunArg>(MultBy<TFunArg, TFunRet>(fun_, fOne<Tensor1s>()), result_);
		//                    return Integrate3D<TFunArg>(fun_,result_);
	}

	static real_t ShapeFun1(const Tensor1s&);
	static real_t ShapeFun2(const Tensor1s&);
	static real_t ShapeFun3(const Tensor1s&);
	static real_t ShapeFun4(const Tensor1s&);
	static real_t ShapeFun5(const Tensor1s&);
	static real_t ShapeFun6(const Tensor1s&);
	static real_t ShapeFun7(const Tensor1s&);
	static real_t ShapeFun8(const Tensor1s&);

	static Tensor1s& ShapeGrad1(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad2(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad3(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad4(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad5(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad6(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad7(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad8(const Tensor1s&, Tensor1s&);

	explicit tIsoParallelepiped8() : tFinitElement(8)
	{
	}
};
#pragma warning(pop)

inline const tNode& tFinitElement::Node(size_t nodeNo_) const
{
#ifdef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_pNodes.size(), "Invalid node # in FE");
#endif
	return *(m_pNodes[--nodeNo_]);
}

inline tFinitElement& tFinitElement::Link(const tMaterial& matlToLink_)
{
#ifdef STRONGCHECK
	Assert(m_pMaterial != &matlToLink_, "Invalid material in FE");
#endif
	// if (pMaterial == &matlToLink_) return *this;
	m_pMaterial = &matlToLink_;
	return *this;
}

inline void tNonIsoparametricFE::DefineCoordShapeFunAndItsGrad(
	size_t number_,
	real_t (*pfun_)(const Tensor1s&),
	Tensor1s& (*pgrad_)(const Tensor1s&, Tensor1s&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= m_coordApproxData.size(),
		"invalid number of shape function in tFinitElement::DefineCoordShapeFunAndItsGrad");
	Assert(
		m_coordApproxData[number_ - 1].m_pFunction == nullptr,
		"shape function is already defined in tFinitElement::DefineCoordShapeFunAndItsGrad");
	Assert(
		m_coordApproxData[number_ - 1].m_pFunctionGrad == nullptr,
		"shape gradient is already defined in tFinitElement::DefineCoordShapeFunAndItsGrad");
#endif
	m_coordApproxData[--number_].m_pFunction = pfun_;
	m_coordApproxData[number_].m_pFunctionGrad = pgrad_;
}

inline void tNonIsoparametricFE::DefineDisplShapeFunAndItsGrad(
	size_t number_,
	real_t (*pfun_)(const Tensor1s&),
	Tensor1s& (*pgrad_)(const Tensor1s&, Tensor1s&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= m_displApproxData.size(),
		"invalid number of shape function in tFinitElement::DefineDisplShapeFunAndItsGrad");
	Assert(
		m_displApproxData[number_ - 1].m_pFunction == nullptr,
		"shape function is already defined in tFinitElement::DefineDisplShapeFunAndItsGrad");
	Assert(
		m_displApproxData[number_ - 1].m_pFunctionGrad == nullptr,
		"shape gradient is already defined in tFinitElement::DefineDisplShapeFunAndItsGrad");
#endif
	m_displApproxData[--number_].m_pFunction = pfun_;
	m_displApproxData[number_].m_pFunctionGrad = pgrad_;
}

inline real_t tNonIsoparametricFE::CoordShapeFun(const tNode& node_, const Tensor1s& lc_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	//       p =
	//       std::find_if(CoordApproxData.begin(),CoordApproxData.end(),tNodeAnd2FunsPtrs::Node_address_is(&node_));
	Assert(
		std::find_if(
			m_coordApproxData.begin(),
			m_coordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_coordApproxData.end(),
		"There are no shape functions corresponding to the node in tFinitElement::CoordShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  m_coordApproxData.begin(),
				  m_coordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->m_pFunction))(lc_);
}

inline Tensor1s& tNonIsoparametricFE::CoordShapeGrad(
	const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_coordApproxData.begin(),
			m_coordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_coordApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement::CoordShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  m_coordApproxData.begin(),
				  m_coordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->m_pFunctionGrad))(lc_, result_);
}

inline real_t tNonIsoparametricFE::DisplShapeFun(const tNode& node_, const Tensor1s& lc_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_displApproxData.begin(),
			m_displApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_displApproxData.end(),
		"There are no shape functions corresponding to the node in tFinitElement::DisplShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  m_displApproxData.begin(),
				  m_displApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->m_pFunction))(lc_);
}

inline Tensor1s& tNonIsoparametricFE::DisplShapeGrad(
	const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_displApproxData.begin(),
			m_displApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_displApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement::DisplShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  m_displApproxData.begin(),
				  m_displApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->m_pFunctionGrad))(lc_, result_);
}
