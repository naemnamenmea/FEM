#pragma once

#include "NumTypes.h"
#include "Tensors.h"
//#include "HybridTensors.h"
#include "NameList.h"
#include "Node.hpp"
#include "GaussIntegr.hpp"
#include <vector>
#include <string>
#include <map>
#include <istream>

using namespace GaussIntegr;

template <typename T>
class Node3d;
template <typename T>
class tMaterial;
class tAnalysis;

template <typename T>
class tFinitElement
{
protected:
	static const tNamesToFuns<tFinitElement<T>> Factory;

	std::vector<const Node3d<T>*> pNodes;
	const tMaterial<T>* pMaterial;

private:
	mutable T StrainLevel;
	//   mutable const Tensor2s<T> *pDeformData;

public:
	//   void LinkDeformData(const Tensor2s<T>* p_) const {pDeformData = p_;}
	//   void UnLinkDeformData() const {pDeformData = nullptr;}

	typename std::vector<const Node3d<T>*>::const_iterator p1stNode() const
	{
		return pNodes.begin();
	}  // for use in T tIsoparametricFE::tArray<T>::operator()

	static tFinitElement<T>* NewFE(const std::string& kindName_)
	{
		return tFinitElement<T>::Factory.CallFunction(kindName_);
	}

protected:
	tFinitElement() : pMaterial(nullptr), StrainLevel(0.)
	{
	}
	tFinitElement(cardinal_t numberOfNodes_)
		: pNodes(numberOfNodes_, nullptr), pMaterial(nullptr), StrainLevel(0.)
	{
	}
	void DefineNextNode(const Node3d<T>&);

public:
	virtual ~tFinitElement()
	{
	}
	virtual const std::string& Kind() const = 0;
	virtual tFinitElement<T>& DefineNextNode(const Node3d<T>&, bool, bool) = 0;
	cardinal_t HowManyNodes() const
	{
		return static_cast<cardinal_t>(pNodes.size());
	}
	virtual cardinal_t HowManyCoordApproxNodes() const = 0;
	virtual cardinal_t HowManyDisplApproxNodes() const = 0;
	virtual cardinal_t SpaceDimension() const = 0;
	const Node3d<T>& Node(cardinal_t) const;
	virtual const Node3d<T>& CoordApproxNode(cardinal_t) const = 0;
	virtual const Node3d<T>& DisplApproxNode(cardinal_t) const = 0;
	bool Has(const Node3d<T>& o_) const
	{
		return std::find(pNodes.begin(), pNodes.end(), &o_) > pNodes.end();
	}
	tFinitElement<T>& Link(const tMaterial<T>&);
	const tMaterial<T>& GetMaterial() const
	{
		return *pMaterial;
	}
	tMaterial<T>& GetMaterial()
	{
		return *pMaterial;
	}
	bool Has(const tMaterial<T>& o_) const
	{
		return pMaterial == &o_;
	}
	virtual const SymmetricTensor4s<T>& MaterialElasTensor(const Tensor1s<T>&) const;

	virtual T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	virtual T DisplShapeFun(cardinal_t, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& DisplShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& CoordShapeGrad(
		const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	virtual T DisplShapeFun(const Node3d<T>&, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& DisplShapeGrad(
		const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const = 0;

protected:
	typedef std::map<std::pair<T, std::pair<T, T>>, Tensor2s<T>> tMapRotToLoc;
	mutable tMapRotToLoc RotToLoc_DB;
	Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&, Tensor2s<T>&)
		const;	// be careful - do not call for 1D or 2D - FEs!!! (because singular)
public:
	virtual const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const = 0;
	virtual T Jacobian_l2g(const Tensor1s<T>&) const = 0;
	class fJacobian_l2g
	{
	private:
		mutable T Result;
		const tFinitElement<T>& FE;

	public:
		fJacobian_l2g(const tFinitElement<T>& fe_) : FE(fe_), Result()
		{
		}
		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return Result = FE.Jacobian_l2g(x_);
		}
	};

	virtual Tensor1s<T>& GlobalCoord(const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual Tensor1s<T>& LocalBasis(
		const Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&) const;
	virtual Tensor2s<T>& TensorOfRotationToLocal(const Tensor1s<T>&, Tensor2s<T>&) const = 0;
	const Tensor2s<T>& TensorOfRotationToLocal(const Tensor1s<T>&) const;

	template <typename TArg, typename TRet>
	struct fIntegrand
	{
		virtual TRet& operator()(const TArg&) const = 0;
	};
	template <typename TArg>
	struct fOne : public fIntegrand<TArg, T>
	{
		mutable T Result;
		fOne() : Result(1.)
		{
		}
		T& operator()(const TArg&) const
		{
			return Result = 1.;
		}
		//         T operator()(const TArg&) const {return 1.;}
		template <typename TArg2>
		T& operator()(const TArg& x_, const TArg2& y_) const
		{
			return Result = 1.;
		}  // need for tParallelepiped1
	};
	template <typename TArg, typename TRet, typename TMult, typename TFac>
	class fMultWithAssgn : public tFinitElement<T>::fIntegrand<TArg, TRet>
	{
	private:
		TMult& Multiplicand;
		const TFac& Factor;

	public:
		fMultWithAssgn(TMult& m_, const TFac& f_) : Multiplicand(m_), Factor(f_)
		{
		}
		TRet& operator()(const TArg& x_) const
		{
			return Multiplicand(x_) *= Factor(x_);
		}
	};
	template <typename TArg, typename TRet, typename TMult, typename TFac>
	static fMultWithAssgn<TArg, TRet, TMult, TFac> MultBy(TMult& m_, const TFac& f_)
	{
		return fMultWithAssgn<TArg, TRet, TMult, TFac>(m_, f_);
	}

	virtual Tensor2a<T>& Integrate_local(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>&, Tensor2a<T>&) const = 0;
	virtual Tensor2a<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>&, Tensor2a<T>&) const = 0;
	virtual T Integrate(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>&) const = 0;

	//   virtual Tensor2s<T>& Integrate (const
	//   tFinitElement<T>::fIntegrand<Tensor1s<T>,Tensor2s<T>>&,Tensor2s<T>&) const =0;//for strain
	virtual SymmetricTensor2s<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, SymmetricTensor2s<T>>&,
		SymmetricTensor2s<T>&) const = 0;  // for strain

	//   tFEsSetOfTensor2& Calculate(tFinitElement<T>::fIntegrand<tFEsSetOfTensor2>& fun_,
	//   tEvalMethod em_, tFEsSetOfTensor2& result_) const {return em_==atCentre?
	//   result_.Swap(fun_(Tensor1s<T>(0.))) : Integrate(fun_,result_) /= Volume();} Tensor1&
	//   Displacement (const Tensor1&, Tensor1&) const; Tensor2s<T>& DisplacementGradient (const
	//   Tensor1s<T>&, Tensor2s<T>&) const; tFEsTensor2& DisplacementGradient (const Tensor1s<T>&,
	//   tFEsTensor2&) const;
	virtual T Volume() const
	{
		const fOne<Tensor1s<T>> one;
		return Integrate(one);
	}
	//   T Mass() const;//tMaterial declared only {return Volume() * Material().Density();}
public:
	Tensor2s<T>& NablaDispl_g(const Tensor1s<T>&, Tensor2s<T>&) const;
	//   Tensor2s<T>& NodalStrain_g(cardinal_t, const Tensor1s<T>&, Tensor2s<T>&) const;
	//   Tensor2& Strain(Tensor2& result_) const {return Strain(Tensor1s<T>(0.),result_);}//
};

template <typename T>
class tNonIsoparametricFE : virtual public tFinitElement<T>
{
private:
	//===Data for providing calls of shape functions and its gradients:
	struct tNodeAnd2FunsPtrs
	{
		const Node3d<T>* pNode;
		T (*pFunction)(const Tensor1s<T>&);
		Tensor1s<T>& (*pFunctionGrad)(const Tensor1s<T>&, Tensor1s<T>&);
		tNodeAnd2FunsPtrs() : pNode(nullptr), pFunction(nullptr), pFunctionGrad(nullptr)
		{
		}
		struct Node_address_is
		{
			const Node3d<T>* AddressToFind;
			Node_address_is(const Node3d<T>* pnode_) : AddressToFind(pnode_)
			{
			}
			bool operator()(const tNodeAnd2FunsPtrs& valueToTest_)
			{
				return valueToTest_.pNode == AddressToFind;
			}
		};
	};
	std::vector<tNodeAnd2FunsPtrs> CoordApproxData, DisplApproxData;

protected:
	void DefineCoordShapeFunAndItsGrad(
		cardinal_t, T (*)(const Tensor1s<T>&), Tensor1s<T>& (*)(const Tensor1s<T>&, Tensor1s<T>&));
	void DefineDisplShapeFunAndItsGrad(
		cardinal_t, T (*)(const Tensor1s<T>&), Tensor1s<T>& (*)(const Tensor1s<T>&, Tensor1s<T>&));
	//===end of data for providing calls of shape functions and its gradients
	tNonIsoparametricFE(cardinal_t coordNod_, cardinal_t displNod_)
		: CoordApproxData(coordNod_), DisplApproxData(displNod_)
	{
	}

public:
	virtual tFinitElement<T>& DefineNextNode(const Node3d<T>&, bool, bool);
	cardinal_t HowManyCoordApproxNodes() const
	{
		return CoordApproxData.size();
	}
	cardinal_t HowManyDisplApproxNodes() const
	{
		return DisplApproxData.size();
	}
	virtual const Node3d<T>& CoordApproxNode(cardinal_t) const;
	virtual const Node3d<T>& DisplApproxNode(cardinal_t) const;
	//   virtual T   CoordShapeFun             (cardinal_t, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	//   virtual T   DisplShapeFun             (cardinal_t, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& DisplShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual T DisplShapeFun(const Node3d<T>&, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& DisplShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const;
};

template <typename T>
class tIsoparametricFE : virtual public tFinitElement<T>
{
	//===Data for providing calls of shape functions and its gradients:
private:
	template <typename tTYPE>
	class tArray : public std::vector<tTYPE>
	{
	public:
		tArray(cardinal_t dim_) : std::vector<tTYPE>(dim_)
		{
		}

		tTYPE operator()(cardinal_t) const;
		tTYPE operator()(const tFinitElement<T>*, const Node3d<T>*) const;
	};
	typedef T (*pShapeFun)(const Tensor1a<T>&);
	typedef Tensor1a<T>& (*pShapeFunGrad)(const Tensor1a<T>&, Tensor1a<T>&);

protected:
	typedef tArray<T (*)(const Tensor1s<T>&)> tArrayOfPtrsToShapeFunctions;
	typedef tArray<Tensor1s<T>& (*)(const Tensor1s<T>&, Tensor1s<T>&)>
		tArrayOfPtrsToShapeFunGradients;
	// friend   T(*)(const Tensor1s<T>&)          tArrayOfPtrsToShapeFunctions::operator()(const
	// tFinitElement<T>*,const Node3d<T>*) const; friend Tensor1&(*)(const Tensor1s<T>&,Tensor1&)
	// tArrayOfPtrsToShapeFunGradients::operator()(const tFinitElement<T>*,const Node3d<T>*) const;

	//**************** ERROR in g++ 4.0.2 **********************
	// friend class tArrayOfPtrsToShapeFunctions;
	// friend class tArrayOfPtrsToShapeFunGradients;
	friend class tArray<T (*)(const Tensor1s<T>&)>;
	friend class tArray<Tensor1s<T>& (*)(const Tensor1s<T>&, Tensor1s<T>&)>;

	//===end of data for providing calls of shape functions and its gradients

	tIsoparametricFE()
	{
	}

public:
	virtual cardinal_t HowManyCoordApproxNodes() const
	{
		return this->HowManyNodes();
	}
	virtual cardinal_t HowManyDisplApproxNodes() const
	{
		return this->HowManyNodes();
	}
	virtual const Node3d<T>& CoordApproxNode(cardinal_t i_) const
	{
		return this->Node(i_);
	}
	virtual const Node3d<T>& DisplApproxNode(cardinal_t i_) const
	{
		return this->Node(i_);
	}
#ifdef STRONGCHECK
	virtual tFinitElement<T>& DefineNextNode(const Node3d<T>&, bool = true, bool = true);
#else	// ifndef STRONGCHECK
	virtual tFinitElement<T>& DefineNextNode(const Node3d<T>& node_, bool = true, bool = true)
	{
		tFinitElement<T>::DefineNextNode(node_);
		return *this;
	}
#endif	// def STRONGCHECK
	virtual T DisplShapeFun(cardinal_t n_, const Tensor1s<T>& c_) const
	{
		return this->CoordShapeFun(n_, c_);
	}
	virtual Tensor1s<T>& DisplShapeGrad(cardinal_t n_, const Tensor1s<T>& c_, Tensor1s<T>& r_) const
	{
		return this->CoordShapeGrad(n_, c_, r_);
	}
	virtual T DisplShapeFun(const Node3d<T>& n_, const Tensor1s<T>& c_) const
	{
		return this->CoordShapeFun(n_, c_);
	}
	virtual Tensor1s<T>& DisplShapeGrad(
		const Node3d<T>& n_, const Tensor1s<T>& c_, Tensor1s<T>& r_) const
	{
		return this->CoordShapeGrad(n_, c_, r_);
	}
};

template <typename T>
class t1D_FE : virtual public tFinitElement<T>
{
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
			this->MultBy<TFunArg, TFunRet>(fun_, tFinitElement<T>::fJacobian_l2g(*this)), result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}

public:
	t1D_FE() : StateIs1dStress()
	{
	}

	virtual Tensor2a<T>& Integrate_local(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}
	virtual Tensor2a<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}
	virtual T Integrate(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		T result(0.);
		return Integrate<Tensor1s<T>>(fun_, result);
	}
	//   virtual Tensor2s<T>& Integrate (const
	//   tFinitElement<T>::fIntegrand<Tensor1s<T>,Tensor2s<T>>& fun_, Tensor2s<T>& result_)
	//   const//for strain
	virtual SymmetricTensor2s<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, SymmetricTensor2s<T>>& fun_,
		SymmetricTensor2s<T>& result_) const  // for strain
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}

protected:
	bool StateIs1dStress;

public:
	virtual cardinal_t SpaceDimension() const
	{
		return 1;
	}
	virtual Tensor1s<T>& GlobalCoord(const Tensor1s<T>&, Tensor1s<T>&) const;

private:
	mutable std::map<T, Tensor2s<T>> JacobyMatrix_DB;
	mutable std::map<T, T> Jacobian_DB;

public:
	virtual const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const;
	virtual T Jacobian_l2g(const Tensor1s<T>&) const;

	virtual Tensor1s<T>& LocalBasis(
		const Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&) const;
	virtual Tensor2s<T>& TensorOfRotationToLocal(const Tensor1s<T>&, Tensor2s<T>&) const;
	virtual Tensor1s<T>& Tangent(const Tensor1s<T>&, Tensor1s<T>&) const;
	t1D_FE& Set1dStress(bool setToStress_)
	{
		StateIs1dStress = setToStress_;
		return *this;
	}
	bool Has1dStressState() const
	{
		return StateIs1dStress;
	}
	virtual t1D_FE& SetCrosSecArea(T) = 0;
	virtual T CrosSecArea(T = 0.) const = 0;
	T CrosSecArea(const Tensor1s<T>& lc_) const
	{
		return CrosSecArea(lc_(1));
	}
	class fCrossSecArea
	{
	private:
		mutable T Result;
		const t1D_FE& FE;

	public:
		fCrossSecArea(const t1D_FE& fe_) : FE(fe_), Result()
		{
		}
		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return Result = FE.CrosSecArea(x_);
		}
	};
	virtual T Length() const
	{
		T result(0.);
		const tFinitElement<T>::template fOne<Tensor1s<T>> one;
		return Integrate1D<Tensor1s<T>>(one, result);
	}

	virtual const SymmetricTensor4s<T>& MaterialElasTensor(const Tensor1s<T>&) const;
};

template <typename T>
class t2D_FE : virtual public tFinitElement<T>
{
protected:
	bool StateIs2dStress;
	t2D_FE() : StateIs2dStress(false)
	{
	}

public:
	virtual Tensor1s<T>& LocalBasis(
		const Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&) const;
	Tensor1s<T>& UnitNormal(const Tensor1s<T>& locoord_, Tensor1s<T>& result_) const
	{
		Tensor1s<T> tmp1, tmp2;
		LocalBasis(locoord_, tmp1, tmp2, result_);
		return result_;
	}
	virtual Tensor1s<T>& GlobalCoord(const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual Tensor2s<T>& TensorOfRotationToLocal(const Tensor1s<T>&, Tensor2s<T>&) const;

private:
	mutable std::map<std::pair<T, T>, Tensor2s<T>> JacobyMatrix_DB;
	mutable std::map<std::pair<T, T>, T> Jacobian_DB;

public:
	virtual const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const;
	virtual T Jacobian_l2g(const Tensor1s<T>& lc_) const;

	virtual T Integrate2D_local(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const = 0;
	virtual T Integrate2D(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const = 0;

	t2D_FE& Set2dStress(bool setToStress_)
	{
		StateIs2dStress = setToStress_;
		return *this;
	}

	bool Has2dStressState() const
	{
		return StateIs2dStress;
	}

	virtual t2D_FE& SetThickness(T) = 0;

	virtual T Thickness(T = 0., T = 0.) const = 0;

	class fThickness
	{
	private:
		mutable T Result;
		const t2D_FE& FE;

	public:
		fThickness(const t2D_FE& fe_) : FE(fe_), Result()
		{
		}
		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return Result = FE.Thickness(x_(1), x_(2));
		}
		template <typename TArg>
		T& operator()(const TArg& x_, const TArg& y_) const
		{
			return Result = FE.Thickness(x_, y_);
		}
	};
	virtual T Area() const
	{
		const tFinitElement<T>::template fOne<Tensor1s<T>> one;
		return Integrate2D(one);
	}

	virtual const SymmetricTensor4s<T>& MaterialElasTensor(const Tensor1s<T>&) const;
};

template <typename T>
class t3D_FE : virtual public tFinitElement<T>
{
private:
	mutable std::map<std::pair<std::pair<T, T>, T>, Tensor2s<T>> JacobyMatrix_DB;
	mutable std::map<std::pair<std::pair<T, T>, T>, T> Jacobian_DB;
	//     template<typename T> T& Integrate3D (const
	//     tFinitElement<T>::tFinitElement<T>::fIntegrand<T>&, T&) const;//2nd parameter must be
	//     zero!!! template<typename T> T& Integrate (const
	//     tFinitElement<T>::tFinitElement<T>::fIntegrand<T>&, T&) const;//2nd parameter must be
	//     zero!!!
public:
	virtual const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const;
	virtual T Jacobian_l2g(const Tensor1s<T>&) const;
	virtual Tensor2s<T>& TensorOfRotationToLocal(const Tensor1s<T>&, Tensor2s<T>&) const;

	/*
	 virtual T&  Integrate3D(const tFinitElement<T>::fIntegrand<T>& fun_, T&  result_)const
		{return Integrate3D<T>(fun_,(result_=0.));}
	 virtual T&  Integrate (const tFinitElement<T>::fIntegrand<T>&  fun_, T&  result_) const
	   {return Integrate3D<T>(fun_,(result_=0.));}
	 virtual Tensor1a<T>& Integrate (const tFinitElement<T>::fIntegrand<Tensor1a<T>>& fun_,
	 Tensor1a<T>& result_)const
	   {return Integrate3D<Tensor1a<T>>(fun_,result_);}
	 virtual Tensor2a<T>& Integrate_local (const
	 tFinitElement<T>::fIntegrand<Tensor1s<T>,Tensor2a<T>>& fun_, Tensor2a<T>& result_)const
	   {return result_;/*Integrate3D<Tensor2a<T>>(fun_,result_);}
	*/

	virtual cardinal_t SpaceDimension() const
	{
		return 3;
	}
};

template <typename T>
class tRectangle : public t2D_FE<T>
{
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
			this->MultBy<TFunArg, TFunRet>(fun_, tFinitElement<T>::fJacobian_l2g(*this)), result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return this->Integrate2D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, t2D_FE<T>::fThickness(*this)), result_);
		//                    return Integrate2D_local<TFunArg>(fun_,result_); //isn't right, for
		//                    experiment only
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return this->Integrate2D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, t2D_FE<T>::fThickness(*this)), result_);
	}

public:
	virtual Tensor2a<T>& Integrate_local(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}
	virtual Tensor2a<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}
	virtual T Integrate(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		T result(0.);
		return Integrate<Tensor1s<T>>(fun_, result);
	}
	//   virtual Tensor2s<T>& Integrate (const
	//   tFinitElement<T>::fIntegrand<Tensor1s<T>,Tensor2s<T>>& fun_, Tensor2s<T>& result_)
	//   const//for strain
	virtual SymmetricTensor2s<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, SymmetricTensor2s<T>>& fun_,
		SymmetricTensor2s<T>& result_) const  // for strain
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}
	virtual T Integrate2D_local(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		fun_;
		T result(0.);
		const tFinitElement<T>::template fOne<Tensor1s<T>> one;
		return Integrate2D_local<Tensor1s<T>>(one, result);
	}
	virtual T Integrate2D(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		fun_;
		T result(0.);
		const tFinitElement<T>::template fOne<Tensor1s<T>> one;
		return Integrate2D<Tensor1s<T>>(one, result);
	}

	virtual cardinal_t SpaceDimension() const
	{
		return 2;
	}
};

#pragma warning(push)
#pragma warning(disable : 4250)
template <typename T>
class tIsoRod2ConstSec : public tIsoparametricFE<T>, public t1D_FE<T>
{
private:
	static const std::string& KindName;

	static T ShapeFun1(const Tensor1s<T>&);
	static T ShapeFun2(const Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad1(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad2(const Tensor1s<T>&, Tensor1s<T>&);
	class tPtrShapeFunArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray() : tIsoparametricFE::tArrayOfPtrsToShapeFunctions(2)
		{
			operator[](0) = ShapeFun1;
			operator[](1) = ShapeFun2;
		}
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;
	class tPtrShapeFunGradArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray() : tIsoparametricFE::tArrayOfPtrsToShapeFunGradients(2)
		{
			operator[](0) = ShapeGrad1;
			operator[](1) = ShapeGrad2;
		}
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	T Area;
	explicit tIsoRod2ConstSec() : tFinitElement<T>(2), Area()
	{
	}

public:
	virtual const std::string& Kind() const
	{
		return KindName;
	}
	virtual T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	//   virtual T   DisplShapeFun             (cardinal_t, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& DisplShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const;
	//   virtual T   DisplShapeFun             (const Node3d<T>&, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& DisplShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&)
	//   const;
	static tFinitElement<T>* NewFE()
	{
		return new tIsoRod2ConstSec();
	}

	//   virtual Tensor1s<T>& Tangent(const Tensor1s<T>&, Tensor1s<T>& result_) const {return
	//   (result_=Node(2).Coord();result_-=Node(1).Coord());} virtual T Length() const {return
	//   (Tensor1s<T> d=Node(2).Coord();d-=Node(1).Coord()).Length();}
	virtual t1D_FE<T>& SetCrosSecArea(T newval_)
	{
		Area = newval_;
		return *this;
	}
	virtual T CrosSecArea(T = 0.) const
	{
		return Area;
	}
	//   virtual T Volume() const {return Length()*CrosSecArea();}
};

template <typename T>
class tIsoQuad4ConsThick : public tIsoparametricFE<T>, public tRectangle<T>
{
private:
	static const std::string& KindName;

	static T ShapeFun1(const Tensor1s<T>&);
	static T ShapeFun2(const Tensor1s<T>&);
	static T ShapeFun3(const Tensor1s<T>&);
	static T ShapeFun4(const Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad1(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad2(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad3(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad4(const Tensor1s<T>&, Tensor1s<T>&);
	class tPtrShapeFunArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;
	class tPtrShapeFunGradArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	T ThicknessValue;
	explicit tIsoQuad4ConsThick() : tFinitElement<T>(4), ThicknessValue()
	{
	}

public:
	virtual T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	//   virtual T   DisplShapeFun             (cardinal_t, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& DisplShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const;
	//   virtual T   DisplShapeFun             (const Node3d<T>&, const Tensor1s<T>&) const;
	//   virtual Tensor1s<T>& DisplShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&)
	//   const;

	static tFinitElement<T>* NewFE()
	{
		return new tIsoQuad4ConsThick();
	}

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	//   virtual T Volume() const {return Area() * ThicknessValue;}
	virtual t2D_FE<T>& SetThickness(T newval_)
	{
		ThicknessValue = newval_;
		return *this;
	}
	virtual T Thickness(T = 0., T = 0.) const
	{
		return ThicknessValue;
	}
};

template <typename T>
class tIsoParallelepiped8 : public tIsoparametricFE<T>, public t3D_FE<T>
{
private:
	static const std::string& KindName;
	static const fIntegrate<3, GAUSS_ORDER_3D> Integrator;

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
			this->MultBy<TFunArg, TFunRet>(fun_, tFinitElement<T>::fJacobian_l2g(*this)), result_);
	}

	/*
	class fOne
		{
		private:
		mutable T Result;
		public:
		template <typename TArg>
			T& operator()(const TArg& x_) const {return Result = 1.;}
		template <typename TArg>
			T& operator()(const TArg& x_, const TArg& y_) const {return Result = 1.;}
		};
	*/

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, tFinitElement<T>::template fOne<Tensor1s<T>>()),
			result_);
		//                    return Integrate3D_local<TFunArg>(fun_,result_);
	}
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, tFinitElement<T>::template fOne<Tensor1s<T>>()),
			result_);
		//                    return Integrate3D<TFunArg>(fun_,result_);
	}

public:
	virtual Tensor2a<T>& Integrate_local(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}
	virtual Tensor2a<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}
	virtual T Integrate(const tFinitElement<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		T result(0.);
		return Integrate<Tensor1s<T>>(fun_, result);
	}
	//   virtual Tensor2s<T>& Integrate (const
	//   tFinitElement<T>::fIntegrand<Tensor1s<T>,Tensor2s<T>>& fun_, Tensor2s<T>& result_)
	//   const//for strain
	virtual SymmetricTensor2s<T>& Integrate(
		const tFinitElement<T>::fIntegrand<Tensor1s<T>, SymmetricTensor2s<T>>& fun_,
		SymmetricTensor2s<T>& result_) const  // for strain
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}

private:
	static T ShapeFun1(const Tensor1s<T>&);
	static T ShapeFun2(const Tensor1s<T>&);
	static T ShapeFun3(const Tensor1s<T>&);
	static T ShapeFun4(const Tensor1s<T>&);
	static T ShapeFun5(const Tensor1s<T>&);
	static T ShapeFun6(const Tensor1s<T>&);
	static T ShapeFun7(const Tensor1s<T>&);
	static T ShapeFun8(const Tensor1s<T>&);

	static Tensor1s<T>& ShapeGrad1(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad2(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad3(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad4(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad5(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad6(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad7(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad8(const Tensor1s<T>&, Tensor1s<T>&);

	class tPtrShapeFunArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;

	class tPtrShapeFunGradArray : public tIsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	explicit tIsoParallelepiped8() : tFinitElement<T>(8)
	{
	}

public:
	virtual T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const;
	virtual Tensor1s<T>& CoordShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const;

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	static tFinitElement<T>* NewFE()
	{
		return new tIsoParallelepiped8();
	}
};
#pragma warning(pop)

- tFinitElement -

	template <typename T>
	inline const Node3d<T>& tFinitElement<T>::Node(cardinal_t nodeNo_)const
{
#ifdef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_pNodes.size(), "Invalid node # in FE");
#endif
	return *(pNodes[--nodeNo_]);
}

template <typename T>
inline tFinitElement<T>& tFinitElement<T>::Link(const tMaterial<T>& matlToLink_)
{
#ifdef STRONGCHECK
	Assert(m_pMaterial != &matlToLink_, "Invalid material in FE");
#endif
	// if (pMaterial == &matlToLink_) return *this;
	pMaterial = &matlToLink_;
	return *this;
}

- tNonIsoparametricFE

	template <typename T>
	inline void tNonIsoparametricFE<T>::DefineCoordShapeFunAndItsGrad(
		cardinal_t number_,
		T (*pfun_)(const Tensor1s<T>&),
		Tensor1s<T>& (*pgrad_)(const Tensor1s<T>&, Tensor1s<T>&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= m_coordApproxData.size(),
		"invalid number of shape function in tFinitElement<T>::DefineCoordShapeFunAndItsGrad");
	Assert(
		m_coordApproxData[number_ - 1].m_pFunction == nullptr,
		"shape function is already defined in tFinitElement<T>::DefineCoordShapeFunAndItsGrad");
	Assert(
		m_coordApproxData[number_ - 1].m_pFunctionGrad == nullptr,
		"shape gradient is already defined in tFinitElement<T>::DefineCoordShapeFunAndItsGrad");
#endif
	CoordApproxData[--number_].pFunction = pfun_;
	CoordApproxData[number_].pFunctionGrad = pgrad_;
}

template <typename T>
inline void tNonIsoparametricFE<T>::DefineDisplShapeFunAndItsGrad(
	cardinal_t number_,
	T (*pfun_)(const Tensor1s<T>&),
	Tensor1s<T>& (*pgrad_)(const Tensor1s<T>&, Tensor1s<T>&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= m_displApproxData.size(),
		"invalid number of shape function in tFinitElement<T>::DefineDisplShapeFunAndItsGrad");
	Assert(
		m_displApproxData[number_ - 1].m_pFunction == nullptr,
		"shape function is already defined in tFinitElement<T>::DefineDisplShapeFunAndItsGrad");
	Assert(
		m_displApproxData[number_ - 1].m_pFunctionGrad == nullptr,
		"shape gradient is already defined in tFinitElement<T>::DefineDisplShapeFunAndItsGrad");
#endif
	DisplApproxData[--number_].pFunction = pfun_;
	DisplApproxData[number_].pFunctionGrad = pgrad_;
}

template <typename T>
inline T tNonIsoparametricFE<T>::CoordShapeFun(const Node3d<T>& node_, const Tensor1s<T>& lc_) const
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
		"There are no shape functions corresponding to the node in "
		"tFinitElement<T>::CoordShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  CoordApproxData.begin(),
				  CoordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunction))(lc_);
}

template <typename T>
inline Tensor1s<T>& tNonIsoparametricFE<T>::CoordShapeGrad(
	const Node3d<T>& node_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_coordApproxData.begin(),
			m_coordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_coordApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement<T>::CoordShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  CoordApproxData.begin(),
				  CoordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunctionGrad))(lc_, result_);
}

template <typename T>
inline T tNonIsoparametricFE<T>::DisplShapeFun(const Node3d<T>& node_, const Tensor1s<T>& lc_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_displApproxData.begin(),
			m_displApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_displApproxData.end(),
		"There are no shape functions corresponding to the node in "
		"tFinitElement<T>::DisplShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  DisplApproxData.begin(),
				  DisplApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunction))(lc_);
}

template <typename T>
inline Tensor1s<T>& tNonIsoparametricFE<T>::DisplShapeGrad(
	const Node3d<T>& node_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			m_displApproxData.begin(),
			m_displApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < m_displApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement<T>::DisplShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  DisplApproxData.begin(),
				  DisplApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunctionGrad))(lc_, result_);
}

template <typename T>
std::shared_ptr<tFinitElement<T>> CreateFiniteElement(size_t dim, size_t nodes);
