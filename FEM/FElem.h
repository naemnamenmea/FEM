// Author: Chekhov Vladimir Valerevich

#ifndef FElemH
#define FElemH
//---------------------------------------------------------------------------
#include <vector>
#include <string>
#include <map>
#include <istream>
//------------------------------------------------------------------------------
#include "NumTypes.h"
#include "Tensors.h"
#include "HybridTensors.h"
#include "NameList.h"
#include "Node.h"
#include "GaussIntegr.hpp"
//------------------------------------------------------------------------------
class tNode;
class tMaterial;
class tAnalysis;
//------------------------------------------------------------------------------
class tFinitElement
{
protected:
	static const tNamesToFuns<tFinitElement> Factory;

	std::vector<const tNode*> pNodes;
	const tMaterial* pMaterial;

private:
	mutable real_t StrainLevel;
	//   mutable const Tensor2s *pDeformData;

public:
	//   void LinkDeformData(const Tensor2s* p_) const {pDeformData = p_;}
	//   void UnLinkDeformData() const {pDeformData = NULL;}

	std::vector<const tNode*>::const_iterator p1stNode() const
	{
		return pNodes.begin();
	}  // for use in T tIsoparametricFE::tArray<T>::operator()

	static tFinitElement* NewFE(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}

protected:
	tFinitElement() : pMaterial(NULL), StrainLevel(0.)
	{
	}
	tFinitElement(cardinal_t numberOfNodes_)
		: pNodes(numberOfNodes_, NULL), pMaterial(NULL), StrainLevel(0.)
	{
	}
	void DefineNextNode(const tNode&);

public:
	virtual ~tFinitElement()
	{
	}
	virtual const std::string& Kind() const = 0;
	virtual tFinitElement& DefineNextNode(const tNode&, bool, bool) = 0;
	cardinal_t HowManyNodes() const
	{
		return pNodes.size();
	}
	virtual cardinal_t HowManyCoordApproxNodes() const = 0;
	virtual cardinal_t HowManyDisplApproxNodes() const = 0;
	virtual cardinal_t SpaceDimension() const = 0;
	const tNode& Node(cardinal_t) const;
	virtual const tNode& CoordApproxNode(cardinal_t) const = 0;
	virtual const tNode& DisplApproxNode(cardinal_t) const = 0;
	bool Has(const tNode& o_) const
	{
		return std::find(pNodes.begin(), pNodes.end(), &o_) > pNodes.end();
	}
	tFinitElement& Link(const tMaterial&);
	const tMaterial& Material() const
	{
		return *pMaterial;
	}
	bool Has(const tMaterial& o_) const
	{
		return pMaterial == &o_;
	}
	virtual const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const;
	//----
	virtual real_t CoordShapeFun(cardinal_t, const Tensor1s&) const = 0;
	virtual Tensor1s& CoordShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t DisplShapeFun(cardinal_t, const Tensor1s&) const = 0;
	virtual Tensor1s& DisplShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const = 0;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t DisplShapeFun(const tNode&, const Tensor1s&) const = 0;
	virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const = 0;
	//----
protected:
	typedef std::map<std::pair<real_t, std::pair<real_t, real_t> >, Tensor2s> tMapRotToLoc;
	mutable tMapRotToLoc RotToLoc_DB;
	Tensor2s& JacobyMatrix_l2g(const Tensor1s&, Tensor2s&)
		const;	// be careful - do not call for 1D or 2D - FEs!!! (because singular)
public:
	virtual const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const = 0;
	virtual real_t Jacobian_l2g(const Tensor1s&) const = 0;
	class fJacobian_l2g
	{
	private:
		mutable real_t Result;
		const tFinitElement& FE;

	public:
		fJacobian_l2g(const tFinitElement& fe_) : FE(fe_)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return Result = FE.Jacobian_l2g(x_);
		}
	};
	//----
	virtual Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const;
	virtual Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const;
	virtual Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const = 0;
	const Tensor2s& TensorOfRotationToLocal(const Tensor1s&) const;
	//----
	template <typename TArg, typename TRet>
	struct fIntegrand
	{
		virtual TRet& operator()(const TArg&) const = 0;
	};
	template <typename TArg>
	struct fOne : public fIntegrand<TArg, real_t>
	{
		mutable real_t Result;
		//         fOne(): Result(1.l) {}
		real_t& operator()(const TArg&) const
		{
			return Result = 1.l;
		}
		//         real_t operator()(const TArg&) const {return 1.l;}
		template <typename TArg2>
		real_t& operator()(const TArg& x_, const TArg2& y_) const
		{
			return Result = 1.;
		}  // need for tParallelepiped1
	};
	template <typename TArg, typename TRet, typename TMult, typename TFac>
	class fMultWithAssgn : public fIntegrand<TArg, TRet>
	{
	private:
		mutable TMult& Multiplicand;
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

	virtual Tensor2a& Integrate_local(const fIntegrand<Tensor1s, Tensor2a>&, Tensor2a&) const = 0;
	virtual Tensor2a& Integrate(const fIntegrand<Tensor1s, Tensor2a>&, Tensor2a&) const = 0;
	virtual real_t Integrate(const fIntegrand<Tensor1s, real_t>&) const = 0;

	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>&,Tensor2s&) const =0;//for
	//   strain
	virtual SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>&,
		SymmetricTensor2s&) const = 0;	// for strain

	//----
	//   tFEsSetOfTensor2& Calculate(fIntegrand<tFEsSetOfTensor2>& fun_, tEvalMethod em_,
	//   tFEsSetOfTensor2& result_) const {return em_==atCentre? result_.Swap(fun_(Tensor1s(0.l))) :
	//   Integrate(fun_,result_) /= Volume();} Tensor1& Displacement (const Tensor1&, Tensor1&)
	//   const; Tensor2s& DisplacementGradient (const Tensor1s&, Tensor2s&) const; tFEsTensor2&
	//   DisplacementGradient (const Tensor1s&, tFEsTensor2&) const;
	virtual real_t Volume() const
	{
		const fOne<Tensor1s> one;
		return Integrate(one);
	}
	//   real_t Mass() const;//tMaterial declared only {return Volume() * Material().Density();}
public:
	Tensor2s& NablaDispl_g(const Tensor1s&, Tensor2s&) const;
	//   Tensor2s& NodalStrain_g(cardinal_t, const Tensor1s&, Tensor2s&) const;
	//   Tensor2& Strain(Tensor2& result_) const {return Strain(Tensor1s(0.l),result_);}//
};
//------------------------------------------------------------------------------
class tNonIsoparametricFE : virtual public tFinitElement
{
private:
	//===Data for providing calls of shape functions and its gradients:
	struct tNodeAnd2FunsPtrs
	{
		const tNode* pNode;
		real_t (*pFunction)(const Tensor1s&);
		Tensor1s& (*pFunctionGrad)(const Tensor1s&, Tensor1s&);
		tNodeAnd2FunsPtrs() : pNode(NULL), pFunction(NULL), pFunctionGrad(NULL)
		{
		}
		struct Node_address_is
		{
			const tNode* AddressToFind;
			Node_address_is(const tNode* pnode_) : AddressToFind(pnode_)
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
		cardinal_t, real_t (*)(const Tensor1s&), Tensor1s& (*)(const Tensor1s&, Tensor1s&));
	void DefineDisplShapeFunAndItsGrad(
		cardinal_t, real_t (*)(const Tensor1s&), Tensor1s& (*)(const Tensor1s&, Tensor1s&));
	//===end of data for providing calls of shape functions and its gradients
	tNonIsoparametricFE(cardinal_t coordNod_, cardinal_t displNod_)
		: CoordApproxData(coordNod_), DisplApproxData(displNod_)
	{
	}

public:
	virtual tFinitElement& DefineNextNode(const tNode&, bool, bool);
	cardinal_t HowManyCoordApproxNodes() const
	{
		return CoordApproxData.size();
	}
	cardinal_t HowManyDisplApproxNodes() const
	{
		return DisplApproxData.size();
	}
	virtual const tNode& CoordApproxNode(cardinal_t) const;
	virtual const tNode& DisplApproxNode(cardinal_t) const;
	//   virtual real_t   CoordShapeFun             (cardinal_t, const Tensor1s&) const;
	//   virtual Tensor1s& CoordShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (cardinal_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;
	virtual real_t DisplShapeFun(const tNode&, const Tensor1s&) const;
	virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;
};
//------------------------------------------------------------------------------
class tIsoparametricFE : virtual public tFinitElement
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
		tTYPE operator()(const tFinitElement*, const tNode*) const;
	};
	typedef real_t (*pShapeFun)(const Tensor1a&);
	typedef Tensor1a& (*pShapeFunGrad)(const Tensor1a&, Tensor1a&);

protected:
	typedef tArray<real_t (*)(const Tensor1s&)> tArrayOfPtrsToShapeFunctions;
	typedef tArray<Tensor1s& (*)(const Tensor1s&, Tensor1s&)> tArrayOfPtrsToShapeFunGradients;
	// friend   real_t(*)(const Tensor1s&)          tArrayOfPtrsToShapeFunctions::operator()(const
	// tFinitElement*,const tNode*) const; friend Tensor1&(*)(const Tensor1s&,Tensor1&)
	// tArrayOfPtrsToShapeFunGradients::operator()(const tFinitElement*,const tNode*) const;

	//**************** ERROR in g++ 4.0.2 **********************
	// friend class tArrayOfPtrsToShapeFunctions;
	// friend class tArrayOfPtrsToShapeFunGradients;
	friend class tArray<real_t (*)(const Tensor1s&)>;
	friend class tArray<Tensor1s& (*)(const Tensor1s&, Tensor1s&)>;

	//===end of data for providing calls of shape functions and its gradients

	tIsoparametricFE()
	{
	}

public:
	virtual cardinal_t HowManyCoordApproxNodes() const
	{
		return HowManyNodes();
	}
	virtual cardinal_t HowManyDisplApproxNodes() const
	{
		return HowManyNodes();
	}
	virtual const tNode& CoordApproxNode(cardinal_t i_) const
	{
		return Node(i_);
	}
	virtual const tNode& DisplApproxNode(cardinal_t i_) const
	{
		return Node(i_);
	}
#ifdef STRONGCHECK
	virtual tFinitElement& DefineNextNode(const tNode&, bool = true, bool = true);
#else	// ifndef STRONGCHECK
	virtual tFinitElement& DefineNextNode(const tNode& node_, bool = true, bool = true)
	{
		tFinitElement::DefineNextNode(node_);
		return *this;
	}
#endif	// def STRONGCHECK
	virtual real_t DisplShapeFun(cardinal_t n_, const Tensor1s& c_) const
	{
		return CoordShapeFun(n_, c_);
	}
	virtual Tensor1s& DisplShapeGrad(cardinal_t n_, const Tensor1s& c_, Tensor1s& r_) const
	{
		return CoordShapeGrad(n_, c_, r_);
	}
	virtual real_t DisplShapeFun(const tNode& n_, const Tensor1s& c_) const
	{
		return CoordShapeFun(n_, c_);
	}
	virtual Tensor1s& DisplShapeGrad(const tNode& n_, const Tensor1s& c_, Tensor1s& r_) const
	{
		return CoordShapeGrad(n_, c_, r_);
	}
};
//------------------------------------------------------------------------------
class t1D_FE : virtual public tFinitElement
{
private:
	static const fIntegrate<1, GAUSS_ORDER_1D> Integrator;
	//--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D_local(TFun fun_, TFunRet& result_) const  // result_ isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}  //--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D(TFun fun_, TFunRet& result_) const	 // result_ isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}  //--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D_local<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}  //--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D<TFunArg>(MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}

public:
	//--------------
	virtual Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}  //--------------
	virtual Tensor2a& Integrate(const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate<Tensor1s>(fun_, result_);
	}  //--------------
	virtual real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const
	{
		real_t result(0.l);
		return Integrate<Tensor1s>(fun_, result);
	}  //--------------
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	virtual SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const  // for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}  //--------------
	//--------------

protected:
	bool StateIs1dStress;

public:
	virtual cardinal_t SpaceDimension() const
	{
		return 1;
	}
	virtual Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const;

private:
	mutable std::map<real_t, Tensor2s> JacobyMatrix_DB;
	mutable std::map<real_t, real_t> Jacobian_DB;

public:
	virtual const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const;
	virtual real_t Jacobian_l2g(const Tensor1s&) const;

	virtual Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const;
	virtual Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const;
	virtual Tensor1s& Tangent(const Tensor1s&, Tensor1s&) const;
	t1D_FE& Set1dStress(bool setToStress_)
	{
		StateIs1dStress = setToStress_;
		return *this;
	}
	bool Has1dStressState() const
	{
		return StateIs1dStress;
	}
	virtual t1D_FE& SetCrosSecArea(real_t) = 0;
	virtual real_t CrosSecArea(real_t = 0.l) const = 0;
	real_t CrosSecArea(const Tensor1s& lc_) const
	{
		return CrosSecArea(lc_(1));
	}
	class fCrossSecArea
	{
	private:
		mutable real_t Result;
		const t1D_FE& FE;

	public:
		fCrossSecArea(const t1D_FE& fe_) : FE(fe_)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return Result = FE.CrosSecArea(x_);
		}
	};
	virtual real_t Length() const
	{
		real_t result(0.l);
		const fOne<Tensor1s> one;
		return Integrate1D<Tensor1s>(one, result);
	}

	virtual const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const;
};
//------------------------------------------------------------------------------
class t2D_FE : virtual public tFinitElement
{
protected:
	bool StateIs2dStress;
	t2D_FE() : StateIs2dStress(false)
	{
	}

public:
	virtual Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const;
	Tensor1s& UnitNormal(const Tensor1s& locoord_, Tensor1s& result_) const
	{
		Tensor1s tmp1, tmp2;
		LocalBasis(locoord_, tmp1, tmp2, result_);
		return result_;
	}
	virtual Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const;
	virtual Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const;
	//----
private:
	mutable std::map<std::pair<real_t, real_t>, Tensor2s> JacobyMatrix_DB;
	mutable std::map<std::pair<real_t, real_t>, real_t> Jacobian_DB;

public:
	virtual const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const;
	virtual real_t Jacobian_l2g(const Tensor1s& lc_) const;
	//----
	virtual real_t Integrate2D_local(const fIntegrand<Tensor1s, real_t>& fun_) const = 0;
	virtual real_t Integrate2D(const fIntegrand<Tensor1s, real_t>& fun_) const = 0;
	//----
	t2D_FE& Set2dStress(bool setToStress_)
	{
		StateIs2dStress = setToStress_;
		return *this;
	}
	bool Has2dStressState() const
	{
		return StateIs2dStress;
	}
	virtual t2D_FE& SetThickness(real_t) = 0;
	virtual real_t Thickness(real_t = 0.l, real_t = 0.l) const = 0;
	class fThickness
	{
	private:
		mutable real_t Result;
		const t2D_FE& FE;

	public:
		fThickness(const t2D_FE& fe_) : FE(fe_)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return Result = FE.Thickness(x_(1), x_(2));
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_, const TArg& y_) const
		{
			return Result = FE.Thickness(x_, y_);
		}
	};
	virtual real_t Area() const
	{
		const fOne<Tensor1s> one;
		return Integrate2D(one);
	}
	//----
	virtual const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const;
};
//------------------------------------------------------------------------------
class t3D_FE : virtual public tFinitElement
{
private:
	mutable std::map<std::pair<std::pair<real_t, real_t>, real_t>, Tensor2s> JacobyMatrix_DB;
	mutable std::map<std::pair<std::pair<real_t, real_t>, real_t>, real_t> Jacobian_DB;
	//     template<typename T> T& Integrate3D (const tFinitElement::fIntegrand<T>&, T&) const;//2nd
	//     parameter must be zero!!! template<typename T> T& Integrate (const
	//     tFinitElement::fIntegrand<T>&, T&) const;//2nd parameter must be zero!!!
public:
	virtual const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const;
	virtual real_t Jacobian_l2g(const Tensor1s&) const;
	virtual Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const;
	/*     virtual real_t&  Integrate3D(const fIntegrand<real_t>& fun_, real_t&  result_)const
							{return Integrate3D<real_t>(fun_,(result_=0.l));}*/
	//     virtual real_t&  Integrate (const fIntegrand<real_t>&  fun_, real_t&  result_) const
	//        {return Integrate3D<real_t>(fun_,(result_=0.l));}
	//     virtual Tensor1a& Integrate (const fIntegrand<Tensor1a>& fun_, Tensor1a& result_)const
	//        {return Integrate3D<Tensor1a>(fun_,result_);}
	//     virtual Tensor2a& Integrate_local (const fIntegrand<Tensor1s,Tensor2a>& fun_, Tensor2a&
	//     result_)const
	//        {return result_;/*Integrate3D<Tensor2a>(fun_,result_);*/}

	virtual cardinal_t SpaceDimension() const
	{
		return 3;
	}
};
//------------------------------------------------------------------------------
class tRectangle : public t2D_FE
{
private:
	static const fIntegrate<2, GAUSS_ORDER_2D> Integrator;
	//--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D_local(
		TFun fun_, TFunRet& result_) const	// 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}
	//--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D(TFun fun_, TFunRet& result_) const	 // 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}
	//--------------
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

public:
	//--------------
	virtual Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}
	virtual Tensor2a& Integrate(const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	virtual real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const
	{
		real_t result(0.l);
		return Integrate<Tensor1s>(fun_, result);
	}
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	virtual SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const  // for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}  //--------------
	virtual real_t Integrate2D_local(const fIntegrand<Tensor1s, real_t>& fun_) const
	{
		real_t result(0.l);
		const fOne<Tensor1s> one;
		return Integrate2D_local<Tensor1s>(one, result);
	}
	virtual real_t Integrate2D(const fIntegrand<Tensor1s, real_t>& fun_) const
	{
		real_t result(0.l);
		const fOne<Tensor1s> one;
		return Integrate2D<Tensor1s>(one, result);
	}
	//--------------
	virtual cardinal_t SpaceDimension() const
	{
		return 2;
	}
};
//------------------------------------------------------------------------------
class tIsoRod2ConstSec : public tIsoparametricFE, public t1D_FE
{
private:
	static const std::string& KindName;
	//----
	static real_t ShapeFun1(const Tensor1s&);
	static real_t ShapeFun2(const Tensor1s&);
	static Tensor1s& ShapeGrad1(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad2(const Tensor1s&, Tensor1s&);
	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
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
	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
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
	//----
	real_t Area;
	explicit tIsoRod2ConstSec() : tFinitElement(2)
	{
	}

public:
	virtual const std::string& Kind() const
	{
		return KindName;
	}
	virtual real_t CoordShapeFun(cardinal_t, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (cardinal_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (const tNode&, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;
	static tFinitElement* NewFE()
	{
		return new tIsoRod2ConstSec();
	}

	//   virtual Tensor1s& Tangent(const Tensor1s&, Tensor1s& result_) const {return
	//   (result_=Node(2).Coord();result_-=Node(1).Coord());} virtual real_t Length() const {return
	//   (Tensor1s d=Node(2).Coord();d-=Node(1).Coord()).Length();}
	virtual t1D_FE& SetCrosSecArea(real_t newval_)
	{
		Area = newval_;
		return *this;
	}
	virtual real_t CrosSecArea(real_t = 0.l) const
	{
		return Area;
	}
	//   virtual real_t Volume() const {return Length()*CrosSecArea();}
};
//------------------------------------------------------------------------------
class tIsoQuad4ConsThick : public tIsoparametricFE, public tRectangle
{
private:
	static const std::string& KindName;
	//----
	static real_t ShapeFun1(const Tensor1s&);
	static real_t ShapeFun2(const Tensor1s&);
	static real_t ShapeFun3(const Tensor1s&);
	static real_t ShapeFun4(const Tensor1s&);
	static Tensor1s& ShapeGrad1(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad2(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad3(const Tensor1s&, Tensor1s&);
	static Tensor1s& ShapeGrad4(const Tensor1s&, Tensor1s&);
	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;
	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;
	//----
	real_t ThicknessValue;
	explicit tIsoQuad4ConsThick() : tFinitElement(4)
	{
	}

public:
	virtual real_t CoordShapeFun(cardinal_t, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (cardinal_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (const tNode&, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;

	static tFinitElement* NewFE()
	{
		return new tIsoQuad4ConsThick();
	}

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	//   virtual real_t Volume() const {return Area() * ThicknessValue;}
	virtual t2D_FE& SetThickness(real_t newval_)
	{
		ThicknessValue = newval_;
		return *this;
	}
	virtual real_t Thickness(real_t = 0.l, real_t = 0.l) const
	{
		return ThicknessValue;
	}
};
//------------------------------------------------------------------------------
class tIsoParallelepiped8 : public tIsoparametricFE, public t3D_FE
{
private:
	static const std::string& KindName;
	static const fIntegrate<3, GAUSS_ORDER_3D> Integrator;
	//--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D_local(
		TFun fun_, TFunRet& result_) const	// 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}
	//--------------
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D(TFun fun_, TFunRet& result_) const	 // 2nd parameter isn't assigned by 0!!!
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			MultBy<TFunArg, TFunRet>(fun_, fJacobian_l2g(*this)), result_);
	}
	//--------------
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

public:
	//--------------
	virtual Tensor2a& Integrate_local(
		const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate_local<Tensor1s>(fun_, result_);
	}
	virtual Tensor2a& Integrate(const fIntegrand<Tensor1s, Tensor2a>& fun_, Tensor2a& result_) const
	{
		return Integrate<Tensor1s>(fun_, result_);
	}
	virtual real_t Integrate(const fIntegrand<Tensor1s, real_t>& fun_) const
	{
		real_t result(0.l);
		return Integrate<Tensor1s>(fun_, result);
	}
	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>& fun_, Tensor2s& result_)
	//   const//for strain
	virtual SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>& fun_,
		SymmetricTensor2s& result_) const  // for strain
	{
		return Integrate<Tensor1s>(fun_, result_);
	}  //--------------
	//--------------

private:
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

	class tPtrShapeFunArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;

	class tPtrShapeFunGradArray : public tIsoparametricFE::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	explicit tIsoParallelepiped8() : tFinitElement(8)
	{
	}

public:
	virtual real_t CoordShapeFun(cardinal_t, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(cardinal_t, const Tensor1s&, Tensor1s&) const;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const;

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	static tFinitElement* NewFE()
	{
		return new tIsoParallelepiped8();
	}
};
//------------------------------------------------------------------------------
//--------------------------tFinitElement---------------------------------------
//------------------------------------------------------------------------------
inline const tNode& tFinitElement::Node(cardinal_t nodeNo_) const
{
#ifdef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= pNodes.size(), "Invalid node # in FE");
#endif
	return *(pNodes[--nodeNo_]);
}
//------------------------------------------------------------------------------
inline tFinitElement& tFinitElement::Link(const tMaterial& matlToLink_)
{
#ifdef STRONGCHECK
	Assert(pMaterial != &matlToLink_, "Invalid material in FE");
#endif
	// if (pMaterial == &matlToLink_) return *this;
	pMaterial = &matlToLink_;
	return *this;
}
//------------------------------------------------------------------------------
//--------------------------tNonIsoparametricFE---------------------------------
//------------------------------------------------------------------------------
inline void tNonIsoparametricFE::DefineCoordShapeFunAndItsGrad(
	cardinal_t number_,
	real_t (*pfun_)(const Tensor1s&),
	Tensor1s& (*pgrad_)(const Tensor1s&, Tensor1s&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= CoordApproxData.size(),
		"invalid number of shape function in tFinitElement::DefineCoordShapeFunAndItsGrad");
	Assert(
		CoordApproxData[number_ - 1].pFunction == NULL,
		"shape function is already defined in tFinitElement::DefineCoordShapeFunAndItsGrad");
	Assert(
		CoordApproxData[number_ - 1].pFunctionGrad == NULL,
		"shape gradient is already defined in tFinitElement::DefineCoordShapeFunAndItsGrad");
#endif
	CoordApproxData[--number_].pFunction = pfun_;
	CoordApproxData[number_].pFunctionGrad = pgrad_;
}
//------------------------------------------------------------------------------
inline void tNonIsoparametricFE::DefineDisplShapeFunAndItsGrad(
	cardinal_t number_,
	real_t (*pfun_)(const Tensor1s&),
	Tensor1s& (*pgrad_)(const Tensor1s&, Tensor1s&))
{
#ifdef STRONGCHECK
	Assert(
		number_ > 0 && number_ <= DisplApproxData.size(),
		"invalid number of shape function in tFinitElement::DefineDisplShapeFunAndItsGrad");
	Assert(
		DisplApproxData[number_ - 1].pFunction == NULL,
		"shape function is already defined in tFinitElement::DefineDisplShapeFunAndItsGrad");
	Assert(
		DisplApproxData[number_ - 1].pFunctionGrad == NULL,
		"shape gradient is already defined in tFinitElement::DefineDisplShapeFunAndItsGrad");
#endif
	DisplApproxData[--number_].pFunction = pfun_;
	DisplApproxData[number_].pFunctionGrad = pgrad_;
}
//------------------------------------------------------------------------------
inline real_t tNonIsoparametricFE::CoordShapeFun(const tNode& node_, const Tensor1s& lc_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	//       p =
	//       std::find_if(CoordApproxData.begin(),CoordApproxData.end(),tNodeAnd2FunsPtrs::Node_address_is(&node_));
	Assert(
		std::find_if(
			CoordApproxData.begin(),
			CoordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < CoordApproxData.end(),
		"There are no shape functions corresponding to the node in tFinitElement::CoordShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  CoordApproxData.begin(),
				  CoordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunction))(lc_);
}
//------------------------------------------------------------------------------
inline Tensor1s& tNonIsoparametricFE::CoordShapeGrad(
	const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			CoordApproxData.begin(),
			CoordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < CoordApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement::CoordShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  CoordApproxData.begin(),
				  CoordApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunctionGrad))(lc_, result_);
}
//------------------------------------------------------------------------------
inline real_t tNonIsoparametricFE::DisplShapeFun(const tNode& node_, const Tensor1s& lc_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			DisplApproxData.begin(),
			DisplApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < DisplApproxData.end(),
		"There are no shape functions corresponding to the node in tFinitElement::DisplShapeFun");
#endif
	// return (*(p->pFunction))(lc_);
	return (*(std::find_if(
				  DisplApproxData.begin(),
				  DisplApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunction))(lc_);
}
//------------------------------------------------------------------------------
inline Tensor1s& tNonIsoparametricFE::DisplShapeGrad(
	const tNode& node_, const Tensor1s& lc_, Tensor1s& result_) const
{
#ifdef STRONGCHECK
	// std::vector<tNodeAnd2FunsPtrs>::const_iterator
	Assert(
		std::find_if(
			DisplApproxData.begin(),
			DisplApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(&node_)) < DisplApproxData.end(),
		"There are no shape function grads corresponding to the node in "
		"tFinitElement::DisplShapeGrad");
#endif
	// return (*(p->pFunctionGrad))(lc_,result_);
	return (*(std::find_if(
				  DisplApproxData.begin(),
				  DisplApproxData.end(),
				  tNodeAnd2FunsPtrs::Node_address_is(&node_))
				  ->pFunctionGrad))(lc_, result_);
}
//------------------------------------------------------------------------------
#endif	// FElemH
