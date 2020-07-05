#pragma once

#include "Tensors.hpp"
#include "NameList.hpp"
#include "Material.hpp"
#include "Node.hpp"

#include <map>

class tFinitElement
{
public:
	template <typename TArg, typename TRet>
	struct fIntegrand
	{
		virtual TRet& operator()(const TArg&) const = 0;
	};

	template <typename TArg, typename TRet, typename TMult, typename TFac>
	class fMultWithAssgn : public fIntegrand<TArg, TRet>
	{
	public:
		fMultWithAssgn(TMult& m_, const TFac& f_) : m_multiplicand(m_), m_factor(f_)
		{
		}
		TRet& operator()(const TArg& x_) const
		{
			return m_multiplicand(x_) *= m_factor(x_);
		}

	private:
		TMult& m_multiplicand;
		const TFac& m_factor;
	};

	static tFinitElement* NewFE(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}
	template <typename TArg, typename TRet, typename TMult, typename TFac>
	static fMultWithAssgn<TArg, TRet, TMult, TFac> MultBy(TMult& m_, const TFac& f_)
	{
		return fMultWithAssgn<TArg, TRet, TMult, TFac>(m_, f_);
	}

	//   void LinkDeformData(const Tensor2s* p_) const {pDeformData = p_;}
	//   void UnLinkDeformData() const {pDeformData = nullptr;}

	std::vector<const tNode*>::const_iterator p1stNode() const
	{
		return m_pNodes.begin();
	}  // for use in T tIsoparametricFE::tArray<T>::operator()
	virtual ~tFinitElement()
	{
	}
	virtual const std::string& Kind() const = 0;
	virtual tFinitElement& DefineNextNode(const tNode&, bool, bool) = 0;
	size_t HowManyNodes() const
	{
		return m_pNodes.size();
	}
	virtual size_t HowManyCoordApproxNodes() const = 0;
	virtual size_t HowManyDisplApproxNodes() const = 0;
	virtual size_t SpaceDimension() const = 0;
	const tNode& Node(size_t) const;
	virtual const tNode& CoordApproxNode(size_t) const = 0;
	virtual const tNode& DisplApproxNode(size_t) const = 0;
	bool Has(const tNode& o_) const
	{
		return std::find(m_pNodes.begin(), m_pNodes.end(), &o_) > m_pNodes.end();
	}
	tFinitElement& Link(const tMaterial&);
	const tMaterial& Material() const
	{
		return *m_pMaterial;
	}
	bool Has(const tMaterial& o_) const
	{
		return m_pMaterial == &o_;
	}
	virtual const SymmetricTensor4s& MaterialElasTensor(const Tensor1s&) const;

	virtual real_t CoordShapeFun(size_t, const Tensor1s&) const = 0;
	virtual Tensor1s& CoordShapeGrad(size_t, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t DisplShapeFun(size_t, const Tensor1s&) const = 0;
	virtual Tensor1s& DisplShapeGrad(size_t, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t CoordShapeFun(const tNode&, const Tensor1s&) const = 0;
	virtual Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const = 0;
	virtual real_t DisplShapeFun(const tNode&, const Tensor1s&) const = 0;
	virtual Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const = 0;

	virtual const Tensor2s& JacobyMatrix_l2g(const Tensor1s&) const = 0;
	virtual real_t Jacobian_l2g(const Tensor1s&) const = 0;

	class fJacobian_l2g
	{
	public:
		fJacobian_l2g(const tFinitElement& fe_) : m_FE(fe_), m_result(0)
		{
		}
		template <typename TArg>
		real_t& operator()(const TArg& x_) const
		{
			return m_result = m_FE.Jacobian_l2g(x_);
		}

	private:
		mutable real_t m_result;
		const tFinitElement& m_FE;
	};

	virtual Tensor1s& GlobalCoord(const Tensor1s&, Tensor1s&) const;
	virtual Tensor1s& LocalBasis(const Tensor1s&, Tensor1s&, Tensor1s&, Tensor1s&) const;
	virtual Tensor2s& TensorOfRotationToLocal(const Tensor1s&, Tensor2s&) const = 0;
	const Tensor2s& TensorOfRotationToLocal(const Tensor1s&) const;

	template <typename TArg>
	class fOne : public fIntegrand<TArg, real_t>
	{
	public:
		fOne() : m_result(0)
		{
		}

		//         fOne(): Result(1.) {}
		real_t& operator()(const TArg&) const
		{
			return m_result = 1.;
		}
		//         real_t operator()(const TArg&) const {return 1.;}
		template <typename TArg2>
		real_t& operator()(const TArg& x_, const TArg2& y_) const
		{
			return m_result = 1.;
		}  // need for tParallelepiped1

	private:
		mutable real_t m_result;
	};

	virtual Tensor2a& Integrate_local(const fIntegrand<Tensor1s, Tensor2a>&, Tensor2a&) const = 0;
	virtual Tensor2a& Integrate(const fIntegrand<Tensor1s, Tensor2a>&, Tensor2a&) const = 0;
	virtual real_t Integrate(const fIntegrand<Tensor1s, real_t>&) const = 0;

	//   virtual Tensor2s& Integrate (const fIntegrand<Tensor1s,Tensor2s>&,Tensor2s&) const =0;//for
	//   strain
	virtual SymmetricTensor2s& Integrate(
		const fIntegrand<Tensor1s, SymmetricTensor2s>&,
		SymmetricTensor2s&) const = 0;	// for strain

	//   tFEsSetOfTensor2& Calculate(fIntegrand<tFEsSetOfTensor2>& fun_, tEvalMethod em_,
	//   tFEsSetOfTensor2& result_) const {return em_==atCentre? result_.Swap(fun_(Tensor1s(0.))) :
	//   Integrate(fun_,result_) /= Volume();} Tensor1& Displacement (const Tensor1&, Tensor1&)
	//   const; Tensor2s& DisplacementGradient (const Tensor1s&, Tensor2s&) const; tFEsTensor2&
	//   DisplacementGradient (const Tensor1s&, tFEsTensor2&) const;
	virtual real_t Volume() const
	{
		const fOne<Tensor1s> one;
		return Integrate(one);
	}
	//   real_t Mass() const;//tMaterial declared only {return Volume() * Material().Density();}
	Tensor2s& NablaDispl_g(const Tensor1s&, Tensor2s&) const;
	//   Tensor2s& NodalStrain_g(size_t, const Tensor1s&, Tensor2s&) const;
	//   Tensor2& Strain(Tensor2& result_) const {return Strain(Tensor1s(0.),result_);}//

protected:
	typedef std::map<std::pair<real_t, std::pair<real_t, real_t> >, Tensor2s> tMapRotToLoc;

	static const tNamesToFuns<tFinitElement> Factory;

	tFinitElement() : m_pMaterial(nullptr), m_strainLevel(0.)
	{
	}
	tFinitElement(size_t numberOfNodes_)
		: m_pNodes(numberOfNodes_, nullptr), m_pMaterial(nullptr), m_strainLevel(0.)
	{
	}
	void DefineNextNode(const tNode&);
	Tensor2s& JacobyMatrix_l2g(const Tensor1s&, Tensor2s&)
		const;	// be careful - do not call for 1D or 2D - FEs!!! (because singular)

	mutable tMapRotToLoc m_rotToLocDB;
	std::vector<const tNode*> m_pNodes;
	const tMaterial* m_pMaterial;

private:
	mutable real_t m_strainLevel;
	//   mutable const Tensor2s *pDeformData;
};

class tNonIsoparametricFE : virtual public tFinitElement
{
public:
	tFinitElement& DefineNextNode(const tNode&, bool, bool) override;
	size_t HowManyCoordApproxNodes() const
	{
		return m_coordApproxData.size();
	}
	size_t HowManyDisplApproxNodes() const
	{
		return m_displApproxData.size();
	}
	const tNode& CoordApproxNode(size_t) const override;
	const tNode& DisplApproxNode(size_t) const override;
	//   virtual real_t   CoordShapeFun             (size_t, const Tensor1s&) const;
	//   virtual Tensor1s& CoordShapeGrad(size_t, const Tensor1s&, Tensor1s&) const;
	//   virtual real_t   DisplShapeFun             (size_t, const Tensor1s&) const;
	//   virtual Tensor1s& DisplShapeGrad(size_t, const Tensor1s&, Tensor1s&) const;
	real_t CoordShapeFun(const tNode&, const Tensor1s&) const override;
	Tensor1s& CoordShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const override;
	real_t DisplShapeFun(const tNode&, const Tensor1s&) const override;
	Tensor1s& DisplShapeGrad(const tNode&, const Tensor1s&, Tensor1s&) const override;

protected:
	tNonIsoparametricFE(size_t coordNod_, size_t displNod_)
		: m_coordApproxData(coordNod_), m_displApproxData(displNod_)
	{
	}
	void DefineCoordShapeFunAndItsGrad(
		size_t, real_t (*)(const Tensor1s&), Tensor1s& (*)(const Tensor1s&, Tensor1s&));
	void DefineDisplShapeFunAndItsGrad(
		size_t, real_t (*)(const Tensor1s&), Tensor1s& (*)(const Tensor1s&, Tensor1s&));
	//===end of data for providing calls of shape functions and its gradients

private:
	//===Data for providing calls of shape functions and its gradients:
	class tNodeAnd2FunsPtrs
	{
	public:
		tNodeAnd2FunsPtrs() : m_pNode(nullptr), m_pFunction(nullptr), m_pFunctionGrad(nullptr)
		{
		}

		class Node_address_is
		{
		public:
			Node_address_is(const tNode* pnode_) : m_addressToFind(pnode_)
			{
			}
			bool operator()(const tNodeAnd2FunsPtrs& valueToTest_)
			{
				return valueToTest_.m_pNode == m_addressToFind;
			}

		private:
			const tNode* m_addressToFind;
		};

		const tNode* m_pNode;
		real_t (*m_pFunction)(const Tensor1s&);
		Tensor1s& (*m_pFunctionGrad)(const Tensor1s&, Tensor1s&);
	};

	std::vector<tNodeAnd2FunsPtrs> m_coordApproxData;
	std::vector<tNodeAnd2FunsPtrs> m_displApproxData;
};

class tIsoparametricFE : virtual public tFinitElement
{
	//===Data for providing calls of shape functions and its gradients:
public:
	size_t HowManyCoordApproxNodes() const override
	{
		return HowManyNodes();
	}
	size_t HowManyDisplApproxNodes() const override
	{
		return HowManyNodes();
	}
	const tNode& CoordApproxNode(size_t i_) const override
	{
		return Node(i_);
	}
	const tNode& DisplApproxNode(size_t i_) const override
	{
		return Node(i_);
	}
#ifdef STRONGCHECK
	tFinitElement& DefineNextNode(const tNode&, bool = true, bool = true) override;
#else	// ifndef STRONGCHECK
	tFinitElement& DefineNextNode(const tNode& node_, bool = true, bool = true) override
	{
		tFinitElement::DefineNextNode(node_);
		return *this;
	}
#endif	// def STRONGCHECK
	real_t DisplShapeFun(size_t n_, const Tensor1s& c_) const override
	{
		return CoordShapeFun(n_, c_);
	}
	Tensor1s& DisplShapeGrad(size_t n_, const Tensor1s& c_, Tensor1s& r_) const override
	{
		return CoordShapeGrad(n_, c_, r_);
	}
	real_t DisplShapeFun(const tNode& n_, const Tensor1s& c_) const override
	{
		return CoordShapeFun(n_, c_);
	}
	Tensor1s& DisplShapeGrad(const tNode& n_, const Tensor1s& c_, Tensor1s& r_) const override
	{
		return CoordShapeGrad(n_, c_, r_);
	}

protected:
	template <typename tTYPE>
	class tArray : public std::vector<tTYPE>
	{
	public:
		tArray(size_t dim_) : std::vector<tTYPE>(dim_)
		{
		}
		tTYPE operator()(size_t) const;
		tTYPE operator()(const tFinitElement*, const tNode*) const;
	};

	typedef real_t (*pShapeFun)(const Tensor1a&);
	typedef Tensor1a& (*pShapeFunGrad)(const Tensor1a&, Tensor1a&);
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
};
