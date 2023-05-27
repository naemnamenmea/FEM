#pragma once

#include "math_constants.hpp"
#include "Tensors.h"
#include "Material.hpp"
#include "GaussIntegr.hpp"
#include "stl_containers_read.hpp"
#include "NumTypes.h"
#include "Node.hpp"
#include <istream>
#include <list>
#include <vector>
#include <map>
#include <unordered_set>

/*
 * 1 dim - 2 nodes
 * 2 dim - 4 nodes
 * 3 dim - 8 nodes
 */

template <typename T>
class Node3d;
class tMaterial;
class tAnalysis;

template <typename T>
class FiniteElementBase
{
public:
	virtual const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const = 0;
	virtual T Jacobian_l2g(const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& LocalBasis(
		const Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&) const;

	class fJacobian_l2g
	{
	private:
		mutable T result;
		const FiniteElementBase& FE;

	public:
		fJacobian_l2g(const FiniteElementBase& fe_) : FE(fe_), result(0)
		{
		}

		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return result = FE.Jacobian_l2g(x_);
		}
	};

	FiniteElementBase();
	FiniteElementBase(size_t nodes);

	const auto& GetType() const
	{
		return m_type;
	}
	const auto& GetNodeNumbers() const
	{
		return m_nodeNumbers;
	}
	const auto& GetMaterial() const
	{
		return *m_material;
	}

	void SetType(std::string type)
	{
		this->m_type = std::move(type);
	}
	void SetNodeNumbers(const std::vector<size_t>& nodeNumbers)
	{
		this->m_nodeNumbers = nodeNumbers;
	}
	void SetMaterial(std::unordered_set<Material>::const_iterator material)
	{
		this->m_material = material;
	}

	virtual void ReadProperties(std::istream& is);
	virtual unsigned int GetDim() const = 0;
	virtual double GetParameter() const
	{
		return mathdef::MAX_DOUBLE;
	}

	cardinal_t HowManyNodes() const
	{
		return pNodes.size();
	}

	const Node3d<T>& Node(cardinal_t nodeNo_) const
	{
		return *(pNodes[--nodeNo_]);
	}
	virtual const Node3d<T>& CoordApproxNode(cardinal_t) const = 0;
	virtual cardinal_t HowManyCoordApproxNodes() const = 0;

	virtual T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	//virtual T DisplShapeFun(cardinal_t, const Tensor1s<T>&) const = 0;
	//virtual Tensor1s<T>& DisplShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	virtual T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const = 0;
	virtual Tensor1s<T>& CoordShapeGrad(
		const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const = 0;
	//virtual T DisplShapeFun(const Node3d<T>&, const Tensor1s<T>&) const = 0;
	//virtual Tensor1s<T>& DisplShapeGrad(
	//	const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const = 0;

	template <typename TArg, typename TRet>
	struct fIntegrand
	{
		virtual TRet& operator()(const TArg&) const = 0;
	};

	template <typename TArg>
	struct fOne : public fIntegrand<TArg, T>
	{
		mutable T result;

		fOne() : result(0)
		{
		}

		T& operator()(const TArg&) const
		{
			return result = 1.;
		}

		template <typename TArg2>
		T& operator()(const TArg& x_, const TArg2& y_) const
		{
			return result = 1.;
		}
	};

	template <typename TArg, typename TRet, typename TMult, typename TFac>
	class fMultWithAssgn : public fIntegrand<TArg, TRet>
	{
	public:
		fMultWithAssgn(TMult& multiplicand, const TFac& factor)
			: m_multiplicand(multiplicand), m_factor(factor)
		{
		}
		TRet& operator()(const TArg& arg) const
		{
			return m_multiplicand(arg) *= m_factor(arg);
		}

	private:
		TMult& m_multiplicand;
		const TFac& m_factor;
	};

	template <typename TArg, typename TRet, typename TMult, typename TFac>
	static fMultWithAssgn<TArg, TRet, TMult, TFac> MultBy(TMult& multiplicand, const TFac& factor)
	{
		return fMultWithAssgn<TArg, TRet, TMult, TFac>(multiplicand, factor);
	}

	virtual Tensor2a<T>& Integrate_local(
		const fIntegrand<Tensor1s<T>, Tensor2a<T>>&, Tensor2a<T>&) const = 0;
	virtual Tensor2a<T>& Integrate(
		const fIntegrand<Tensor1s<T>, Tensor2a<T>>&, Tensor2a<T>&) const = 0;
	virtual T Integrate(const fIntegrand<Tensor1s<T>, T>&) const = 0;

	virtual T Volume() const
	{
		const fOne<Tensor1s<T>> one;
		return Integrate(one);
	}

protected:
	Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&, Tensor2s<T>&)
		const;	// be careful - do not call for 1D or 2D - FEs!!! (because singular)

	std::string m_type;
	std::vector<const Node3d<T>*> pNodes;
	std::vector<size_t> m_nodeNumbers;
	std::unordered_set<Material>::const_iterator m_material;
};

template <typename T>
class IsoparametricFE : virtual public FiniteElementBase<T>
{
private:
	template <typename tTYPE>
	class tArray : public std::vector<tTYPE>
	{
	public:
		tArray(cardinal_t dim_) : std::vector<tTYPE>(dim_)
		{
		}
		tTYPE operator()(cardinal_t i_) const
		{
			return this->operator[](--i_);
		}
		tTYPE operator()(const FiniteElementBase<T>*, const Node3d<T>*) const;
	};

protected:
	typedef tArray<T(*)(const Tensor1s<T>&)> tArrayOfPtrsToShapeFunctions;
	typedef tArray<Tensor1s<T> & (*)(const Tensor1s<T>&, Tensor1s<T>&)> tArrayOfPtrsToShapeFunGradients;

public:
	const Node3d<T>& CoordApproxNode(cardinal_t i_) const override
	{
		return FiniteElementBase<T>::Node(i_);
	}

	cardinal_t HowManyCoordApproxNodes() const override
	{
		return FiniteElementBase<T>::HowManyNodes();
	}
};

template <typename T>
class NonIsoparametricFE : virtual public FiniteElementBase<T>
{
private:
	//===Data for providing calls of shape functions and its gradients:
	struct tNodeAnd2FunsPtrs
	{
		tNodeAnd2FunsPtrs() : pNode(NULL), pFunction(NULL), pFunctionGrad(NULL)
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

		T (*pFunction)(const Tensor1s<T>&);
		Tensor1s<T>& (*pFunctionGrad)(const Tensor1s<T>&, Tensor1s<T>&);

		const Node3d<T>* pNode;
	};

public:
	cardinal_t HowManyCoordApproxNodes() const
	{
		return CoordApproxData.size();
	}
	const Node3d<T>& CoordApproxNode(cardinal_t) const override;

protected:
	NonIsoparametricFE(cardinal_t coordNod_, cardinal_t displNod_)
		: CoordApproxData(coordNod_), DisplApproxData(displNod_)
	{
	}

private:
	std::vector<tNodeAnd2FunsPtrs> CoordApproxData;
	std::vector<tNodeAnd2FunsPtrs> DisplApproxData;
};

template <typename T>
class FiniteElement1d : virtual public FiniteElementBase<T>
{
public:
	virtual Tensor1s<T>& Tangent(const Tensor1s<T>&, Tensor1s<T>&) const;

	const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const override;
	T Jacobian_l2g(const Tensor1s<T>&) const override;

	Tensor1s<T>& LocalBasis(
		const Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&, Tensor1s<T>&) const override;
	virtual T CrosSecArea(T = 0.l) const = 0;
	T CrosSecArea(const Tensor1s<T>& lc_) const
	{
		return CrosSecArea(lc_(1));
	}

	class fCrossSecArea
	{
	private:
		mutable T result;
		const FiniteElement1d& FE;

	public:
		fCrossSecArea(const FiniteElement1d& fe_) : FE(fe_), result(0)
		{
		}
		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return result = FE.CrosSecArea(x_);
		}
	};

	Tensor2a<T>& Integrate_local(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}

	Tensor2a<T>& Integrate(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}

	T Integrate(const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, T>& fun_) const
	{
		T result(0.l);
		return Integrate<Tensor1s<T>>(fun_, result);
	}

private:
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, FiniteElementBase<T>::fJacobian_l2g(*this)),
			result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate1D_local(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate1D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, fCrossSecArea(*this)), result_);
	}

	static const GaussIntegr::fIntegrate<1, GAUSS_ORDER_1D> Integrator;

	mutable std::map<T, T> Jacobian_DB;
	mutable std::map<T, Tensor2s<T>> JacobyMatrix_DB;
};

template <typename T>
class FiniteElement2d : virtual public FiniteElementBase<T>
{
private:
	mutable std::map<std::pair<T, T>, T> Jacobian_DB;
	mutable std::map<std::pair<T, T>, Tensor2s<T>> JacobyMatrix_DB;

public:
	Tensor1s<T>& UnitNormal(const Tensor1s<T>& locoord_, Tensor1s<T>& result_) const
	{
		Tensor1s<T> tmp1, tmp2;
		LocalBasis(locoord_, tmp1, tmp2, result_);
		return result_;
	}

	const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const override;
	T Jacobian_l2g(const Tensor1s<T>& lc_) const override;

	virtual T Thickness(T = 0.l, T = 0.l) const = 0;

	class fThickness
	{
	private:
		mutable T result;
		const FiniteElement2d<T>& FE;

	public:
		fThickness(const FiniteElement2d<T>& fe_) : FE(fe_), result(0)
		{
		}

		template <typename TArg>
		T& operator()(const TArg& x_) const
		{
			return result = FE.Thickness(x_(1), x_(2));
		}

		template <typename TArg>
		T& operator()(const TArg& x_, const TArg& y_) const
		{
			return result = FE.Thickness(x_, y_);
		}
	};
};

template <typename T>
class Rectangle : public FiniteElement2d<T>
{
public:
	Tensor2a<T>& Integrate_local(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}

	Tensor2a<T>& Integrate(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}

	T Integrate(const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, T>& fun_) const override
	{
		T result(0.l);
		return Integrate<Tensor1s<T>>(fun_, result);
	}

private:
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, FiniteElementBase<T>::fJacobian_l2g(*this)),
			result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate2D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, FiniteElement2d<T>::fThickness(*this)), result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate2D_local(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate2D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, FiniteElement2d<T>::fThickness(*this)), result_);
	}

	static const GaussIntegr::fIntegrate<2, GAUSS_ORDER_2D> Integrator;
};

template <typename T>
class FiniteElement3d : virtual public FiniteElementBase<T>
{
public:
	const Tensor2s<T>& JacobyMatrix_l2g(const Tensor1s<T>&) const override;
	T Jacobian_l2g(const Tensor1s<T>& lc_) const override;

private:
	mutable std::map<std::pair<std::pair<T, T>, T>, T> Jacobian_DB;
	mutable std::map<std::pair<std::pair<T, T>, T>, Tensor2s<T>> JacobyMatrix_DB;
};

#pragma warning(push)
#pragma warning(disable : 4250)
template <typename T>
class IsoRod2ConstSec : public IsoparametricFE<T>, public FiniteElement1d<T>
{
public:
	explicit IsoRod2ConstSec() : FiniteElementBase<T>(2), m_crossSectionalArea(0)
	{
	}
	IsoRod2ConstSec(size_t nodes);

	void ReadProperties(std::istream& is) override;

	double GetParameter() const override
	{
		return m_crossSectionalArea;
	}

	unsigned int GetDim() const override
	{
		return 1;
	}

	Tensor1s<T>& CoordShapeGrad(
		cardinal_t no_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const override
	{
		return (*(pShapeFunGrads(no_)))(lc_, result_);
	}

	T CrosSecArea(T = 0.l) const override
	{
		return m_crossSectionalArea;
	}

private:
	static const std::string& KindName;

	static T ShapeFun1(const Tensor1s<T>&);
	static T ShapeFun2(const Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad1(const Tensor1s<T>&, Tensor1s<T>&);
	static Tensor1s<T>& ShapeGrad2(const Tensor1s<T>&, Tensor1s<T>&);
	class tPtrShapeFunArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray() : IsoparametricFE<T>::tArrayOfPtrsToShapeFunctions(2)
		{
			operator[](0) = ShapeFun1;
			operator[](1) = ShapeFun2;
		}
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;
	class tPtrShapeFunGradArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray() : IsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients(2)
		{
			operator[](0) = ShapeGrad1;
			operator[](1) = ShapeGrad2;
		}
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	double m_crossSectionalArea;
};

template <typename T>
class IsoQuad4ConsThick : public IsoparametricFE<T>, public Rectangle<T>
{
public:
	explicit IsoQuad4ConsThick() : FiniteElementBase<T>(4), m_thickness(0)
	{
	}
	IsoQuad4ConsThick(size_t nodes);

	void ReadProperties(std::istream& is) override;

	double GetParameter() const override
	{
		return m_thickness;
	}

	unsigned int GetDim() const override
	{
		return 2;
	}

	T Thickness(T = 0.l, T = 0.l) const override
	{
		return m_thickness;
	}

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
	class tPtrShapeFunArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;
	class tPtrShapeFunGradArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

	double m_thickness;
};

template <typename T>
class IsoParallelepiped8 : public IsoparametricFE<T>, public FiniteElement3d<T>
{
public:
	explicit IsoParallelepiped8() : FiniteElementBase<T>(8)
	{
	}
	IsoParallelepiped8(size_t nodes);

	unsigned int GetDim() const override
	{
		return 3;
	}

	Tensor2a<T>& Integrate_local(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate_local<Tensor1s<T>>(fun_, result_);
	}

	Tensor2a<T>& Integrate(
		const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, Tensor2a<T>>& fun_,
		Tensor2a<T>& result_) const override
	{
		return Integrate<Tensor1s<T>>(fun_, result_);
	}

	T Integrate(const FiniteElementBase<T>::fIntegrand<Tensor1s<T>, T>& fun_) const override
	{
		T result(0.l);
		return Integrate<Tensor1s<T>>(fun_, result);
	}

private:
	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D_local(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(fun_, result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate3D(TFun fun_, TFunRet& result_) const
	{
		return Integrator.ByPlusAssgn<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(fun_, FiniteElementBase<T>::fJacobian_l2g(*this)),
			result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate_local(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D_local<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(
				fun_, FiniteElementBase<T>::template fOne<Tensor1s<T>>()),
			result_);
	}

	template <typename TFunArg, typename TFunRet, typename TFun>
	TFunRet& Integrate(TFun& fun_, TFunRet& result_) const
	{
		return Integrate3D<TFunArg>(
			this->MultBy<TFunArg, TFunRet>(
				fun_, FiniteElementBase<T>::template fOne<Tensor1s<T>>()),
			result_);
	}

	static const GaussIntegr::fIntegrate<3, GAUSS_ORDER_3D> Integrator;

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

	class tPtrShapeFunArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunctions
	{
	public:
		tPtrShapeFunArray();
	};
	static const tPtrShapeFunArray pShapeFunctions;
	friend class tPtrShapeFunArray;

	class tPtrShapeFunGradArray : public IsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients
	{
	public:
		tPtrShapeFunGradArray();
	};
	static const tPtrShapeFunGradArray pShapeFunGrads;
	friend class tPtrShapeFunGradArray;

public:
	T CoordShapeFun(cardinal_t, const Tensor1s<T>&) const override;
	Tensor1s<T>& CoordShapeGrad(cardinal_t, const Tensor1s<T>&, Tensor1s<T>&) const override;
	T CoordShapeFun(const Node3d<T>&, const Tensor1s<T>&) const override;
	Tensor1s<T>& CoordShapeGrad(const Node3d<T>&, const Tensor1s<T>&, Tensor1s<T>&) const override;

	// virtual const std::string& Kind() const { return KindName; }
	static FiniteElementBase<T>* NewFE()
	{
		return new IsoParallelepiped8<T>();
	}
};
#pragma warning(pop)

template <typename T>
std::shared_ptr<FiniteElementBase<T>> CreateFiniteElement(size_t dim, size_t nodes);
