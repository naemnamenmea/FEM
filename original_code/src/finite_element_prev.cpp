#include "stdafx.hpp"
#include "finite_element.hpp"
#include <memory>

using namespace std;

template <typename T>
const GaussIntegr::fIntegrate<1, GAUSS_ORDER_1D> FiniteElement1d<T>::Integrator;

template <typename T>
const GaussIntegr::fIntegrate<2, GAUSS_ORDER_2D> Rectangle<T>::Integrator;

template <typename T>
const GaussIntegr::fIntegrate<3, GAUSS_ORDER_3D> IsoParallelepiped8<T>::Integrator;

template <typename T>
FiniteElementBase<T>::FiniteElementBase()
{
}

template <typename T>
FiniteElementBase<T>::FiniteElementBase(size_t nodes) : m_nodeNumbers(nodes)
{
}

template <typename T>
void FiniteElementBase<T>::ReadProperties(std::istream& /*is*/)
{
}

template <typename T>
const Node3d<T>& NonIsoparametricFE<T>::CoordApproxNode(cardinal_t nodeNo_) const
{
	return *(CoordApproxData[--nodeNo_].pNode);
}

template <typename T>
IsoRod2ConstSec<T>::IsoRod2ConstSec(size_t nodes)
	: FiniteElementBase<T>(nodes), m_crossSectionalArea(0)
{
}

template <typename T>
void IsoRod2ConstSec<T>::ReadProperties(std::istream& is)
{
	is >> m_crossSectionalArea;
}

template <typename T>
IsoQuad4ConsThick<T>::IsoQuad4ConsThick(size_t nodes) : FiniteElementBase<T>(nodes), m_thickness(0)
{
}

template <typename T>
void IsoQuad4ConsThick<T>::ReadProperties(std::istream& is)
{
	is >> m_thickness;
}

template <typename T>
IsoParallelepiped8<T>::IsoParallelepiped8(size_t nodes) : FiniteElementBase<T>(nodes)
{
}

template <typename T>
std::shared_ptr<FiniteElementBase<T>> CreateFiniteElement(size_t dim, size_t nodes)
{
	std::shared_ptr<FiniteElementBase<T>> res;
	switch (dim)
	{
		case 1:
			res = std::make_shared<IsoRod2ConstSec<T>>(nodes);
			break;
		case 2:
			res = std::make_shared<IsoQuad4ConsThick<T>>(nodes);
			break;
		case 3:
			res = std::make_shared<IsoParallelepiped8<T>>(nodes);
			break;
		default:
			break;
	}
	return res;
}

template <typename T>
T FiniteElement1d<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<T, T>::iterator p = Jacobian_DB.find(lc_(1));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<T, T>::iterator, bool> ins =
			Jacobian_DB.insert(map<T, T>::value_type(lc_(1), JacobyMatrix_l2g(lc_).Det()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t1D_FE::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB.insert(typename map<T, T>::value_type(lc_(1), JacobyMatrix_l2g(lc_).Det()))
				.first;
#endif
	}
	return p->second;
}

template <typename T>
T FiniteElement2d<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<T, T>, T>::iterator p = Jacobian_DB.find(make_pair(lc_(1), lc_(2)));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<T, T>, T>::iterator, bool> ins = Jacobian_DB.insert(
			map<pair<T, T>, T>::value_type(make_pair(lc_(1), lc_(2)), JacobyMatrix_l2g(lc_).Det()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t2D_FE::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB
				.insert(typename map<pair<T, T>, T>::value_type(
					make_pair(lc_(1), lc_(2)), JacobyMatrix_l2g(lc_).Det()))
				.first;
#endif
	}
	return p->second;
}

template <typename T>
T FiniteElement3d<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<pair<T, T>, T>, T>::iterator p =
		Jacobian_DB.find(make_pair(make_pair(lc_(1), lc_(2)), lc_(3)));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<pair<T, T>, T>, T>::iterator, bool> ins =
			Jacobian_DB.insert(map<pair<pair<T, T>, T>, T>::value_type(
				make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), JacobyMatrix_l2g(lc_).Det()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t3D_FE::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB
				.insert(typename map<pair<pair<T, T>, T>, T>::value_type(
					make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), JacobyMatrix_l2g(lc_).Det()))
				.first;
#endif
	}
	return p->second;
}

template <typename T>
Tensor2s<T>& FiniteElementBase<T>::JacobyMatrix_l2g(
	const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	result_.Assign0();
	Tensor1s<T> tmp1;
	Tensor2s<T> tmp2;
	const Node3d<T>* pnode;
	for (cardinal_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		pnode = &CoordApproxNode(i);
		result_ += CoordShapeGrad(*pnode, lc_, tmp1).DirectProduct(pnode->Coord(), tmp2);
	}
	return result_;
}

template <typename T>
Tensor1s<T>& FiniteElement1d<T>::Tangent(const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	result_.Assign0();
	Tensor1s<T> nodalRadiusVector, shapeGrad;
	for (cardinal_t i = 1; i <= this->HowManyCoordApproxNodes(); ++i)
	{
		CoordShapeGrad(i, lc_, shapeGrad);
		nodalRadiusVector = FiniteElementBase<T>::Node(i).Coord();
		result_ += (nodalRadiusVector *= shapeGrad(1));
	}
	return result_;
}

template <typename T>
Tensor1s<T>& FiniteElement1d<T>::LocalBasis(
	const Tensor1s<T>& lc_,
	Tensor1s<T>& tangent_,
	Tensor1s<T>& normal1_,
	Tensor1s<T>& normal2_) const
{
	Tangent(lc_, tangent_);
	cardinal_t i = fabs(tangent_(1)) < fabs(tangent_(2)) ? 2 : 1;
	if (fabs(tangent_(i)) < fabs(tangent_(3)))
		i = 1;
	else
		++i;
	tangent_.CrossProductWithBasisVector(i, normal2_).Normalize();
	normal2_.CrossProduct(tangent_, normal1_).Normalize();
	return tangent_;
}

template <typename T>
const Tensor2s<T>& FiniteElement1d<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<T, Tensor2s<T>>::iterator p = JacobyMatrix_DB.find(lc_(1));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<T, Tensor2s>::iterator, bool> ins =
			JacobyMatrix_DB.insert(map<T, Tensor2s>::value_type(lc_(1), Tensor2s()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t1D_FE::JacobyMatrix_l2g(const Tensor1s&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB.insert(map<T, Tensor2s<T>>::value_type(lc_(1), Tensor2s<T>())).first;
#endif
		FiniteElementBase<T>::JacobyMatrix_l2g(lc_, p->second);
		Tensor1s<T> normal1, normal2, tmp;
		LocalBasis(lc_, tmp, normal1, normal2);
		// p->second = Second_vector_of_basis.DiadeProdact(normal1):
		p->second(2, 1) = normal1(1);
		p->second(2, 2) = normal1(2);
		p->second(2, 3) = normal1(3);
		// p->second = Third_vector_of_basis.DiadeProdact(normal2):
		p->second(3, 1) = normal2(1);
		p->second(3, 2) = normal2(2);
		p->second(3, 3) = normal2(3);
	}
	return p->second;
}

template <typename T>
const Tensor2s<T>& FiniteElement2d<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<T, T>, Tensor2s<T>>::iterator p =
		JacobyMatrix_DB.find(make_pair(lc_(1), lc_(2)));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<T, T>, Tensor2s>::iterator, bool> ins = JacobyMatrix_DB.insert(
			map<pair<T, T>, Tensor2s>::value_type(make_pair(lc_(1), lc_(2)), Tensor2s()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t2D_FE::JacobyMatrix_l2g(const Tensor1s&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB
				.insert(map<pair<T, T>, Tensor2s<T>>::value_type(
					make_pair(lc_(1), lc_(2)), Tensor2s<T>()))
				.first;
#endif
		FiniteElementBase<T>::JacobyMatrix_l2g(lc_, p->second);
		Tensor1s<T> normal1;
		UnitNormal(lc_, normal1);
		/* p->second += ThirdBasicVector.DiadeProdact(normal1):*/
		p->second(3, 1) += normal1(1);
		p->second(3, 2) += normal1(2);
		p->second(3, 3) += normal1(3);
		// if (!lc_.is0(3)) {/*must be defined for use of this function beyond the median surface of
		// FE*/}
	}
	return p->second;
}

template <typename T>
const Tensor2s<T>& FiniteElement3d<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<pair<T, T>, T>, Tensor2s<T>>::iterator p =
		JacobyMatrix_DB.find(make_pair(make_pair(lc_(1), lc_(2)), lc_(3)));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<pair<T, T>, T>, Tensor2s>::iterator, bool> ins =
			JacobyMatrix_DB.insert(map<pair<pair<T, T>, T>, Tensor2s>::value_type(
				make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), Tensor2s()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in t3D_FE::JacobyMatrix_l2g(const Tensor1s&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB
				.insert(map<pair<pair<T, T>, T>, Tensor2s<T>>::value_type(
					make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), Tensor2s<T>()))
				.first;
#endif
		FiniteElementBase<T>::JacobyMatrix_l2g(lc_, p->second);
	}
	return p->second;
}

template std::shared_ptr<FiniteElementBase<double>> CreateFiniteElement(size_t dim, size_t nodes);

template class IsoRod2ConstSec<double>;
template class IsoQuad4ConsThick<double>;
template class IsoParallelepiped8<double>;
