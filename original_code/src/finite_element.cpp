//------------------------------------------------------------------------------
#include "finite_element.hpp"
#include "test_runner.h"
#include "Material.hpp"
#include "NumTypes.h"
#ifdef STRONGCHECK
#include "test_runner.h"
#endif
//------------------------------------------------------------------------------
#include <algorithm>
//------------------------------------------------------------------------------
using namespace std;
using namespace GaussIntegr;
//------------------------------------------------------------------------------
//--------------------------tFinitElement<T>---------------------------------------
//------------------------------------------------------------------------------
static const cardinal_t number_of_FE_kinds = 3;
static const char* initNames[] = {
	"Rod2Const", "MembraneQuad4Const", "Parallelepiped8"};	// Add new kind names here

static tFinitElement<__real_t>* (*initFunctions[])() = {
	tIsoRod2ConstSec<__real_t>::NewFE,
	tIsoQuad4ConsThick<__real_t>::NewFE,
	tIsoParallelepiped8<__real_t>::NewFE};	// Add new pointers to creaters here

template <typename T>
const tNamesToFuns<tFinitElement<T>> tFinitElement<T>::Factory(
	initNames, initFunctions, number_of_FE_kinds);

template <typename T>
const string& tIsoRod2ConstSec<T>::KindName =
	tFinitElement<T>::Factory.FindName(tIsoRod2ConstSec<T>::NewFE);
template <typename T>
const string& tIsoQuad4ConsThick<T>::KindName =
	tFinitElement<T>::Factory.FindName(tIsoQuad4ConsThick<T>::NewFE);
template <typename T>
const string& tIsoParallelepiped8<T>::KindName =
	tFinitElement<T>::Factory.FindName(tIsoParallelepiped8<T>::NewFE);
//------------------------------------------------------------------------------
template <typename T>
void tFinitElement<T>::DefineNextNode(const Node3d<T>& node_)
{
#ifdef STRONGCHECK
	std::vector<const Node3d<T>*>::iterator ppnode =
		std::find(pNodes.begin(), pNodes.end(), (const Node3d<T>*)NULL);
	ASSERT(ppnode < pNodes.end(), "Attempt to link redundant node to FE");
#endif
	*std::find(pNodes.begin(), pNodes.end(), (const Node3d<T>*)NULL) = &node_;
}
//------------------------------------------------------------------------------
template <typename T>
const SymmetricTensor4s<T>& tFinitElement<T>::MaterialElasTensor(const Tensor1s<T>& lc_) const
{
	return this->pMaterial->ElasTensor(TensorOfRotationToLocal(lc_), StrainLevel);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tFinitElement<T>::GlobalCoord(const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	result_.Assign0();
	Tensor1s<T> radiusvector;
	const Node3d<T>* pnode;
	for (cardinal_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		pnode = &CoordApproxNode(i);
		radiusvector = pnode->GetCoord();
		result_ += (radiusvector *= CoordShapeFun(*pnode, lc_));
	}
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor2s<T>& tFinitElement<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	result_.Assign0();
	Tensor1s<T> tmp1;
	Tensor2s<T> tmp2;
	const Node3d<T>* pnode;
	for (cardinal_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		pnode = &CoordApproxNode(i);
		result_ += CoordShapeGrad(*pnode, lc_, tmp1).DirectProduct(pnode->GetCoord(), tmp2);
	}
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tFinitElement<T>::LocalBasis(
	const Tensor1s<T>& lc_, Tensor1s<T>& e1_, Tensor1s<T>& e2_, Tensor1s<T>& e3_) const
{
	e1_.Assign0();
	e2_.Assign0();
	e3_.Assign0();
	Tensor1s<T> nablaNi, r1, r2, r3;
	for (cardinal_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		CoordShapeGrad(i, lc_, nablaNi);
		r1 = r2 = r3 = CoordApproxNode(i).GetCoord();
		e1_ += (r1 *= nablaNi(1));
		e2_ += (r2 *= nablaNi(2));
		e3_ += (r3 *= nablaNi(3));
	}
	return e1_;
}
//------------------------------------------------------------------------------
template <typename T>
const Tensor2s<T>& tFinitElement<T>::TensorOfRotationToLocal(const Tensor1s<T>& lc_) const
{
	typename tMapRotToLoc::iterator p = RotToLoc_DB.find(
		make_pair(lc_(1), make_pair(lc_(2), lc_(3))));	// modify to find with given precision!!!
	if (p == RotToLoc_DB.end())
	{
#ifdef STRONGCHECK
		pair<tMapRotToLoc::iterator, bool> ins = RotToLoc_DB.insert(
			tMapRotToLoc::value_type(make_pair(lc_(1), make_pair(lc_(2), lc_(3))), Tensor2s<T>()));
		ASSERT(ins.second, "failed insertion in tFinitElement<T>::TensorOfRotationToLocal");
		p = ins.first;
#else
		p = RotToLoc_DB
				.insert(tMapRotToLoc::value_type(
					make_pair(lc_(1), make_pair(lc_(2), lc_(3))), Tensor2s<T>()))
				.first;
#endif
		TensorOfRotationToLocal(lc_, p->second);
	}
	return p->second;
}
//------------------------------------------------------------------------------
/*Tensor2s<T>& tFinitElement<T>::DisplacementGradient (const Tensor1s<T>& lc_, Tensor2s<T>& result_)
const
{
 result_.Assign0();   Tensor1s<T> tmp1;   Tensor2s<T> tmp2;   const Node3d<T>* pnode;
 cardinal_t j=0;
 for (cardinal_t i=1; i<=HowManyDisplApproxNodes(); ++i)
   {
	pnode = &Node(i);
	if (!pnode->DisplacementIs0())
	  result_ += (DisplShapeGrad(i,lc_,tmp1).DiadeProduct(pnode->Displacement(),tmp2));
   }
 if (!result_.is0())
   result_.ScalarMultiply_left(JacobyMatrix_l2g(lc_,tmp2).Invert());
 return result_;
}*/
//------------------------------------------------------------------------------
// T tFinitElement<T>::Mass() const
//{
// return /*Volume() */ Material().Density();
//}
//------------------------------------------------------------------------------
template <typename T>
Tensor2s<T>& tFinitElement<T>::NablaDispl_g(const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	Tensor1s<T> tmp1;
	Tensor2s<T> tmp2;
	bool not_added(true);
	for (cardinal_t i = 1; i <= HowManyDisplApproxNodes(); ++i)
		if (!DisplApproxNode(i).DisplacementIs0())
		{
			DisplShapeGrad(DisplApproxNode(i), lc_, tmp1);
			//      tmp1.DotMultiply_left(JacobyMatrix_l2g(lc_,tmp2).Invert()); -- singular matrix
			//      because JacobyMatrix_l2g(lc_,tmp2) is non-virtual
			tmp1.DotMultiply_left((tmp2 = JacobyMatrix_l2g(lc_)).Invert());
			tmp1.DirectProduct(DisplApproxNode(i).Displacement(), tmp2);
			if (not_added)
			{
				result_ = tmp2;
				not_added = false;
			}
			else
				result_ += tmp2;
		}
	if (not_added)
		return result_.Assign0();

	return result_;
}
//------------------------------------------------------------------------------
//--------------------------tNonIsoparametricFE<T>---------------------------------
//------------------------------------------------------------------------------
template <typename T>
tFinitElement<T>& tNonIsoparametricFE<T>::DefineNextNode(
	const Node3d<T>& node_, bool coordApproxing_, bool displApproxing_)
{
#ifdef STRONGCHECK
	ASSERT(coordApproxing_ || displApproxing_, "Node approxes neither coord nor displ");
#endif
	tFinitElement<T>::DefineNextNode(node_);
	const Node3d<T>* null = NULL;
	typename vector<tNodeAnd2FunsPtrs>::iterator ppapproxdata;
	if (coordApproxing_)
	{
		ppapproxdata = find_if(
			CoordApproxData.begin(),
			CoordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(null));
#ifdef STRONGCHECK
		ASSERT(
			ppapproxdata < CoordApproxData.end(),
			"Attempt to link redundant coord approxing node to FE");
#endif
		ppapproxdata->pNode = &node_;
	}
	if (displApproxing_)
	{
		ppapproxdata = find_if(
			DisplApproxData.begin(),
			DisplApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(null));
#ifdef STRONGCHECK
		ASSERT(
			ppapproxdata < DisplApproxData.end(),
			"Attempt to link redundant coord approxing node to FE");
#endif
		ppapproxdata->pNode = &node_;
	}
	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
const Node3d<T>& tNonIsoparametricFE<T>::CoordApproxNode(cardinal_t nodeNo_) const
{
#ifndef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= CoordApproxData.size(), "Invalid coord approx node # in FE");
#endif
	return *(CoordApproxData[--nodeNo_].pNode);
}
//------------------------------------------------------------------------------
template <typename T>
const Node3d<T>& tNonIsoparametricFE<T>::DisplApproxNode(cardinal_t nodeNo_) const
{
#ifndef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= DisplApproxData.size(), "Invalid displ approx node # in FE");
#endif
	return *(DisplApproxData[--nodeNo_].pNode);
}
//------------------------------------------------------------------------------
//--------------------------tIsoparametricFE<T>------------------------------------
//------------------------------------------------------------------------------
template <typename T>
template <typename R>
R tIsoparametricFE<T>::tArray<R>::operator()(cardinal_t i_) const
{
#ifdef STRONGCHECK
	ASSERT(
		i_ > 0 && i_ <= this->size(), "invalid index of shape function or its grad in operator()");
#endif

	return this->operator[](--i_);
}
//------------------------------------------------------------------------------
template <typename T>
template <typename R>
R tIsoparametricFE<T>::tArray<R>::operator()(
	const tFinitElement<T>* pFE_, const Node3d<T>* pnode_) const
{
#ifdef STRONGCHECK
	ASSERT(
		pFE_->HowManyNodes() == this->size(),
		"number of nodes in FE is not equal to number of shape funs(grads)");
#endif

	typename vector<R>::const_iterator ptr = this->begin();
	typename vector<const Node3d<T>*>::const_iterator ppnode = pFE_->p1stNode();
	while (ptr != this->end() && *ppnode != pnode_)
	{
		++ptr;
		++ppnode;
	}
#ifdef STRONGCHECK
	ASSERT(*ppnode == pnode_, "the node is not belongs to the FE");
#endif
	return *ptr;
}
//------------------------------------------------------------------------------
#ifdef STRONGCHECK
tFinitElement<T>& tIsoparametricFE<T>::DefineNextNode(
	const Node3d<T>& node_, bool coordApproxing_, bool displApproxing_)
{
	ASSERT(
		coordApproxing_ && displApproxing_,
		"Node  does not approx coord or displ in isoparametric FE");
	tFinitElement<T>::DefineNextNode(node_);
	return *this;
}
#endif
//------------------------------------------------------------------------------
//-----------------------------t1D_FE<T>-------------------------------------------
//------------------------------------------------------------------------------
template <typename T>
const fIntegrate<1, GAUSS_ORDER_1D> t1D_FE<T>::Integrator;
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& t1D_FE<T>::GlobalCoord(const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	tFinitElement<T>::GlobalCoord(lc_, result_);
	// if (!lc_.is0(2)  ||  !lc_.is0(3))
	//   {
	Tensor1s<T> normal1, normal2, tmp;
	LocalBasis(lc_, tmp, normal1, normal2);
	(result_ += (normal1 *= lc_(2))) += (normal2 *= lc_(3));
	//   }
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
const Tensor2s<T>& t1D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<T, Tensor2s<T>>::iterator p = JacobyMatrix_DB.find(lc_(1));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<T, Tensor2s<T>>::iterator, bool> ins =
			JacobyMatrix_DB.insert(map<T, Tensor2s<T>>::value_type(lc_(1), Tensor2s<T>()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t1D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB.insert(map<T, Tensor2s<T>>::value_type(lc_(1), Tensor2s<T>())).first;
#endif
		tFinitElement<T>::JacobyMatrix_l2g(lc_, p->second);
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
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& t1D_FE<T>::Tangent(const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	result_.Assign0();
	Tensor1s<T> nodalRadiusVector, shapeGrad;
	for (cardinal_t i = 1; i <= this->HowManyCoordApproxNodes(); ++i)
	{
		this->CoordShapeGrad(i, lc_, shapeGrad);
		nodalRadiusVector = this->Node(i).GetCoord();
		result_ += (nodalRadiusVector *= shapeGrad(1));
	}
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& t1D_FE<T>::LocalBasis(
	const Tensor1s<T>& lc_,
	Tensor1s<T>& tangent_,
	Tensor1s<T>& normal1_,
	Tensor1s<T>& normal2_) const
{
	Tangent(lc_, tangent_);
	cardinal_t i = fabs(tangent_(1)) < fabs(tangent_(2)) ? 2 : 1;
	if (fabs(tangent_(static_cast<small_t>(i))) < fabs(tangent_(3)))
		i = 1;
	else
		++i;
	tangent_.CrossProductWithBasisVector(static_cast<small_t>(i), normal2_).Normalize();
	normal2_.CrossProduct(tangent_, normal1_).Normalize();
	return tangent_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor2s<T>& t1D_FE<T>::TensorOfRotationToLocal(const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	Tensor1s<T> rotatedBasis1, rotatedBasis2, rotatedBasis3;
	// Assigning: rB1 - tangent to the FE direction, rB2,rB3 - unit normals to the direction:
	LocalBasis(lc_, rotatedBasis1, rotatedBasis2, rotatedBasis3)
		.Normalize();  // tangent is normalised
	// now rB123 correspond to the basis of global coord system, rotated to the FE surface
	return result_.AssignRows(rotatedBasis1, rotatedBasis2, rotatedBasis3);
}
//------------------------------------------------------------------------------
template <typename T>
const SymmetricTensor4s<T>& t1D_FE<T>::MaterialElasTensor(const Tensor1s<T>& lc_) const
{
	if (StateIs1dStress)
		return this->pMaterial->ElasTensor_1D(tFinitElement<T>::TensorOfRotationToLocal(lc_));
	else
		return tFinitElement<T>::MaterialElasTensor(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
T t1D_FE<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<T, T>::iterator p = Jacobian_DB.find(lc_(1));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<T, T>::iterator, bool> ins =
			Jacobian_DB.insert(map<T, T>::value_type(lc_(1), JacobyMatrix_l2g(lc_).Det()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t1D_FE<T>::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB.insert(map<T, T>::value_type(lc_(1), JacobyMatrix_l2g(lc_).Det())).first;
#endif
	}
	return p->second;
}
//------------------------------------------------------------------------------
template <typename T>
T t2D_FE<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<T, T>, T>::iterator p = Jacobian_DB.find(make_pair(lc_(1), lc_(2)));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<T, T>, T>::iterator, bool> ins = Jacobian_DB.insert(
			map<pair<T, T>, T>::value_type(make_pair(lc_(1), lc_(2)), JacobyMatrix_l2g(lc_).Det()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t2D_FE<T>::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB
				.insert(map<pair<T, T>, T>::value_type(
					make_pair(lc_(1), lc_(2)), JacobyMatrix_l2g(lc_).Det()))
				.first;
#endif
	}
	return p->second;
}
//------------------------------------------------------------------------------
template <typename T>
T t3D_FE<T>::Jacobian_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<pair<T, T>, T>, T>::iterator p =
		Jacobian_DB.find(make_pair(make_pair(lc_(1), lc_(2)), lc_(3)));
	if (p == Jacobian_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<pair<T, T>, T>, T>::iterator, bool> ins =
			Jacobian_DB.insert(map<pair<pair<T, T>, T>, T>::value_type(
				make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), JacobyMatrix_l2g(lc_).Det()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t3D_FE<T>::Jacobian_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = Jacobian_DB
				.insert(map<pair<pair<T, T>, T>, T>::value_type(
					make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), JacobyMatrix_l2g(lc_).Det()))
				.first;
#endif
	}
	return p->second;
}
//------------------------------------------------------------------------------
//-----------------------------t2D_FE<T>-------------------------------------------
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& t2D_FE<T>::LocalBasis(
	const Tensor1s<T>& lc_,
	Tensor1s<T>& tangent1_,
	Tensor1s<T>& tangent2_,
	Tensor1s<T>& normal_) const
{
	tangent1_.Assign0();
	tangent2_.Assign0();
	Tensor1s<T> nablaNi, radiusVectOfNodei;

	for (cardinal_t i = 1; i <= this->HowManyCoordApproxNodes(); ++i)
	{
		this->CoordShapeGrad(i, lc_, nablaNi);
		radiusVectOfNodei = normal_ =
			this->CoordApproxNode(i).GetCoord();  // here normal_ is used as tmp storage
		tangent1_ += (radiusVectOfNodei *= nablaNi(1));
		tangent2_ += (normal_ *= nablaNi(2));
	}
	tangent1_.CrossProduct(tangent2_, normal_).Normalize();
	return tangent1_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& t2D_FE<T>::GlobalCoord(const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	tFinitElement<T>::GlobalCoord(lc_, result_);
	Tensor1s<T> tmp;
	// if (!lc_.is0(3))
	result_ += (UnitNormal(lc_, tmp) *= lc_(3));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
const Tensor2s<T>& t2D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<T, T>, Tensor2s<T>>::iterator p =
		JacobyMatrix_DB.find(make_pair(lc_(1), lc_(2)));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<T, T>, Tensor2s<T>>::iterator, bool> ins = JacobyMatrix_DB.insert(
			map<pair<T, T>, Tensor2s<T>>::value_type(make_pair(lc_(1), lc_(2)), Tensor2s<T>()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t2D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB
				.insert(map<pair<T, T>, Tensor2s<T>>::value_type(
					make_pair(lc_(1), lc_(2)), Tensor2s<T>()))
				.first;
#endif
		tFinitElement<T>::JacobyMatrix_l2g(lc_, p->second);
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
//------------------------------------------------------------------------------
template <typename T>
Tensor2s<T>& t2D_FE<T>::TensorOfRotationToLocal(const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	Tensor1s<T> rotatedBasis1, rotatedBasis2, rotatedBasis3;
	// Assigning: rB1, rB2 - coord.vectors of FE surface, rB3 - normal to the surface:
	LocalBasis(lc_, rotatedBasis1, rotatedBasis2, rotatedBasis3);
	// Because rB1, rB2 are not orthogonal, rB2 is replaced by vector orthogonal to rB1,rB3:
	rotatedBasis1.Normalize();
	rotatedBasis3.CrossProduct(rotatedBasis1, rotatedBasis2);
	// now rB123 corresponds to the basis of global coord system, rotated to FE surface
	return result_.AssignRows(rotatedBasis1, rotatedBasis2, rotatedBasis3);
}
//------------------------------------------------------------------------------
template <typename T>
const SymmetricTensor4s<T>& t2D_FE<T>::MaterialElasTensor(const Tensor1s<T>& lc_) const
{
	if (StateIs2dStress)
		return this->pMaterial->ElasTensor_2D(tFinitElement<T>::TensorOfRotationToLocal(lc_));
	else
		return tFinitElement<T>::MaterialElasTensor(lc_);
}
//------------------------------------------------------------------------------
//-----------------------------t3D_FE<T>-------------------------------------------
//------------------------------------------------------------------------------
template <typename T>
const Tensor2s<T>& t3D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>& lc_) const
{
	typename map<pair<pair<T, T>, T>, Tensor2s<T>>::iterator p =
		JacobyMatrix_DB.find(make_pair(make_pair(lc_(1), lc_(2)), lc_(3)));
	if (p == JacobyMatrix_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<pair<pair<T, T>, T>, Tensor2s<T>>::iterator, bool> ins =
			JacobyMatrix_DB.insert(map<pair<pair<T, T>, T>, Tensor2s<T>>::value_type(
				make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), Tensor2s<T>()));
		ASSERT(
			ins.second,
			"failed insertion of ElasTensor in t3D_FE<T>::JacobyMatrix_l2g(const Tensor1s<T>&)");
		p = ins.first;
#else
		p = JacobyMatrix_DB
				.insert(map<pair<pair<T, T>, T>, Tensor2s<T>>::value_type(
					make_pair(make_pair(lc_(1), lc_(2)), lc_(3)), Tensor2s<T>()))
				.first;
#endif
		tFinitElement<T>::JacobyMatrix_l2g(lc_, p->second);
	}
	return p->second;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor2s<T>& t3D_FE<T>::TensorOfRotationToLocal(const Tensor1s<T>& lc_, Tensor2s<T>& result_) const
{
	Tensor1s<T> rotatedBasis1, rotatedBasis2, rotatedBasis3;
	this->LocalBasis(lc_, rotatedBasis1, rotatedBasis2, rotatedBasis3);
	rotatedBasis1.Normalize().CrossProduct(rotatedBasis2, rotatedBasis3);
	rotatedBasis3.Normalize().CrossProduct(rotatedBasis1, rotatedBasis2);
	return result_.AssignRows(rotatedBasis1, rotatedBasis2, rotatedBasis3);
}
//------------------------------------------------------------------------------
//-----------------------------tRectangle<T>---------------------------------------
//------------------------------------------------------------------------------
template <typename T>
const fIntegrate<2, GAUSS_ORDER_2D> tRectangle<T>::Integrator;
//------------------------------------------------------------------------------
//-----------------------------tIsoRod2ConstSec<T>--------------------------
//------------------------------------------------------------------------------
template <typename T>
T tIsoRod2ConstSec<T>::ShapeFun1(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(lc_(1) >= -1. && lc_(1) <= 1., "value of local coord must be between -1 & 1");
#endif
	return 0.5 * (1. - lc_(1));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoRod2ConstSec<T>::ShapeFun2(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(lc_(1) >= -1. && lc_(1) <= 1., "value of local coord must be between -1 & 1");
#endif
	return 0.5 * (1. + lc_(1));
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoRod2ConstSec<T>::ShapeGrad1(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(lc_(1) >= -1. && lc_(1) <= 1., "value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, -0.5);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoRod2ConstSec<T>::ShapeGrad2(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(lc_(1) >= -1. && lc_(1) <= 1., "value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.5);
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
const typename tIsoRod2ConstSec<T>::tPtrShapeFunArray pShapeFunctions;
template <typename T>
const typename tIsoRod2ConstSec<T>::tPtrShapeFunGradArray pShapeFunGrads;
//------------------------------------------------------------------------------
template <typename T>
T tIsoRod2ConstSec<T>::CoordShapeFun(cardinal_t no_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(no_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoRod2ConstSec<T>::CoordShapeGrad(
	cardinal_t no_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(no_)))(lc_, result_);
}
//------------------------------------------------------------------------------
/*T   tIsoRod2ConstSec<T>::DisplShapeFun (cardinal_t no_, const Tensor1s<T>& lc_) const
{
 return (*(pShapeFunctions(no_)))(lc_);
}
//------------------------------------------------------------------------------
Tensor1s<T>& tIsoRod2ConstSec<T>::DisplShapeGrad(cardinal_t no_, const Tensor1s<T>& lc_,
Tensor1s<T>& result_) const
{
 return (*(pShapeFunGrads(no_)))(lc_,result_);
}*/
//------------------------------------------------------------------------------
template <typename T>
T tIsoRod2ConstSec<T>::CoordShapeFun(const Node3d<T>& node_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(this, &node_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoRod2ConstSec<T>::CoordShapeGrad(
	const Node3d<T>& node_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(this, &node_)))(lc_, result_);
}
//------------------------------------------------------------------------------
/*T tIsoRod2ConstSec<T>::DisplShapeFun (const Node3d<T>& node_, const Tensor1s<T>& lc_) const
{
 return (*(pShapeFunctions(this,&node_)))(lc_);
}
//------------------------------------------------------------------------------
Tensor1& tIsoRod2ConstSec<T>::DisplShapeGrad (const Node3d<T>& node_, const Tensor1s<T>& lc_,
Tensor1s<T>& result_) const
{
 return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
}*/
//------------------------------------------------------------------------------
//-----------------------------tIsoQuad4ConsThick<T>-------------------------------
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::ShapeFun1(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.25 * (1. + lc_(1)) * (1. - lc_(2));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::ShapeFun2(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.25 * (1. - lc_(1)) * (1. - lc_(2));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::ShapeFun3(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.25 * (1. - lc_(1)) * (1. + lc_(2));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::ShapeFun4(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.25 * (1. + lc_(1)) * (1. + lc_(2));
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::ShapeGrad1(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();	// to do assign only last component
	result_.Assign(1, 0.25 * (1. - lc_(2)));
	result_.Assign(2, 0.25 * (-1. - lc_(1)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::ShapeGrad2(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.25 * (-1. + lc_(2)));
	result_.Assign(2, 0.25 * (-1. + lc_(1)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::ShapeGrad3(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		tMessage("value of local coord must be between -1 & 1\nvaluse are:\n") << lc_(1) << '\n'
																			   << lc_(2));
#endif
	result_.Assign0();
	result_.Assign(1, 0.25 * (-1. - lc_(2)));
	result_.Assign(2, 0.25 * (1. - lc_(1)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::ShapeGrad4(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.25 * (1. + lc_(2)));
	result_.Assign(2, 0.25 * (1. + lc_(1)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
tIsoQuad4ConsThick<T>::tPtrShapeFunArray::tPtrShapeFunArray()
	: tIsoparametricFE<T>::tArrayOfPtrsToShapeFunctions(4)
{
	this->operator[](0) = ShapeFun1;
	this->operator[](1) = ShapeFun2;
	this->operator[](2) = ShapeFun3;
	this->operator[](3) = ShapeFun4;
}
//------------------------------------------------------------------------------
template <typename T>
const typename tIsoQuad4ConsThick<T>::tPtrShapeFunArray tIsoQuad4ConsThick<T>::pShapeFunctions;
//------------------------------------------------------------------------------
template <typename T>
tIsoQuad4ConsThick<T>::tPtrShapeFunGradArray::tPtrShapeFunGradArray()
	: tIsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients(4)
{
	this->operator[](0) = ShapeGrad1;
	this->operator[](1) = ShapeGrad2;
	this->operator[](2) = ShapeGrad3;
	this->operator[](3) = ShapeGrad4;
}
//------------------------------------------------------------------------------
template <typename T>
const typename tIsoQuad4ConsThick<T>::tPtrShapeFunGradArray tIsoQuad4ConsThick<T>::pShapeFunGrads;
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::CoordShapeFun(cardinal_t no_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(no_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::CoordShapeGrad(
	cardinal_t no_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(no_)))(lc_, result_);
}
//------------------------------------------------------------------------------
/*T   tIsoQuad4ConsThick<T>::DisplShapeFun (cardinal_t no_, const Tensor1s<T>& lc_) const
{
 return (*(pShapeFunctions(no_)))(lc_);
}
//------------------------------------------------------------------------------
Tensor1s<T>& tIsoQuad4ConsThick<T>::DisplShapeGrad(cardinal_t no_, const Tensor1s<T>& lc_,
Tensor1s<T>& result_) const
{
 return (*(pShapeFunGrads(no_)))(lc_,result_);
}*/
//------------------------------------------------------------------------------
template <typename T>
T tIsoQuad4ConsThick<T>::CoordShapeFun(const Node3d<T>& node_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(this, &node_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoQuad4ConsThick<T>::CoordShapeGrad(
	const Node3d<T>& node_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(this, &node_)))(lc_, result_);
}
//------------------------------------------------------------------------------
/*T tIsoQuad4ConsThick<T>::DisplShapeFun (const Node3d<T>& node_, const Tensor1& lc_) const
{
 return (*(pShapeFunctions(this,&node_)))(lc_);
}
//------------------------------------------------------------------------------
Tensor1s<T>& tIsoQuad4ConsThick<T>::DisplShapeGrad (const Node3d<T>& node_, const Tensor1s<T>& lc_,
Tensor1s<T>& result_) const
{
 return (*(pShapeFunGrads(this,&node_)))(lc_,result_);
}*/
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template <typename T>
const fIntegrate<3, GAUSS_ORDER_3D> tIsoParallelepiped8<T>::Integrator;
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun1(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. + lc_(1)) * (1. - lc_(2)) * (1. - lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun2(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. - lc_(1)) * (1. - lc_(2)) * (1. - lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun3(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. - lc_(3)) * (1. + lc_(2)) * (1. - lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun4(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. + lc_(1)) * (1. + lc_(2)) * (1. - lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun5(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. + lc_(1)) * (1. - lc_(2)) * (1. + lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun6(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. - lc_(1)) * (1. - lc_(2)) * (1. + lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun7(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. - lc_(1)) * (1. + lc_(2)) * (1. + lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::ShapeFun8(const Tensor1s<T>& lc_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	return 0.125 * (1. + lc_(1)) * (1. + lc_(2)) * (1. + lc_(3));
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad1(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.125 * (1. - lc_(2)) * (1. - lc_(3)));
	result_.Assign(2, -0.125 * (1. + lc_(1)) * (1. - lc_(3)));
	result_.Assign(3, -0.125 * (1. + lc_(1)) * (1. - lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad2(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, -0.125 * (1. - lc_(2)) * (1. - lc_(3)));
	result_.Assign(2, -0.125 * (1. - lc_(1)) * (1. - lc_(3)));
	result_.Assign(3, -0.125 * (1. - lc_(1)) * (1. - lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad3(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, -0.125 * (1. + lc_(2)) * (1. - lc_(3)));
	result_.Assign(2, 0.125 * (1. - lc_(1)) * (1. - lc_(3)));
	result_.Assign(3, -0.125 * (1. - lc_(1)) * (1. + lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad4(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.125 * (1. + lc_(2)) * (1. - lc_(3)));
	result_.Assign(2, 0.125 * (1. + lc_(1)) * (1. - lc_(3)));
	result_.Assign(3, -0.125 * (1. + lc_(1)) * (1. + lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad5(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.125 * (1. - lc_(2)) * (1. + lc_(3)));
	result_.Assign(2, -0.125 * (1. + lc_(1)) * (1. + lc_(3)));
	result_.Assign(3, 0.125 * (1. + lc_(1)) * (1. - lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad6(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, -0.125 * (1. - lc_(2)) * (1. + lc_(3)));
	result_.Assign(2, -0.125 * (1. - lc_(1)) * (1. + lc_(3)));
	result_.Assign(3, 0.125 * (1. - lc_(1)) * (1. - lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad7(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, -0.125 * (1. + lc_(2)) * (1. + lc_(3)));
	result_.Assign(2, 0.125 * (1. - lc_(1)) * (1. + lc_(3)));
	result_.Assign(3, 0.125 * (1. - lc_(1)) * (1. + lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::ShapeGrad8(const Tensor1s<T>& lc_, Tensor1s<T>& result_)
{
#ifdef STRONGCHECK
	ASSERT(
		lc_(1) >= -1. && lc_(1) <= 1. && lc_(2) >= -1. && lc_(2) <= 1. && lc_(3) >= -1. &&
			lc_(3) <= 1.,
		"value of local coord must be between -1 & 1");
#endif
	result_.Assign0();
	result_.Assign(1, 0.125 * (1. + lc_(2)) * (1. + lc_(3)));
	result_.Assign(2, 0.125 * (1. + lc_(1)) * (1. + lc_(3)));
	result_.Assign(3, 0.125 * (1. + lc_(1)) * (1. + lc_(2)));
	return result_;
}
//------------------------------------------------------------------------------
template <typename T>
tIsoParallelepiped8<T>::tPtrShapeFunArray::tPtrShapeFunArray()
	: tIsoparametricFE<T>::tArrayOfPtrsToShapeFunctions(8)
{
	this->operator[](0) = ShapeFun1;
	this->operator[](1) = ShapeFun2;
	this->operator[](2) = ShapeFun3;
	this->operator[](3) = ShapeFun4;
	this->operator[](4) = ShapeFun5;
	this->operator[](5) = ShapeFun6;
	this->operator[](6) = ShapeFun7;
	this->operator[](7) = ShapeFun8;
}
//------------------------------------------------------------------------------
template <typename T>
const typename tIsoParallelepiped8<T>::tPtrShapeFunArray tIsoParallelepiped8<T>::pShapeFunctions;
//------------------------------------------------------------------------------
template <typename T>
tIsoParallelepiped8<T>::tPtrShapeFunGradArray::tPtrShapeFunGradArray()
	: tIsoparametricFE<T>::tArrayOfPtrsToShapeFunGradients(8)
{
	this->operator[](0) = ShapeGrad1;
	this->operator[](1) = ShapeGrad2;
	this->operator[](2) = ShapeGrad3;
	this->operator[](3) = ShapeGrad4;
	this->operator[](4) = ShapeGrad5;
	this->operator[](5) = ShapeGrad6;
	this->operator[](6) = ShapeGrad7;
	this->operator[](7) = ShapeGrad8;
}
//------------------------------------------------------------------------------
template <typename T>
const typename tIsoParallelepiped8<T>::tPtrShapeFunGradArray tIsoParallelepiped8<T>::pShapeFunGrads;
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::CoordShapeFun(cardinal_t no_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(no_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::CoordShapeGrad(
	cardinal_t no_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(no_)))(lc_, result_);
}
//------------------------------------------------------------------------------
template <typename T>
T tIsoParallelepiped8<T>::CoordShapeFun(const Node3d<T>& node_, const Tensor1s<T>& lc_) const
{
	return (*(pShapeFunctions(this, &node_)))(lc_);
}
//------------------------------------------------------------------------------
template <typename T>
Tensor1s<T>& tIsoParallelepiped8<T>::CoordShapeGrad(
	const Node3d<T>& node_, const Tensor1s<T>& lc_, Tensor1s<T>& result_) const
{
	return (*(pShapeFunGrads(this, &node_)))(lc_, result_);
}
//------------------------------------------------------------------------------
// template<typename T>
// T& t3D_FE<T>::Integrate3D (const tFinitElement<T>::fIntegrand<T>& function_, T& result_) const
//{
//  Tensor1s<T> localCoord(0.);
//  T gaussWeightX, gaussWeightY, gaussWeightZ;   cardinal_t i=0, j, k;
//  const vector<tFinitElement<T>::tGaussCoeff> &gausscoeffsX = *pGaussCoeffsX,
//                                           &gausscoeffsY = *pGaussCoeffsY,
//                                           &gausscoeffsZ = *pGaussCoeffsZ;
//  for (; i<gausscoeffsX.size(); ++i)
//    {
//     gausscoeffsX[i].WhatCoeffs(localCoord(1),gaussWeightX);
//     for (j=0; j<gausscoeffsY.size(); ++j)
//       {
//        gausscoeffsY[j].WhatCoeffs(localCoord(2),gaussWeightY);
//        for (k=0; k<gausscoeffsZ.size(); ++k)
//          {
//           gausscoeffsZ[k].WhatCoeffs(localCoord(3),gaussWeightZ);
//           result_ += (function_(localCoord) *=
//                       (gaussWeightZ * gaussWeightY * gaussWeightX *
//                        Jacobian_l2g(localCoord)));
//          }
//       }
//    }
//  return result_;
//}

template <typename T>
std::shared_ptr<tFinitElement<T>> CreateFiniteElement(size_t dim, size_t nodes)
{
	std::shared_ptr<tFinitElement<T>> res;
	switch (dim)
	{
		case 1:
			res = std::make_shared<tIsoRod2ConstSec<T>>(nodes);
			break;
		case 2:
			res = std::make_shared<tIsoQuad4ConsThick<T>>(nodes);
			break;
		case 3:
			res = std::make_shared<tIsoParallelepiped8<T>>(nodes);
			break;
		default:
			break;
	}
	return res;
}