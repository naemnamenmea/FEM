#include "FiniteElement.hpp"

#include <vector>

using namespace std;

void tFinitElement::DefineNextNode(const tNode& node_)
{
#ifdef STRONGCHECK
	vector<const tNode*>::iterator ppnode =
		find(m_pNodes.begin(), m_pNodes.end(), (const tNode*)nullptr);
	Assert(ppnode < m_pNodes.end(), "Attempt to link redundant node to FE");
#endif
	*find(m_pNodes.begin(), m_pNodes.end(), (const tNode*)nullptr) = &node_;
}

const SymmetricTensor4s& tFinitElement::MaterialElasTensor(const Tensor1s& lc_) const
{
	return m_pMaterial->ElasTensor(TensorOfRotationToLocal(lc_), m_strainLevel);
}

Tensor1s& tFinitElement::GlobalCoord(const Tensor1s& lc_, Tensor1s& result_) const
{
	result_.Assign0();
	Tensor1s radiusvector;
	const tNode* pnode;
	for (size_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		pnode = &CoordApproxNode(i);
		radiusvector = pnode->Coord();
		result_ += (radiusvector *= CoordShapeFun(*pnode, lc_));
	}
	return result_;
}

Tensor2s& tFinitElement::JacobyMatrix_l2g(const Tensor1s& lc_, Tensor2s& result_) const
{
	result_.Assign0();
	Tensor1s tmp1;
	Tensor2s tmp2;
	const tNode* pnode;
	for (size_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		pnode = &CoordApproxNode(i);
		result_ += CoordShapeGrad(*pnode, lc_, tmp1).DirectProduct(pnode->Coord(), tmp2);
	}
	return result_;
}

Tensor1s& tFinitElement::LocalBasis(
	const Tensor1s& lc_, Tensor1s& e1_, Tensor1s& e2_, Tensor1s& e3_) const
{
	e1_.Assign0();
	e2_.Assign0();
	e3_.Assign0();
	Tensor1s nablaNi, r1, r2, r3;
	for (size_t i = 1; i <= HowManyCoordApproxNodes(); ++i)
	{
		CoordShapeGrad(i, lc_, nablaNi);
		r1 = r2 = r3 = CoordApproxNode(i).Coord();
		e1_ += (r1 *= nablaNi(1));
		e2_ += (r2 *= nablaNi(2));
		e3_ += (r3 *= nablaNi(3));
	}
	return e1_;
}

const Tensor2s& tFinitElement::TensorOfRotationToLocal(const Tensor1s& lc_) const
{
	tMapRotToLoc::iterator p = m_rotToLocDB.find(
		make_pair(lc_(1), make_pair(lc_(2), lc_(3))));	// modify to find with given precision!!!
	if (p == m_rotToLocDB.end())
	{
#ifdef STRONGCHECK
		pair<tMapRotToLoc::iterator, bool> ins = m_rotToLocDB.insert(
			tMapRotToLoc::value_type(make_pair(lc_(1), make_pair(lc_(2), lc_(3))), Tensor2s()));
		Assert(ins.second, "failed insertion in tFinitElement::TensorOfRotationToLocal");
		p = ins.first;
#else
		p = m_rotToLocDB
				.insert(tMapRotToLoc::value_type(
					make_pair(lc_(1), make_pair(lc_(2), lc_(3))), Tensor2s()))
				.first;
#endif
		TensorOfRotationToLocal(lc_, p->second);
	}
	return p->second;
}

Tensor2s& tFinitElement::NablaDispl_g(const Tensor1s& lc_, Tensor2s& result_) const
{
	Tensor1s tmp1;
	Tensor2s tmp2;
	bool not_added(true);
	for (size_t i = 1; i <= HowManyDisplApproxNodes(); ++i)
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

/*Tensor2s& tFinitElement::DisplacementGradient (const Tensor1s& lc_, Tensor2s& result_) const
{
 result_.Assign0();   Tensor1s tmp1;   Tensor2s tmp2;   const tNode* pnode;
 size_t j=0;
 for (size_t i=1; i<=HowManyDisplApproxNodes(); ++i)
   {
  pnode = &Node(i);
  if (!pnode->DisplacementIs0())
	result_ += (DisplShapeGrad(i,lc_,tmp1).DiadeProduct(pnode->Displacement(),tmp2));
   }
 if (!result_.is0())
   result_.ScalarMultiply_left(JacobyMatrix_l2g(lc_,tmp2).Invert());
 return result_;
}*/

// real_t tFinitElement::Mass() const
//{
// return /*Volume() */ Material().Density();
//}

tFinitElement& tNonIsoparametricFE::DefineNextNode(
	const tNode& node_, bool coordApproxing_, bool displApproxing_)
{
#ifdef STRONGCHECK
	Assert(coordApproxing_ || displApproxing_, "Node approxes neither coord nor displ");
#endif
	tFinitElement::DefineNextNode(node_);
	const tNode* null = nullptr;
	vector<tNodeAnd2FunsPtrs>::iterator ppapproxdata;
	if (coordApproxing_)
	{
		ppapproxdata = find_if(
			m_coordApproxData.begin(),
			m_coordApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(null));
#ifdef STRONGCHECK
		Assert(
			ppapproxdata < m_coordApproxData.end(),
			"Attempt to link redundant coord approxing node to FE");
#endif
		ppapproxdata->m_pNode = &node_;
	}
	if (displApproxing_)
	{
		ppapproxdata = find_if(
			m_displApproxData.begin(),
			m_displApproxData.end(),
			tNodeAnd2FunsPtrs::Node_address_is(null));
#ifdef STRONGCHECK
		Assert(
			ppapproxdata < m_displApproxData.end(),
			"Attempt to link redundant coord approxing node to FE");
#endif
		ppapproxdata->m_pNode = &node_;
	}
	return *this;
}

const tNode& tNonIsoparametricFE::CoordApproxNode(size_t nodeNo_) const
{
#ifndef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_coordApproxData.size(), "Invalid coord approx node # in FE");
#endif
	return *(m_coordApproxData[--nodeNo_].m_pNode);
}

const tNode& tNonIsoparametricFE::DisplApproxNode(size_t nodeNo_) const
{
#ifndef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_displApproxData.size(), "Invalid displ approx node # in FE");
#endif
	return *(m_displApproxData[--nodeNo_].m_pNode);
}

template <typename T>
T tIsoparametricFE::tArray<T>::operator()(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && i_ <= this->size(), "invalid index of shape function or its grad in operator()");
#endif

	return this->operator[](--i_);
}

template <typename T>
T tIsoparametricFE::tArray<T>::operator()(const tFinitElement* pFE_, const tNode* pnode_) const
{
#ifdef STRONGCHECK
	Assert(
		pFE_->HowManyNodes() == this->size(),
		"number of nodes in FE is not equal to number of shape funs(grads)");
#endif

	typename vector<T>::const_iterator ptr = this->begin();
	vector<const tNode*>::const_iterator ppnode = pFE_->p1stNode();
	while (ptr != this->end() && *ppnode != pnode_)
	{
		++ptr;
		++ppnode;
	}
#ifdef STRONGCHECK
	Assert(*ppnode == pnode_, "the node is not belongs to the FE");
#endif
	return *ptr;
}

#ifdef STRONGCHECK
tFinitElement& tIsoparametricFE::DefineNextNode(
	const tNode& node_, bool coordApproxing_, bool displApproxing_)
{
	Assert(
		coordApproxing_ && displApproxing_,
		"Node  does not approx coord or displ in isoparametric FE");
	tFinitElement::DefineNextNode(node_);
	return *this;
}
#endif

template class tIsoparametricFE::tArray<real_t (*)(const Tensor1s&)>;
template class tIsoparametricFE::tArray<Tensor1s& (*)(const Tensor1s&, Tensor1s&)>;
