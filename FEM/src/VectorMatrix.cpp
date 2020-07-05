
#include "VectorMatrix.hpp"
#include "Node.hpp"
#include "FiniteElementModel.hpp"
#include "HybridTensors.hpp"
#include "GraphTheory.hpp"
#include "ThrowMessage.hpp"
#include "ProcessDisplay.hpp"

using namespace std;

#include <numeric>

tNodalTensor1Column::tNodalTensor1Column(const tFE_model& model_, const string& name_)
	: m_caption(name_),
	  m_components(model_.HowManyNodes()),
	  m_constNodeCursor(m_components.begin()),
	  m_nodeCursor(m_components.begin())
{
	vector<tNodalTensor1>::iterator p = m_components.begin();
	for (size_t i = 1; p < m_components.end(); ++i, ++p) p->Link(model_.GetNode(i));
}

tNodalTensor1Column& tNodalTensor1Column::Link(const tFE_model& model_)
{
#ifdef STRONGCHECK
	Assert(m_components.empty(), "data already allocated in tNodalTensor1Column::Link");
#endif
	m_components.resize(model_.HowManyNodes());
	vector<tNodalTensor1>::iterator p = m_components.begin();
	for (size_t i = 1; p < m_components.end(); ++i, ++p) p->Link(model_.GetNode(i));
	m_constNodeCursor = m_components.begin();
	m_nodeCursor = m_components.begin();
	return *this;
}

struct fAddDOF
{
	size_t operator()(size_t prev_sum_, const tNodalTensor1& o_) const
	{
		return prev_sum_ + o_.Node().NumberOfDOFs();
	}
};

size_t tNodalTensor1Column::NumberOfDOFs() const
{
	return accumulate(m_components.begin(), m_components.end(), size_t(0), fAddDOF());
}

struct tNode_addr_is
{
	tNode_addr_is(const tNode* address_) : m_address(address_)
	{
	}
	bool operator()(const tNodalTensor1& v_) const
	{
		return &v_.Node() == m_address;
	}

	const tNode* const m_address;
};

const tNodalTensor1& tNodalTensor1Column::operator()(const tNode& node_) const
{
#ifdef STRONGCHECK
	Assert(
		find_if(m_components.begin(), m_components.end(), tNode_addr_is(&node_)) <
			m_components.end(),
		"Looking for absent node in tNodalTensor1Column::operator(tNode&) const");
#endif
	if (m_constNodeCursor == m_components.end())
		m_constNodeCursor = m_components.begin();
	m_constNodeCursor = find_if(m_constNodeCursor, m_components.end(), tNode_addr_is(&node_));
	if (m_constNodeCursor == m_components.end())
		m_constNodeCursor =
			find_if(m_components.begin(), m_components.end(), tNode_addr_is(&node_));
#ifdef STRONGCHECK
	Assert(
		m_constNodeCursor >= m_components.begin() && m_constNodeCursor < m_components.end(),
		"implementation error in tNodalTensor1Column::operator(tNode&) const");
#endif
	return *(m_constNodeCursor++);
	// return *find_if(Components.begin(),Components.end(),tNode_addr_is(&node_));
}

tNodalTensor1& tNodalTensor1Column::operator()(const tNode& node_)
{
#ifdef STRONGCHECK
	Assert(
		find_if(m_components.begin(), m_components.end(), tNode_addr_is(&node_)) <
			m_components.end(),
		"Looking for absent node in tNodalTensor1Column::operator(tNode&)");
#endif
	if (m_nodeCursor == m_components.end())
		m_nodeCursor = m_components.begin();
	m_nodeCursor = find_if(m_nodeCursor, m_components.end(), tNode_addr_is(&node_));
	if (m_nodeCursor == m_components.end())
		m_nodeCursor = find_if(m_components.begin(), m_components.end(), tNode_addr_is(&node_));
#ifdef STRONGCHECK
	Assert(
		m_nodeCursor >= m_components.begin() && m_nodeCursor < m_components.end(),
		"implementation error in tNodalTensor1Column::operator(tNode&)");
#endif
	return *(m_nodeCursor++);
	// return *find_if(Components.begin(),Components.end(),tNode_addr_is(&node_));
}

tNodalTensor1Column& tNodalTensor1Column::operator+=(const tNodalTensor1Column& o_)
{
	vector<tNodalTensor1>::iterator p = m_components.begin();
	vector<tNodalTensor1>::const_iterator po = o_.m_components.begin();
	for (; p < m_components.end(); ++p, ++po)
	{
#ifdef STRONGCHECK
		Assert(
			&p->Node() == &po->Node(),
			"Non-corresponding NodalTensors in tNodalTensor1Column::operator+=");
#endif
		*p += *po;
	}
	return *this;
}

tNodalTensor1Column& tNodalTensor1Column::operator-=(const tNodalTensor1Column& o_)
{
	vector<tNodalTensor1>::iterator p = m_components.begin();
	vector<tNodalTensor1>::const_iterator po = o_.m_components.begin();
	for (; p < m_components.end(); ++p, ++po)
	{
#ifdef STRONGCHECK
		Assert(
			&p->Node() == &po->Node(),
			"Non-corresponding NodalTensors in tNodalTensor1Column::operator-=");
#endif
		*p -= *po;
	}
	return *this;
}

/*tFEsTensor2Column::tFEsTensor2Column (const tFE_model& model_, const string& name_):
	   Caption(name_), Components(model_.HowManyFinEls())
{
 vector<tFEsTensor2>::iterator p=Components.begin();
 for (size_t i=1; p<Components.end(); ++i,++p)
		   p->Link(model_.FE(i));
}

tFEsTensor2Column& tFEsTensor2Column::Link(const tFE_model& model_)
{
 Components.resize(model_.HowManyFinEls());
 vector<tFEsTensor2>::iterator p = Components.begin();
 for (size_t i=1; p<Components.end(); ++i,++p)
		   p->Link(model_.FE(i));
 return *this;
}

tFEsTensor2Column& tFEsTensor2Column::Transform (const tFEsTensor2Column::fTransformer& function_)
{
// for_each(Components.begin(),Components.end(),function_);
 for (vector<tFEsTensor2>::iterator pc=Components.begin(); pc<Components.end(); ++pc)
		   function_(*pc);
 return *this;
}

struct tFE_addr_is
{
  const tFinitElement* const Address;
  tFE_addr_is (const tFinitElement* address_): Address(address_) {}
  bool operator() (const tFEsTensor2& v_) const
	  {return &v_.FE() == Address;}
};

const tFEsTensor2& tFEsTensor2Column::operator() (const tFinitElement& fe_) const
{
#ifdef STRONGCHECK
 Assert(find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_)) < Components.end(), "Looking
for absent FE in tFEsTensor2Column::operator(tFinitElement&) const"); #endif static
vector<tFEsTensor2>::const_iterator cursor = Components.begin(); cursor =
find_if(cursor,Components.end(),tFE_addr_is(&fe_)); if (cursor == Components.end()) cursor =
find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_)); #ifdef STRONGCHECK Assert(cursor <
Components.end(), "implementation error in tFEsTensor2Column::operator(tFinitElement&) const");
#endif
 return *(cursor++);
// return *find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_));
}

tFEsTensor2& tFEsTensor2Column::operator() (const tFinitElement& fe_)
{
#ifdef STRONGCHECK
 Assert(find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_)) < Components.end(), "Looking
for absent FE in tFEsTensor2Column::operator(tFinitElement&)"); #endif static
vector<tFEsTensor2>::iterator cursor = Components.begin(); cursor =
find_if(cursor,Components.end(),tFE_addr_is(&fe_)); if (cursor == Components.end()) cursor =
find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_)); #ifdef STRONGCHECK Assert(cursor <
Components.end(), "implementation error in tFEsTensor2Column::operator(tFinitElement&)"); #endif
 return *(cursor++);
// return *find_if(Components.begin(),Components.end(),tFE_addr_is(&fe_));
}*/

/*tFEsSetOfTensor2Column::tFEsSetOfTensor2Column(const tFE_model& model_,size_t reserv_):
	   Kinds(), Components(model_.HowManyFinEls())
{
 vector<tFEsSetOfTensor2>::iterator pc=Components.begin();
 for (size_t i=1; pc<Components.end(); ++i,++pc)
		  pc->Link(model_.FE(i)).Reserve(reserv_);
}

tFEsSetOfTensor2Column& tFEsSetOfTensor2Column::Include(tTensor2Kind kind_)
{
 for (vector<tFEsSetOfTensor2>::iterator pc=Components.begin(); pc<Components.end(); ++pc)
	pc->AddElement();
 Kinds.push_back(kind_);
 return *this;
}

size_t tFEsSetOfTensor2Column::ColNo (tTensor2Kind kind_) const
{
 size_t result = 1;
 for (vector<tTensor2Kind>::const_iterator pk=Kinds.begin(); pk<Kinds.end(); ++result,++pk)
   if ((*pk) == kind_) return result;
 return 0;
}

tFEsSetOfTensor2Column& tFEsSetOfTensor2Column::Transform (const
tFEsSetOfTensor2Column::fTransformer& function_)
{
 for (vector<tFEsSetOfTensor2>::iterator pc=Components.begin(); pc<Components.end(); ++pc)
	function_(*pc);
 return *this;
}

tFEsSetOfTensor2Column& tFEsSetOfTensor2Column::Swap (tFEsSetOfTensor2Column::tTensor2Kind kind_,
tFEsTensor2Column& subst_)
{
#ifdef STRONGCHECK
 Assert(Components.size() == subst_.Dimension(), "different dimensions in
tFEsSetOfTensor2Column::Swap"); #endif size_t numberOfColContainingThisKind = ColNo(kind_);
#ifdef STRONGCHECK
 Assert(numberOfColContainingThisKind > 0, "absent kind in tFEsSetOfTensor2Column::Swap");
#endif
 vector<tFEsSetOfTensor2>::iterator pc=Components.begin();
 for (size_t i=1; pc<Components.end(); ++i,++pc)
	   pc->Swap(subst_(i),numberOfColContainingThisKind);
 return *this;
}*/

tNodalTensor2SymMatrix& tNodalTensor2SymMatrix::Link(const tFE_model& model_)
{
#ifdef STRONGCHECK
	Assert(m_pNodes.empty(), "data already allocated in tNodalTensor2SymMatrix::Link");
#endif
	m_pNodes.resize(model_.HowManyNodes());
	m_components.resize(m_pNodes.size());
	vector<const tNode*>::iterator ppnode = m_pNodes.begin();
	vector<vector<Tensor2a> >::iterator p_row = m_components.begin();
	for (size_t i = 1; i <= m_pNodes.size(); ++i, ++ppnode, ++p_row)
	{
		*ppnode = &model_.GetNode(i);
		p_row->resize(i);
	}
	return *this;
}

struct fAddDOF2
{
	size_t operator()(size_t prev_sum_, const tNode* p_) const
	{
#ifdef STRONGCHECK
		Assert(p_ != nullptr, "null pointer to tNode in fAddDOF2::operator()");
#endif
		return prev_sum_ + p_->NumberOfDOFs();
	}
};

size_t tNodalTensor2SymMatrix::NumberOfDOFs() const
{
	return accumulate(m_pNodes.begin(), m_pNodes.end(), size_t(0), fAddDOF2());
	/* size_t result = 0;
	 for (vector<const tNode*>::const_iterator ppnode=pNodes.begin(); ppnode<pNodes.end(); ++ppnode)
	   {
	#ifdef STRONGCHECK
		Assert(*ppnode != nullptr, "non-allocated node in tNodalTensor2SymMatrix::NumberOfDOFs");
	#endif
		result += (*ppnode)->NumberOfDOFs();
	   }
	 return result;*/
}

struct assign0_2D
{
	void operator()(vector<Tensor2a>& o_) const
	{
		for_each(o_.begin(), o_.end(), mem_fn(&Tensor2a::Assign0));
	}
};

tNodalTensor2SymMatrix& tNodalTensor2SymMatrix::Assign0()
{
	for_each(m_components.begin(), m_components.end(), assign0_2D());
	/* size_t j;
	 vector<Tensor2s**>::iterator ppptensor = pComponents.begin();
	 for (size_t i=0; i<pComponents.size(); ++i)
	   {
		for (j=0; j<=i; ++j)
		  if ((*ppptensor)[j] != nullptr)
			{delete (*ppptensor)[j]; (*ppptensor)[j] = nullptr;}
		++ppptensor;
	   }*/
	return *this;
}

struct not_is_0
{
	bool operator()(const Tensor2a& o_) const
	{
		return !o_.is0();
	}
};

size_t tNodalTensor2SymMatrix::BandWidth(size_t i_) const
{
	--i_;
#ifdef STRONGCHECK
	Assert(
		i_ >= 0 && i_ < m_components.size(), "Invalid index in tNodalTensor2SymMatrix::BandWidth");
	Assert(
		static_cast<size_t>(m_components[i_].end() - m_components[i_].begin()) == i_ + 1,
		"Invalid bandwidth in tNodalTensor2SymMatrix::BandWidth");
	Assert(
		find_if(m_components[i_].begin(), m_components[i_].end(), not_is_0()) <
			m_components[i_].end(),
		"Invalid bandwidth in tNodalTensor2SymMatrix::BandWidth");
#endif
	return m_components[i_].end() -
		   find_if(m_components[i_].begin(), m_components[i_].end(), not_is_0()) - 1;
	/* size_t result = (--rowNo_);
	 Tensor2s** pptensor = pComponents[rowNo_];
	 for (size_t i=0; i<rowNo_; ++i)
	   {
		if (*pptensor != nullptr) break;
		++pptensor;
		--result;
	   }
	 return result;*/
}

struct fmax_band_width
{
	size_t operator()(size_t prev_max_, const vector<Tensor2a>& o_) const
	{
		return max(prev_max_, size_t(o_.end() - find_if(o_.begin(), o_.end(), not_is_0())));
	}
};

size_t tNodalTensor2SymMatrix::MaxBandWidth() const
{
	return accumulate(m_components.begin(), m_components.end(), size_t(0), fmax_band_width()) - 1;
	/*
	 size_t currentWidth, result = 0;
	 Tensor2s** pptensor;
	 vector<Tensor2s**>::const_iterator prow=pComponents.begin();
	 for (size_t i=1; i<pComponents.size(); ++i)
	   {
		++prow;
		pptensor = (*prow);
		currentWidth = i;
		for (size_t j=0; j<i; ++j)
		  {
		   if (*pptensor != nullptr) break;
		   ++pptensor;
		   --currentWidth;
		  }
		if (currentWidth > result)  result = currentWidth;
	   }
	 return result;*/
}

tNodalTensor2SymMatrix& tNodalTensor2SymMatrix::operator+=(tNodalTensor2& o_)
{
#ifdef STRONGCHECK
	Assert(!o_.is0(), "zero addend in tNodalTensor2SymMatrix::operator+= (const tNodalTensor2&)");
#endif
	// if (o_.is0())   return *this;

	size_t i = find(m_pNodes.begin(), m_pNodes.end(), &o_.Node(1)) - m_pNodes.begin(),
		   j = &o_.Node(1) == &o_.Node(2)
				   ? i
				   : find(m_pNodes.begin(), m_pNodes.end(), &o_.Node(2)) - m_pNodes.begin();
#ifdef STRONGCHECK
	Assert(
		i < m_pNodes.size() && j < m_pNodes.size(),
		"non-existing node in tNodalTensor2SymMatrix::operator+=(const tNodalTensor2&)");
#endif
	if (i < j)
		m_components[j][i].AddTransposed(Tensor2CaptureDataAdapter(o_));
	else
		m_components[i][j] += Tensor2CaptureDataAdapter(o_);
	return *this;
}

tSparseColumn::tSparseColumn(const tConstColumnVector& original_, bool copyContents_)
	: m_pComponents(original_.Dimension())
{
	vector<real_t*>::iterator ppcurrent = m_pComponents.begin();
	if (copyContents_)
		for (size_t i = 1; i <= m_pComponents.size(); ++i, ++ppcurrent)
			*ppcurrent = original_.is0(i) ? nullptr : new real_t(original_(i));
	else
		for (; ppcurrent < m_pComponents.end(); ++ppcurrent) *ppcurrent = nullptr;
}

tSparseColumn::~tSparseColumn()
{
	for_each(m_pComponents.begin(), m_pComponents.end(), fDelete_object<real_t>());
	// for (vector<real_t*>::iterator ppcurrent = pComponents.begin(); ppcurrent <
	// pComponents.end(); ++ppcurrent)
	//   /*if (*ppcurrent != nullptr)*/ delete *ppcurrent;
}

tColumnVector& tSparseColumn::Assign0()
{
	for_each(m_pComponents.begin(), m_pComponents.end(), fDelete_object<real_t>());
	// for (vector<real_t*>::iterator ppcurrent = pComponents.begin(); ppcurrent <
	// pComponents.end(); ++ppcurrent)
	//   if (*ppcurrent != nullptr) {delete *ppcurrent; *ppcurrent = nullptr;}
	return *this;
}

tColumnVector& tSparseColumn::Assign(size_t index_, real_t value_)
{
#ifdef STRONGCHECK
	Assert(index_ > 0 && index_ <= m_pComponents.size(), "invalid index in tSparseColumn::Assign");
#endif
	if (m_pComponents[--index_] == nullptr)
		m_pComponents[index_] = new real_t(value_);
	else
		(*m_pComponents[index_]) = value_;
	return *this;
}

tColumnVector& tSparseColumn::Add(size_t index_, real_t addend_)
{
#ifdef STRONGCHECK
	Assert(index_ > 0 && index_ <= m_pComponents.size(), "invalid index in tSparseColumn::Add");
#endif
	if (m_pComponents[--index_] == nullptr)
		m_pComponents[index_] = new real_t(addend_);
	else
		(*m_pComponents[index_]) += addend_;
	return *this;
}

tConstRefColumn::tConstRefColumn(
	const tNodalTensor1Column& dataSource_, const tMarkedGraph& dataSourceGraph_)
	: m_size(dataSourceGraph_.NumberOfNodes()), m_pComponents(new prtToConstReal[m_size])
{
	prtToConstReal* ppcurrent = m_pComponents;
	size_t i = 1, j;
	small_t k;
	for (; i <= m_size; ++i, ++ppcurrent)
	{
		dataSourceGraph_.SourceRowColNo(i, j, k);
		if (dataSource_.is0(j) || dataSource_(j).is0(k))
			*ppcurrent = nullptr;
		else
			*ppcurrent = &(dataSource_(j)(k));
	}
}

tRefTensor1ComponentsColumn::tRefTensor1ComponentsColumn(
	tNodalTensor1Column& dataSource_, const tMarkedGraph& dataSourceGraph_)
	: m_size(dataSourceGraph_.NumberOfNodes()), m_components(new tElement[m_size])
{
	tElement* pcurrent = m_components;
	size_t i = 1, j;
	small_t k;
	for (; i <= m_size; ++i, ++pcurrent)
	{
		dataSourceGraph_.SourceRowColNo(i, j, k);
		pcurrent->m_pTensor = &dataSource_(j);
		pcurrent->m_index = k;
	}
	dataSource_.Assign0();
}

template <typename T>
struct fDelete_Element
{
	void operator()(T& o_) const
	{
		delete o_.m_pValue;
		o_.m_pValue = nullptr;
	}
};

tColumnVector& tRefTensor1ComponentsColumn::Assign0()
{
	for_each(m_components, m_components + m_size, fDelete_Element<tElement>());
	// for (tElement* pcurrent=Components; pcurrent<Components+Size; ++pcurrent)
	//   if (!pcurrent->is0())   delete pcurrent->pValue;
	return *this;
}

tColumnVector& tRefTensor1ComponentsColumn::Assign(size_t index_, real_t value_)
{
#ifdef STRONGCHECK
	Assert(index_ > 0 && index_ <= m_size, "invalid index in tRefTensor1ComponentsColumn::Assign");
#endif
	if (m_components[--index_].is0())
		m_components[index_].m_pValue = new real_t(value_);
	else
		m_components[index_]() = value_;
	return *this;
}

tColumnVector& tRefTensor1ComponentsColumn::Add(size_t index_, real_t addend_)
{
#ifdef STRONGCHECK
	Assert(index_ > 0 && index_ <= m_size, "invalid index in tRefTensor1ComponentsColumn::Add");
#endif
	if (m_components[--index_].is0())
		m_components[index_].m_pValue = new real_t(addend_);
	else
		m_components[index_]() += addend_;
	return *this;
}

tRefTensor1ComponentsColumn& tRefTensor1ComponentsColumn::UpdateSource()
{
	for (tElement *pcurrent = m_components, *pend = m_components + m_size; pcurrent < pend;
		 ++pcurrent)
		if (!pcurrent->is0())
			pcurrent->m_pTensor->Add(pcurrent->m_index, *(pcurrent->m_pValue));
	return *this;
}

tSymConstRefMatrix::tSymConstRefMatrix(
	const tNodalTensor2SymMatrix& dataSource_, const tMarkedGraph& dataSourceGraph_)
	: m_size(dataSourceGraph_.NumberOfNodes()), m_pComponents(new vector<const real_t*>[m_size])
{
	tProcessDisplay output("Sparse matrix creating", static_cast<data_t>(m_size));
	size_t i = 1, j, m, n;
	small_t l, k;
	vector<const real_t*>::iterator ppcurrent;
	vector<const real_t*>* prow = m_pComponents;
	for (; i <= m_size; ++prow, ++i, ++output)
	{
		prow->resize(dataSourceGraph_.BandWidth(i) + 1);
		dataSourceGraph_.SourceRowColNo(i, m, k);
		ppcurrent = prow->begin();
		*ppcurrent = &dataSource_(m, m)(k, k);
		for (j = i - 1; (++ppcurrent) < prow->end(); --j)
		{
			dataSourceGraph_.SourceRowColNo(j, n, l);
			if (dataSourceGraph_.HasRib(m, k, n, l))
				*ppcurrent = &(m < n ? dataSource_(n, m)(l, k) : dataSource_(m, n)(k, l));
			else
				*ppcurrent = nullptr;
		}
	}
}

size_t tSymConstRefMatrix::BandWidth(size_t rowNo_) const
{
#ifdef STRONGCHECK
	Assert(rowNo_ > 0 && rowNo_ <= m_size, "Invalid index in tSymConstRefMatrix::BandWidth");
#endif
	vector<const real_t*>* prow = m_pComponents + (--rowNo_);
	size_t result = prow->size();
	--result;
	vector<const real_t*>::const_iterator ppcurrent = prow->end(), ppbegin = prow->begin();
	for (; --ppcurrent > ppbegin; --result)
		if (*ppcurrent != nullptr)
			break;
	return result;
}

const tSymConstRefMatrix& tSymConstRefMatrix::ProfileAndMaxBandWidth(
	size_t& profile_, size_t& max_) const
{
	size_t currentWidth;
	profile_ = 0;
	max_ = 0;
	vector<const real_t*>* prow = m_pComponents;
	vector<const real_t*>* const ppastLastRow = prow + m_size;
	vector<const real_t*>::const_iterator ppcurrent;
	for (++prow; prow < ppastLastRow; ++prow)
	{
		currentWidth = prow->size();
		--currentWidth;
		for (ppcurrent = prow->end(); --ppcurrent > prow->begin(); --currentWidth)
			if (*ppcurrent != nullptr)
				break;
		profile_ += currentWidth;
		if (currentWidth > max_)
			max_ = currentWidth;
	}
	return *this;
}

tSymConstRefMatrix& tSymConstRefMatrix::Pack0()
{
	vector<const real_t*>::iterator pcurrent;
	vector<const real_t*>* prow = m_pComponents;
	vector<const real_t*>* const ppastLastRow = prow + m_size;
	for (; prow < ppastLastRow; ++prow)
		for (pcurrent = prow->begin(); pcurrent < prow->end(); ++pcurrent)
			if (*pcurrent != nullptr && **pcurrent == 0.)
				*pcurrent = nullptr;
	return *this;
}

tSymConstRefMatrix& tSymConstRefMatrix::Swap(size_t i_, size_t j_)
{
	if (i_ == j_)
		return *this;
	size_t k;
	if (i_ < j_)
	{
		k = i_;
		i_ = j_;
		j_ = k;
	}
#ifdef STRONGCHECK
	Assert(j_ > 0 && i_ <= m_size, "Invalid index in tSymConstRefMatrix::Swap");
#endif
	--i_;
	--j_;
	// const real_t* tmp;
	vector<const real_t*>* prow = m_pComponents + m_size;
	--prow;
	for (k = m_size - i_ - 1; k > 0; --k, --prow)
		if (prow->size() > k)
		{
			if (prow->size() > k + i_ - j_)
			{
				//         tmp = prow->at(k);   prow->at(k) = prow->at(k+i_-j_);   prow->at(k+i_-j_)
				//         = tmp;
				swap(prow->operator[](k), prow->operator[](k + i_ - j_));
			}
			else if (prow->operator[](k) != nullptr)
			{
				prow->resize(k + 1 + i_ - j_, nullptr);
				prow->back() = prow->operator[](k);
				prow->operator[](k) = nullptr;
			}
		}
	// vector<vector<const real_t*> >::iterator prow_i = prow,  prow_j = prow_i;
	vector<const real_t*>*prow_i = prow, *prow_j = prow_i;
	prow_j -= (i_ - j_);

	--prow;
	if (prow_j->size() + i_ > prow_i->size() + j_)
	{
		if (prow_j->size() > 1)
			prow_i->resize(prow_j->size() + i_ - j_, nullptr);
		else
			prow_i->reserve(1 + i_ - j_);
	}
	else if (prow_j->size() + i_ < prow_i->size() + j_)
		prow_j->resize(prow_i->size() + j_ - i_, nullptr);
	vector<const real_t*>::iterator p_ik = prow_i->begin(), p_jk = prow_j->begin();
	// tmp = *p_ik;   *p_ik = *p_jk;   *p_jk = tmp;   ++p_ik; ++p_jk;
	iter_swap(p_ik, p_jk);
	++p_ik;
	++p_jk;
	for (k = i_ - j_ - 1; k > 0; --k, --prow, ++p_ik)
	{
		if (p_ik < prow_i->end())
		{
			if (k >= prow->size())
			{
				if (*p_ik != nullptr)
				{
					prow->resize(k + 1, nullptr);
					prow->back() = *p_ik;
					*p_ik = nullptr;
				}
			}
			else if (*p_ik != prow->operator[](k))	// at(k))
			{
				//               tmp = *p_ik;   *p_ik = prow->at(k);   prow->at(k) = tmp;
				swap(*p_ik, prow->operator[](k));
			}
		}
		else if (prow->operator[](k) != nullptr)
		{
			prow_i->resize(i_ + 1 - j_ - k, nullptr);
			p_ik = prow_i->end();
			--p_ik;
			(*p_ik) = prow->operator[](k);
			prow->operator[](k) = nullptr;
		}
	}
	if (p_ik < prow_i->end())
		++p_ik;
	for (; p_ik < prow_i->end(); ++p_jk, ++p_ik)
	{
		//    tmp  = *p_ik;    *p_ik = *p_jk;    *p_jk = tmp;
		iter_swap(p_ik, p_jk);
	}
	return *this;
}

tColumnVector& tSymConstRefMatrix::SolveLinAlgSystem(
	const tConstColumnVector& rightSide_, tColumnVector& solution_) const
{
	tTriangularMatrix lowerTriangular(*this);
	tSparseColumn auxilSolution(rightSide_, false);
	lowerTriangular.SolveLinAlgSystem(rightSide_, auxilSolution);
	lowerTriangular.Transpose().SolveLinAlgSystem(auxilSolution, solution_);
	return solution_;
}

tTriangularMatrix::tTriangularMatrix(const tSymConstRefMatrix& original_)
	: m_size(original_.Dimension()), m_components(new vector<real_t>[m_size]), m_isLower(true)
{
	tProcessDisplay output("Triangular matrix creating", static_cast<data_t>(m_size));
	size_t j;
	vector<real_t>::iterator p_ij, p_ik, p_jk;
	vector<real_t>*prow_i = m_components, *prow_j;
	for (size_t i = 1; i <= m_size; ++i, ++prow_i, ++output)
	{
#ifdef STRONGCHECK
		Assert(
			original_.DiagonalComponent(i) > 0.,
			"L*Lt expansion cannot be used for non-positive definite matrices");
#endif
		prow_i->resize(original_.BandWidth(i) + 1);
		p_ij = prow_i->begin();
		for (j = i - original_.BandWidth(i); j <= i; ++j, ++p_ij)
		{
			*p_ij = original_.is0(i, j) ? 0. : original_(i, j);
			prow_j = &m_components[j - 1];
			p_jk = prow_j->end();
			--p_jk;
			p_ik = p_ij;
			while (p_ik > prow_i->begin() && p_jk > prow_j->begin())
				*p_ij -= (*(--p_ik)) * (*(--p_jk));
			if (i == j)
				*p_ij = sqrt(*p_ij);
			else
				*p_ij /= prow_j->back();
		}
	}
}

tColumnVector& tTriangularMatrix::SolveLowerSystem(
	const tConstColumnVector& rightSide_, tColumnVector& solution_) const
{
	tProcessDisplay output("Lower system solving", static_cast<data_t>(m_size));
#ifdef STRONGCHECK
	Assert(
		rightSide_.Dimension() == Dimension() && solution_.Dimension() == Dimension(),
		"unequal dimensions in tTriangularMatrix::SolveLowerSystem");
#endif
	solution_.Assign0();
	size_t k;
	vector<real_t>::const_iterator p_ik;
	const vector<real_t>* prow_i = m_components;
	for (size_t i = 1; i <= m_size; ++i, ++prow_i, ++output)
	{
		if (!rightSide_.is0(i))
			solution_.Assign(i, rightSide_(i));
		p_ik = prow_i->end();
		--p_ik;
		k = i;
		while (p_ik > prow_i->begin())
		{
			--p_ik;
			--k;
			if (!solution_.is0(k))
				solution_.Add(i, -(*p_ik) * solution_(k));
		}
#ifdef STRONGCHECK
		Assert(
			mathdef::is_not_zero(prow_i->back()),
			"principal element =0 in tTriangularMatrix::SolveLowerSystem");
#endif
		if (!solution_.is0(i))
			solution_(i) /= prow_i->back();
	}
	return solution_;
}

tColumnVector& tTriangularMatrix::SolveUpperSystem(
	const tConstColumnVector& rightSide_, tColumnVector& solution_) const
{
	tProcessDisplay output("Upper system solving", static_cast<data_t>(m_size));
#ifdef STRONGCHECK
	Assert(
		rightSide_.Dimension() == Dimension() && solution_.Dimension() == Dimension(),
		"unequal dimensions in tTriangularMatrix::SolveUpperSystem");
#endif
	tSparseColumn auxRightSide(rightSide_, true);
	solution_.Assign0();
	size_t i;
	const real_t* psol_k;
	vector<real_t>::const_iterator p_ik;
	const vector<real_t>* pcol_k = m_components + m_size;
	for (size_t k = m_size; k > 0; --k, ++output)
	{
		--pcol_k;
#ifdef STRONGCHECK
		Assert(
			mathdef::is_not_zero(pcol_k->back()),
			"principal element =0 in tTriangularMatrix::SolveLowerSystem");
#endif
		if (!auxRightSide.is0(k))
			solution_.Assign(k, auxRightSide(k) / pcol_k->back());
		if (!solution_.is0(k))
		{
			psol_k = &solution_(k);
			p_ik = pcol_k->end();
			--p_ik;
			i = k;
			while (p_ik > pcol_k->begin())
			{
				--p_ik;
				--i;
				auxRightSide.Add(i, -(*p_ik) * (*psol_k));
			}
		}
	}
	return solution_;
}
