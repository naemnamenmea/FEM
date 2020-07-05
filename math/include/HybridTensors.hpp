#pragma once

#include "Tensors.hpp"
#include "Node.hpp"

#include <vector>
#include <string>
//#include "Duo.h"

// class tNode;
class tFinitElement;
class tFE_model;

class tNodalTensor1 : public Tensor1a
{
public:
	tNodalTensor1() : Tensor1a(), m_pNode(nullptr)
	{
	}
	~tNodalTensor1()
	{
	}
	tNodalTensor1(const tNodalTensor1& o_) : Tensor1a(o_), m_pNode(o_.m_pNode)
	{
	}
	Tensor1a& operator=(const Tensor1a& o_)
	{
		return Tensor1a::operator=(o_);
	}
	tNodalTensor1& operator=(const tNodalTensor1& o_)
	{
		Tensor1a::operator=(o_);
		m_pNode = o_.m_pNode;
		return *this;
	}
	tNodalTensor1& Link(const tNode& node_)
	{
		m_pNode = &node_;
		return *this;
	}
	const tNodalTensor1& LinkAsDispl() const;
	const tNode& Node() const;

protected:
	const tNode* m_pNode;
};

class tNodalTensor2 : public Tensor2a
{
public:
	tNodalTensor2() : Tensor2a(), m_pNode1(nullptr), m_pNode2(nullptr)
	{
	}
	tNodalTensor2(const tNodalTensor2& o_)
		: Tensor2a(o_), m_pNode1(&o_.Node(1)), m_pNode2(&o_.Node(2))
	{
	}
	tNodalTensor2(const tNodalTensor2& o_, bool)
		: Tensor2a(), m_pNode1(&o_.Node(1)), m_pNode2(&o_.Node(2))
	{
	}
	~tNodalTensor2()
	{
	}
	tNodalTensor2& operator=(const Tensor2s& o_)
	{
		Tensor2a::operator=(o_);
		return *this;
	}
	tNodalTensor2& operator=(const Tensor2a& o_)
	{
		Tensor2a::operator=(o_);
		return *this;
	}
	tNodalTensor2& operator=(const tNodalTensor2& o_)
	{
		Tensor2a::operator=(o_);
		m_pNode1 = o_.m_pNode1;
		m_pNode2 = o_.m_pNode2;
		return *this;
	}
	tNodalTensor2& Link(const tNode& n1_, const tNode& n2_)
	{
		m_pNode1 = &n1_;
		m_pNode2 = &n2_;
		return *this;
	}
	const tNode& Node(small_t) const;

private:
	const tNode* m_pNode1;
	const tNode* m_pNode2;
};

// class tFEsTensor2r;
/*class tFEsTensor2: public Tensor2a
{
 friend class tFEsSetOfTensor2;
 protected:
   const tFinitElement* pFE;
//   tFEsTensor2(Tensor2s* addr_, const tFinitElement& fe_): Tensor2a(addr_), pFE(&fe_) {}
 public:
   tFEsTensor2(): Tensor2a(), pFE(nullptr) {}
   tFEsTensor2(const Tensor2a& o_): Tensor2a(o_), pFE(nullptr) {}
   tFEsTensor2(const Tensor2a& o_, const tFinitElement& fe_): Tensor2a(o_), pFE(&fe_) {}
   tFEsTensor2(const tFEsTensor2& o_): Tensor2a(o_), pFE(o_.pFE) {}
//   tFEsTensor2(const tFEsTensor2& o_, bool): Tensor2a(), pFE(o_.pFE) {}
//   tFEsTensor2(tFEsTensor2r&);
   ~tFEsTensor2() {}
   tFEsTensor2& Link (const tFinitElement& fe_) {pFE=&fe_; return *this;}
   const tFinitElement& FE() const;
};*/

/* class tFEsTensor2r: public tFEsTensor2
   {
	friend class tFEsTensor2;
	private:
	  typedef Tensor2s* tpTensor2s;
//      tpTensor2s& rpData;
	public:
//      tFEsTensor2r(tpTensor2s& rpdata_, const tFinitElement& fe_): tFEsTensor2(rpdata_,fe_),
rpData(rpdata_) {rpData=nullptr;}
//      tFEsTensor2r(tFEsTensor2r& o_): tFEsTensor2(o_.pData,o_.FE()), rpData(o_.rpData)
{o_.pData=nullptr;}
//      ~tFEsTensor2r() {rpData = pData; pData=nullptr;}
   };

inline tFEsTensor2::tFEsTensor2(tFEsTensor2r& o_): Tensor2a(), pFE(o_.pFE) {pData = o_.pData;
o_.pData = nullptr;}
*/

// class tFEsSetOfTensor2
//{
// friend class tFEsSetOfTensor2;
/* class tFEsTensor2r: public tFEsTensor2
   {
	friend class tFEsTensor2;
	private:
	  typedef Tensor2s* tpTensor2s;
	  tpTensor2s& rpData;
	public:
	  tFEsTensor2r(tpTensor2s& rpdata_, const tFinitElement& fe_): tFEsTensor2(rpdata_,fe_),
   rpData(rpdata_) {rpData=nullptr;} tFEsTensor2r(tFEsTensor2r& o_): tFEsTensor2(o_.pData,o_.FE()),
   rpData(o_.rpData) {o_.pData=nullptr;} ~tFEsTensor2r() {rpData = pData; pData=nullptr;}
   };*/
/* private:
   const tFinitElement* pFE;
   std::vector<Tensor2s*> pData;
 public:
   tFEsSetOfTensor2(): pFE(nullptr), pData() {}
   tFEsSetOfTensor2(const tFinitElement& fe_, size_t size_): pFE(&fe_), pData(size_,nullptr) {}
   ~tFEsSetOfTensor2() {Assign0();}
   tFEsSetOfTensor2& Link(const tFinitElement& fe_) {pFE=&fe_; return *this;}
   tFEsSetOfTensor2& Assign0();
   bool isEmpty() const {return pData.empty();}
   const tFinitElement& FE() const  {return *pFE;}
   tFEsSetOfTensor2& Reserve(size_t howmany_) {pData.reserve(howmany_); return *this;}
   tFEsSetOfTensor2& AddElement(){pData.push_back(nullptr); return *this;}
//   tFEsTensor2r operator()(size_t);
   tFEsSetOfTensor2& Swap(tFEsSetOfTensor2&);
   tFEsSetOfTensor2& operator+=(const tFEsSetOfTensor2& add_);
   tFEsSetOfTensor2& operator*=(real_t);
   tFEsSetOfTensor2& operator/=(real_t);
   tFEsSetOfTensor2& Swap(tFEsTensor2&,size_t);
};*/

inline const tNode& tNodalTensor1::Node() const
{
#ifdef STRONGCHECK
	Assert(m_pNode != nullptr, "attempt to return non-existing node in tNodalTensor1::Node");
#endif
	return *m_pNode;
}

inline const tNodalTensor1& tNodalTensor1::LinkAsDispl() const
{
#ifdef STRONGCHECK
	Assert(m_pNode != nullptr, "non-existing node in tNodalTensor1::LinkAsDispl");
#endif
	m_pNode->LinkDisplacement(m_pData);
	return *this;
}

inline const tNode& tNodalTensor2::Node(const small_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ == 1 || i_ == 2, "invalid index in tNodalTensor2::Node");
	Assert(
		m_pNode1 != nullptr && m_pNode2 != nullptr,
		"attempt to return non-existing node in tNodalTensor2::Node");
#endif
	return *(i_ == 1 ? m_pNode1 : m_pNode2);
}

/*inline const tFinitElement& tFEsTensor2::FE() const
{
#ifdef STRONGCHECK
 Assert(pFE != nullptr, "attempt to return non-existing FE in tFEsTensor2::FE");
#endif
 return *pFE;
}*/

// inline /*tFEsSetOfTensor2::*/tFEsTensor2r tFEsSetOfTensor2::operator()(size_t no_)
/*{
#ifdef STRONGCHECK
 Assert(no_>0 && no_<=pData.size(), "different dimensions in tFEsSetOfTensor2::operator+=");
#endif
 return tFEsTensor2r(pData[--no_],*pFE);
} */

/*inline tFEsSetOfTensor2& tFEsSetOfTensor2::Swap (tFEsSetOfTensor2& subst_)
{
#ifdef STRONGCHECK
 Assert(pFE == subst_.pFE, "unequal FEs in tFEsSetOfTensor2::Swap");
#endif
 pData.swap(subst_.pData);
 return *this;
}

inline tFEsSetOfTensor2& tFEsSetOfTensor2::Swap (tFEsTensor2& subst_, size_t no_)
{
#ifdef STRONGCHECK
 Assert(no_>0 && no_<=pData.size(), "invalid number in tFEsSetOfTensor2::Swap");
 Assert(pFE == subst_.pFE, "unequal FEs in tFEsSetOfTensor2::Swap");
#endif//def STRONGCHECK
 std::swap(pData[--no_],subst_.pData);
 return *this;
}*/
