#pragma once

#include "NumTypes.hpp"	 // size_t, real_t
#include "femDllDefinitions.hpp"

#include <vector>  // - to improve effectiveness must be replaced by <deque> !!!
#include <set>
#include <string>

class tNodalTensor2SymMatrix;

class FEM_API tMarkedGraph
{
public:
	enum tMarkMethod
	{
		no_permutations,
		reverse_Cuthill_McKee,
		minimal_degree,
		minimal_degree2
	};

	tMarkedGraph()
	{
	}
	tMarkedGraph(const tNodalTensor2SymMatrix& o_)
	{
		Link(o_);
	}
	tMarkedGraph(const tMarkedGraph& o_) : m_nodes(o_.m_nodes)
	{
		operator=(o_);
	}
	tMarkedGraph(const tMarkedGraph&, const std::string&);

	bool isEmpty() const
	{
		return m_nodes.empty();
	}
	tMarkedGraph& Link(const tNodalTensor2SymMatrix&);
	const tMarkedGraph& SourceRowColNo(size_t, size_t&, small_t&) const;
	size_t BandWidth(size_t) const;
	bool HasRib(size_t, size_t, size_t, size_t) const;
	size_t NumberOfNodes() const
	{
		return m_nodes.size();
	}
	tMarkedGraph& ReMark(tMarkMethod method_)
	{
		return method_ == reverse_Cuthill_McKee
				   ? RCM_Algorithm()
				   : (method_ == minimal_degree
						  ? MinDegreeAlgorithm()
						  : (method_ == minimal_degree2 ? MinDegreeAlgorithm(false) : *this));
	}
	tMarkedGraph& Prolog_ReMark(unsigned int);
	tMarkedGraph& operator=(const tMarkedGraph&);

private:
#pragma warning(push)
#pragma warning(disable : 4251)
	class tNode
	{
	public:
		static bool fLessDegree(const tNode* pnode1_, const tNode* pnode2_)
		{
			return pnode1_->Degree() < pnode2_->Degree();
		}

		tNode(size_t i_, small_t j_)
			: m_rowColNo(i_), m_componentNo(j_), m_pLinkedNodes(), m_previousLocation(nullptr)
		{
			m_previousLocation = this;
		}
		tNode() : m_rowColNo(0), m_componentNo(0), m_pLinkedNodes(), m_previousLocation(nullptr)
		{
			m_previousLocation = this;
		}
		tNode(const tNode& orig_)
			: m_rowColNo(orig_.SourceRowColNo()),
			  m_componentNo(orig_.SourceComponentNo()),
			  m_pLinkedNodes(orig_.m_pLinkedNodes),
			  m_previousLocation(&orig_)
		{
		}
		tNode& operator=(const tNode& orig_)
		{
			m_rowColNo = orig_.SourceRowColNo();
			m_componentNo = orig_.SourceComponentNo();
			m_pLinkedNodes = orig_.m_pLinkedNodes;
			m_previousLocation = &orig_;
			return *this;
		}
		bool operator==(const tNode& o_) const
		{
			return m_rowColNo == o_.SourceRowColNo() && m_componentNo == o_.SourceComponentNo();
		}
		bool operator!=(const tNode& o_) const
		{
			return !operator==(o_);
		}
		bool operator<(const tNode& o_) const
		{
			return m_rowColNo < o_.SourceRowColNo() ? true : m_componentNo < o_.SourceComponentNo();
		}
		tNode& Link(size_t i_, small_t j_)
		{
			m_rowColNo = i_;
			m_componentNo = j_;
			return *this;
		}
#if defined STRONGCHECK && (!defined __BORLANDC__ || __BORLANDC__ >= 0x550)
		bool Link(const tNode* pnodeToLink_) const
		{
			return pnodeToLink_ == this ? false : m_pLinkedNodes.insert(pnodeToLink_).second;
		}
#else	// ifndef STRONGCHECK
		bool Link(const tNode* pnodeToLink_) const
		{
			return m_pLinkedNodes.insert(pnodeToLink_).second;
		}
#endif	// def STRONGCHECK
		bool HasLink(const tNode* pnodeToCheck_) const
		{
			return m_pLinkedNodes.find(pnodeToCheck_) != m_pLinkedNodes.end();
		}
		bool UnLink(const tNode* pnodeToUnLink_) const
		{
			return m_pLinkedNodes.erase(pnodeToUnLink_) == 1;
		}
		const tNode& UnLinkAll() const;
		const tNode& LinkAdjToEachOther() const;
		const tNode& PointAbsentAdj(
			std::vector<const tNode*>&, std::vector<const tNode*>::iterator&, bool) const;
		bool UpdateLinks() const;
		bool UpdateLink(const tNode&) const;
		const size_t& SourceRowColNo() const
		{
			return m_rowColNo;
		}
		const small_t& SourceComponentNo() const
		{
			return m_componentNo;
		}
		size_t Degree() const
		{
			return m_pLinkedNodes.size();
		}

	private:
		friend class tMarkedGraph;
		// friend class tNode;

		size_t m_rowColNo;		// coord of diagonal component in symmetric matrix of Tensor2
		small_t m_componentNo;	// coord of diagonal Tensor2 component in diagonal component in
								// symmetric matrix of Tensor2
		mutable std::set<const tNode*, std::less<const tNode*> > m_pLinkedNodes;
		mutable const tNode* m_previousLocation;
	};

	const tNode& FindPseudoPeriferalNode(const tNode&) const;
	const tNode& FindEccentricityAndOppositeNode(const tNode&, size_t&) const;
	tMarkedGraph& SwapNodes(const tNode&, const tNode&);
	tMarkedGraph& SwapNodes(size_t, size_t);
	tMarkedGraph& ExcludeNode(size_t, bool = true);
	size_t MinDegNodeNo() const;
	tMarkedGraph& RCM_Algorithm();
	tMarkedGraph& MinDegreeAlgorithm(bool = true);

	std::vector<tNode> m_nodes;
#pragma warning(pop)
};

inline bool tMarkedGraph::tNode::UpdateLink(const tMarkedGraph::tNode& nodeToUpdate_) const
{
#ifdef STRONGCHECK
	if (nodeToUpdate_.UnLink(m_previousLocation))
	{
		Assert(
			nodeToUpdate_.Link(this),
			"error with addition of this to set in tMarkedGraph::tNode::UpdateLink");
		return true;
	}
	return false;
#else
	return nodeToUpdate_.UnLink(m_previousLocation) ? nodeToUpdate_.Link(this) : false;
#endif
}

inline const tMarkedGraph& tMarkedGraph::SourceRowColNo(
	size_t nodeNo_, size_t& rowColNo_, small_t& componentNo_) const
{
#ifdef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_nodes.size(), "invalid index in tMarkedGraph::RowColNo");
#endif
	rowColNo_ = m_nodes[--nodeNo_].SourceRowColNo();
	componentNo_ = m_nodes[nodeNo_].SourceComponentNo();
	return *this;
}
