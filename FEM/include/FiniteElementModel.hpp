#pragma once

#include "NumTypes.hpp"
#include "Node.hpp"
#include "Material.hpp"
#include "VectorMatrix.hpp"
#include "femDllDefinitions.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

class tFinitElement;
class tNodalTensor1Column;

class FEM_API tFE_model
{
public:
	~tFE_model();
	const std::string& Name() const
	{
		return m_caption;
	}
	tFE_model& SetName(const char* newname_)
	{
		m_caption = newname_;
		return *this;
	}
	tFE_model& SetName(const std::string& newname_)
	{
		m_caption = newname_;
		return *this;
	}
	tFinitElement& AddFE(const std::string&, size_t);
	tMaterial& AddMaterial(const std::string&, const std::string&);
	tNode& GetNode(size_t);
	const tNode& GetNode(size_t) const;
	tFinitElement& GetFE(size_t);
	const tFinitElement& GetFE(size_t) const;
	std::vector<tFinitElement*>::const_iterator pFE_begin() const
	{
		return m_pFinEls.begin();
	}
	std::vector<tFinitElement*>::const_iterator pFE_end() const
	{
		return m_pFinEls.end();
	}
	std::vector<tFinitElement*>::iterator pFE_begin()
	{
		return m_pFinEls.begin();
	}
	std::vector<tFinitElement*>::iterator pFE_end()
	{
		return m_pFinEls.end();
	}
	size_t GetNodeNum(const tNode&) const;
	size_t GetFENum(const tFinitElement&) const;
	tMaterial& Material(const std::string&);
	const tMaterial& Material(const std::string&) const;
	const tMaterial& GetMaterial(size_t) const;
	const std::string& Name(const tMaterial&) const;
	bool IsEmpty() const
	{
		return m_nodes.empty() && m_pFinEls.empty() && m_pMaterials.empty();
	}
	bool Contains(const tNode& node_) const
	{
		return std::find_if(m_nodes.begin(), m_nodes.end(), address_is<tNode>(&node_)) <
			   m_nodes.end();
	}
	bool Contains(const tFinitElement& FE_) const
	{
		return std::find(m_pFinEls.begin(), m_pFinEls.end(), &FE_) < m_pFinEls.end();
	}
	bool Contains(const tMaterial& matl_) const
	{
		return std::find(m_pMaterials.begin(), m_pMaterials.end(), &matl_) < m_pMaterials.end();
	}
	const tFE_model& LinkDisplacements(const tNodalTensor1Column*) const;
	size_t HowManyNodes() const
	{
		return m_nodes.size();
	}
	size_t HowManyFinEls() const
	{
		return m_pFinEls.size();
	}
	size_t HowManyMaterials() const
	{
		return m_pMaterials.size();
	}
	size_t MaxNodesPerFE() const;
	size_t MaxFinElsPerNode() const;
	size_t MaxFinElsPerMatl() const;
	real_t Volume() const;
	//   real_t Mass() const;
	tFE_model& ImportNodesFromFile(ifstream_XML&);
	tFE_model& ImportMatlsFromFile(ifstream_XML&);
	tFE_model& ImportFinElsFromFile(ifstream_XML&);
	tFE_model& ImportFromFile(const char* fileName_);
	void Clear();
	void Dump(std::ostream& out) const;

private:
#pragma warning(push)
#pragma warning(disable : 4251)
	std::string m_caption;
	std::vector<tNode> m_nodes;
	std::vector<tFinitElement*> m_pFinEls;
	std::vector<tMaterial*> m_pMaterials;
#pragma warning(pop)
};

inline tNode& tFE_model::GetNode(size_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_nodes.size(), tMessage("Invalid node # = ") << i_ << " in FE-model");
#endif
	return m_nodes[--i_];
}

inline const tNode& tFE_model::GetNode(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_nodes.size(), tMessage("Invalid node # = ") << i_ << " in FE-model");
#endif
	return m_nodes[--i_];
}

inline tFinitElement& tFE_model::GetFE(size_t i_)
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_pFinEls.size(), tMessage("Invalid FE # = ") << i_ << " in FE-model");
#endif
	return *m_pFinEls[--i_];
}

inline const tFinitElement& tFE_model::GetFE(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_pFinEls.size(), tMessage("Invalid FE # = ") << i_ << " in FE-model");
#endif
	return *m_pFinEls[--i_];
}

inline const tMaterial& tFE_model::GetMaterial(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(
		i_ > 0 && i_ <= m_pMaterials.size(),
		"invalid material number in tFE_model::Material(i)const");
#endif
	return *m_pMaterials[--i_];
}

inline size_t tFE_model::GetNodeNum(const tNode& node_) const
{
#ifdef STRONGCHECK
	Assert(&node_ != nullptr, "Attempt to search # of absent node");
	Assert(
		std::find_if(m_nodes.begin(), m_nodes.end(), address_is<tNode>(&node_)) < m_nodes.end(),
		"Node is absent in FE-model");
#endif
	// return &node_ - Nodes.begin() + 1;
	return std::find_if(m_nodes.begin(), m_nodes.end(), address_is<tNode>(&node_)) -
		   m_nodes.begin() + 1;
}

inline size_t tFE_model::GetFENum(const tFinitElement& fe_) const
{
#ifdef STRONGCHECK
	Assert(&fe_ != nullptr, "No # of absent FE");
	Assert(
		std::find(m_pFinEls.begin(), m_pFinEls.end(), &fe_) < m_pFinEls.end(),
		"FE is absent in FE-model");
#endif
	return std::find(m_pFinEls.begin(), m_pFinEls.end(), &fe_) - m_pFinEls.begin() + 1;
}

inline const std::string& tFE_model::Name(const tMaterial& matl_) const
{
#ifdef STRONGCHECK
	Assert(&matl_ != nullptr, "No # of absent material");
	Assert(
		std::find(m_pMaterials.begin(), m_pMaterials.end(), &matl_) < m_pMaterials.end(),
		"Material is absent in FE-model");
#endif
	return (*std::find(m_pMaterials.begin(), m_pMaterials.end(), &matl_))->GetName();
}

inline const tFE_model& tFE_model::LinkDisplacements(const tNodalTensor1Column* pdispls_) const
{
#ifdef STRONGCHECK
// Assert(pdispls_!=nullptr, "non-existing displacement tFE_model::LinkDisplacements");
#endif
	if (pdispls_ != nullptr)
		pdispls_->LinkAsDisplacements();
	return *this;
}
