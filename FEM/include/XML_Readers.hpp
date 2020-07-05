#pragma once

#include "FileInOut.hpp"
#include "Material.hpp"
#include "FiniteElementPool.hpp"
#include "Node.hpp"
#include "Load.hpp"
#include "Analysis.hpp"
#include "FiniteElementModel.hpp"

#include <sstream>
#include <string>

class tNode_XML_Reader : public tXML_Reader
{
public:
	void ReadData(const std::string&, ifstream_XML&) override;
	void SetValue(const std::string& attr_name_, const std::string& value_) override;
	size_t NodeNum() const
	{
		return m_num;
	}
	void SetValsTo(tNode&) const;

private:
	size_t m_num;
	real_t m_x, m_y, m_z;
	bool m_readCoord;
	bool m_fixedX, m_fixedY, m_fixedZ;
};

class tMaterial_XML_Reader : public tXML_Reader
{
public:
	void ReadData(const std::string&, ifstream_XML&) override;
	void SetValue(const std::string& attr_name_, const std::string& value_) override;
	tMaterial* CreateMaterial();

private:
	std::string m_type;
	std::string m_name;
	std::string m_data;
};

class tFinEl_XML_Reader : public tXML_Reader
{
public:
	tFinEl_XML_Reader(const tFE_model& model_)
		: m_model(model_), m_num(), m_planeOr1dStressState(), m_crossValue()
	{
	}
	void ReadData(const std::string& tag_, ifstream_XML& xmlfile_) override;
	void SetValue(const std::string& attr_name_, const std::string& value_) override;
	tFinitElement* CreateFE();
	size_t FE_Num() const
	{
		return m_num;
	}

private:
	const tFE_model& m_model;
	size_t m_num;
	std::string m_type;
	std::string m_materialName;
	std::string m_nodesNums;
	real_t m_crossValue;
	bool m_planeOr1dStressState;
};

class tLoadCase_XML_Reader : public tXML_Reader
{
public:
	tLoadCase_XML_Reader() : m_nodeNum(0)
	{
	}

	void SetValue(const std::string& attr_name_, const std::string& value_) override;
	void SetLoadCaseName(tLoadCase& loadcase_)
	{
		loadcase_.SetName(m_loadCaseName);
	}
	void CreateForce(tLoadCase&, const tFE_model&);

private:
	std::string m_loadCaseName;
	std::string m_forceType;
	std::string m_components;
	size_t m_nodeNum;
};

class tAnalysis_XML_Reader : public tXML_Reader
{
public:
	tAnalysis_XML_Reader(const tFE_model& model_, const tLoads& loads_)
		: m_model(model_), m_loads(loads_), m_evalPerFEMethod(), m_pCase()
	{
	}
	void ReadData(const std::string&, ifstream_XML&) override;
	void SetValue(const std::string& attr_name_, const std::string& value_) override;
	tAnalysis* CreateAnalysis();

private:
	const tFE_model& m_model;
	const tLoads& m_loads;
	const tLoadCase* m_pCase;
	std::string m_type;
	std::string m_name;
	std::string m_oFileName;
	tAnalysis::tEvalPerFEMethod m_evalPerFEMethod;
};

inline tMaterial* tMaterial_XML_Reader::CreateMaterial()
{
	tMaterial* pMaterial = tMaterial::CreateMaterial(m_type);
	pMaterial->SetName(m_name);
	std::istringstream data_in(m_data);
	pMaterial->ReadData(data_in);
	m_data.clear();
	return pMaterial;
}

inline void tNode_XML_Reader::SetValsTo(tNode& node_) const
{
	node_.AssignCoord(X, m_x).AssignCoord(Y, m_y).AssignCoord(Z, m_z);
	if (m_fixedX)
		node_.Fix(X);
	else
		node_.Free(X);
	if (m_fixedY)
		node_.Fix(Y);
	else
		node_.Free(Y);
	if (m_fixedZ)
		node_.Fix(Z);
	else
		node_.Free(Z);
}

inline void tLoadCase_XML_Reader::CreateForce(tLoadCase& loadcase_, const tFE_model& modelToLoad_)
{
	std::istringstream data_in(m_components);
	loadcase_.AddNodalForce(m_forceType, modelToLoad_.GetNode(m_nodeNum)).ReadData(data_in);
	m_components.clear();
}

inline tAnalysis* tAnalysis_XML_Reader::CreateAnalysis()
{
	tAnalysis* panalysis = tAnalysis::NewAnalysis(m_type);
	panalysis->SetName(m_name.c_str());
	panalysis->EvalPerFEMethod = m_evalPerFEMethod;
	panalysis->SetMinimizeBandMethod(tMarkedGraph::reverse_Cuthill_McKee);	// to change!!
	panalysis->Link(m_model);
	panalysis->Link(*m_pCase);
	if (m_oFileName.empty())
		(((((m_oFileName = m_model.Name()) += ".") += m_pCase->Name()) += ".") += m_name) +=
			".result.xml";
	panalysis->SetOutputFileName(m_oFileName.c_str());
	return panalysis;
}
