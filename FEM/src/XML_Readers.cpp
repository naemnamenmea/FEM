
#include "XML_Readers.hpp"
#include "Analysis.hpp"

using namespace std;

void tNode_XML_Reader::ReadData(const string& tag_, ifstream_XML& xmlfile_)
{
	if (tag_ == "node")
	{
		streampos beg_pos = xmlfile_.tellg();
		xmlfile_.FindChar('<');
		string first_subtag;
		xmlfile_.ReadWord(first_subtag);
		xmlfile_.seekg(beg_pos);
		ReadData(first_subtag, xmlfile_);
	}
	else
	{
		m_readCoord = (tag_ == "coord");
		xmlfile_.ReadObject(*this, tag_);
	}
}

void tNode_XML_Reader::SetValue(const string& attr_name_, const string& value_)
{
	if (attr_name_ == "id")
		istringstream(value_) >> m_num;
	else if (m_readCoord)
	{
		istringstream val_in(value_);
		if (attr_name_ == "x")
			val_in >> m_x;
		else if (attr_name_ == "y")
			val_in >> m_y;
		else if (attr_name_ == "z")
			val_in >> m_z;
	}
	else
	{
		if (attr_name_ == "x")
			m_fixedX = value_ == "fixed";
		else if (attr_name_ == "y")
			m_fixedY = value_ == "fixed";
		else if (attr_name_ == "z")
			m_fixedZ = value_ == "fixed";
	}
}

void tMaterial_XML_Reader::ReadData(const string& tag_, ifstream_XML& xmlfile_)
{
	if (tag_ == "material")
		xmlfile_.ReadObject(*this, "SigmaEpsilonDiagram");
	else if (tag_ == "SigmaEpsilonDiagram")
	{
		//      ifstream_XML::pos_type pos = xmlfile_.tellg();
		const size_t points_num = xmlfile_.WordCountBefore("<point", "/SigmaEpsilonDiagram");
		//      xmlfile_.seekg(pos);
		xmlfile_.ReadObject(*this, "function");
		for (size_t i = 0; i < points_num; ++i) xmlfile_.ReadObject(*this, "point");
	}
}

void tMaterial_XML_Reader::SetValue(const string& attr_name_, const string& value_)
{
	if (attr_name_ == "type")
		m_type = value_;
	else if (attr_name_ == "id" || attr_name_ == "name")
		m_name = value_;
	else
		(((m_data += '\t') += attr_name_) += " = ") +=
			value_;	 // " = " but not '=' because of BCB6's gluck in corresponding ReadData
}

void tFinEl_XML_Reader::ReadData(const string& tag_, ifstream_XML& xmlfile_)
{
	if (tag_ == "fe")
	{
		if (m_nodesNums.empty())
			xmlfile_.ReadObject(*this, "");
	}
	else if (tag_ == "nodes")
	{
		GETLINE(xmlfile_, m_nodesNums, '<');
		xmlfile_.unget();
	}
}

void tFinEl_XML_Reader::SetValue(const string& attr_name_, const string& value_)
{
	if (attr_name_ == "id")
		istringstream(value_) >> m_num;
	else if (attr_name_ == "type")
		m_type = value_;
	else if (attr_name_ == "thickness" || attr_name_ == "cross-sec_area")
		istringstream(value_) >> m_crossValue;
	else if (attr_name_ == "state")
		m_planeOr1dStressState = value_.substr(value_.length() - 6, 6) == "stress";
	else if (attr_name_ == "material")
		m_materialName = value_;
	else if (attr_name_ == "nodes")
		m_nodesNums = value_;
}

tFinitElement* tFinEl_XML_Reader::CreateFE()
{
	tFinitElement* pFE = tFinitElement::NewFE(m_type);
	switch (pFE->SpaceDimension())
	{
		case 1:
			dynamic_cast<t1D_FE*>(pFE)
				->SetCrosSecArea(m_crossValue)
				.Set1dStress(m_planeOr1dStressState);
			break;
		case 2:
			dynamic_cast<t2D_FE*>(pFE)
				->SetThickness(m_crossValue)
				.Set2dStress(m_planeOr1dStressState);
			break;
		default:
			THROW_MESSAGE("Unsupported space dimension");
	}
	pFE->Link(m_model.Material(m_materialName));
	istringstream nums_in(m_nodesNums);
	for (size_t i = 0, num; i < pFE->HowManyNodes(); ++i)
	{
		nums_in >> num;
		pFE->DefineNextNode(m_model.GetNode(num), true, true);
	}
	m_nodesNums.clear();
	return pFE;
}

void tLoadCase_XML_Reader::SetValue(const string& attr_name_, const string& value_)
{
	if (attr_name_ == "id" || attr_name_ == "name")
		m_loadCaseName = value_;
	else if (attr_name_ == "type")
		m_forceType = value_;
	else if (attr_name_ == "node")
		istringstream(value_) >> m_nodeNum;
	else if (attr_name_ == "x" || attr_name_ == "y" || attr_name_ == "z")
		(((m_components += '\t') += attr_name_) += " = ") +=
			value_;	 // " = " but not '=' because of BCB6's gluck in corresponding ReadData

	/* else if (attr_name_ == "x")
			istringstream(value_) >> X;
	//      ComponentsStr.insert(0,1,'\t').insert(0,value_);
	 else if (attr_name_ == "y")
			istringstream(value_) >> Y;
	 else if (attr_name_ == "z")
			istringstream(value_) >> Z;*/
}

void tAnalysis_XML_Reader::ReadData(const string& tag_, ifstream_XML& xmlfile_)
{
	if (tag_ == "analysis")
		ReadData(string("minprofile"), xmlfile_);
	else
		xmlfile_.ReadObject(*this, tag_);
}

void tAnalysis_XML_Reader::SetValue(const string& attr_name_, const string& value_)
{
	if (attr_name_ == "type")
		m_type = value_;
	else if (attr_name_ == "id" || attr_name_ == "name")
		m_name = value_;
	else if (attr_name_ == "eval_per_fe")
		m_evalPerFEMethod = value_ == "average" ? tAnalysis::eMeanOverFE : tAnalysis::eAtFECentre;
	else if (attr_name_ == "loadcase")
		m_pCase = m_loads.HasCase(value_) ? &m_loads.Case(value_) : &m_loads.Case(1);
	else if (attr_name_ == "outputfile")
		m_oFileName = value_;
}
