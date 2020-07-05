
#include "FiniteElementModel.hpp"
#include "FiniteElementPool.hpp"
#include "Node.hpp"
#include "FileInOut.hpp"
#include "XML_Readers.hpp"
#include "VectorMatrix.hpp"
#include "ProcessDisplay.hpp"
#include "ThrowMessage.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>

using namespace std;

tFE_model::~tFE_model()
{
	for_each(m_pFinEls.begin(), m_pFinEls.end(), fDelete_object<tFinitElement>());
	for_each(m_pMaterials.begin(), m_pMaterials.end(), fDelete_object<tMaterial>());
}

tMaterial& tFE_model::AddMaterial(const string& kindName_, const string& identifier_)
{
	vector<tMaterial*>::iterator pfreePosition =
		find(m_pMaterials.begin(), m_pMaterials.end(), (tMaterial*)nullptr);
#ifdef STRONGCHECK
	Assert(pfreePosition < m_pMaterials.end(), "there are no positions for new materials");
#endif
	*pfreePosition = tMaterial::CreateMaterial(kindName_);
	(*pfreePosition)->SetName(identifier_);
	return **pfreePosition;
}

tFinitElement& tFE_model::AddFE(const string& kindName_, size_t number_)
{
	Assert(number_ > 0 && number_ <= m_pFinEls.size(), "invalid # of FE");
	tFinitElement*& pcurrentFE = m_pFinEls[--number_];
#ifdef STRONGCHECK
	Assert(pcurrentFE == nullptr, "requested position for new FE is already engaged");
#endif
	pcurrentFE = tFinitElement::NewFE(kindName_);
	return *pcurrentFE;
}

struct material_has_name
{
	material_has_name(const string& name_) : m_name(name_)
	{
	}
	bool operator()(const tMaterial* pmaterial_) const
#ifdef STRONGCHECK
	{
		return pmaterial_ != nullptr && pmaterial_->GetName() == m_name;
	}
#else
	{
		return pmaterial_->Name() == Name;
	}
#endif

	const string& m_name;
};

tMaterial& tFE_model::Material(const string& name_)
{
#ifdef STRONGCHECK
	Assert(
		find_if(m_pMaterials.begin(), m_pMaterials.end(), material_has_name(name_)) <
			m_pMaterials.end(),
		tMessage("Material not found with name ") << name_);
#endif
	return **find_if(m_pMaterials.begin(), m_pMaterials.end(), material_has_name(name_));
}

const tMaterial& tFE_model::Material(const string& name_) const
{
#ifdef STRONGCHECK
	Assert(
		find_if(m_pMaterials.begin(), m_pMaterials.end(), material_has_name(name_)) <
			m_pMaterials.end(),
		tMessage("Material not found with name ") << name_);
#endif
	return **find_if(m_pMaterials.begin(), m_pMaterials.end(), material_has_name(name_));
}

struct fmax_nodes_count
{
	size_t operator()(size_t prev_sum_, const tFinitElement* pfe_) const
	{
		return max(prev_sum_, pfe_->HowManyNodes());
	}
};

size_t tFE_model::MaxNodesPerFE() const
{
	return accumulate(m_pFinEls.begin(), m_pFinEls.end(), size_t(0), fmax_nodes_count());
	/* size_t result = 0;
	 for (vector<tFinitElement*>::const_iterator ppFE=pFinEls.begin(); ppFE<pFinEls.end(); ++ppFE)
	   if   (*ppFE != nullptr)
		 if ((*ppFE)->HowManyNodes() > result)
		   result = (*ppFE)->HowManyNodes();
	 return result;*/
}

template <typename T>
struct fadd1_if_FE_has
{
	fadd1_if_FE_has(const T& o_) : m_ptr(&o_)
	{
	}
	size_t operator()(size_t prev_sum_, const tFinitElement* pfe_) const
	{
		return pfe_->Has(*m_ptr) ? prev_sum_ + 1 : prev_sum_;
	}

	const T* m_ptr;
};

template <typename T>
struct fmax_FEs_count_in
{
	fmax_FEs_count_in(
		vector<tFinitElement*>::const_iterator p1_, vector<tFinitElement*>::const_iterator p2_)
		: m_pBegin(p1_), m_pEnd(p2_)
	{
	}
	size_t operator()(size_t prev_max_, const T& o_) const
	{
		return max(
			prev_max_, size_t(accumulate(m_pBegin, m_pEnd, size_t(0), fadd1_if_FE_has<T>(o_))));
	}
	size_t operator()(size_t prev_max_, const T* o_) const
	{
		return max(
			prev_max_, size_t(accumulate(m_pBegin, m_pEnd, size_t(0), fadd1_if_FE_has<T>(*o_))));
	}

	vector<tFinitElement*>::const_iterator m_pBegin, m_pEnd;
};

size_t tFE_model::MaxFinElsPerNode() const
{
	return accumulate(
		m_nodes.begin(),
		m_nodes.end(),
		size_t(0),
		fmax_FEs_count_in<tNode>(m_pFinEls.begin(), m_pFinEls.end()));
	/* size_t result = 0;
	 for (vector<tNode>::const_iterator pnode=Nodes.begin(); pnode<Nodes.end(); ++pnode)
	   if (pnode->NumberOfLinkedFinEls() > result)
		 result = pnode->NumberOfLinkedFinEls();
	 return result;*/
}

size_t tFE_model::MaxFinElsPerMatl() const
{
	return accumulate(
		m_pMaterials.begin(),
		m_pMaterials.end(),
		size_t(0),
		fmax_FEs_count_in<tMaterial>(m_pFinEls.begin(), m_pFinEls.end()));
	/* size_t result = 0;
	 for (vector<tMaterial*>::const_iterator ppmatl = pMaterials.begin(); ppmatl < pMaterials.end();
	 ++ppmatl) if ((*ppmatl)->NumberOfLinkedFinEls() > result) result =
	 (*ppmatl)->NumberOfLinkedFinEls(); return result;*/
}

struct AddVolumes
{
	real_t operator()(real_t prev_vol_, const tFinitElement* pfe_) const
	{
		return prev_vol_ + pfe_->Volume();
	}
};

real_t tFE_model::Volume() const
{
	return accumulate(m_pFinEls.begin(), m_pFinEls.end(), 0., AddVolumes());
	/* real_t result = 0.;
	  for (vector<tFinitElement*>::const_iterator ppFE=pFinEls.begin(); ppFE<pFinEls.end(); ++ppFE)
	   result += (*ppFE)->Volume();
	 return result;*/
}

// real_t tFE_model::Mass() const
//{
// real_t result = 0.;
// for (vector<tFinitElement*>::const_iterator ppFE=pFinEls.begin(); ppFE<pFinEls.end(); ++ppFE)
//   result += (*ppFE)->Mass();
// return result;
//}

tFE_model& tFE_model::ImportFromFile(const char* fileName_)
{
	ifstream_XML xml_in(fileName_);
	Assert(!!xml_in, tMessage("Cannot open file ") << fileName_);
	xml_in.clear();

	if (xml_in.FindHighLevelTag("model") && xml_in.get() != '>')
	{
		xml_in.FindChar('"');
		GETLINE(xml_in, m_caption, '"');
	}
	else
	{
		m_caption = fileName_;
		if (m_caption.substr(m_caption.length() - 4, 4)[0] == '.')
			m_caption.resize(m_caption.length() - 4);
	}
	streampos beg_pos = xml_in.tellg();
	ImportNodesFromFile(xml_in);

	// xml_in.close();
	// xml_in.open(fileName_);
	xml_in.seekg(beg_pos);
	// xml_in.clear();
	ImportMatlsFromFile(xml_in);

	// xml_in.close();
	// xml_in.open(fileName_);
	xml_in.seekg(beg_pos);
	// xml_in.clear();
	ImportFinElsFromFile(xml_in);

	return *this;
}

tFE_model& tFE_model::ImportNodesFromFile(ifstream_XML& xml_in_)
{
	const size_t numberOfNodes = xml_in_.CountArray("nodes");
	Assert(numberOfNodes > 1, "No nodes has been found in file");
	m_nodes.resize(numberOfNodes);
	tProcessDisplay display("Importing nodes", static_cast<data_t>(numberOfNodes));
	tNode_XML_Reader reader;
	for (size_t i = 0; i < numberOfNodes; ++i, ++display)
	{
		xml_in_.ReadObject(reader);
		reader.SetValsTo(m_nodes.at(reader.NodeNum() - 1));
	}
	return *this;
}

tFE_model& tFE_model::ImportMatlsFromFile(ifstream_XML& xml_in_)
{
	size_t numberOfMatls = xml_in_.CountArray("materials");
	Assert(numberOfMatls > 0, "No materials has been found in file");
	tProcessDisplay display("Importing Materials", static_cast<data_t>(numberOfMatls));

	for_each(m_pMaterials.begin(), m_pMaterials.end(), fDelete_object<tMaterial>());
	m_pMaterials.resize(numberOfMatls, nullptr);

	tMaterial_XML_Reader reader;
	for (size_t i = 0; i < numberOfMatls; ++i, ++display)
	{
		xml_in_.ReadObject(reader);
		m_pMaterials[i] = reader.CreateMaterial();
	}
	return *this;
}

tFE_model& tFE_model::ImportFinElsFromFile(ifstream_XML& xml_in_)
{
	tFinEl_XML_Reader reader(*this);
	size_t numberOfFinEls = xml_in_.CountArray("finite_elements");
	Assert(numberOfFinEls > 0, "No finite elements has been found in file");
	tProcessDisplay display("Importing Finite Elements", static_cast<data_t>(numberOfFinEls));
	// vector<tFinitElement*>::iterator  pptr = pFinEls.begin();
	for_each(m_pFinEls.begin(), m_pFinEls.end(), fDelete_object<tFinitElement>());
	m_pFinEls.resize(numberOfFinEls, nullptr);
	for (size_t i = 0; i < numberOfFinEls; ++i, ++display)
	{
		xml_in_.ReadObject(reader);
		m_pFinEls[reader.FE_Num() - 1] = reader.CreateFE();
	}
	return *this;
}

void tFE_model::Clear()
{
	if (IsEmpty())
		return;
	m_nodes.clear();
	for_each(m_pFinEls.begin(), m_pFinEls.end(), fDelete_object<tFinitElement>());
	for_each(m_pMaterials.begin(), m_pMaterials.end(), fDelete_object<tMaterial>());
}

void tFE_model::Dump(ostream& out) const
{
	out << "Nodes:\n#\tX\tY\tZ\t(X\tY\tZ) - fixed DOFs\n";
	for (size_t i = 0; i < m_nodes.size(); ++i)
		out << i + 1 << '\t' << m_nodes[i].Coord()(1) << '\t' << m_nodes[i].Coord()(2) << '\t'
			<< m_nodes[i].Coord()(3) << '\t' << m_nodes[i].HasFixed(X) << '\t'
			<< m_nodes[i].HasFixed(Y) << '\t' << m_nodes[i].HasFixed(Z) << '\n';

	out << "\nMaterials:\n#\tType\tName\tParams\n";
	for (size_t i = 0; i < m_pMaterials.size(); ++i)
		out << i + 1 << ' ' << m_pMaterials[i]->Kind().c_str() << '\t'
			<< m_pMaterials[i]->GetName().c_str() << "\tDensity = " << m_pMaterials[i]->Density()
			<< "\tE = " << dynamic_cast<tIsotropicMaterial*>(m_pMaterials[i])->YoungModule() << '\t'
			<< "\tNu = " << dynamic_cast<tIsotropicMaterial*>(m_pMaterials[i])->PoissonRatio()
			<< '\n';

	out << "\nFinite elements:\n#\tType\tMaterial\tNodes\tParams\n";
	for (size_t i = 0; i < m_pFinEls.size(); ++i)
	{
		out << i + 1 << "\t";
		out << m_pFinEls[i]->Material().Kind().c_str() << '\t'
			<< m_pFinEls[i]->Material().GetName().c_str() << '\t';
		for (size_t j = 1; j <= m_pFinEls[i]->HowManyNodes(); ++j)
			out << (&m_pFinEls[i]->Node(j) - &m_nodes.front()) + 1 << ' ';
		switch (m_pFinEls[i]->SpaceDimension())
		{
			case 1:
				out << "\tArea = " << dynamic_cast<t1D_FE*>(m_pFinEls[i])->CrosSecArea()
					<< "\t1d-Str"
					<< (dynamic_cast<t1D_FE*>(m_pFinEls[i])->Has1dStressState() ? "ess" : "ain");
				break;
			case 2:
				out << "\tThick = " << dynamic_cast<t2D_FE*>(m_pFinEls[i])->Thickness()
					<< "\t2d-Str"
					<< (dynamic_cast<t2D_FE*>(m_pFinEls[i])->Has2dStressState() ? "ess" : "ain");
				break;
			default:
				THROW_MESSAGE("Unsupported space dimension");
		}
		out << '\n';
	}
}
