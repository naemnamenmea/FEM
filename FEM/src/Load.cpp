
#include "Load.hpp"
#include "VectorMatrix.hpp"
#include "XML_Readers.hpp"

using namespace std;

static const size_t nodal_load_kinds_number = 1;
static const char* initNames[] = {"Static"};  // Add new kind names here
static tNodalLoad* (*initFunctions[])() = {
	tStatictNodalLoad::NewLoad};  // Add new pointers to creaters here

const tNamesToFuns<tNodalLoad> tNodalLoad::Factory(
	initNames, initFunctions, nodal_load_kinds_number);

const string& tStatictNodalLoad::KindName =
	tNodalLoad::Factory.FindName(tStatictNodalLoad::NewLoad);

tNodalLoad& tStatictNodalLoad::ReadData(istream& str_in_)
{
	real_t value;
	/* str_in_ >> value;   Assign(X,value);
 str_in_ >> value;   Assign(Y,value);
 str_in_ >> value;
	 * Assign(Z,value);*/
	string buffer;
	for (int i = 0; i < 3; ++i)
	{
		getline(str_in_, buffer, '=');
		str_in_.get();
		while (isspace(buffer[buffer.size() - 1])) buffer.resize(buffer.size() - 1);  // pop_back();
		while (isspace(buffer[0])) buffer.erase(0, 1);
		str_in_ >> value;
		if (buffer == "x")
			Assign(X, value);
		else if (buffer == "y")
			Assign(Y, value);
		else if (buffer == "z")
			Assign(Z, value);
		else
			throw tMessage("invalid data for components of nodal force");
		//         {str_in_ >> buffer;
		//          --i;}
	}
	return *this;
}

tNodalTensor1& tStatictNodalLoad::Level(tNodalTensor1& result_, real_t part_) const
{
#ifdef STRONGCHECK
	Assert(part_ >= 0., "negative time in tStaticNodalLoad::Level");
	Assert(mathdef::is_lte(part_, 1.), "time > 1 in tStaticNodalLoad::Level");
	Assert(m_pNode != nullptr, "load with no linked nodes in tStaticNodalLoad::Level");
#endif
	result_ = *this;
	if (!is0())
		result_ *= part_;
	// if (EQ0(part_))
	//          Assign0();
	return result_;
}

tLoadCase::~tLoadCase()
{
	for_each(m_pNodalLoads.begin(), m_pNodalLoads.end(), fDelete_object<tNodalLoad>());
	// for (vector<tNodalLoad*>::iterator p=pNodalLoads.begin(); p<pNodalLoads.end(); ++p)
	//    delete *p;
}

tNodalLoad& tLoadCase::AddNodalForce(const string& kindName_, const tNode& nodeToLink_)
{
	vector<tNodalLoad*>::iterator pfreePosition =
		find(m_pNodalLoads.begin(), m_pNodalLoads.end(), (tNodalLoad*)nullptr);
#ifdef STRONGCHECK
	Assert(pfreePosition < m_pNodalLoads.end(), "there are no positions for new load forces");
#endif
	*pfreePosition = tNodalLoad::NewNodalLoad(kindName_);
	(*pfreePosition)->Link(nodeToLink_);
	// if (kindName_ == tStatictNodalLoad::KindName) *plastEngagedPosition = new
	// tStatictNodalLoad(nodeToLink_); else throw tMessage("Invalid load kind name");
	return **pfreePosition;
}

struct fAddLoadLevelFrom1node
{
	fAddLoadLevelFrom1node(tNodalTensor1Column& result_, real_t level_, tNodalTensor1& tmp_)
		: m_result(result_), m_level(level_), m_tmp(tmp_)
	{
	}
	void operator()(const tNodalLoad* p_) const
	{
		//                                              result += p_->Level(tmp,level);
		m_result.AssignComponent(p_->Level(m_tmp, m_level));
	}

	tNodalTensor1Column& m_result;
	real_t m_level;
	tNodalTensor1& m_tmp;
};

tNodalTensor1Column& tLoadCase::LoadLevel(tNodalTensor1Column& result_, real_t level_) const
{
#ifdef STRONGCHECK
	Assert(!m_pNodalLoads.empty(), "no loads in tLoadCase::LoadLevel");
	Assert(level_ >= 0., "negative time in tLoadCase::LoadLevel");
#endif
	tNodalTensor1 tmp;
	for_each(
		m_pNodalLoads.begin(), m_pNodalLoads.end(), fAddLoadLevelFrom1node(result_, level_, tmp));
	return result_;
}

tLoadCase& tLoadCase::ImportFromFile(ifstream_XML& xml_in, const tFE_model& modelToLoad_)
{
	xml_in.FindWord("<loadcase");
	tLoadCase_XML_Reader reader;
	xml_in.ReadAttribs(reader);
	reader.SetLoadCaseName(*this);
	size_t numOfForces = xml_in.CountItems("loadcase", "force");
	for_each(m_pNodalLoads.begin(), m_pNodalLoads.end(), fDelete_object<tNodalLoad>());
	m_pNodalLoads.resize(numOfForces, nullptr);
	for (size_t i = 0; i < numOfForces; ++i)
	{
		xml_in.ReadObject(reader, "force");
		reader.CreateForce(*this, modelToLoad_);
	}
	return *this;
}

tLoads& tLoads::ImportFromFile(const char* fileName_, const tFE_model& modelToLoad_)
{
	ifstream_XML xml_in(fileName_);
	Assert(!!xml_in, tMessage("Cannot open file ") << fileName_);
	size_t numOfLoadCases = xml_in.CountArray("loads");
	Assert(numOfLoadCases > 0, "No loads has been found in file");
	m_loadCases.resize(numOfLoadCases);
	for (size_t i = 0; i < numOfLoadCases; ++i) m_loadCases[i].ImportFromFile(xml_in, modelToLoad_);
	return *this;
}

tLoadCase& tLoads::Case(const string& name_)
{
	vector<tLoadCase>::iterator pcase = m_loadCases.begin();
	for (; pcase < m_loadCases.end(); ++pcase)
		if (pcase->Name() == name_)
			return *pcase;
	throw tMessage("invalid name in tLoads::Case");
}

const tLoadCase& tLoads::Case(const string& name_) const
{
	vector<tLoadCase>::const_iterator pcase = m_loadCases.begin();
	for (; pcase < m_loadCases.end(); ++pcase)
		if (pcase->Name() == name_)
			return *pcase;
	throw tMessage("invalid name in tLoads::Case const");
}

void tLoads::Dump(ostream& out, const tFE_model& model_) const
{
	tNodalTensor1Column loadlevel;
	loadlevel.Link(model_);
	for (size_t i = 0; i < m_loadCases.size(); ++i)
	{
		out << "Case " << i + 1 << "\t\"" << m_loadCases[i].Name().c_str()
			<< "\"\n\tNode #\tForce: X Y Z\n";
		loadlevel.Assign0();
		m_loadCases[i].LoadLevel(loadlevel);
		for (size_t j = 1; j <= loadlevel.Dimension(); ++j)
			if (!loadlevel(j).is0())
				out << '\t' << j << '\t' << loadlevel(j)(1) << ' ' << loadlevel(j)(2) << ' '
					<< loadlevel(j)(3) << '\n';
	}
}
