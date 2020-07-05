
#include "Analysis.hpp"
#include "LinearElasticAnalysis.hpp"
#include "NonLAnalys.hpp"
#include "XML_Readers.hpp"

#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif
#include <iterator>
#include <fstream>	//for debug
#include <ctime>	//for debug
#include <cstdlib>	//exit() only

using namespace std;

static const size_t analysis_kinds_number = 2;	// Set the number of analyses here!
static const char* initNames[] = {"LinearElastic", "PhysNonlinear"};  // Add new kind names here!
static tAnalysis* (*initFunctions[])() = {
	tLinearElasticAnalysis::NewAnalysis,
	PhysNonlinearAnalysis::NewAnalysis};  // Add new pointers to creators here!

const tNamesToFuns<tAnalysis> tAnalysis::Factory(initNames, initFunctions, analysis_kinds_number);

const string& PhysNonlinearAnalysis::KindName =
	tAnalysis::Factory.FindName(PhysNonlinearAnalysis::NewAnalysis);

tAnalysis::tAnalysis()
	: m_caption(),
	  m_pFEModel(nullptr),
	  m_pLoadCase(nullptr),
	  m_minimizeBandMethod(tMarkedGraph::no_permutations),
	  m_pInitialDisplLevel(nullptr),
	  m_pInitialLoadLevel(nullptr),
	  /*InitialTime(0.), */ m_finalTime(1.),
	  m_pOrigGraph(nullptr),
	  m_pActiveGraph(nullptr),
	  m_pMinimizedGraph(nullptr),
	  m_pPrologMinGraph(nullptr),
	  m_solutionFound(false),
	  EvalPerFEMethod(eAtFECentre),
	  m_searchDepth(0)
{
}

tAnalysis& tAnalysis::Run()
{
	if (m_pMinimizedGraph != nullptr)
	{
		if (m_pMinimizedGraph->isEmpty())
			*m_pMinimizedGraph = *m_pOrigGraph;

		m_pMinimizedGraph->ReMark(m_minimizeBandMethod);
		m_pActiveGraph = m_pMinimizedGraph;
	}

	if (m_pPrologMinGraph != nullptr)
	{
		*m_pPrologMinGraph = m_pMinimizedGraph == nullptr ? *m_pOrigGraph : *m_pMinimizedGraph;

		/*    if(pMinimizedGraph==nullptr) *pPrologMinGraph = *pOrigGraph;
			else *pPrologMinGraph = *pMinimizedGraph;*/

		m_pPrologMinGraph->Prolog_ReMark(m_searchDepth);
		m_pActiveGraph = m_pPrologMinGraph;
	}
	return *this;
}

/*void dump_all_fes(const tFE_model* const pModel_)
{
  size_t fes_num = pModel_->HowManyFinEls();
  ofstream out("__All_FEs_Data.txt");
  for(size_t i = 1, nodes_num=0; i <= fes_num; ++i)
  {
	nodes_num = pModel_->FE(i).HowManyNodes();
	out << "FE " << i << " Number of nodes " << nodes_num << endl;
	for(size_t j=1; j <= nodes_num; ++j)
	  out << "\tNode " << j << "   " << pModel_->FE(i).Node(j).Coord()(1) << '\t'
		  << pModel_->FE(i).Node(j).Coord()(2) << '\t'
	  << pModel_->FE(i).Node(j).Coord()(3) << endl;
  }
} */

struct fDelete_analysis
{
	void operator()(const AnalysesArray::PtrToAnalysis& o_) const
	{
		delete o_.m_pAnalysis;
	}
};

AnalysesArray::~AnalysesArray()
{
	for_each(m_analyses.begin(), m_analyses.end(), fDelete_analysis());
}

void AnalysesArray::Delete(const char* NameToFind)
{
	set<PtrToAnalysis>::iterator p = m_analyses.begin();
	for (; p != m_analyses.end(); ++p)
		if (p->m_pAnalysis->Name() == NameToFind)
		{
			if (p->m_pAnalysis == m_pCurrent)
				m_pCurrent = nullptr;
			delete p->m_pAnalysis;
			m_analyses.erase(p);
			return;
		}
#ifdef STRONGCHECK
	throw tMessage("Can't find analysis with name ") << NameToFind;
#endif
}

void AnalysesArray::ImportFromFile(
	const char* fileName_, const tFE_model& modelToLoad_, const tLoads& loadCases_)
{
	ifstream_XML xml_in(fileName_);
	Assert(!!xml_in, tMessage("Cannot open file ") << fileName_);

	tAnalysis_XML_Reader reader(modelToLoad_, loadCases_);

	size_t numberOfAnalyses = xml_in.CountArray("analyses");
	Assert(numberOfAnalyses > 0, "No analyses has been found in file");

	for (size_t i = 0; i < numberOfAnalyses; ++i)
	{
		xml_in.ReadObject(reader);
#ifdef STRONGCHECK
		Assert(
			m_analyses.insert(reader.CreateAnalysis()).second,
			"Can't insert new analisys Analyses::Add");
#else
		m_analyses.insert(reader.CreateAnalysis());
#endif
	}
	m_pCurrent = m_analyses.begin()->m_pAnalysis;
}

tAnalysis& AnalysesArray::operator()(const char* name_) const
{
	set<PtrToAnalysis>::const_iterator p = m_analyses.begin();
	for (; p != m_analyses.end(); ++p)
		if (p->m_pAnalysis->Name() == name_)
			break;
	return *(p->m_pAnalysis);
}

tAnalysis& AnalysesArray::operator()(int n_) const
{
	set<PtrToAnalysis>::const_iterator p = m_analyses.begin();
	advance(p, --n_);
	return *(p->m_pAnalysis);
}

void AnalysesArray::DeleteCurrent()
{
	set<PtrToAnalysis>::iterator iter = m_analyses.find(m_pCurrent);
#ifdef STRONGCHECK
	Assert(iter != m_analyses.end(), "Current analysis not found in Analyses");
#endif
	delete iter->m_pAnalysis;
	m_analyses.erase(iter);
	m_pCurrent = nullptr;
}

struct has_assembled_matrix
{
	bool operator()(const AnalysesArray::PtrToAnalysis& p_) const
	{
		return p_.m_pAnalysis->Completed();
	}
};

void AnalysesArray::SetInitialGraph()
{
	// ofstream log("SenInitGr.txt");

	// cout << "orig graph is empty " << boolalpha << !OrigGraph.isEmpty() << endl;

	if (!m_origGraph.isEmpty())
		return;
	set<PtrToAnalysis>::const_iterator p =
		find_if(m_analyses.begin(), m_analyses.end(), has_assembled_matrix());

	const tAnalysis* p_anls = p == m_analyses.end() ? m_pCurrent : p->m_pAnalysis;

	const t1StepAnalysis* p_1st_anls = &p_anls->Step(1);

	if (p == m_analyses.end() && !p_1st_anls->Completed())
	{  // cout << "Assemble all matrices" << endl;
		const_cast<t1StepAnalysis*>(p_1st_anls)->AssembleAllMatrices();
		// cout << "done" << endl;
	}
	// cout << "OrigGraph.Link(p_1st_anls->Matrix(1));" << endl;
	m_origGraph.Link(p_1st_anls->Matrix(1));
	// cout << "done" <<endl;
}

void AnalysesArray::RunCurrent()
{
#ifdef STRONGCHECK
	Assert(!isEmpty(), "There are no analyses");
#endif

	//  ofstream log1("analysis_log.txt");
	//  cout << "setting initial graph... ";

	SetInitialGraph();

	//  cout << "done\n";

	switch (m_graphkind)
	{
		case Minimize:

			//       cout << "use prev minimized graph\n";

			m_pCurrent->Link(m_origGraph, &m_minimizedGraph);
			break;
		case PrologMin:

			//       cout << "use prev min and prolog\n";

			m_pCurrent->Link(m_origGraph, &m_minimizedGraph, &m_prologMinGraph);
			break;
		case OnlyPrologMin:

			//       cout << "prolog only" << endl;

			m_pCurrent->Link(m_origGraph, nullptr, &m_prologMinGraph);
			break;
		default:

			//       cout << "use original graph" << endl;

			m_pCurrent->Link(m_origGraph);
	}

	//  cout << "runing analysis... " << endl;

	m_pCurrent->Run();

	//  cout << "analysis complate\n";
}

void AnalysesArray::Clear()
{
	if (!isEmpty())
	{
		for_each(m_analyses.begin(), m_analyses.end(), fDelete_analysis());
		m_analyses.clear();
	}
	m_pCurrent = nullptr;
	m_graphkind = Original;
}
