
#include "GraphTheory.hpp"
#include "Node.hpp"
#include "VectorMatrix.hpp"
#include "ThrowMessage.hpp"
#include "ProcessDisplay.hpp"

#include <cstdlib>	// sprintf
#include <algorithm>
#include <functional>
#include <fstream>
#include <iterator>

//#define PL_KERNEL
//#include "SWI-cpp.h"

using namespace std;

const tMarkedGraph::tNode& tMarkedGraph::tNode::UnLinkAll() const
{
	for (set<const tNode*>::const_iterator pnode = m_pLinkedNodes.begin();
		 pnode != m_pLinkedNodes.end();
		 ++pnode)
#ifdef STRONGCHECK
		Assert((*pnode)->UnLink(this), "unlink error in tMarkedGraph::tNode::UnLinkAll");
#else
		(*pnode)->UnLink(this);
#endif
	m_pLinkedNodes.clear();
	return *this;
}

const tMarkedGraph::tNode& tMarkedGraph::tNode::LinkAdjToEachOther() const
{
	set<const tNode*>::const_iterator pnode1 = m_pLinkedNodes.begin(), pnode2;
	for (; pnode1 != m_pLinkedNodes.end(); ++pnode1)
		for (pnode2 = pnode1; (++pnode2) != m_pLinkedNodes.end();)
			if ((*pnode1)->Link(*pnode2))
#ifdef STRONGCHECK
				Assert(
					(*pnode2)->Link(*pnode1),
					"link error in tMarkedGraph::tNode::LinkAdjToEachOther");
#else
				(*pnode2)->Link(*pnode1);
#endif
	return *this;
}

const tMarkedGraph::tNode& tMarkedGraph::tNode::PointAbsentAdj(
	vector<const tMarkedGraph::tNode*>& pointed_,
	vector<const tMarkedGraph::tNode*>::iterator& plastOfPreviouslyPointed_,
	bool sortIsRequired_) const
// Checks Agj for presence pointers to these in [pointed_.begin(),plastOfPreviouslyPointed_],
// allocates pointers to those Adj elements which are absent there,
// sorts (if sortIsRequired_) all allocated by order of degrees increase,
// and moves plastOfPreviouslyPointed_ to the last of allocated
{
#ifdef STRONGCHECK
	Assert(
		find(pointed_.begin(), pointed_.end(), *plastOfPreviouslyPointed_) != pointed_.end(),
		"invalid arguments in tMarkedGraph::PointAbsentAdj");
#endif
	set<const tNode*>::const_iterator padj = m_pLinkedNodes.begin();
	vector<const tNode*>::iterator padjcopy = plastOfPreviouslyPointed_;
	for (++plastOfPreviouslyPointed_; padj != m_pLinkedNodes.end(); ++padj)
		if (find(pointed_.begin(), plastOfPreviouslyPointed_, *padj) == plastOfPreviouslyPointed_)
			*(++padjcopy) = *padj;
	if (sortIsRequired_)
		sort(plastOfPreviouslyPointed_, padjcopy + 1, fLessDegree);
	plastOfPreviouslyPointed_ = padjcopy;
	return *this;
}

bool tMarkedGraph::tNode::UpdateLinks() const
{
	if (m_previousLocation != this)
	{
		for (set<const tNode*>::const_iterator ppnode = m_pLinkedNodes.begin();
			 ppnode != m_pLinkedNodes.end();
			 ++ppnode)
		{
#ifdef STRONGCHECK
			Assert(
				(*ppnode)->m_pLinkedNodes.erase(m_previousLocation) > 0 &&
					(*ppnode)->m_pLinkedNodes.insert(this).second,
				"bad link in tMarkedGraph::tNode::UpdateLinks");
#else
			(*ppnode)->pLinkedNodes.erase(PreviousLocation);
			(*ppnode)->pLinkedNodes.insert(this);
#endif
		}
		m_previousLocation = this;
		return true;
	}
	return false;
}

tMarkedGraph& tMarkedGraph::operator=(const tMarkedGraph& original_)
{
	m_nodes = original_.m_nodes;
	vector<tNode>::iterator pnode1 = m_nodes.begin(), pnode2;
	for (; pnode1 < m_nodes.end(); ++pnode1)
		for (pnode2 = pnode1; (++pnode2) < m_nodes.end();)
			if (pnode2->UpdateLink(*pnode1))
#ifdef STRONGCHECK
				Assert(
					pnode1->UpdateLink(*pnode2),
					"non-mutual link in tMarkedGraph::tMarkedGraph(tMarkedGraph&)");
#else
				pnode1->UpdateLink(*pnode2);
#endif
	return *this;
}

tMarkedGraph::tMarkedGraph(const tMarkedGraph& original_, const string& NodesList_)
	: m_nodes(original_.m_nodes.size())
{
	unsigned int j = 0;
	if (NodesList_[j] != '[' && NodesList_[NodesList_.size() - 1] != ']')
		throw tMessage("Not valid nodes list");
	string Node_num;
	//   Node_num.reserve(10);
	for (size_t i = original_.m_nodes.size(); i > 0; --i)
	{
		size_t i_ = i - 1;
		++j;
		while (!(NodesList_[j] == ',' || NodesList_[j] == ']' || NodesList_[j] == '\0'))
		{
			Node_num += NodesList_[j];
			++j;
		}
		m_nodes[i_] = original_.m_nodes[static_cast<size_t>(atoi(Node_num.c_str()) - 1)];
		Node_num.clear();
	}

	vector<tNode>::iterator pnode1 = m_nodes.begin(), pnode2;
	for (; pnode1 < m_nodes.end(); ++pnode1)
		for (pnode2 = pnode1; (++pnode2) < m_nodes.end();)
			if (pnode2->UpdateLink(*pnode1))
#ifdef STRONGCHECK
				Assert(
					pnode1->UpdateLink(*pnode2),
					"non-mutual link in tMarkedGraph::tMarkedGraph(tMarkedGraph&)");
#else
				pnode1->UpdateLink(*pnode2);
#endif
}

tMarkedGraph& tMarkedGraph::Link(const tNodalTensor2SymMatrix& SourceMatrix_)
{
	m_nodes.resize(SourceMatrix_.NumberOfDOFs(), tNode());
	vector<tNode>::iterator pnode1 = m_nodes.begin();
	size_t i = 1, j;
	const size_t sourceDim = SourceMatrix_.Dimension();
	small_t k, l;
	tProcessDisplay output("Matrix renumbering", static_cast<data_t>(sourceDim));
	const ::tNode *pgridNodei, *pgridNodej;
	for (; i <= sourceDim; ++i, ++output)
	{
		pgridNodei = &SourceMatrix_.Node(i);
		for (k = 1; k <= 3; ++k)
			if (!pgridNodei->HasFixed(k))
			{
				Assert(
					!(SourceMatrix_.is0(i, i) || SourceMatrix_(i, i).is0(k, k) ||
					  mathdef::is_lte(SourceMatrix_(i, i)(k, k), 0.)),
					tMessage("Pivot corresponding node ")
						<< i << " component " << k << " equals "
						<< (SourceMatrix_.is0(i, i) || SourceMatrix_(i, i).is0(k, k)
								? real_t(0)
								: SourceMatrix_(i, i)(k, k))
						<< "!!!\nIt must be >0 to the matrix to be positive definite");
				pnode1->Link(i, k);
				++pnode1;
			}
	}
	const Tensor2a* ptensor;
	vector<tNode>::iterator pnode2;
	// output = 0;
	output.Reset("Matrix renumbering", static_cast<data_t>(sourceDim));
	for (i = 1; i <= sourceDim; ++i, ++output)
	{
		pgridNodei = &SourceMatrix_.Node(i);
		for (j = i - SourceMatrix_.BandWidth(i); j <= i; ++j)
			if (!SourceMatrix_.is0(i, j))
			{
				pgridNodej = &SourceMatrix_.Node(j);
				ptensor = &SourceMatrix_(i, j);
				for (k = 1; k <= 3; ++k)
					if (!pgridNodei->HasFixed(k))
					{
						pnode1 = m_nodes.end();
						for (l = 1; l < (i == j ? k : 4); ++l)
							if ((!pgridNodej->HasFixed(l)) &&
								(!ptensor->is0(k, l)) /* && (*ptensor)(k,l)!=0.*/)
							{
								if (pnode1 == m_nodes.end())
									pnode1 = find(m_nodes.begin(), m_nodes.end(), tNode(i, k));
								pnode2 = find(m_nodes.begin(), m_nodes.end(), tNode(j, l));
#ifdef STRONGCHECK
								Assert(
									pnode1 < m_nodes.end(),
									"node1 not found in "
									"tMarkedGraph::Link(tNodalTensor2SymMatrix&)");
								Assert(
									pnode2 < m_nodes.end(),
									"node2 not found in "
									"tMarkedGraph::Link(tNodalTensor2SymMatrix&)");
								Assert(
									pnode1->Link(&*pnode2),
									"link error in tMarkedGraph::Link(tNodalTensor2SymMatrix&)");
								Assert(
									pnode2->Link(&*pnode1),
									"link error found in "
									"tMarkedGraph::Link(tNodalTensor2SymMatrix&)");
#else
								pnode1->Link(&*pnode2);
								pnode2->Link(&*pnode1);
#endif
							}
					}
			}
	}
	return *this;
}

tMarkedGraph& tMarkedGraph::SwapNodes(size_t node1No_, size_t node2No_)
{
#ifdef STRONGCHECK
	Assert(
		node1No_ > 0 && node2No_ > 0 && node1No_ <= m_nodes.size() && node2No_ <= m_nodes.size(),
		"invalid index in tMarkedGraph::SwapNodes");
#endif
	if ((node1No_--) == (node2No_--))
		return *this;
	const tNode tmp = m_nodes[node1No_];
	tmp.UpdateLinks();
	m_nodes[node1No_] = m_nodes[node2No_];
	m_nodes[node1No_].UpdateLinks();
	m_nodes[node2No_] = tmp;
	m_nodes[node2No_].UpdateLinks();
	return *this;
}

tMarkedGraph& tMarkedGraph::SwapNodes(
	const tMarkedGraph::tNode& node1_, const tMarkedGraph::tNode& node2_)
{
#ifdef STRONGCHECK
	Assert(
		&*find(m_nodes.begin(), m_nodes.end(), node1_) == &node1_,
		"absent nodes in tMarkedGraph::SwapNodes");
	Assert(
		&*find(m_nodes.begin(), m_nodes.end(), node2_) == &node2_,
		"absent nodes in tMarkedGraph::SwapNodes");
#endif
	if (&node1_ == &node2_)
		return *this;
	const tNode tmp = node1_;
	tmp.UpdateLinks();
	const_cast<tNode&>(node1_) = node2_;
	node1_.UpdateLinks();
	const_cast<tNode&>(node2_) = tmp;
	node2_.UpdateLinks();
	return *this;
}

tMarkedGraph& tMarkedGraph::ExcludeNode(size_t nodeNo_, bool linkAdj_)
{
#ifdef STRONGCHECK
	Assert(nodeNo_ > 0 && nodeNo_ <= m_nodes.size(), "invalid index in tMarkedGraph::ExcludeNode");
#endif
	// tNode* pnodeToExclude = &Nodes[--nodeNo_];
	vector<tNode>::iterator pnodeToExclude = m_nodes.begin() + --nodeNo_;
	if (linkAdj_)
		pnodeToExclude->LinkAdjToEachOther();
	pnodeToExclude->UnLinkAll();
	// Nodes.erase(pnodeToExclude); - doesn't work properly (because it uses tNode::operator=
	// without UpdateLinks())
	if (nodeNo_ + 1 < m_nodes.size())
		for (vector<tNode>::iterator pnode = pnodeToExclude, pnextNode = pnode;
			 (++pnextNode) < m_nodes.end();
			 ++pnode)
			(*pnode = *(pnextNode)).UpdateLinks();
	m_nodes.pop_back();
	return *this;
}

size_t tMarkedGraph::MinDegNodeNo() const
{
	vector<tNode>::const_iterator pnode = m_nodes.begin();
	size_t minDeg = pnode->Degree(), currentNo = 2, result = 1;
	for (; (++pnode) < m_nodes.end(); ++currentNo)
		if (pnode->Degree() < minDeg)
		{
			minDeg = pnode->Degree();
			result = currentNo;
		}
	return result;
}

size_t tMarkedGraph::BandWidth(size_t i_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && i_ <= m_nodes.size(), "invalid index in tMarkedGraph::BandWidth");
#endif
	const tNode* pnodei = &m_nodes[--i_];
	size_t result = 0, currentResult;
	set<const tNode*>::const_iterator pnode = pnodei->m_pLinkedNodes.begin(),
									  end = pnodei->m_pLinkedNodes.end(), pnode1;
	for (;;)
	{
		pnode1 = find_if(pnode, end, bind(less<const tNode*>(), placeholders::_1, pnodei));
		if (pnode1 != end)
		{
			currentResult = pnodei - *pnode1;
			if (currentResult > result)
				result = currentResult;
			pnode = (++pnode1);
		}
		else
			break;
	}
	return result;
}

bool tMarkedGraph::HasRib(size_t r1_, size_t c1_, size_t r2_, size_t c2_) const
{
#ifdef STRONGCHECK
	Assert(
		r1_ > 0 && r2_ > 0 && c1_ > 0 && c2_ > 0 && c1_ <= 3 && c2_ <= 3,
		"invalid index in tMarkedGraph::HasRib");
#endif
	vector<tNode>::const_iterator
		pnode1 = find(m_nodes.begin(), m_nodes.end(), tNode(r1_, static_cast<small_t>(c1_))),
		pnode2 = find(m_nodes.begin(), m_nodes.end(), tNode(r2_, static_cast<small_t>(c2_)));
#ifdef STRONGCHECK
	Assert(pnode1 < m_nodes.end() && pnode2 < m_nodes.end(), "absent node in tMarkedGraph::HasRib");
	Assert(pnode1 != pnode2, "identical nodes in tMarkedGraph::HasRib");
#endif
	return pnode1->HasLink(&*pnode2);
}

const tMarkedGraph::tNode& tMarkedGraph::FindPseudoPeriferalNode(
	const tMarkedGraph::tNode& initApprox_) const
{
	size_t node1eccentricity, node2eccentricity;
	const tNode *pnode1 = &initApprox_,
				*pnode2 = &FindEccentricityAndOppositeNode(*pnode1, node1eccentricity),
				*pnextOppositeNode = &FindEccentricityAndOppositeNode(*pnode2, node2eccentricity);
	while (node2eccentricity > node1eccentricity)
	{
		node1eccentricity = node2eccentricity;
		pnode1 = pnode2;
		pnode2 = pnextOppositeNode;
		pnextOppositeNode = &FindEccentricityAndOppositeNode(*pnode2, node2eccentricity);
	}
	return *pnode2;
}

const tMarkedGraph::tNode& tMarkedGraph::FindEccentricityAndOppositeNode(
	const tMarkedGraph::tNode& node_, size_t& eccentricity_) const
{
#ifdef STRONGCHECK
	Assert(
		&*find(m_nodes.begin(), m_nodes.end(), node_) == &node_,
		"absent node in tMarkedGraph::Define1stForRCM");
#endif
	eccentricity_ = 0;
	vector<const tNode*> prootedLevelStructure(m_nodes.size());
	vector<const tNode*>::iterator p1stOnCurrentLevel,
		pcurrentOnCurrentLevel = prootedLevelStructure.begin(),
		plastInTheNextLevel = pcurrentOnCurrentLevel, p1stBeyondCurrentLevel;
	(*pcurrentOnCurrentLevel) = &node_;
	do
	{
		++eccentricity_;
		p1stBeyondCurrentLevel = plastInTheNextLevel;
		++p1stBeyondCurrentLevel;
		p1stOnCurrentLevel = pcurrentOnCurrentLevel;
		while (pcurrentOnCurrentLevel < p1stBeyondCurrentLevel)
			(*(pcurrentOnCurrentLevel++))
				->PointAbsentAdj(prootedLevelStructure, plastInTheNextLevel, false);
	} while (p1stBeyondCurrentLevel <= plastInTheNextLevel);
	--eccentricity_;
	return **min_element(p1stOnCurrentLevel, p1stBeyondCurrentLevel, tNode::fLessDegree);
}

tMarkedGraph& tMarkedGraph::RCM_Algorithm()
{
	tProcessDisplay output("RCM renumbering", 1);
	vector<tNode>::iterator pfirstOfOrdered = m_nodes.end();
	--pfirstOfOrdered;
	vector<const tNode*> psortedNodes(m_nodes.size());
	vector<const tNode*>::iterator pplastOfSorted = psortedNodes.begin(),
								   ppcurrentSorted = pplastOfSorted, ppfirstOfSorted;
	while (pfirstOfOrdered > m_nodes.begin())
	{
		(*ppcurrentSorted) = &FindPseudoPeriferalNode(*pfirstOfOrdered);
		ppfirstOfSorted = ppcurrentSorted;
		do
			(*(ppcurrentSorted++))->PointAbsentAdj(psortedNodes, pplastOfSorted, true);
		while (ppcurrentSorted <= pplastOfSorted);
		for (pplastOfSorted = ppfirstOfSorted; pplastOfSorted < ppcurrentSorted;
			 ++pplastOfSorted, --pfirstOfOrdered)
			SwapNodes(*pfirstOfOrdered, **pplastOfSorted);
	}
	++output;
	return *this;
}

tMarkedGraph& tMarkedGraph::MinDegreeAlgorithm(bool excludeWithLinkAdj_)
{
	tProcessDisplay output("\"Min degree\" renumbering", static_cast<data_t>(m_nodes.size()));
	size_t minDegNo, beginOfRestNo = 0;
	tMarkedGraph copyOfThis(*this);
	while ((++beginOfRestNo) < m_nodes.size())
	{
		minDegNo = copyOfThis.MinDegNodeNo();
		copyOfThis.SwapNodes(1, minDegNo);
		minDegNo += (beginOfRestNo - 1);
		SwapNodes(beginOfRestNo, minDegNo);
		copyOfThis.ExcludeNode(1, excludeWithLinkAdj_);
		++output;
	}
	return *this;
}

tMarkedGraph& tMarkedGraph::Prolog_ReMark(unsigned int depth_)
{
#if defined __WIN32__ || defined _WIN32
	string RibsList("[");
	for (size_t i = 0; i < m_nodes.size(); ++i)
	{
		char CurrNode[10], SecondNode[10];
		//     itoa(i+1, CurrNode, 10);
#pragma warning(push)
#pragma warning(disable : 4996)
		sprintf(CurrNode, "%zu", i + 1);
#pragma warning(pop)

		for (set<const tNode*, less<const tNode*> >::const_iterator ppLinkedNode =
				 m_nodes[i].m_pLinkedNodes.begin();
			 ppLinkedNode != m_nodes[i].m_pLinkedNodes.end();
			 ++ppLinkedNode)
		{
			vector<tNode>::const_iterator pNode =
				find(m_nodes.begin(), m_nodes.end(), **ppLinkedNode);
#ifdef STRONGCHECK
			Assert(pNode < m_nodes.end(), "unknown link in tMarkedGraph::tNode::pLinkedNodes");
#endif
			size_t SecondNum = static_cast<size_t>(pNode - m_nodes.begin());
			++SecondNum;
			if (SecondNum > i + 1)
			{
				//          itoa(SecondNum, SecondNode, 10);
#pragma warning(push)
#pragma warning(disable : 4996)
				sprintf(SecondNode, "%zu", SecondNum);
#pragma warning(pop)
				RibsList += CurrNode;
				RibsList += '-';
				RibsList += SecondNode;
				RibsList += ',';
			}
		}
	}
	RibsList[RibsList.size() - 1] = ']';
#ifdef STRONGCHECK
/*   AllocConsole();
   HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
   SetConsoleScreenBufferSize(hOut,GetLargestConsoleWindowSize(hOut));
   DWORD tmp;
   WriteConsole(hOut, RibsList.c_str(), RibsList.size(), &tmp, nullptr);*/
#endif
	ofstream RibsListData(/*(ExtractFilePath(Application->ExeName)+*/ "ribslist.txt" /*).c_str()*/);
	Assert(!!RibsListData, "Can't create ribslist.txt");
	copy(RibsList.begin(), RibsList.end(), ostream_iterator<char>(RibsListData));
	//  RibsListData << RibsList;
	RibsListData << '.';
	RibsListData << endl;
	RibsListData << depth_;
	RibsListData << '.';
	RibsListData.close();
	// Âûחמג ןנמכמדא
	STARTUPINFO si;	 // ZeroMemory( &si, sizeof(STARTUPINFO));
	si.cb = sizeof(STARTUPINFO);
	si.lpReserved = nullptr;
	// si.lpDesktop = nullptr;
	si.lpTitle = nullptr;  //"Created process";
	si.dwFlags = STARTF_USESHOWWINDOW;
	si.wShowWindow = SW_SHOWNOACTIVATE;
	si.cbReserved2 = 0;
	si.lpReserved2 = nullptr;
	PROCESS_INFORMATION pi;
	if (!CreateProcess(
			/*(ExtractFilePath(Application->ExeName)+*/ "MatrixProfile.exe" /*).c_str()*/,
			nullptr,
			nullptr,
			nullptr,
			FALSE,
			NORMAL_PRIORITY_CLASS,
			nullptr,
			/*ExtractFilePath(Application->ExeName).c_str()*/ nullptr,
			&si,
			&pi))
	{
		CloseHandle(pi.hProcess);
		CloseHandle(pi.hThread);
		throw tMessage("Cant launch state space search");
	}
	WaitForSingleObject(pi.hProcess, INFINITE);
	/*  char* Result(nullptr);
	try {
	  PlEngine Prolog(Application->ExeName.c_str());
	  PlTermv Parametrs(3);
	  const char* temp = RibsList.c_str();
	  Parametrs[0] = temp;
	  WriteConsole(hOut, (char*)Parametrs[0], RibsList.size(), &tmp, nullptr);
	  Parametrs[1] = 6;
	  int a = PlCall("hill_climbing", Parametrs);
	  temp = (char*)Parametrs[2];
	  unsigned int len = strlen(temp);
	  Result = new char[len];
	  strcpy(Result, temp);
  #ifdef STRONGCHECK
	  WriteConsole(hOut, Result, len, &tmp, nullptr);
	  FreeConsole();
  #endif
	}
	catch (PlException& Msg)
	{
	  delete[] Result;
	  throw tMessage((char*)Msg);
	}*/
	//  string Result;
	ifstream NewNodesList(/*(ExtractFilePath(Application->ExeName)+*/ "result.txt" /*).c_str()*/);
	Assert(!!NewNodesList, "Can't open file with state space search results");
	string Result((istreambuf_iterator<char>(NewNodesList)), istreambuf_iterator<char>());
	//  copy(istream_iterator<char>(NewNodesList), istream_iterator<char>(), back_inserter(Result));
	/*#ifdef STRONGCHECK
	  WriteConsole(hOut, Result.c_str(), Result.size(), &tmp, nullptr);
	  FreeConsole();
	#endif*/
	tMarkedGraph NewGraph(*this, Result);
	*this = NewGraph;
#endif	// def __WIN32__  ||  _WIN32
	return *this;
}
