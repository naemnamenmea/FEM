#pragma once

#include "VectorMatrix.hpp"
#include "GraphTheory.hpp"

class tAlgEqSystem
{
};

class tLinAlgEqSystem : public tAlgEqSystem
{
public:
	tLinAlgEqSystem(
		const tMarkedGraph*,
		const tNodalTensor2SymMatrix&,
		const tNodalTensor1Column& /*, tMarkedGraph::tMarkMethod=tMarkedGraph::minimal_degree*/);
	tNodalTensor1Column& Solve(tNodalTensor1Column&);
	size_t Dimension()
	{
		return m_leftSide.Dimension();
	}
	void ProfileAndMaxBandWidth(size_t& p_, size_t& m_)
	{
		m_leftSide.ProfileAndMaxBandWidth(p_, m_);
	}

private:
	const tMarkedGraph* m_pGraph;
	tSymConstRefMatrix m_leftSide;
	tConstRefColumn m_rightSide;
};

inline tLinAlgEqSystem::tLinAlgEqSystem(
	const tMarkedGraph* pgraph_,
	const tNodalTensor2SymMatrix& left_,
	const tNodalTensor1Column& right_ /*, tMarkedGraph::tMarkMethod minimizeBand_*/)
	: m_pGraph(pgraph_),
	  m_leftSide(left_, *m_pGraph /*.ReMark(minimizeBand_)*/),
	  m_rightSide(right_, *m_pGraph)
{
}

inline tNodalTensor1Column& tLinAlgEqSystem::Solve(tNodalTensor1Column& tensorSolution_)
{
	tRefTensor1ComponentsColumn solution(tensorSolution_, *m_pGraph);
	m_leftSide.SolveLinAlgSystem(m_rightSide, solution);
	return tensorSolution_;
}
