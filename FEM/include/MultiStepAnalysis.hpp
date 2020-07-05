#pragma once

#include "Analysis.hpp"
#include <vector>

class MultiStepAnalysis : public tAnalysis
{
public:
	//   MultiStepAnalysis(const tFE_model& model_, const tLoadCase& loadCase_,
	//   tFinitElement::tEvalMethod evalMethod_, tMarkedGraph::tMarkMethod
	//   minBand_=tMarkedGraph::minimal_degree, real_t stepsnum_, const tNodalTensor1Column*
	//   initApprox_=nullptr);
	MultiStepAnalysis(size_t iternum_ = 0) : tAnalysis(), m_pSteps(iternum_, nullptr)
	{
	}
	virtual ~MultiStepAnalysis()
	{
	}
	tAnalysis& Link(const tFE_model&) override;
	tAnalysis& Link(const tLoadCase&) override;
	virtual void SetStepsNum(size_t);
	// Solving equation system and calculation of results:
	size_t HowManySteps() const override
	{
		return m_pSteps.size();
	}
	const t1StepAnalysis& Step(size_t) const override;
	tAnalysis& Run() override;
	// Show results and other values:
	const tNodalTensor1Column& Load() const override
	{
		return m_pSteps.back()->Load();
	}
	const tNodalTensor1Column& Displacements() const override
	{
		return m_pSteps.back()->Displacements();
	}
	const tAnalysis& WriteResult(tOutput_XML&) const override
	{
		return *this;
	}  //заглушка!!

protected:
	std::vector<t1StepAnalysis*> m_pSteps;

private:
	//   std::vector<tAnalysis::tStepResult> Results;
	//   std::vector<tAnalysis::tStepResult>::iterator pCurrentStepResult;
};
