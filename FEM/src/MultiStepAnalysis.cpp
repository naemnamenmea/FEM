#include "1StepAnalysis.hpp"
#include "MultiStepAnalysis.hpp"
#include "ProcessDisplay.hpp"

using namespace std;

template <typename T>
struct tLink
{
	tLink(const T& o_) : m_what(o_)
	{
	}
	void operator()(t1StepAnalysis* p_) const
	{
		p_->Link(m_what);
	}

	const T& m_what;
};

tAnalysis& MultiStepAnalysis::Link(const tFE_model& model_)
{
	tAnalysis::Link(model_);
	for_each(m_pSteps.begin(), m_pSteps.end(), tLink<tFE_model>(model_));
	/* for (vector<t1StepAnalysis*>::iterator pp=pSteps.begin(); pp!=pSteps.end(); ++pp)
		(*pp)->Link(model_);*/
	return *this;
}

tAnalysis& MultiStepAnalysis::Link(const tLoadCase& loadCase_)
{
	tAnalysis::Link(loadCase_);
	for_each(m_pSteps.begin(), m_pSteps.end(), tLink<tLoadCase>(loadCase_));
	/* for (vector<t1StepAnalysis*>::iterator pp=pSteps.begin(); pp!=pSteps.end(); ++pp)
		(*pp)->Link(loadCase_);*/
	return *this;
}

void MultiStepAnalysis::SetStepsNum(size_t steps_num_)
{
	for_each(m_pSteps.begin(), m_pSteps.end(), fDelete_object<t1StepAnalysis>());
	m_pSteps.resize(steps_num_, nullptr);
}

tAnalysis& MultiStepAnalysis::Run()
{
	tAnalysis::Run();
	real_t time_step = 1. / HowManySteps(), cur_time(time_step);
	tProcessDisplay output("Load steps execution", static_cast<data_t>(HowManySteps()));
	vector<t1StepAnalysis*>::iterator ppcurStep = m_pSteps.begin();
	(*ppcurStep)->SetFinalTime(cur_time);
	// pFE_model->LinkDisplacements(pInitialDisplLevel);
	(*ppcurStep)->AssembleAllMatrices();
	(*ppcurStep)->Link(*m_pActiveGraph).Run();
	for (const tNodalTensor1Column *pPrevDispl, *pPrevloads;;)
	{
		cur_time += time_step;
		pPrevDispl = &(*ppcurStep)->Displacements();
		pPrevloads = &(*ppcurStep)->Load();
		++ppcurStep;
		++output;
		if (ppcurStep == m_pSteps.end())
			break;
		(*ppcurStep)->Link(*m_pActiveGraph).SetFinalTime(cur_time);
		(*ppcurStep)->SetInitialLoad(pPrevloads).SetInitialDisplacement(pPrevDispl).Run();
	}
	m_solutionFound = true;
	return *this;
}

const t1StepAnalysis& MultiStepAnalysis::Step(size_t stepnum_) const
{
#ifdef STRONGCHECK
	Assert(
		stepnum_ > 0 && stepnum_ <= HowManySteps(), tMessage("Invalid step number = ") << stepnum_);
#endif	// def STRONGCHECK
	return *m_pSteps[--stepnum_];
}
