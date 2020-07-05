
#include "NonLAnalys.hpp"
#include "LinearElasticAnalysis.hpp"

struct fCreateNewAnalysis
{
	void operator()(t1StepAnalysis*& p_) const
	{
		p_ = new tLinearElasticAnalysis;
	}
};

PhysNonlinearAnalysis::PhysNonlinearAnalysis(size_t steps_num_) : MultiStepAnalysis(steps_num_)
{
	for_each(m_pSteps.begin(), m_pSteps.end(), fCreateNewAnalysis());
	/* vector<t1StepAnalysis*>::iterator ppCurrentIteration = pSteps.begin();
	 for(;ppCurrentIteration != pSteps.end(); ++ppCurrentIteration)
		*ppCurrentIteration = new tLinearElasticAnalysis;*/
}

PhysNonlinearAnalysis::~PhysNonlinearAnalysis()
{
	for_each(m_pSteps.begin(), m_pSteps.end(), fDelete_object<t1StepAnalysis>());
	/*  vector<t1StepAnalysis*>::iterator ppCurrentIteration = pSteps.begin();
	  for(;ppCurrentIteration != pSteps.end(); ++ppCurrentIteration)
		 delete *ppCurrentIteration;*/
}

void PhysNonlinearAnalysis::SetStepsNum(size_t steps_num_)
{
	MultiStepAnalysis::SetStepsNum(steps_num_);
	for_each(m_pSteps.begin(), m_pSteps.end(), fCreateNewAnalysis());
	/* vector<t1StepAnalysis*>::iterator ppCurrentIteration = pSteps.begin();
	 for(;ppCurrentIteration != pSteps.end(); ++ppCurrentIteration)
		*ppCurrentIteration = new tLinearElasticAnalysis;*/
	m_pSteps.front()
		->SetInitialDisplacement(m_pInitialDisplLevel)
		.SetInitialLoad(m_pInitialLoadLevel);
}
