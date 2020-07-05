#pragma once

#include "Analysis.hpp"
#include "1StepAnalysis.hpp"
#include "MultiStepAnalysis.hpp"

class PhysNonlinearAnalysis : public MultiStepAnalysis
{
public:
	static tAnalysis* NewAnalysis(/*unsigned int ItersNum*/)
	{
		return new PhysNonlinearAnalysis(/*ItersNum*/);
	}

	// PhysNonlinearAnalysis(unsigned int=0);
	virtual ~PhysNonlinearAnalysis();
	void SetStepsNum(size_t stepnum_) override;

	//   virtual tFEsTensor2& Strain(const Tensor1s& localCoord_, tFEsTensor2&
	//   res_/*=tFEsTensor2()*/) const
	//                               {return pSteps.front()->Strain(localCoord_, res_);}

	//   virtual tFEsTensor2& Stress(const Tensor1s& localCoord_, tFEsTensor2& strain_, tFEsTensor2&
	//   res_/*=tFEsTensor2()*/) const
	//                               {return pSteps.front()->Stress(localCoord_, strain_, res_);}

	//   virtual const tFEsTensor2Column& Strains(tFEsSetOfTensor2Column::tTensor2Kind& kind_) const
	//                               {return pSteps.front()->Strains(kind_);}

	//   virtual const tFEsTensor2Column& Stresses(tFEsSetOfTensor2Column::tTensor2Kind& kind_)
	//   const
	//                               {return pSteps.front()->Stresses(kind_);}

	const tAnalysis& ProfileAndMaxBandWidth(size_t& pr_, size_t& max_) const override
	{
		return m_pSteps.front()->ProfileAndMaxBandWidth(pr_, max_);
	}

	size_t MatrixDimension() const override
	{
		return m_pSteps.front()->MatrixDimension();
	}

	virtual const tNodalTensor2SymMatrix& Matrix(size_t = 1) const
	{
		return m_pSteps.front()->Matrix();
	}

private:
	static const std::string& KindName;

	explicit PhysNonlinearAnalysis(size_t = 0);
};
