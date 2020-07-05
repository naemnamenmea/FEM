#pragma once

#include "Analysis.hpp"

struct tStateInfo
{  // state of structure at some time
	tNodalTensor1Column m_displs;
	std::vector<SymmetricTensor2s> m_linearStrains_g;

	//                       Tuple<tNodalTensor1Column,std::vector<SymmetricTensor2s>,std::vector<real_t>
	//                       > mutable tFEsSetOfTensor2Column Tensor2Values; tStepResult(const
	//                       tFE_model& model_): Displs(model_), Tensor2Values(model_) {}
	//                       tStepResult(const tStepResult& o_): Displs(o_.Displs),
	//                       Tensor2Values(o_.Tensor2Values) {} mutable
	//                       std::map<tFEsSetOfTensor2Column::tTensor2Kind,tFEsTensor2Column>
	//                       Tensor2Values; tStepResult(): Tensor2Values() {} tStepResult(const
	//                       tFE_model& model_): Displs(model_), Tensor2Values() {}
	//                       tStepResult(const tStepResult& orig_):
	//                       Tensor2Values(orig_.Tensor2Values) {}
	//                     tStepResult& Link(const tFE_model& model_) {Displs.Link(model_); return
	//                     *this;}
};

class t1StepAnalysis : public tAnalysis
{
public:
	//   t1StepAnalysis(const tFE_model& model_, const tLoadCase& loadCase_,
	//   tFinitElement::tEvalMethod evalMethod_, tMarkedGraph::tMarkMethod
	//   minBand_=tMarkedGraph::minimal_degree,tAnalysis::tStepResult* pres_=nullptr):
	//   tAnalysis(model_,loadCase_,evalMethod_,minBand_),pResult(pres_){}
	t1StepAnalysis() : m_matricesAssembled(false)
	{
	}
	virtual ~t1StepAnalysis()
	{
	}
	// Assembling equation system:
	/*   class fIntegrand: public tFinitElement::fIntegrand<Tensor2a>
			{
			 private:
			 mutable tNodalTensor2 Result;
			 const t1StepAnalysis& tAnalysis;
			 const tFinitElement& FE;
			 public:
			 fIntegrand(const t1StepAnalysis& anal_, const tFinitElement& fe_): tAnalysis(anal_),
		 FE(fe_) {} fIntegrand& LinkNodes(const tNodalTensor2& o_)
		 {Result.Link(o_.Node(1),o_.Node(2)); return *this;} Tensor2a& operator()(const Tensor1s&
	   lc_) const {return tAnalysis.Intergrand(FE,lc_,Result);}
			};*/
	virtual t1StepAnalysis& AssembleAllMatrices();
	//   virtual tNodalTensor2& Intergrand(const tFinitElement&, const Tensor1s&, tNodalTensor2&)
	//   const =0;
	// Solving equation system and calculation of results:
	// Show results and other values:
	size_t HowManySteps() const override
	{
		return 1;
	}
	const t1StepAnalysis& Step(size_t) const override
	{
		return *this;
	}
	//   virtual const tNodalTensor1Column& Displacements() const {return Result.Displs;}
	tAnalysis& Run() override
	{
		tAnalysis::Run();
		if (!m_matricesAssembled)
			AssembleAllMatrices();
		return *this;
	}
	virtual const tNodalTensor2SymMatrix& Matrix(size_t = 1) const = 0;
	const tNodalTensor1Column& Load() const override
	{
		return m_finalLoadLevel;
	}
	virtual t1StepAnalysis& EvalState() = 0;

protected:
	tNodalTensor1Column m_finalLoadLevel;
	//   tStateInfo Result;
	bool m_matricesAssembled;
	//   tFEsTensor2Column& Tensor2Result(tFEsSetOfTensor2Column::tTensor2Kind) const;

private:
	//   void CalculateStrainAndStress
	//   (tFEsSetOfTensor2Column::tTensor2Kind,tFEsSetOfTensor2Column::tTensor2Kind,tFEsTensor2Column&,tFEsTensor2Column&)
	//   const;
};
