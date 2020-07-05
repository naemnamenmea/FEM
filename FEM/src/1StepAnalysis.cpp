#include "1StepAnalysis.hpp"
#include "FiniteElementModel.hpp"
#include "Load.hpp"

/*void t1StepAnalysis::CalculateStrainAndStress (tFEsSetOfTensor2Column::tTensor2Kind strainKind_,
												 tFEsSetOfTensor2Column::tTensor2Kind stressKind_,
												 tFEsTensor2Column& strains_,
												 tFEsTensor2Column& stresses_) const
{
 tFEsSetOfTensor2Column strainAndStress(*pFE_model,2);
 strainAndStress.Include(strainKind_).Include(stressKind_);
 strainAndStress.Transform(fFEsStrainStress(*this,strainAndStress.ColNo(strainKind_),strainAndStress.ColNo(stressKind_)));
 strainAndStress.Swap(strainKind_,strains_);
 strainAndStress.Swap(stressKind_,stresses_);
}*/

/*tFEsTensor2Column& t1StepAnalysis::Tensor2Result(tFEsSetOfTensor2Column::tTensor2Kind kind_) const
{
 const tFEsTensor2Column tmp;
 pair<map<tFEsSetOfTensor2Column::tTensor2Kind,tFEsTensor2Column>::iterator,bool>
	 ppres = pResult->Tensor2Values.insert(make_pair(kind_,tmp));
 tFEsTensor2Column* presult = &((ppres.first)->second);
 if (ppres.second)// i.e. if insertion has been successfull (i.e. this column did not exist before)
	 {
	if (kind_.BelongsToGlobalCoordSystem())
		{
		 presult->Link(*pFE_model);
//       if (pFE_model->??????)   pFE_model->LinkDisplacements(pResult->Displs);
		 tFEsSetOfTensor2Column::tTensor2Kind kind2 = kind_.ConjByMatLow();
		 ppres = pResult->Tensor2Values.insert(make_pair(kind2,tmp));
		 tFEsTensor2Column& conjugate = (ppres.first)->second;
		 if (ppres.second) conjugate.Link(*pFE_model);
		 if (kind_.isStrain()) CalculateStrainAndStress(kind_,kind2,*presult,conjugate);
		 else                CalculateStrainAndStress(kind2,kind_,conjugate,*presult);
		}
	else
		{
		 presult = &(pResult->Tensor2Values[kind_] = Tensor2Result(kind_.ConjByRotation()));
		 const Tensor1s centreOfFE(0.);
		 presult->Transform(RotateToFELocalAxes(centreOfFE));
		}
	 }
 return *presult;
}*/

t1StepAnalysis& t1StepAnalysis::AssembleAllMatrices()
{
	m_pLoadCase->LoadLevel(m_finalLoadLevel, m_finalTime);
	m_pFEModel->LinkDisplacements(m_pInitialDisplLevel);
	return *this;
}
