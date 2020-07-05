#include "LinearElasticAnalysis.hpp"
#include "FiniteElementModel.hpp"
#include "Analysis.hpp"
#include "ProcessDisplay.hpp"
#include "LinAlgEq.hpp"
#include <string>
#include <ctime>

using namespace std;

const string& tLinearElasticAnalysis::KindName =
	tAnalysis::Factory.FindName(tLinearElasticAnalysis::NewAnalysis);

t1StepAnalysis& tLinearElasticAnalysis::AssembleAllMatrices()
{
	t1StepAnalysis::AssembleAllMatrices();

	tProcessDisplay output(
		"Stiffness matrix assembling", static_cast<data_t>(m_pFEModel->HowManyFinEls()));
	const tFinitElement* pcurrentFE;
	const tNode* nodei;
	tNodalTensor2 currentFEstiffnessMatrixComponent;

	m_stiffnessMatrix.Assign0();

	const Tensor1s centreOfFE(0.);
	const SymmetricTensor4s* pcurrentElasTensor = 0;

	ofstream out("time_check.txt");	 //=========================!!!!!!!!!!!!!!!!!
	time_t time0 = time(nullptr);
	clock_t clock0 = clock();
	for (size_t numOfDisplNodes, i, j, k = 1; k <= m_pFEModel->HowManyFinEls(); ++k)
	{
		pcurrentFE = &m_pFEModel->GetFE(k);

		if (pcurrentFE->Material().isIsotropic())
		{
			pcurrentElasTensor = &pcurrentFE->MaterialElasTensor(centreOfFE);
		}

		numOfDisplNodes = pcurrentFE->HowManyDisplApproxNodes();
		fIntegrand subintegralFunction(*this, *pcurrentFE);
		for (i = 1; i <= numOfDisplNodes; ++i)
		{
			nodei = &pcurrentFE->DisplApproxNode(i);
			for (j = 1; j <= i;
				 ++j)  // j<=i (!!! not j<=numOfDisplNodes) because StiffnessMatrix is symmetric
			{
				currentFEstiffnessMatrixComponent.Link(*nodei, pcurrentFE->DisplApproxNode(j));
				subintegralFunction.LinkNodes(currentFEstiffnessMatrixComponent);
				if (!currentFEstiffnessMatrixComponent.is0())
					currentFEstiffnessMatrixComponent.Assign0toData();
				pcurrentFE->Integrate_local(subintegralFunction, currentFEstiffnessMatrixComponent);
#ifdef STRONGCHECK
				Assert(!currentFEstiffnessMatrixComponent.is0(), "zero stiffness matrix component");
#endif

				if (pcurrentFE->Material().isIsotropic())
				{
					currentFEstiffnessMatrixComponent.Dot2Multiply(*pcurrentElasTensor, 2, 4);
				}

				m_stiffnessMatrix += currentFEstiffnessMatrixComponent;
			}
		}
		++output;
	}
	out << time(nullptr) - time0 << " = ";
	out << double(clock() - clock0) / CLK_TCK << " secs" << endl;
	// out << double(clock()-clock0)/CLOCKS_PER_SEC << endl;
	out.close();  //=========================!!!!!!!!!!!!!!!!!
	// exit(EXIT_SUCCESS);		 //=========================!!!!!!!!!!!!!!!!!
	// throw tMessage("forced program termination after matrix
	// assemblage");//=========================!!!!!!!!!!!!!!!!!
	out.open("matrix.txt");	 //=========================!!!!!!!!!!!!!!!!!
	out << m_stiffnessMatrix.NumberOfDOFs() << ' ' << m_stiffnessMatrix.Dimension() << '\n';
	/*
	 output.Reset("Stiffness matrix output to file",pFE_model->HowManyNodes());
	 for (size_t i=1, j, k, l; i<=StiffnessMatrix.Dimension(); ++i,++output)
	 {
		const tNode &Nodei = pFE_model->Node(i), *pNodej;
		for (k=1; k<=3; ++k)
		 if (!Nodei.HasFixed(k))
		 {
		for (j=1; j<i; ++j)
		 {
			pNodej = &pFE_model->Node(j);
			if (StiffnessMatrix.is0(i,j))
			 {
				 for (l=1; l<=3; ++l)
				if (!pNodej->HasFixed(l))
				 out << "           0\t";
			 }
			 else
				 for (l=1; l<=3; ++l)
				if (!pNodej->HasFixed(l))
				 out << (StiffnessMatrix(i,j).is0(k,l)? 0.:
						 StiffnessMatrix(i,j)(k,l)) << '\t';
		 }
			// j==i :
			if (StiffnessMatrix.is0(i,i))
			 {
				 for (l=1; l<=k; ++l)
				if (!Nodei.HasFixed(l))
				 out << "           0\t";
			 }
			 else
				 for (l=1; l<=k; ++l)
				if (!Nodei.HasFixed(l))
				 out << (StiffnessMatrix(i,i).is0(k,l)? 0.:
						 StiffnessMatrix(i,i)(k,l)) << '\t';
		out << '\n';
		 }
	 }
	*/
	out.close();  //=========================!!!!!!!!!!!!!!!!!
	output.Reset("tLoads column setting", 2);

	// pLoadCase->LoadLevel(LoadsColumn, FinalTime);

	m_loadsColumn = m_finalLoadLevel;
	++output;
	if (m_pInitialLoadLevel != nullptr)
		m_loadsColumn -= *m_pInitialLoadLevel;
	m_matricesAssembled = true;
	++output;
	return *this;
}

tAnalysis& tLinearElasticAnalysis::Run()
{
	t1StepAnalysis::Run();
	// tProcessDisplay output("Equation system forming",2);
	tLinAlgEqSystem eqSystem(
		m_pActiveGraph, m_stiffnessMatrix, m_loadsColumn /*,MinimizeBandMethod*/);
	// output.Reset("Linear equation system solving",2) = 2;
	tProcessDisplay output("Linear equation system solving", 2);
	tNodalTensor1Column& displs_column(m_stateInfo.m_displacements);
	eqSystem.Solve(displs_column);
	++output;
	if (m_pInitialDisplLevel != nullptr)
		displs_column += *m_pInitialDisplLevel;
	m_dimension = eqSystem.Dimension();
	eqSystem.ProfileAndMaxBandWidth(m_profile, m_maxbandwidth);
	m_pFEModel->LinkDisplacements(&displs_column);
	output.Reset("Stress-strain state evaluating", 1);
	EvalState();
	m_solutionFound = true;
	++output;
	return *this;
}

t1StepAnalysis& tLinearElasticAnalysis::EvalState()
{
	tProcessDisplay output(
		"Evaluation of strains", static_cast<data_t>(m_pFEModel->HowManyFinEls()));
	m_stateInfo.m_linearStrains.resize(m_pFEModel->HowManyFinEls());
	vector<tFinitElement*>::const_iterator ppf = m_pFEModel->pFE_begin();
	vector<SymmetricTensor2s>::iterator ps = m_stateInfo.m_linearStrains.begin();
	switch (EvalPerFEMethod)
	{
		case eAtFECentre:
			for (const Tensor1s fe_centre(0.); ppf < m_pFEModel->pFE_end(); ++ppf, ++ps, ++output)
				Strain_g(fe_centre, **ppf, *ps);
			break;

		case eMeanOverFE:
			for (; ppf < m_pFEModel->pFE_end(); ++ppf, ++ps, ++output) Strain_g_mean(**ppf, *ps);
			break;
		default:
			THROW_MESSAGE("Unknown tEvalPerFEMethod");
	}
	return *this;
}

tLinearElasticAnalysis::tLinearElasticAnalysis()
	: t1StepAnalysis(),
	  m_stiffnessMatrix(),
	  m_loadsColumn(),
	  m_dimension(0),
	  m_profile(0),
	  m_maxbandwidth(0)
{
	m_stiffnessMatrix.SetName("Stiffness matrix");
	m_loadsColumn.SetName("tLoads");
}

const tAnalysis& tLinearElasticAnalysis::WriteResult(tOutput_XML& xml_out_) const
{
	if (!m_solutionFound)
		return *this;
	tProcessDisplay output(
		"Writing displacements", static_cast<data_t>(Displacements().Dimension()));
	if (Name().empty())
		xml_out_.OpenTag("results", false);
	else
	{
		xml_out_.OpenTag("results", true);
		xml_out_.AddAtr("id", Name().c_str(), true);
	}
	xml_out_.OpenTag("displacements", true);
	xml_out_.AddAtr("type", "nodal", true);
	for (size_t i = 1; i <= Displacements().Dimension(); ++i, ++output)
	{
		xml_out_.OpenTag("node", true);
		xml_out_.AddAtr("id", i);
		if (Displacements().is0(i))
			xml_out_.AddAtr("x", 0.).AddAtr("y", 0.).AddAtr("z", 0.);
		else
			xml_out_.AddAtr("x", Displacements()(i).is0(1) ? 0. : Displacements()(i)(1))
				.AddAtr("y", Displacements()(i).is0(2) ? 0. : Displacements()(i)(2))
				.AddAtr("z", Displacements()(i).is0(3) ? 0. : Displacements()(i)(3));
		xml_out_.CloseTag();
	}
	xml_out_.CloseTag("displacements");

	output.Reset("Writing strains", static_cast<data_t>(m_pFEModel->HowManyFinEls()));
	xml_out_.OpenTag("strains", true);
	xml_out_.AddAtr("type", "fe-wise");
	xml_out_.AddAtr("id", "linear strain");
	xml_out_.AddAtr("coord", "global", true);
	//         vector<SymmetricTensor2s>::const_iterator pstrain = Result.LinearStrains_g.begin();
	vector<SymmetricTensor2s>::const_iterator pstrain = m_stateInfo.m_linearStrains.begin();
	for (size_t i = 1; pstrain < m_stateInfo.m_linearStrains.end(); ++i, ++pstrain, ++output)
	{
		xml_out_.OpenTag("fe", true);
		xml_out_.AddAtr("id", i);
		xml_out_.AddAtr("x", (*pstrain)(1, 1));
		xml_out_.AddAtr("y", (*pstrain)(2, 2));
		xml_out_.AddAtr("z", (*pstrain)(3, 3));
		xml_out_.AddAtr("xy", (*pstrain)(1, 2));
		xml_out_.AddAtr("yz", (*pstrain)(2, 3));
		xml_out_.AddAtr("zx", (*pstrain)(3, 1));
		xml_out_.CloseTag();
	}
	xml_out_.CloseTag("strains");
	xml_out_.CloseTag("results");
	return *this;
}

/*tFEsTensor2& tLinearElasticAnalysis::Strain (const Tensor1s& lc_,
											 tFEsTensor2& result_) const
{
 result_.FE().DisplacementGradient(lc_,result_);
 tFEsTensor2 gradT = result_;
 gradT.Transpose();      result_ += gradT;     result_ *= 0.5;
 return result_;
}*/

/*tFEsTensor2& tLinearElasticAnalysis::Stress (const Tensor1s& localCoord_,
												 tFEsTensor2& strain_,
												 tFEsTensor2& result_) const
{
 if (strain_.is0())
		{
		 result_.Assign0();
		 return result_;
		}
 Tensor4 elastensor;    result_.FE().MaterialElasTensor(strain_,localCoord_,elastensor);
 elastensor.DoubleScalarProduct(strain_,result_);
 return result_;
}*/

/*const tFEsTensor2Column& tLinearElasticAnalysis::Strains (tFEsSetOfTensor2Column::tTensor2Kind&
kind_) const
{
 kind_.Assign(tFEsSetOfTensor2Column::linear_strain);
 tFEsTensor2Column& result = Tensor2Result(kind_);
 if (result.Name().empty())
	 {
	string name = "Linear strains";
	if (!kind_.BelongsToGlobalCoordSystem())
					 name += ", local";
	result.SetName(name);
	 }
 return result;
}

const tFEsTensor2Column& tLinearElasticAnalysis::Stresses (tFEsSetOfTensor2Column::tTensor2Kind&
kind_) const
{
 kind_.Assign(tFEsSetOfTensor2Column::Cauchy_stress);
 tFEsTensor2Column& result = Tensor2Result(kind_);
 if (result.Name().empty())
	 {
	string name = "Cauchy stresses";
	if (!kind_.BelongsToGlobalCoordSystem())
					 name += ", local";
	result.SetName(name);
	 }
 return result;
}*/
