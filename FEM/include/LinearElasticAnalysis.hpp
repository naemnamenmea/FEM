#pragma once

#include "1StepAnalysis.hpp"

class tLinearElasticAnalysis : public t1StepAnalysis
{
public:
	static tAnalysis* NewAnalysis()
	{
		return new tLinearElasticAnalysis();
	}

	class fIntegrand : public tFinitElement::fIntegrand<Tensor1s, Tensor2a>
	{
	public:
		fIntegrand(const tLinearElasticAnalysis& anal_, const tFinitElement& fe_)
			: m_analysis(anal_), m_FE(fe_)
		{
		}
		fIntegrand& LinkNodes(const tNodalTensor2& o_)
		{
			m_result.Link(o_.Node(1), o_.Node(2));
			return *this;
		}
		virtual ~fIntegrand()
		{
		}

		Tensor2a& operator()(const Tensor1s& lc_) const
		{
			//          FE.JacobyMatrix_l2g(lc_,jacobyMatrix).Invert();
			(m_jacobyMatrix = m_FE.JacobyMatrix_l2g(lc_)).Invert();
			m_FE.DisplShapeGrad(m_result.Node(1), lc_, m_grad1);
			m_FE.DisplShapeGrad(m_result.Node(2), lc_, m_grad2) *= m_FE.Jacobian_l2g(lc_);
			m_grad1.DotMultiply_left(m_jacobyMatrix);
			m_grad2.DotMultiply_left(m_jacobyMatrix);
			m_grad1.DirectProduct(m_grad2, m_result);
			if (!m_FE.Material().isIsotropic())
				m_result.Dot2Multiply(
					m_FE.MaterialElasTensor(lc_),
					2,
					4);	 // додати StrainLevel  // нема повороту до локальних восей
			return m_result;
		}

	private:
		const tLinearElasticAnalysis& m_analysis;
		const tFinitElement& m_FE;
		mutable tNodalTensor2 m_result;
		mutable Tensor2s m_jacobyMatrix;
		mutable Tensor1s m_grad1, m_grad2;
	};

	class fStrain_g_twofold : public tFinitElement::fIntegrand<Tensor1s, SymmetricTensor2s>
	{
	public:
		fStrain_g_twofold(const tFinitElement& fe_) : m_FE(fe_)
		{
		}
		SymmetricTensor2s& operator()(const Tensor1s& lc_) const
		{
			return m_FE.NablaDispl_g(lc_, m_nablaDispl).Symmetrization_twofold(m_result);
		}

	private:
		mutable Tensor2s m_nablaDispl;
		mutable SymmetricTensor2s m_result;
		const tFinitElement& m_FE;
	};

	struct tStateInfo
	{
		tNodalTensor1Column m_displacements;
		std::vector<SymmetricTensor2s> m_linearStrains;
		std::vector<real_t> m_linearStrainsVonMises;
	};

	tLinearElasticAnalysis();
	tAnalysis& Link(const tFE_model& model_) override
	{
		m_pFEModel = &model_;
		m_stiffnessMatrix.Link(model_);
		m_loadsColumn.Link(model_);
		m_stateInfo.m_displacements.Link(model_);
		m_finalLoadLevel.Link(model_);
		return *this;
	}
	tAnalysis& Link(const tLoadCase& loadCase_) override
	{
		return tAnalysis::Link(loadCase_);
	}
	//   virtual tFEsTensor2& Strain(const Tensor1s&, tFEsTensor2&/* result_=tFEsTensor2()*/) const;
	//   virtual tFEsTensor2& Stress(const Tensor1s&, tFEsTensor2&, tFEsTensor2&/*
	//   result_=tFEsTensor2()*/) const; virtual const tFEsTensor2Column&
	//   Strains(tFEsSetOfTensor2Column::tTensor2Kind&) const; virtual const tFEsTensor2Column&
	//   Stresses(tFEsSetOfTensor2Column::tTensor2Kind&) const;
	t1StepAnalysis& AssembleAllMatrices() override;
	size_t MatrixDimension() const override
	{
		return m_dimension;
	}
	const tAnalysis& ProfileAndMaxBandWidth(size_t& pr_, size_t& max_) const override
	{
		pr_ = m_profile;
		max_ = m_maxbandwidth;
		return *this;
	}
	tAnalysis& Run() override;
	const tNodalTensor2SymMatrix& Matrix(size_t = 1) const override
	{
		return m_stiffnessMatrix;
	}

	SymmetricTensor2s& Strain_g(
		const Tensor1s& lc_, const tFinitElement& fe_, SymmetricTensor2s& result_) const
	{
		return (result_ = fStrain_g_twofold(fe_)(lc_)) *= 0.5;
	}

	SymmetricTensor2s& Strain_g_mean(const tFinitElement& fe_, SymmetricTensor2s& result_) const
	{
		fe_.Integrate(fStrain_g_twofold(fe_), result_.Assign0());
		return result_ *= (0.5 / fe_.Volume());
	}

	t1StepAnalysis& EvalState() override;

	const tAnalysis& WriteResult(tOutput_XML&) const override;
	const tNodalTensor1Column& Displacements() const override
	{
		return m_stateInfo.m_displacements;
	}

	tStateInfo m_stateInfo;

private:
	static const std::string& KindName;

	tNodalTensor2SymMatrix m_stiffnessMatrix;
	tNodalTensor1Column m_loadsColumn;
	size_t m_dimension;
	size_t m_profile;
	size_t m_maxbandwidth;
};
