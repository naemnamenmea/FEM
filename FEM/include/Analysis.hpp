#pragma once

#include "Material.hpp"
#include "FiniteElementPool.hpp"
#include "VectorMatrix.hpp"
#include "FileInOut.hpp"
#include "NameList.hpp"
#include "GraphTheory.hpp"
#include "femDllDefinitions.hpp"

#include <vector>
#include <map>
#include <tuple>

class t1StepAnalysis;
class tLoadCase;
class tLoads;

class tAnalysis
{
public:
	enum tEvalPerFEMethod
	{
		eAtFECentre,
		eMeanOverFE
	} EvalPerFEMethod;

	static tAnalysis* NewAnalysis(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}
	//   tAnalysis(const tFE_model&, const tLoadCase&, tFinitElement::tEvalMethod,
	//   tMarkedGraph::tMarkMethod=tMarkedGraph::minimal_degree);
	tAnalysis();
	virtual ~tAnalysis()
	{
	}

	// Change options:
	tAnalysis& SetName(const char* name_)
	{
		m_caption = name_;
		return *this;
	}
	const std::string& Name() const
	{
		return m_caption;
	}
	virtual tAnalysis& Link(const tFE_model& model_)
	{
		m_pFEModel = &model_;
		return *this;
	}
	virtual tAnalysis& Link(const tLoadCase& loadCase_)
	{
		m_pLoadCase = &loadCase_;
		return *this;
	}
	tAnalysis& Link(
		const tMarkedGraph& og_, tMarkedGraph* pmg_ = nullptr, tMarkedGraph* ppg_ = nullptr)
	{
		m_pOrigGraph = m_pActiveGraph = &og_;
		m_pMinimizedGraph = pmg_;
		m_pPrologMinGraph = ppg_;
		return *this;
	}
	tAnalysis& SetMinimizeBandMethod(tMarkedGraph::tMarkMethod method_)
	{
		m_minimizeBandMethod = method_;
		return *this;
	}
	tAnalysis& SetInitialDisplacement(const tNodalTensor1Column* pdispls_)
	{
		m_pInitialDisplLevel = pdispls_;
		return *this;
	}
	tAnalysis& SetInitialLoad(const tNodalTensor1Column* ploads_)
	{
		m_pInitialLoadLevel = ploads_;
		return *this;
	}
	//   tAnalysis& SetInitialTime(real_t newTime_) {InitialTime = newTime_; return *this;}
	tAnalysis& SetFinalTime(real_t newTime_)
	{
		m_finalTime = newTime_;
		return *this;
	}

	// Solving equation system and calculation of results:
	virtual tAnalysis& Run();
	bool Completed() const
	{
		return m_solutionFound;
	}
	virtual size_t HowManySteps() const = 0;
	virtual const t1StepAnalysis& Step(size_t) const = 0;

	virtual const tAnalysis& WriteResult(tOutput_XML&) const = 0;
	const tAnalysis& WriteResultXML() const
	{
		tOutput_XML result_file(m_outputFileName.c_str());
		WriteResult(result_file);
		return *this;
	}
	tAnalysis& SetOutputFileName(const char* filename_)
	{
		m_outputFileName = filename_;
		return *this;
	}
	//   virtual tFEsTensor2& Strain(const Tensor1s&,tFEsTensor2&) const =0;
	//   virtual tFEsTensor2& Stress(const Tensor1s&,tFEsTensor2&,tFEsTensor2&) const =0;
	//   virtual const tFEsTensor2Column& Strains (tFEsSetOfTensor2Column::tTensor2Kind&) const =0;
	//   virtual const tFEsTensor2Column& Stresses(tFEsSetOfTensor2Column::tTensor2Kind&) const =0;

	/*   class fLocalStrainStress: public tFinitElement::fIntegrand<tFEsSetOfTensor2>
	   {
		private:
		  mutable tFEsSetOfTensor2 Result;
		  const tAnalysis& tAnalysis;
		  size_t StrainNo, StressNo;
		public:
		  fLocalStrainStress(const tFinitElement& fe_,const tAnalysis& anal_, size_t k1_,
	   size_t k2_): Result(fe_,2),tAnalysis(anal_),StrainNo(k1_),StressNo(k2_) {}
		  tFEsSetOfTensor2& operator()(const Tensor1s& lc_) const
				 {*/
	// /*!!!*/              tFEsTensor2/*&*/ strain = Result(StrainNo),  stress = Result(StressNo);
	//              tAnalysis.Strain(lc_,strain);
	//              tAnalysis.Stress(lc_,strain,/*Result(StressNo)*/stress);
	/*              return Result;
				 }
	   };
	   class fFEsStrainStress: public tFEsSetOfTensor2Column::fTransformer
	   {
		private:
		  const tAnalysis& tAnalysis;
		  size_t Tensor1No, Tensor2No;
		public:
		  fFEsStrainStress(const tAnalysis& anal_, size_t k1_, size_t k2_):
	 tAnalysis(anal_),Tensor1No(k1_),Tensor2No(k2_){} void operator()(tFEsSetOfTensor2& set_) const
				 {set_.FE().Calculate(fLocalStrainStress(set_.FE(),tAnalysis,Tensor1No,Tensor2No),tAnalysis.ValuesEvalMethod,set_);}
	   };
	 friend class fFEsStrainStress;*/
	/*class RotateToFELocalAxes: public tFEsTensor2Column::fTransformer
	{
	 private:
	   const Tensor1s& LocalCoord;
	 public:
	   RotateToFELocalAxes(const Tensor1s& lc_): LocalCoord(lc_){}
	   void operator()(tFEsTensor2& festensor_) const
		   {
			if (festensor_.is0()) return;
			Tensor2s rotationTensor;
			festensor_.FE().TensorOfRotationToLocal(LocalCoord,rotationTensor);
			festensor_.ScalarMultiply(rotationTensor).ScalarMultiply_left(rotationTensor.Transpose());
		   }
	};*/
	// Show results and other values:
	virtual const tAnalysis& ProfileAndMaxBandWidth(size_t&, size_t&) const = 0;
	virtual size_t MatrixDimension() const = 0;
	virtual const tNodalTensor1Column& Load() const = 0;
	//   virtual const tNodalTensor2SymMatrix& Matrix(size_t=1) const =0;//changes are required
	virtual const tNodalTensor1Column& Displacements() const = 0;
	//   /*virtual*/const tNodalTensor1Column& Displacements(const std::string&) const {return
	//   Displacements();}
	//   /*virtual*/const tNodalTensor1Column& Displacements(size_t) const {return
	//   Displacements();}

	unsigned int m_searchDepth;

protected:
	static const tNamesToFuns<tAnalysis> Factory;

	const tFE_model* m_pFEModel;
	const tLoadCase* m_pLoadCase;
	tMarkedGraph::tMarkMethod m_minimizeBandMethod;
	const tNodalTensor1Column *m_pInitialDisplLevel, *m_pInitialLoadLevel;
	/*const*/ real_t /*InitialTime,*/ m_finalTime;
	const tMarkedGraph *m_pOrigGraph, *m_pActiveGraph;
	tMarkedGraph *m_pMinimizedGraph, *m_pPrologMinGraph;
	bool m_solutionFound;

private:
	std::string m_caption;
	std::string m_outputFileName;
};

class FEM_API AnalysesArray
{
public:
	enum whatgraph
	{
		Original,
		Minimize,
		PrologMin,
		OnlyPrologMin
	};

	//   AnalysesArray(): Analyses(), pCurrent(nullptr) {};
	AnalysesArray() : m_pCurrent(nullptr), m_graphkind(Original), m_indexOfCurrGraph(0)
	{
	}
	~AnalysesArray();
	size_t HowManyAnalyses() const
	{
		return m_analyses.size();
	}
	//   void Push_back(tAnalysis* pvalue_) {Analyses.insert(make_pair(&(pvalue_->Name()),
	//   pvalue_));}
	void ImportFromFile(const char*, const tFE_model&, const tLoads&);
	void Delete(const char*);
	void DeleteCurrent();
	bool isEmpty() const
	{
		return m_analyses.empty();
	}
	tAnalysis& operator()(const char*) const;
	tAnalysis& operator()(int) const;
	tAnalysis& SetCurrent(const char* name_)
	{
		m_pCurrent = &operator()(name_);
		return *m_pCurrent;
	}
	tAnalysis& GetCurrent() const
	{
		return *m_pCurrent;
	}
	void RunCurrent();
	const tMarkedGraph* GetPtrToOrigGraph() const
	{
		return &m_origGraph;
	}
	const tMarkedGraph& GetBestGraph()
	{
		return m_prologMinGraph.isEmpty() ? m_minimizedGraph : m_prologMinGraph;
	}
	void SetGraphKind(whatgraph kind_)
	{
		m_graphkind = kind_;
	}
	//   void SetCurrGraph(tMarkedGraph* pgraph_, unsigned int index_) {pGraphs[index_] = pgraph_;}
	//   tAnalysis& SetTypeCurrent(const char*);
	void Clear();

private:
#pragma warning(push)
#pragma warning(disable : 4251)
	friend struct has_assembled_matrix;
	friend struct fDelete_analysis;

	class PtrToAnalysis
	{
	public:
		PtrToAnalysis(tAnalysis* panalis_) : m_pAnalysis(panalis_)
		{
		}
		PtrToAnalysis(const AnalysesArray::PtrToAnalysis& orig_) : m_pAnalysis(orig_.m_pAnalysis)
		{
		}
		// ~PtrToAnalysis() {delete pAnalysis;}
		bool operator<(const AnalysesArray::PtrToAnalysis& anltocmp_) const
		{
			return m_pAnalysis->Name() < anltocmp_.m_pAnalysis->Name();
		}
		bool operator==(const AnalysesArray::PtrToAnalysis& anltocmp_) const
		{
			return m_pAnalysis->Name() == anltocmp_.m_pAnalysis->Name();
		}
		PtrToAnalysis& operator=(const PtrToAnalysis& orig_)
		{
			m_pAnalysis = orig_.m_pAnalysis;
			return *this;
		}

		tAnalysis* m_pAnalysis;
	};

	void SetInitialGraph();

	std::set<PtrToAnalysis> m_analyses;
	//   std::vector<tMarkedGraph*> pGraphs;
	tMarkedGraph m_origGraph, m_minimizedGraph, m_prologMinGraph;
	tAnalysis* m_pCurrent;
	unsigned int m_indexOfCurrGraph;
	whatgraph m_graphkind;
#pragma warning(pop)
};
