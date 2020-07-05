#pragma once

#include "FiniteElementModel.hpp"
#include "Load.hpp"
#include "Analysis.hpp"
#include <fstream>

struct tStructure
{
	tStructure() : m_FEModel(), m_loads(), m_analyses()
	{
	}
	tStructure(const char* fileName_)
	{
		ImportFromFile(fileName_);
	}
	//   ~tStructure() {}

	tStructure& ImportFromFile(const char*);
	tStructure& Clear();
	void Dump(const char*) const;

	tFE_model m_FEModel;
	tLoads m_loads;
	AnalysesArray m_analyses;
};

inline tStructure& tStructure::ImportFromFile(const char* fileName_)
{
	m_FEModel.ImportFromFile(fileName_);
	m_loads.ImportFromFile(fileName_, m_FEModel);
	m_analyses.ImportFromFile(fileName_, m_FEModel, m_loads);
	return *this;
}

inline tStructure& tStructure::Clear()
{
	m_FEModel.Clear();
	m_loads.Clear();
	m_analyses.Clear();
	return *this;
}

inline void tStructure::Dump(const char* filename_) const
{
	std::ofstream dump_out(filename_);
	dump_out << "FE-model\n";
	m_FEModel.Dump(dump_out);
	dump_out << "\nLoads\n";
	m_loads.Dump(dump_out, m_FEModel);
}
