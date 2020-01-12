#pragma once

#include <filesystem>
#include <fstream>
#include <list>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>
#include "app_constants.hpp"
#include "auxiliary.hpp"
#include "finite_element.hpp"
#include "material.hpp"
#include "numeric_math.hpp"
#include "stl_containers_read.hpp"

namespace
{
namespace fs = std::filesystem;
}

template <size_t DIM>
class FiniteElementModel
{
public:
	FiniteElementModel();

	const auto& GetNodes() const
	{
		return nodes;
	}

	const auto& GetMaterials() const
	{
		return materialCollection;
	}

	const auto& GetFiniteElements() const
	{
		return finiteElements;
	}

	void SetParams(unsigned int argc, const char* argv[]);

	void ReadData(fs::path path, FEM::IO_FORMAT inputFormat = FEM::IO_FORMAT::CHEKHOV);

	void WriteData(fs::path path, FEM::IO_FORMAT outputFormat = FEM::IO_FORMAT::CHEKHOV) const
	{
	}

private:
	void ReadFiniteElement(std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement)
	{
		ReadGenericFEEntry(is, finiteElement);
		ReadSpecFEEntry(is, finiteElement);
	}

	void ReadGenericFEEntry(std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement);

	void ReadSpecFEEntry(std::istream& is, std::shared_ptr<FiniteElementBase> finiteElement);

	std::unordered_map<size_t, Node<DIM>> nodes;
	std::unordered_set<Material> materialCollection;
	std::unordered_map<size_t, shared_ptr<FiniteElementBase>> finiteElements;
};
