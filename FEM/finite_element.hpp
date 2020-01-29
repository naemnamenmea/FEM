#pragma once

#include "math_constants.hpp"
#include "material.hpp"
#include "stl_containers_read.hpp"
#include <istream>
#include <list>
#include <vector>

template <size_t DIM>
class Node
{
public:
	Node() : coord(DIM), freedomDegrees(DIM)
	{
	}
	Node(const std::vector<double>& coord, std::vector<FEM::DEGREES_OF_FREEDOM> fd)
		: coord(DIM), freedomDegrees(fd)
	{
	}

	const auto& GetCoord() const
	{
		return coord;
	}

	const auto& GetFreedomDegrees() const
	{
		return freedomDegrees;
	}

  /*
  * example of input:
  * 1 0. −10. 0. fixed fixed fixed
  */
  friend std::istream& operator>>(std::istream& is, Node& node)
  {
    std::vector<std::string> _freedomDegrees(node.freedomDegrees.size());
    is >> node.coord >> _freedomDegrees;
    for (size_t i = 0; i < node.freedomDegrees.size(); ++i)
    {
      node.freedomDegrees[i] = _freedomDegrees[i] == "fixed" ? FEM::DEGREES_OF_FREEDOM::FIXED
        : FEM::DEGREES_OF_FREEDOM::FREE;
    }
    return is;
  }

private:
	std::vector<double> coord;
	std::vector<FEM::DEGREES_OF_FREEDOM> freedomDegrees;
};

/*
 * 1 dim - 2 nodes
 * 2 dim - 4 nodes
 * 3 dim - 8 nodes
 */
class FiniteElementBase
{
public:
	FiniteElementBase();
	FiniteElementBase(size_t nodes);

	const auto& GetType() const
	{
		return type;
	}
	const auto& GetNodeNumbers() const
	{
		return nodeNumbers;
	}
	const auto& GetMaterial() const
	{
		return *material;
	}

	void SetType(std::string type)
	{
		this->type = std::move(type);
	}
	void SetNodeNumbers(const std::vector<size_t>& nodeNumbers)
	{
		this->nodeNumbers = nodeNumbers;
	}
	void SetMaterial(std::unordered_set<Material>::const_iterator material)
	{
		this->material = material;
	}

	virtual void ReadProperties(std::istream& is);
	virtual unsigned int GetDim() const = 0;
	virtual double GetParameter() const
	{
		return mathdef::MAX_DOUBLE;
	}

protected:
	std::string type;
	std::vector<size_t> nodeNumbers;
  std::unordered_set<Material>::const_iterator material;
};

class FiniteElement1d : public FiniteElementBase
{
public:
	FiniteElement1d(size_t nodes);

	void ReadProperties(std::istream& is) override;

	double GetParameter() const override
	{
		return crossSectionalArea;
	}

	unsigned int GetDim() const override
	{
		return 1;
	}

private:
	double crossSectionalArea;
};

class FiniteElement2d : public FiniteElementBase
{
public:
	FiniteElement2d(size_t nodes);

	void ReadProperties(std::istream& is) override;

	double GetParameter() const override
	{
		return thickness;
	}

	unsigned int GetDim() const override
	{
		return 2;
	}

private:
	double thickness;
};

class FiniteElement3d : public FiniteElementBase
{
public:
	FiniteElement3d(size_t nodes);

	unsigned int GetDim() const override
	{
		return 3;
	}
};

std::shared_ptr<FiniteElementBase> CreateFiniteElement(size_t dim, size_t nodes);
