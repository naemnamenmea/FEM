#pragma once

#include "app_constants.hpp"
#include "NumTypes.h"  // cardinal_t, T
#include "Tensors.h"   //Tensor1
#include <vector>

template <typename T>
class tFinitElement;

template <typename T>
class Node3d
{
private:
	Tensor1s<T> m_RadiusVector;	 // m_coord
	mutable const Tensor1s<T>* m_pDisplacement;
	std::vector<FEM::DEGREES_OF_FREEDOM> m_DegreesOfFreedom;

public:
	Node3d() : m_RadiusVector(), m_pDisplacement(NULL), m_DegreesOfFreedom()
	{
	}
	//   explicit Node(const Node& o_):  RadiusVector(o_.RadiusVector),
	//   m_pDisplacement(o_.m_pDisplacement), m_DegreesOfFreedom(o_.m_DegreesOfFreedom) {} Node&
	//   AssignCoord(const Tensor1& newCoord_) {RadiusVector = newCoord_; return *this;}

	Node3d<T>& AssignCoord(const Tensor1s<T>& newCoord_)
	{
		m_RadiusVector = newCoord_;
		return *this;
	}

	bool DisplacementIs0() const
	{
		return m_pDisplacement == NULL /*|| m_pDisplacement->is0()*/;
	}

	const Tensor1s<T>& Displacement() const
	{
		return *m_pDisplacement;
	}

	const Tensor1s<T>& GetCoord() const
	{
		return m_RadiusVector;
	}

	const std::vector<FEM::DEGREES_OF_FREEDOM>& GetDegreesOfFreedom() const
	{
		return m_DegreesOfFreedom;
	}

	//   Tensor1 DeformedCoord()      const {return m_pDisplacement==NULL? Coord() : Coord() +
	//   *m_pDisplacement;}

	T Coord(COMPONENT i_) const
	{
		return m_RadiusVector(i_);
	}

	Node3d<T>& AssignCoord(COMPONENT i_, const T val_)
	{
		m_RadiusVector(i_) = val_;
		return *this;
	}

	bool HasFixed(COMPONENT i_) const
	{
		return !(
			i_ == X ? m_DegreesOfFreedom.FreeX
					: (i_ == Y ? m_DegreesOfFreedom.FreeY : m_DegreesOfFreedom.FreeZ));
	}

	bool HasFixed(cardinal_t i_) const
	{
		return !(
			i_ == 1 ? m_DegreesOfFreedom.FreeX
					: (i_ == 2 ? m_DegreesOfFreedom.FreeY : m_DegreesOfFreedom.FreeZ));
	}

	cardinal_t NumberOfDOFs() const
	{
		return (m_DegreesOfFreedom.FreeX ? 1 : 0) + (m_DegreesOfFreedom.FreeY ? 1 : 0) +
			   (m_DegreesOfFreedom.FreeZ ? 1 : 0);
	}

	Node3d<T>& Fix(COMPONENT i_)
	{
		(i_ == X ? m_DegreesOfFreedom.FreeX
				 : (i_ == Y ? m_DegreesOfFreedom.FreeY : m_DegreesOfFreedom.FreeZ)) = false;
		return *this;
	}

	Node3d<T>& Free(COMPONENT i_)
	{
		(i_ == X ? m_DegreesOfFreedom.FreeX
				 : (i_ == Y ? m_DegreesOfFreedom.FreeY : m_DegreesOfFreedom.FreeZ)) = true;
		return *this;
	}

	void LinkDisplacement(const Tensor1s<T>* p_) const
	{
		m_pDisplacement = p_;
	}

	void UnLinkDisplacement() const
	{
		m_pDisplacement = NULL;
	}

	/*
	 * example of input:
	 * 1 0. −10. 0. fixed fixed fixed
	 */
	template <typename T>
	friend std::istream& operator>>(std::istream& is, Node3d<T>& node)
	{
		std::vector<std::string> _freedomDegrees(node.m_DegreesOfFreedom.size());
		is >> node.m_RadiusVector >> _freedomDegrees;
		for (size_t i = 0; i < node.m_DegreesOfFreedom.size(); ++i)
		{
			node.m_DegreesOfFreedom[i] = _freedomDegrees[i] == "free"
											 ? FEM::DEGREES_OF_FREEDOM::FREE
											 : FEM::DEGREES_OF_FREEDOM::FIXED;
		}
		return is;
	}
};
