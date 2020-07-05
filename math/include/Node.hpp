#pragma once

#include "NumTypes.hpp"	 // size_t, real_t
#include "Tensors.hpp"	 //Tensor1
#include <vector>

class tFinitElement;

class tNode
{
public:
	struct DegreesOfFreedom
	{
		DegreesOfFreedom() : m_freeX(false), m_freeY(false), m_freeZ(false){};
		DegreesOfFreedom(bool tX_, bool tY_, bool tZ_ = false)
			: m_freeX(tX_), m_freeY(tY_), m_freeZ(tZ_)
		{
		}

		bool m_freeX, m_freeY, m_freeZ;	 // true <=> free, false <=> fixed
	};

	tNode() : m_radiusVector(), m_pDisplacement(nullptr), m_DOF()
	{
	}
	//   explicit tNode(const tNode& o_):  RadiusVector(o_.RadiusVector),
	//   pDisplacement(o_.pDisplacement), DOF(o_.DOF) {} tNode& AssignCoord(const Tensor1&
	//   newCoord_) {RadiusVector = newCoord_; return *this;}
	tNode& AssignCoord(const Tensor1s& newCoord_)
	{
		m_radiusVector = newCoord_;
		return *this;
	}
	const Tensor1s& Coord() const
	{
		return m_radiusVector;
	}
	bool DisplacementIs0() const
	{
		return m_pDisplacement == nullptr /*|| pDisplacement->is0()*/;
	}
	const Tensor1s& Displacement() const
	{
		return *m_pDisplacement;
	}
	//   Tensor1 DeformedCoord()      const {return pDisplacement==nullptr? Coord() : Coord() +
	//   *pDisplacement;}
	real_t Coord(component_t i_) const
	{
		return m_radiusVector(static_cast<small_t>(i_));
	}
	tNode& AssignCoord(component_t i_, const real_t val_)
	{
		m_radiusVector(static_cast<small_t>(i_)) = val_;
		return *this;
	}
	bool HasFixed(component_t i_) const
	{
		return !(i_ == X ? m_DOF.m_freeX : (i_ == Y ? m_DOF.m_freeY : m_DOF.m_freeZ));
	}
	bool HasFixed(size_t i_) const
	{
		return !(i_ == 1 ? m_DOF.m_freeX : (i_ == 2 ? m_DOF.m_freeY : m_DOF.m_freeZ));
	}
	size_t NumberOfDOFs() const
	{
		return static_cast<size_t>(
			(m_DOF.m_freeX ? 1 : 0) + (m_DOF.m_freeY ? 1 : 0) + (m_DOF.m_freeZ ? 1 : 0));
	}
	tNode& Fix(component_t i_)
	{
		(i_ == X ? m_DOF.m_freeX : (i_ == Y ? m_DOF.m_freeY : m_DOF.m_freeZ)) = false;
		return *this;
	}
	tNode& Free(component_t i_)
	{
		(i_ == X ? m_DOF.m_freeX : (i_ == Y ? m_DOF.m_freeY : m_DOF.m_freeZ)) = true;
		return *this;
	}
	void LinkDisplacement(const Tensor1s* p_) const
	{
		m_pDisplacement = p_;
	}
	void UnLinkDisplacement() const
	{
		m_pDisplacement = nullptr;
	}

private:
	Tensor1s m_radiusVector;
	mutable const Tensor1s* m_pDisplacement;
	DegreesOfFreedom m_DOF;
};
