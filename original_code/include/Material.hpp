#pragma once

#include "NumTypes.h"
#include "Tensors.h"
#include "NameList.h"
//#include "FileInOut.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>

template<typename T>
class tFinitElement;

template <typename T>
class tMaterial
{
protected:
	static const tNamesToFuns<tMaterial<T>> Factory;
	std::string Caption;
	T SpecificWeight;

private:
public:
	tMaterial() : SpecificWeight()
	{
	}

	static tMaterial<T>* NewMaterial(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}
	virtual ~tMaterial()
	{
	}
	virtual const std::string& Kind() const = 0;
	const std::string& GetName() const
	{
		return Caption;
	}
	virtual bool isIsotropic() const = 0;
	virtual bool isStrainDependent() const
	{
		return true;
	}
	tMaterial<T>& SetName(const std::string& newName_)
	{
		Caption = newName_;
		return *this;
	}

	virtual tMaterial<T>& ReadData(std::istream&) = 0;

	virtual const SymmetricTensor4s<T>& ElasTensor(
		const Tensor2s<T>& rot_, T strain_ = 0.) const = 0;
	virtual const SymmetricTensor4s<T>& ElasTensor_1D(
		const Tensor2s<T>& rot_, T strain_ = 0.) const = 0;	 // for 1D-stress state
	virtual const SymmetricTensor4s<T>& ElasTensor_2D(
		const Tensor2s<T>& rot_, T strain_ = 0.) const = 0;	 // for 2D-stress state
	T Density() const
	{
		return SpecificWeight;
	}
	//   static T StressIntensity(const Tensor2& stress_) {Tensor2a tmp; return
	//   sqrt(3.l*stress_.Dev(tmp).I2());}
	/*   static*/ T StrainIntensity(const Tensor2s<T>& strain_) const
	{
		Tensor2s<T> tmp;
		return 2. * sqrt(-strain_.Dev(tmp).I2() / 3.l);
	}
};
//------------------------------------------------------------------------------
template <typename T>
class tIsotropicMaterial : virtual public tMaterial<T>
{
public:
	virtual bool isIsotropic() const
	{
		return true;
	}
	virtual T Mu(T = 0.) const = 0;
	virtual T YoungModule(T = 0.) const = 0;
	virtual T ShearModule(T = 0.) const = 0;
	virtual T PoissonRatio(T = 0.) const = 0;
	virtual T Lambda(T = 0.) const = 0;
};
//------------------------------------------------------------------------------
template <typename T>
class tStrainIndependentMaterial : virtual public tMaterial<T>
{
protected:
	SymmetricTensor4s<T> ElasTensor_prinAx, ElasTensor1D_prinAx,
		ElasTensor2D_prinAx;  // elastic constants tensors: usual, for 1d-stress, for 2d-stress
private:
	typedef std::map<std::pair<T, std::pair<T, T>>, SymmetricTensor4s<T>> tMapRotSymT4;
	mutable tMapRotSymT4 ElasTensor_DB, ElasTensor1D_DB, ElasTensor2D_DB;
	static const SymmetricTensor4s<T>& __RotateElasTensor(
		const SymmetricTensor4s<T>& src_, const Tensor2s<T>& rot_, tMapRotSymT4& db_);

public:
	tStrainIndependentMaterial() : ElasTensor_prinAx(), ElasTensor1D_prinAx(), ElasTensor2D_prinAx()
	{
	}

	virtual bool isStrainDependent() const
	{
		return false;
	}

	//   virtual SymmetricTensor4s<T>& ElasTensor   (const Tensor2s<T>& rot_, T strain_=0.,
	//   SymmetricTensor4s<T>&) const =0; virtual SymmetricTensor4s<T>& ElasTensor_1D(const
	//   Tensor2s<T>& rot_, T strain_=0., SymmetricTensor4s<T>&) const =0;//for 1D-stress state
	//   virtual SymmetricTensor4s<T>& ElasTensor_2D(const Tensor2s<T>& rot_, T strain_=0.,
	//   SymmetricTensor4s<T>&) const =0;//for 2D-stress state

	virtual const SymmetricTensor4s<T>& ElasTensor(const Tensor2s<T>& rot_, T = 0.) const
	{
		return __RotateElasTensor(ElasTensor_prinAx, rot_, ElasTensor_DB);
	}
	virtual const SymmetricTensor4s<T>& ElasTensor_1D(const Tensor2s<T>& rot_, T = 0.) const
	{
		return __RotateElasTensor(ElasTensor1D_prinAx, rot_, ElasTensor1D_DB);
	}
	virtual const SymmetricTensor4s<T>& ElasTensor_2D(const Tensor2s<T>& rot_, T = 0.) const
	{
		return __RotateElasTensor(ElasTensor2D_prinAx, rot_, ElasTensor2D_DB);
	}
};
//------------------------------------------------------------------------------
#pragma warning(push)
#pragma warning(disable : 4250)
template <typename T>
class tIsotropicLinearElasticMaterial : public tIsotropicMaterial<T>,
										public tStrainIndependentMaterial<T>
{
private:
	static const std::string& KindName;
	tIsotropicLinearElasticMaterial()
	{
	}

	T E, Nu, lambda, mu;

public:
	static tMaterial<T>* NewMaterial()
	{
		return new tIsotropicLinearElasticMaterial<T>();
	}

	virtual tMaterial<T>& ReadData(std::istream&);

	virtual const SymmetricTensor4s<T>& ElasTensor(const Tensor2s<T>&, T = 0.) const
	{
		return this->ElasTensor_prinAx;
	}

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	virtual T Mu(T = 0.) const
	{
		return mu;
	}
	virtual T YoungModule(T = 0.) const
	{
		return E;
	}
	virtual T ShearModule(T = 0.) const
	{
		return Mu();
	}
	virtual T PoissonRatio(T = 0.) const
	{
		return Nu;
	}
	virtual T Lambda(T = 0.) const
	{
		return lambda;
	}
};
#pragma warning(pop)

template <typename T>
class tIsotropicPlasticMaterial : public tIsotropicMaterial<T>
{
private:
	static const std::string& KindName;
	explicit tIsotropicPlasticMaterial() : tmp_to_kill()
	{
	}
	T E, Nu /*, lambda, mu*/;
	std::vector<T> Sigma, Epsilon;
	mutable std::map<T, SymmetricTensor4s<T>> ElasTensor_DB;

public:

	static tMaterial<T>* NewMaterial()
	{
		return new tIsotropicPlasticMaterial<T>();
	}

	SymmetricTensor4s<T> tmp_to_kill;
	virtual SymmetricTensor4s<T>& ElasTensor_1D(
		const Tensor2s<T>&, SymmetricTensor4s<T>& res_) const
	{
		return res_;
	}
	virtual SymmetricTensor4s<T>& ElasTensor_2D(
		const Tensor2s<T>&, SymmetricTensor4s<T>& res_) const
	{
		return res_;
	}  //��������!!

	virtual const SymmetricTensor4s<T>& ElasTensor(const Tensor2s<T>& rot_, T strain_ = 0.) const
	{
		rot_;
		strain_;
		return tmp_to_kill;
	}  //��������!!
	virtual const SymmetricTensor4s<T>& ElasTensor_1D(const Tensor2s<T>& rot_, T strain_ = 0.) const
	{
		rot_;
		strain_;
		return tmp_to_kill;
	}  //��������!!
	virtual const SymmetricTensor4s<T>& ElasTensor_2D(const Tensor2s<T>& rot_, T strain_ = 0.) const
	{
		rot_;
		strain_;
		return tmp_to_kill;
	}  //��������!!

	virtual tMaterial<T>& ReadData(std::istream&);

	virtual const std::string& Kind() const
	{
		return KindName;
	}
	virtual bool isIsotropic() const
	{
		return true;
	}
	virtual bool isStrainDependent() const
	{
		return true;
	}
	virtual const SymmetricTensor4s<T>& ElasTensor(T = 0) const;
	virtual T YoungModule(T = 0.) const;
	virtual T ShearModule(T = 0.) const
	{
		return Mu();
	}
	virtual T PoissonRatio(T = 0.) const
	{
		return Nu;
	}
	virtual T Lambda(T eps_ = 0.) const
	{
		return YoungModule(eps_) * Nu / (1. + Nu) * (1. - 2. * Nu);
	}
	virtual T Mu(T eps_ = 0.) const
	{
		return YoungModule(eps_) / 2. * (1. + Nu);
	}
	T SigmaEpsilonDiagram(T eps_);
};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

class Material
{
public:
	Material();
	Material(std::string name);

	class Properties
	{
	public:
		double GetDensity() const
		{
			return m_density;
		}
		void SetDensity(double newDensity)
		{
			m_density = newDensity;
		}

		friend std::istream& operator>>(std::istream& is, Properties& prop);

	private:
		double m_density;
	};

	const std::string& GetName() const
	{
		return m_name;
	}

	const Properties& GetProperties() const
	{
		return m_properties;
	}
	void SetProperties(const Properties& prop)
	{
		m_properties = prop;
	}

	friend std::istream& operator>>(std::istream& is, Material& material);

private:
	std::string m_name;
	Properties m_properties;
};

inline bool operator==(const Material& lhs, const Material& rhs)
{
	return lhs.GetName() == rhs.GetName();
}

namespace std
{
template <>
struct hash<Material>
{
	std::size_t operator()(Material const& material) const noexcept
	{
		std::size_t h = std::hash<std::string>{}(material.GetName());
		return h;
	}
};
}  // namespace std
