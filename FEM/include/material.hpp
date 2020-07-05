#pragma once

#include "NumTypes.hpp"
#include "Tensors.hpp"
#include "NameList.hpp"
#include "FileInOut.hpp"

#include <string>
#include <map>
#include <vector>
#include <istream>

class tFinitElement;

class tMaterial
{
public:
	static tMaterial* CreateMaterial(const std::string& kindName_)
	{
		return Factory.CallFunction(kindName_);
	}

	tMaterial() : m_specificWeight()
	{
	}
	virtual ~tMaterial()
	{
	}
	virtual const std::string& Kind() const = 0;
	const std::string& GetName() const
	{
		return m_name;
	}
	virtual bool isIsotropic() const = 0;
	virtual bool isStrainDependent() const
	{
		return true;
	}
	tMaterial& SetName(const std::string& name)
	{
		m_name = name;
		return *this;
	}

	virtual tMaterial& ReadData(std::istream&) = 0;

	virtual const SymmetricTensor4s& ElasTensor(
		const Tensor2s& rot_, real_t strain_ = 0.) const = 0;
	virtual const SymmetricTensor4s& ElasTensor_1D(
		const Tensor2s& rot_, real_t strain_ = 0.) const = 0;  // for 1D-stress state
	virtual const SymmetricTensor4s& ElasTensor_2D(
		const Tensor2s& rot_, real_t strain_ = 0.) const = 0;  // for 2D-stress state
	real_t Density() const
	{
		return m_specificWeight;
	}
	//   static real_t StressIntensity(const Tensor2& stress_) {Tensor2a tmp; return
	//   sqrt(3.*stress_.Dev(tmp).I2());}
	/*   static*/ real_t StrainIntensity(const Tensor2s& strain_) const
	{
		Tensor2s tmp;
		return 2. * sqrt(-strain_.Dev(tmp).I2() / 3.);
	}

protected:
	static const tNamesToFuns<tMaterial> Factory;

	std::string m_name;
	real_t m_specificWeight;
};

class tIsotropicMaterial : virtual public tMaterial
{
public:
	virtual ~tIsotropicMaterial()
	{
	}

	bool isIsotropic() const override
	{
		return true;
	}
	virtual real_t Mu(real_t = 0.) const = 0;
	virtual real_t YoungModule(real_t = 0.) const = 0;
	virtual real_t ShearModule(real_t = 0.) const = 0;
	virtual real_t PoissonRatio(real_t = 0.) const = 0;
	virtual real_t Lambda(real_t = 0.) const = 0;
};

class tStrainIndependentMaterial : virtual public tMaterial
{
public:
	tStrainIndependentMaterial()
		: m_elasTensor_prinAx(), m_elasTensor1D_prinAx(), m_elasTensor2D_prinAx()
	{
	}

	bool isStrainDependent() const override
	{
		return false;
	}

	//   virtual SymmetricTensor4s& ElasTensor   (const Tensor2s& rot_, real_t strain_=0.,
	//   SymmetricTensor4s&) const =0; virtual SymmetricTensor4s& ElasTensor_1D(const Tensor2s&
	//   rot_, real_t strain_=0., SymmetricTensor4s&) const =0;//for 1D-stress state virtual
	//   SymmetricTensor4s& ElasTensor_2D(const Tensor2s& rot_, real_t strain_=0.,
	//   SymmetricTensor4s&) const =0;//for 2D-stress state

	const SymmetricTensor4s& ElasTensor(const Tensor2s& rot_, real_t = 0.) const override
	{
		return __RotateElasTensor(m_elasTensor_prinAx, rot_, m_elasTensor_DB);
	}
	const SymmetricTensor4s& ElasTensor_1D(const Tensor2s& rot_, real_t = 0.) const override
	{
		return __RotateElasTensor(m_elasTensor1D_prinAx, rot_, m_elasTensor1D_DB);
	}
	const SymmetricTensor4s& ElasTensor_2D(const Tensor2s& rot_, real_t = 0.) const override
	{
		return __RotateElasTensor(m_elasTensor2D_prinAx, rot_, m_elasTensor2D_DB);
	}

protected:
	// elastic constants tensors: usual, for 1d-stress, for 2d-stress
	SymmetricTensor4s m_elasTensor_prinAx;
	SymmetricTensor4s m_elasTensor1D_prinAx;
	SymmetricTensor4s m_elasTensor2D_prinAx;

private:
	typedef std::map<std::pair<real_t, std::pair<real_t, real_t> >, SymmetricTensor4s> tMapRotSymT4;

	static const SymmetricTensor4s& __RotateElasTensor(
		const SymmetricTensor4s& src_, const Tensor2s& rot_, tMapRotSymT4& db_);

	mutable tMapRotSymT4 m_elasTensor_DB, m_elasTensor1D_DB, m_elasTensor2D_DB;
};

#pragma warning(push)
#pragma warning(disable : 4250)
class tIsotropicLinearElasticMaterial : public tIsotropicMaterial, public tStrainIndependentMaterial
{
public:
	static tMaterial* NewMaterial()
	{
		return new tIsotropicLinearElasticMaterial();
	}

	tMaterial& ReadData(std::istream&) override;

	const SymmetricTensor4s& ElasTensor(const Tensor2s&, real_t = 0.) const override
	{
		return m_elasTensor_prinAx;
	}

	const std::string& Kind() const override
	{
		return KindName;
	}
	real_t Mu(real_t = 0.) const override
	{
		return m_Mu;
	}
	real_t YoungModule(real_t = 0.) const override
	{
		return m_E;
	}
	real_t ShearModule(real_t = 0.) const override
	{
		return Mu();
	}
	real_t PoissonRatio(real_t = 0.) const override
	{
		return m_Nu;
	}
	real_t Lambda(real_t = 0.) const override
	{
		return m_lambda;
	}

private:
	static const std::string& KindName;

	tIsotropicLinearElasticMaterial() : m_E(), m_Nu(), m_lambda(), m_Mu()
	{
	}

	real_t m_E, m_Nu, m_lambda, m_Mu;
};
#pragma warning(pop)

class tIsotropicPlasticMaterial : public tIsotropicMaterial
{
public:
	virtual ~tIsotropicPlasticMaterial()
	{
	}

	static tMaterial* NewMaterial()
	{
		return new tIsotropicPlasticMaterial();
	}

	virtual SymmetricTensor4s& ElasTensor_1D(const Tensor2s&, SymmetricTensor4s& res_) const
	{
		return res_;
	}  //заглушка!!
	virtual SymmetricTensor4s& ElasTensor_2D(const Tensor2s&, SymmetricTensor4s& res_) const
	{
		return res_;
	}  //заглушка!!

	const SymmetricTensor4s& ElasTensor(
		const Tensor2s& /*rot_*/, real_t /*strain_ = 0.*/) const override
	{
		return m_tmpToKill;
	}  //заглушка!!
	const SymmetricTensor4s& ElasTensor_1D(
		const Tensor2s& /*rot_*/, real_t /*strain_ = 0.*/) const override
	{
		return m_tmpToKill;
	}  //заглушка!!
	const SymmetricTensor4s& ElasTensor_2D(
		const Tensor2s& /*rot_*/, real_t /*strain_ = 0.*/) const override
	{
		return m_tmpToKill;
	}  //заглушка!!

	tMaterial& ReadData(std::istream&) override;

	const std::string& Kind() const override
	{
		return KindName;
	}
	bool isIsotropic() const override
	{
		return true;
	}
	bool isStrainDependent() const override
	{
		return true;
	}
	virtual const SymmetricTensor4s& ElasTensor(real_t = 0) const;
	real_t YoungModule(real_t = 0.) const override;
	real_t ShearModule(real_t = 0.) const override
	{
		return Mu();
	}
	real_t PoissonRatio(real_t = 0.) const override
	{
		return m_Nu;
	}
	real_t Lambda(real_t eps_ = 0.) const override
	{
		return YoungModule(eps_) * m_Nu / (1. + m_Nu) * (1. - 2. * m_Nu);
	}
	real_t Mu(real_t eps_ = 0.) const override
	{
		return YoungModule(eps_) / 2. * (1. + m_Nu);
	}
	real_t SigmaEpsilonDiagram(real_t eps_);

	SymmetricTensor4s m_tmpToKill;

private:
	static const std::string& KindName;

	explicit tIsotropicPlasticMaterial() : m_tmpToKill(), m_E(), m_Nu()
	{
	}

	real_t m_E, m_Nu /*, lambda, mu*/;
	std::vector<real_t> m_sigma, m_epsilon;
	mutable std::map<real_t, SymmetricTensor4s> m_elasTensorDB;
};
