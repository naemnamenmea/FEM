#include "stdafx.hpp"

#include <functional>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
//------------------------------------------------------------------------------
#include "Material.hpp"
#ifdef STRONGCHECK
#include "ThrowMessage.h"
#endif
//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
//----------------------------tMaterial<T>-----------------------------------------
//------------------------------------------------------------------------------
static const cardinal_t material_kinds_number = 2;
static const char* initNames[] = {"IsotropicLinearElastic",
								  "IsotropicPlastic"};	// Add new kind names here

static tMaterial<__real_t>* (*initFunctions[])() = {
	tIsotropicLinearElasticMaterial<__real_t>::NewMaterial,
	tIsotropicPlasticMaterial<__real_t>::NewMaterial};	 // Add new pointers to creaters here

template <typename T>
const tNamesToFuns<tMaterial<T>> tMaterial<T>::Factory(
	initNames, initFunctions, material_kinds_number);

template <typename T>
const string& tIsotropicLinearElasticMaterial<T>::KindName =
	tMaterial<T>::Factory.FindName(tIsotropicLinearElasticMaterial<T>::NewMaterial);
template <typename T>
const string& tIsotropicPlasticMaterial<T>::KindName =
	tMaterial<T>::Factory.FindName(tIsotropicPlasticMaterial<T>::NewMaterial);
//------------------------------------------------------------------------------
//---------------------------tStrainIndependentMaterial<T>-------------------------
//------------------------------------------------------------------------------
template <typename T>
const SymmetricTensor4s<T>& tStrainIndependentMaterial<T>::__RotateElasTensor(
	const SymmetricTensor4s<T>& init_,
	const Tensor2s<T>& rot_,
	tMapRotSymT4& db_)	//, SymmetricTensor4s<T>& result_)
{
	static Tensor1s<T> vec;
	rot_.AttachedVector(vec);
	typename tMapRotSymT4::iterator p = db_.find(
		make_pair(vec(1), make_pair(vec(2), vec(3))));	// modify to find with given precision!!!
	if (p == db_.end())
	{
#ifdef STRONGCHECK
		pair<tMapRotSymT4::iterator, bool> ins = db_.insert(tMapRotSymT4::value_type(
			make_pair(vec(1), make_pair(vec(2), vec(3))), SymmetricTensor4s<T>(init_)));
		Assert(ins.second, "failed insertion in tStrainIndependentMaterial<T>::RotateElasTensor");
		p = ins.first;
#else
		p = db_
				.insert(tMapRotSymT4::value_type(
					make_pair(vec(1), make_pair(vec(2), vec(3))), SymmetricTensor4s<T>(init_)))
				.first;
#endif
		p->second.TransformAll(rot_);
	}
	return p->second;
	// return result_=p->second;
}
//------------------------------------------------------------------------------
//---------------------------tIsotropicLinearElasticMaterial<T>--------------------
//------------------------------------------------------------------------------
template <typename T>
tMaterial<T>& tIsotropicLinearElasticMaterial<T>::ReadData(istream& str_in_)
{
	string buffer;
	for (int i = 0; i < 3; ++i)
	{
		getline(str_in_, buffer, '=');
		str_in_.get();
		while (/*std::*/ isspace(buffer[buffer.size() - 1]))
			buffer.resize(buffer.size() - 1);  // pop_back();
		while (/*std::*/ isspace(buffer[0])) buffer.erase(0, 1);
		if (buffer == "E")
			str_in_ >> E;
		else if (buffer == "rho")
			str_in_ >> this->SpecificWeight;
		else if (buffer == "nu")
			str_in_ >> Nu;
		else
		{
			str_in_ >> buffer;
			--i;
		}
	}
	mu = E / (1. + Nu);
	lambda = mu * Nu / (1. - 2. * Nu);
	mu *= 0.5;

	this->ElasTensor_prinAx.Assign0();
	this->ElasTensor_prinAx(1, 1) = this->ElasTensor_prinAx(2, 2) = this->ElasTensor_prinAx(3, 3) =
		lambda + 2. * mu;
	this->ElasTensor_prinAx(1, 2) = this->ElasTensor_prinAx(1, 3) = this->ElasTensor_prinAx(2, 3) =
		lambda;
	this->ElasTensor_prinAx(4, 4) = this->ElasTensor_prinAx(5, 5) = this->ElasTensor_prinAx(6, 6) =
		mu;

	this->ElasTensor1D_prinAx = this->ElasTensor_prinAx;
	this->ElasTensor1D_prinAx(1, 1) = E;

	this->ElasTensor2D_prinAx = this->ElasTensor_prinAx;
	this->ElasTensor2D_prinAx(1, 1) = this->ElasTensor2D_prinAx(2, 2) =
		mu * (lambda / (.5 * lambda + mu) + 2.);
	this->ElasTensor2D_prinAx(1, 2) = lambda * mu / (.5 * lambda + mu);

	return *this;
}
//------------------------------------------------------------------------------
//-------------------------tIsotropicPlasticMaterial<T>----------------------------
//------------------------------------------------------------------------------
// tMaterial<T>& tIsotropicPlasticMaterial<T>::ReadData (InputFileStream& stream_)
// tMaterial<T>& tIsotropicPlasticMaterial<T>::ReadData (InputFileStreamXML& stream_)
//{

template <typename T>
tMaterial<T>& tIsotropicPlasticMaterial<T>::ReadData(istream& str_in_)
{
	str_in_ >> this->SpecificWeight >> E >> Nu;

	string DiagramType;
	str_in_ >> DiagramType;
	if (DiagramType == "Pointwise")
	{
		str_in_.peek();

		while (!str_in_.eof())
		{
			T tmp;
			str_in_ >> tmp;
			Sigma.push_back(tmp);
			str_in_ >> tmp;
			Epsilon.push_back(tmp);
			str_in_.peek();
		}
	}

	return *this;
}
//------------------------------------------------------------------------------
template <typename T>
T tIsotropicPlasticMaterial<T>::YoungModule(T Epsi_) const
{
	using namespace std::placeholders;
#ifdef STRONGCHECK
	Assert(Epsi_ >= Epsilon.front(), "Invalid sigma/epsilon diagram");
#endif
	typename vector<T>::const_iterator p =
		find_if(Epsilon.begin(), Epsilon.end(), bind(greater<T>(), _1, Epsi_));
	if (p == Epsilon.end())
		return (Sigma.back() - *(Sigma.end() - 2)) / (Epsilon.back() - *(Epsilon.end() - 2));
	const size_t dist = distance(Epsilon.begin(), p);
	return Sigma[dist] - Sigma[dist - 1] / (*p - *(p - 1));
}
//------------------------------------------------------------------------------
template <typename T>
const SymmetricTensor4s<T>& tIsotropicPlasticMaterial<T>::ElasTensor(T strain_level_) const
{
	// T epsilon = StrainIntensity(strain_);
	typename map<T, SymmetricTensor4s<T>>::iterator p = ElasTensor_DB.find(strain_level_);
	if (p == ElasTensor_DB.end())
	{
#ifdef STRONGCHECK
		pair<map<T, SymmetricTensor4s<T>>::iterator, bool> ins = ElasTensor_DB.insert(
			map<T, SymmetricTensor4s<T>>::value_type(strain_level_, SymmetricTensor4s<T>()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in tIsotropicPlasticMaterial<T>::ElasTensor(const "
			"Tensor2s<T>&)");
		p = ins.first;
#else
		p = ElasTensor_DB
				.insert(
					map<T, SymmetricTensor4s<T>>::value_type(strain_level_, SymmetricTensor4s<T>()))
				.first;
#endif
		SymmetricTensor4s<T>& result = p->second;
		result.Assign0();

		T nu = PoissonRatio(strain_level_), tmp = YoungModule(strain_level_) / (1. + nu);
		result(4, 4) = result(5, 5) = result(6, 6) = .5 * tmp;
		tmp /= (1. - 2. * nu);
		result(1, 1) = result(2, 2) = result(3, 3) = (1. - nu) * tmp;
		result(1, 2) = result(1, 3) = result(2, 3) = nu * tmp;
	}
	return p->second;
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

std::istream& operator>>(std::istream& is, Material::Properties& prop)
{
	is >> prop.m_density;
	return is;
}

/*
 * example of input:
 * Al_Alloy 0.0028
 */
std::istream& operator>>(std::istream& is, Material& material)
{
	is >> material.m_name >> material.m_properties;
	return is;
}

Material::Material() : m_properties()
{
}

Material::Material(std::string name) : m_properties(), m_name(std::move(name))
{
}
