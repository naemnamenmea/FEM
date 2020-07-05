#include "Material.hpp"
#ifdef STRONGCHECK
#include "ThrowMessage.hpp"
#endif

#include <functional>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>

using namespace std;

static const size_t material_kinds_number = 2;
static const char* initNames[] = {"IsotropicLinearElastic",
								  "IsotropicPlastic"};	// Add new kind names here
static tMaterial* (*initFunctions[])() = {
	tIsotropicLinearElasticMaterial::NewMaterial,
	tIsotropicPlasticMaterial::NewMaterial};  // Add new pointers to creaters here

const tNamesToFuns<tMaterial> tMaterial::Factory(initNames, initFunctions, material_kinds_number);

const string& tIsotropicLinearElasticMaterial::KindName =
	tMaterial::Factory.FindName(tIsotropicLinearElasticMaterial::NewMaterial);
const string& tIsotropicPlasticMaterial::KindName =
	tMaterial::Factory.FindName(tIsotropicPlasticMaterial::NewMaterial);

//отображение индексов: 11->1, 22->2, 33->3, 12->4, 23->5, 13->6
/*inline unsigned short int map1(unsigned short int i_)
{ switch (i_)
  {
   case 1: return 1;
   case 2: return 2;
   case 3: return 3;
   case 4: return 1;
   case 5: return 2;
  }
   return 1;
}
inline unsigned short int map2(unsigned short int i_)
{ switch (i_)
  {
   case 1: return 1;
   case 2: return 2;
   case 3: return 3;
   case 4: return 2;
   case 5: return 3;
  }
   return 3;
}*/

/*SymmetricTensor4s& tStrainIndependentMaterial::ElasTensor_1D (const Tensor2s& rot_,
SymmetricTensor4s& result_) const
{
// ofstream out("test.txt");
// out << "ElasTensor (matrix view):\n";
// for (int j=1; j<=6; ++j)
// {
//  for (int i=1; i<=j; ++i)
//   out << i << j << ": " << ElasTensor_prinAx(i,j) << ' ';
//  out << '\n';
// }
// out << "ElasTensor1d (matrix view):\n";
// for (int j=1; j<=6; ++j)
// {
//  for (int i=1; i<=j; ++i)
//   out << i << j << ": " << ElasTensor1D_prinAx(i,j) << ' ';
//  out << '\n';
// }

// static Tensor4s tmp; tmp = ElasTensor1D_prinAx;
 // tmp.TransformAll(rot_);

 result_ = ElasTensor1D_prinAx;
// Tensor2s rotT(rot_); rotT.Transpose();
 result_.TransformAll(rot_);
// result_.TransformAll(rotT);

// out << "ElasTensor1d (index view):\n";
// for (int i=1; i<=3; ++i)
// for (int j=1; j<=3; ++j)
// {
// for (int k=1; k<=3; ++k)
// for (int l=1; l<=3; ++l)
//   out << i << j << k << l << ": " << ElasTensor1D_prinAx(i,j,k,l) << " = " << tmp(i,j,k,l) << " =
" << result_(i,j,k,l) << "\n";
//  out << '\n';
// }

// out << "result (matrix view):\n";
// for (int j=1; j<=6; ++j)
// {
//  for (int i=1; i<=j; ++i)
//   out << i << j << ": " << result_(i,j) << ' ';
//  out << '\n';
// }

// out << "rotation:\n";
// for (int i=1; i<=3; ++i)
// {
//  for (int j=1; j<=3; ++j)
//    out << rot_(i,j) << '\t';
//  out << '\n';
// }
// out << "det = " << rot_.Det() << '\n';
// out.close();


// for (int j=1; j<=6; ++j)
// for (int i=1; i<=j; ++i) //   i.e. result_ = tmp;
// result_(i,j) = tmp(map1(i),map2(i),map1(j),map2(j));

 return result_;
}*/

const SymmetricTensor4s& tStrainIndependentMaterial::__RotateElasTensor(
	const SymmetricTensor4s& init_,
	const Tensor2s& rot_,
	tMapRotSymT4& db_)	//, SymmetricTensor4s& result_)
{
	static Tensor1s vec;
	rot_.AttachedVector(vec);
	tMapRotSymT4::iterator p = db_.find(
		make_pair(vec(1), make_pair(vec(2), vec(3))));	// modify to find with given precision!!!
	if (p == db_.end())
	{
#ifdef STRONGCHECK
		pair<tMapRotSymT4::iterator, bool> ins = db_.insert(tMapRotSymT4::value_type(
			make_pair(vec(1), make_pair(vec(2), vec(3))), SymmetricTensor4s(init_)));
		Assert(ins.second, "failed insertion in tStrainIndependentMaterial::RotateElasTensor");
		p = ins.first;
#else
		p = db_.insert(tMapRotSymT4::value_type(
						   make_pair(vec(1), make_pair(vec(2), vec(3))), SymmetricTensor4s(init_)))
				.first;
#endif
		p->second.TransformAll(rot_);
	}
	return p->second;
	// return result_=p->second;
}

tMaterial& tIsotropicLinearElasticMaterial::ReadData(istream& str_in_)
{
	string buffer;
	for (int i = 0; i < 3; ++i)
	{
		getline(str_in_, buffer, '=');
		str_in_.get();
		while (isspace(buffer[buffer.size() - 1])) buffer.resize(buffer.size() - 1);  // pop_back();
		while (isspace(buffer[0])) buffer.erase(0, 1);
		if (buffer == "E")
			str_in_ >> m_E;
		else if (buffer == "rho")
			str_in_ >> m_specificWeight;
		else if (buffer == "nu")
			str_in_ >> m_Nu;
		else
		{
			str_in_ >> buffer;
			--i;
		}
	}
	m_Mu = m_E / (1. + m_Nu);
	m_lambda = m_Mu * m_Nu / (1. - 2. * m_Nu);
	m_Mu *= 0.5;

	m_elasTensor_prinAx.Assign0();
	m_elasTensor_prinAx(1, 1) = m_elasTensor_prinAx(2, 2) = m_elasTensor_prinAx(3, 3) =
		m_lambda + 2. * m_Mu;
	m_elasTensor_prinAx(1, 2) = m_elasTensor_prinAx(1, 3) = m_elasTensor_prinAx(2, 3) = m_lambda;
	m_elasTensor_prinAx(4, 4) = m_elasTensor_prinAx(5, 5) = m_elasTensor_prinAx(6, 6) = m_Mu;

	m_elasTensor1D_prinAx = m_elasTensor_prinAx;
	m_elasTensor1D_prinAx(1, 1) = m_E;

	m_elasTensor2D_prinAx = m_elasTensor_prinAx;
	m_elasTensor2D_prinAx(1, 1) = m_elasTensor2D_prinAx(2, 2) =
		m_Mu * (m_lambda / (.5 * m_lambda + m_Mu) + 2.);
	m_elasTensor2D_prinAx(1, 2) = m_lambda * m_Mu / (.5 * m_lambda + m_Mu);

	return *this;
}

// tMaterial& tIsotropicPlasticMaterial::ReadData (InputFileStream& stream_)
// tMaterial& tIsotropicPlasticMaterial::ReadData (InputFileStreamXML& stream_)
//{
/* string str;
 SpecificWeight = stream_.ReadValue<real_t>(str);
 E  = stream_.ReadValue<real_t>(str);
 Nu  = stream_.ReadValue<real_t>(str);*/

/* SpecificWeight = atof(stream_.ReadElement("Rho").c_str());
 E  = atof(stream_.ReadElement("E").c_str());
 Nu  = atof(stream_.ReadElement("Nu").c_str());*/

/* mu = E / (1. + Nu),
 lambda = mu*Nu / (1. - 2.*Nu);
 mu *= 0.5;*/

/* streampos position;
 int n = stream_.HowManyBetween("{", "}", "/", position);
 stream_.seekg(position);
 Sigma.resize(n, 0.);
 Epsilon.resize(n, 0.);
 for(int i = 0; i < n; ++i)
   stream_.ReadPair<real_t>('/',Sigma[i],Epsilon[i]);*/

/* size_t n = stream_.ChildElemNum("SIGMA_EPSILON","Sigma");
 Sigma.resize(n, 0.);
 Epsilon.resize(n, 0.);
 for(size_t i = 0; i < n; ++i)
 {
   Sigma[i] = atof(stream_.ReadElement("Sigma").c_str());
   Epsilon[i] = atof(stream_.ReadElement("Epsilon").c_str());
 }

 return *this;
}*/

tMaterial& tIsotropicPlasticMaterial::ReadData(istream& str_in_)
{
	str_in_ >> m_specificWeight >> m_E >> m_Nu;

	string DiagramType;
	str_in_ >> DiagramType;
	if (DiagramType == "Pointwise")
	{
		str_in_.peek();

		while (!str_in_.eof())
		{
			real_t tmp;
			str_in_ >> tmp;
			m_sigma.push_back(tmp);
			str_in_ >> tmp;
			m_epsilon.push_back(tmp);
			str_in_.peek();
		}
	}

	return *this;
}

real_t tIsotropicPlasticMaterial::YoungModule(real_t Epsi_) const
{
#ifdef STRONGCHECK
	Assert(Epsi_ >= m_epsilon.front(), "Invalid sigma/epsilon diagram");
#endif
	vector<real_t>::const_iterator p = find_if(
		m_epsilon.begin(), m_epsilon.end(), bind(greater<real_t>(), placeholders::_1, Epsi_));
	if (p == m_epsilon.end())
		return (m_sigma.back() - *(m_sigma.end() - 2)) /
			   (m_epsilon.back() - *(m_epsilon.end() - 2));
	const int dist = static_cast<int>(distance(m_epsilon.begin(), p));
	return m_sigma[dist] - m_sigma[dist - 1] / (*p - *(p - 1));
}

const SymmetricTensor4s& tIsotropicPlasticMaterial::ElasTensor(real_t strain_level_) const
{
	// real_t epsilon = StrainIntensity(strain_);
	map<real_t, SymmetricTensor4s>::iterator p = m_elasTensorDB.find(strain_level_);
	if (p == m_elasTensorDB.end())
	{
#ifdef STRONGCHECK
		pair<map<real_t, SymmetricTensor4s>::iterator, bool> ins = m_elasTensorDB.insert(
			map<real_t, SymmetricTensor4s>::value_type(strain_level_, SymmetricTensor4s()));
		Assert(
			ins.second,
			"failed insertion of ElasTensor in tIsotropicPlasticMaterial::ElasTensor(const "
			"Tensor2s&)");
		p = ins.first;
#else
		p = ElasTensor_DB
				.insert(
					map<real_t, SymmetricTensor4s>::value_type(strain_level_, SymmetricTensor4s()))
				.first;
#endif
		SymmetricTensor4s& result = p->second;
		result.Assign0();

		real_t nu = PoissonRatio(strain_level_), tmp = YoungModule(strain_level_) / (1. + nu);
		result(4, 4) = result(5, 5) = result(6, 6) = .5 * tmp;
		tmp /= (1. - 2. * nu);
		result(1, 1) = result(2, 2) = result(3, 3) = (1. - nu) * tmp;
		result(1, 2) = result(1, 3) = result(2, 3) = nu * tmp;
		//    for (small_t i=1, j; i<=3; ++i)
		//     for (j=1; j<=3; ++j)
		//       {
		//        result(i,i,j,j) += lambda;
		//        result(i,j,i,j) += mu;
		//        result(i,j,j,i) += mu;
		//       }
	}
	return p->second;
}
