// Author: Chekhov Vladimir Valerevich

//------------------------------------------------------------------------------
#include "Tensors.h"
//#include "ThrowMessage.h"
#include "math_constants.hpp"
#include "test_runner.h"

//------------------------------------------------------------------------------
using namespace std;
//------------------------------------------------------------------------------
//==========================Tensor1s<T>============================================
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor1s<T>::DirectProduct(const Tensor1s<T>& o_, Tensor2s<T>& result_) const
{
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = Data[i] * o_[j];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor1s<T>::DotMultiply(const Tensor2s<T>& o_)
{
	T data_0 = Data[0], data_1 = Data[1];
	for (small_t i = 0; i < 3; ++i)
		Data[i] = data_0 * o_[0][i] + data_1 * o_[1][i] + Data[2] * o_[2][i];
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor1s<T>::DotMultiply_left(const Tensor2s<T>& o_)
{
	T data_0 = Data[0], data_1 = Data[1];
	for (small_t i = 0; i < 3; ++i)
		Data[i] = o_[i][0] * data_0 + o_[i][1] * data_1 + o_[i][2] * Data[2];
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor1s<T>::DirectProduct(const Tensor2s<T>& o_, Tensor3s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = Data[i] * o_[j][k];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor1s<T>::DotProduct(const Tensor3s<T>& o_, Tensor2s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k) result_[i][j] += Data[k] * o_[k][i][j];
		}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor1s<T>::DirectProduct(const Tensor3s<T>& o_, Tensor4s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = Data[i] * o_[j][k][l];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor1s<T>::DotProduct(const Tensor4s<T>& o_, Tensor3s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l) result_[i][j][k] += Data[l] * o_[l][i][j][k];
			}
	return result_;
}
//------------------------------------------------------------------------------
//==========================Tensor1a<T>=============================================
//------------------------------------------------------------------------------
template<typename T>
Tensor1a<T>& Tensor1a<T>::operator-=(const Tensor1a<T>& o_)
{
#ifdef STRONGCHECK
	Assert(this != &o_, "the same parameter in Tensor1a<T>::operator-=");
#endif
	// if (&o_ == this) return Assign0();
	if (__AdjustAdd(o_))
		pData->operator-=(*o_.pData);
	else if (!is0())
		pData->Negate();
	return *this;
}
//------------------------------------------------------------------------------
//==========================Tensor2s<T>=======================================
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>::Tensor2s(const Tensor1s<T>& col1_, const Tensor1s<T>& col2_, const Tensor1s<T>& col3_)
	: Data{}, Rows{}
{
	__SET_ROWS;
	for (small_t i = 0; i < 3; ++i)
	{
		Rows[i][0] = col1_.Data[i];
		Rows[i][1] = col2_.Data[i];
		Rows[i][2] = col3_.Data[i];
	}
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>::Tensor2s(const Tensor1s<T>& col1_, const Tensor1s<T>& col2_, const Tensor1s<T>& col3_, bool)
	: Data{}, Rows{}
{
	__SET_ROWS;
	for (small_t i = 0; i < 3; ++i)
	{
		Rows[0][i] = col1_.Data[i];
		Rows[1][i] = col2_.Data[i];
		Rows[2][i] = col3_.Data[i];
	}
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::operator-=(const Tensor2s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 9; ++pd, ++po) *pd -= *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::operator+=(const Tensor2s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 9; ++pd, ++po) *pd += *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::AssignRows(const Tensor1s<T>& row1_, const Tensor1s<T>& row2_, const Tensor1s<T>& row3_)
{
	T *p0 = *Rows, *p1 = Rows[1], *p2 = Rows[2];
	const T *psrc0 = row1_.Data, *psrc1 = row2_.Data, *psrc2 = row3_.Data;
	for (; p0 < Rows[1]; ++p0, ++p1, ++p2, ++psrc0, ++psrc1, ++psrc2)
	{
		*p0 = *psrc0;
		*p1 = *psrc1;
		*p2 = *psrc2;
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::AssignCols(const Tensor1s<T>& col1_, const Tensor1s<T>& col2_, const Tensor1s<T>& col3_)
{
	T *p0 = *Rows, *p1 = Rows[1], *p2 = Rows[2];
	const T *psrc0 = col1_.Data, *psrc1 = col2_.Data, *psrc2 = col3_.Data;
	for (; p0 <= Rows[2]; p0 += 3, p1 += 3, p2 += 3, ++psrc0, ++psrc1, ++psrc2)
	{
		*p0 = *psrc0;
		*p1 = *psrc1;
		*p2 = *psrc2;
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
T Tensor2s<T>::Minor(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3, "invalid index in Tensor2s<T>::Minor");
#endif
	small_t i1 = 4 - i_, j1 = i1 == 1 ? 1 : 2, i2 = 4 - j_, j2 = i2 == 1 ? 1 : 2;
	i1 -= j1;
	i2 -= j2;
	return Rows[i1][i2] * Rows[j1][j2] - Rows[i1][j2] * Rows[j1][i2];
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::Invert()
{
	T det1 = Det();
#ifdef STRONGCHECK
	Assert(NE0(det1), tMessage("singular tensor in Tensor2s<T>::Invert\ndet = ") << det1);
#endif
	det1 = 1. / det1;
	T tmp[9] = {(Data[4] * Data[8] - Data[5] * Data[7]) * det1,
					 (Data[2] * Data[7] - Data[1] * Data[8]) * det1,
					 (Data[1] * Data[5] - Data[2] * Data[4]) * det1,
					 (Data[5] * Data[6] - Data[3] * Data[8]) * det1,
					 (Data[0] * Data[8] - Data[2] * Data[6]) * det1,
					 (Data[2] * Data[3] - Data[0] * Data[5]) * det1,
					 (Data[3] * Data[7] - Data[4] * Data[6]) * det1,
					 (Data[1] * Data[6] - Data[0] * Data[7]) * det1,
					 (Data[0] * Data[4] - Data[1] * Data[3]) * det1};
	copy(tmp, tmp + 9, Data);
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
 void Tensor2s<T>::EigenValues(T eigVals_[3]) const
{
	// Characteristic equation eigVal^3 - I1*eigVal^2 + I2*eigVal - I3 =0
	// reduced to x^3 + p*x + q =0  by substitution {x = eigVal-I1/3}:
	const T i1 = I1(), i2 = I2(), i3 = I3(), i1_3 = i1 / 3.l, p = i2 - i1 * i1_3,
				 q = i1_3 * (i2 - 2. * i1_3 * i1_3) - i3;	// coefs of reduced eqn
	Assert(mathdef::is_lte(p, 0.), "Eigen values with Im!=0 in Tensor2s<T>::EigenValues(T[])");
	const T ro = sqrt(-p / 3.l);  // auxiliary param
	if (mathdef::is_not_zero(ro))
	{
		const T fi = acos(-q / (2. * ro * ro * ro)) / 3.l,  // auxiliary param
			pi23 = acos(-0.5);									  //=2Pi/3
		eigVals_[0] = 2. * ro * cos(fi) + i1_3;
		eigVals_[1] = 2. * ro * cos(fi + pi23) + i1_3;
		eigVals_[2] = 2. * ro * cos(fi + pi23 * 2.) + i1_3;
	}
	else
		eigVals_[0] = eigVals_[1] = eigVals_[2] =
			(sign(pow(fabs(q), T(1. / 3.l)), -q) + i1_3);
}
//------------------------------------------------------------------------------
template<typename T>
struct less_in_abs
{
	bool operator()(const T& a, const T& b) const
	{
		return fabs(a) < fabs(b);
	}
};
//------------------------------------------------------------------------------
template<typename T>
void Tensor2s<T>::MaxAbsIndex(small_t& row_num_, small_t& col_num_) const
{
	const small_t lin_num = static_cast<small_t>(max_element(Data, Data + 9, less_in_abs<T>()) - Data);
	row_num_ = lin_num / 3;
	col_num_ = lin_num % 3;
}
//------------------------------------------------------------------------------
template<typename T>
small_t Tensor2s<T>::EigenVectors(T eigVal_, Tensor2s<T>& eigVects_, small_t begPos_) const
{
#ifdef STRONGCHECK
	Assert(
		begPos_ > 0 && begPos_ <= 3, "invalid position of eigenvector in Tensor2s<T>::EigenVectors");
#endif
	Tensor2s<T> singular(*this);
	T tmp;
	small_t i = 0, j, i1, i2, j1, j2, multiplicity = 1;
	const small_t posOfCurrent = begPos_ - 1;
	for (; i < 3; ++i) singular[i][i] -= eigVal_;  // i.e. singular = *this - eigval*E
#ifdef STRONGCHECK
	Assert(
		EQ0(singular.Det(), PRECISION * 50.l),
		tMessage(
			"value expected to be eigen is not eigen in Tensor2s<T>::EigenVectors\ndet(A-lambdaE) = ")
			<< singular.Det() << " (abs must be <= " << PRECISION * 20.l << ')');
#endif
	// Index index;
	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			tmp = singular.Minor(i + 1, j + 1);
			if (mathdef::is_not_zero(tmp))	// If exist nonzero minor of 2nd order
			{  // then rang of singular equals 2 and current eigen value is 1-fold:
				i1 = i == 0 ? 1 : 0;
				i2 = 3 - i - i1;
				j1 = j == 0 ? 1 : 0;
				j2 = 3 - j - j1;
				eigVects_[j][posOfCurrent] = tmp;
				eigVects_[j1][posOfCurrent] =
					singular[i1][j2] * singular[i2][j] - singular[i2][j2] * singular[i1][j];
				eigVects_[j2][posOfCurrent] =
					singular[i2][j1] * singular[i1][j] - singular[i1][j1] * singular[i2][j];
				goto normalizing;
			}
		}
#ifdef STRONGCHECK
	Assert(begPos_ <= 2, "invalid position of eigenvector in Tensor2s<T>::EigenVectors");
#endif
	// If all minors of 2nd order are zero then rang=1 and eigen value is 2-fold (or 3-fold):
	MaxAbsIndex(i, j);
	j1 = j == 0 ? 1 : 0;
	j2 = 3 - j - j1;
	eigVects_[j1][posOfCurrent] = eigVects_[j2][posOfCurrent + 1] = 1.;
	eigVects_[j2][posOfCurrent] = eigVects_[j1][posOfCurrent + 1] = 0.;
	if (mathdef::is_not_zero(singular[i][j]))
	{
		eigVects_[j][posOfCurrent] = -singular[i][j1] / singular[i][j];
		eigVects_[j][posOfCurrent + 1] = -singular[i][j2] / singular[i][j];
		multiplicity = 2;
	}
	else
	{  // if tensor is spherical then eigen value is 3-fold and any vector is eigen one:
#ifdef STRONGCHECK
		Assert(begPos_ == 1, "invalid position of eigenvector in Tensor2s<T>::EigenVectors");
#endif
		eigVects_[j][0] = eigVects_[j][1] = eigVects_[j1][2] = eigVects_[j2][2] = 0.;
		eigVects_[j][2] = 1.;
		return 3;
	}
normalizing:  // make all eigen vectors to have unit length:
	for (i = posOfCurrent; i < posOfCurrent + multiplicity; ++i)
	{
		tmp = sqrt(pow(eigVects_[0][i], 2) + pow(eigVects_[1][i], 2) + pow(eigVects_[2][i], 2));
		for (j = 0; j < 3; ++j) eigVects_[j][i] /= tmp;
	}
	return multiplicity;
}
//------------------------------------------------------------------------------
template<typename T>
void Tensor2s<T>::Eigen(T eigvals_[3], Tensor2s<T>& eigvects_) const
{
	EigenValues(eigvals_);
	small_t eig_vects_count = EigenVectors(eigvals_[0], eigvects_);
	if (eig_vects_count == 3)
		return;
	small_t next_eig_val_no = mathdef::is_eq(eigvals_[0], eigvals_[1]) ? 2 : 1;
	EigenVectors(eigvals_[next_eig_val_no], eigvects_, eig_vects_count + 1);
#ifdef STRONGCHECK
	eig_vects_count += EigenVectors(eigvals_[next_eig_val_no], eigvects_, eig_vects_count + 1);
	Assert(
		next_eig_val_no < 2 || eig_vects_count == 3,
		"eigen vectors are not correspond to eigen values in Tensor2s<T>::Eigen");
#endif
	if (eig_vects_count < 3)
		EigenVectors(eigvals_[2], eigvects_, 3);
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::RotateToEigenSystem(Tensor2s<T>& rotation_)
{  // assigns columns of arg to eigenvectors and assigns *this to diagonal form
	T eigenvals[3];
	Eigen(eigenvals, rotation_);
	Assign0();
	for (small_t i = 0; i < 3; ++i) Data[3 * i] = eigenvals[i];
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::DotMultiply(const Tensor2s<T>& o_)	 //(T*F)[ij]=T[ip]*F[pj]
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "mult by itself in Tensor2s<T>::DotMultiply(const Tensor2s<T>&)");
#endif
	// if (&o_ == this)
	//       return DotMultiply(Tensor2s<T>(*this));
	T rows_i0, rows_i1;
	for (small_t i = 0; i < 3; ++i)
	{
		rows_i0 = Rows[i][0];
		rows_i1 = Rows[i][1];
		Rows[i][0] = rows_i0 * o_[0][0] + rows_i1 * o_[1][0] + Rows[i][2] * o_[2][0];
		Rows[i][1] = rows_i0 * o_[0][1] + rows_i1 * o_[1][1] + Rows[i][2] * o_[2][1];
		Rows[i][2] = rows_i0 * o_[0][2] + rows_i1 * o_[1][2] + Rows[i][2] * o_[2][2];
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::DotMultiply_left(const Tensor2s<T>& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "mult by itself in Tensor2s<T>::DotMultiply_left(const Tensor2s<T>&)");
#endif
	// if (&o_ == this)
	//        return DotMultiply_left(Tensor2s<T>(Data));
	T rows_0j, rows_1j;
	for (small_t j = 0; j < 3; ++j)
	{
		rows_0j = Rows[0][j];
		rows_1j = Rows[1][j];
		Rows[0][j] = o_[0][0] * rows_0j + o_[0][1] * rows_1j + o_[0][2] * Rows[2][j];
		Rows[1][j] = o_[1][0] * rows_0j + o_[1][1] * rows_1j + o_[1][2] * Rows[2][j];
		Rows[2][j] = o_[2][0] * rows_0j + o_[2][1] * rows_1j + o_[2][2] * Rows[2][j];
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::CrossMultiply(const Tensor1s<T>& o_)  //(T*f)[ij]=T[ip]*f[q]*LevyChivita[pqj]
{
	T rows_i0, rows_i1;
	for (small_t i = 0; i < 3; ++i)
	{
		rows_i0 = Rows[i][0];
		rows_i1 = Rows[i][1];
		Rows[i][0] = rows_i1 * o_[2] - Rows[i][2] * o_[1];
		Rows[i][1] = Rows[i][2] * o_[0] - rows_i0 * o_[2];
		Rows[i][2] = rows_i0 * o_[1] - rows_i1 * o_[0];
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::CrossMultiply_left(const Tensor1s<T>& o_)	//(t*F)[ij]=t[p]*F[qj]*LevyChivita[pqi]
{
	T rows_0j, rows_1j;
	for (small_t j = 0; j < 3; ++j)
	{
		rows_0j = Rows[0][j];
		rows_1j = Rows[1][j];
		Rows[0][j] = o_[1] * Rows[2][j] - o_[2] * rows_1j;
		Rows[1][j] = o_[2] * rows_0j - o_[0] * Rows[2][j];
		Rows[2][j] = o_[0] * rows_1j - o_[1] * rows_0j;
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor2s<T>::DirectProduct(const Tensor1s<T>& o_, Tensor3s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = Rows[i][j] * o_[k];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
T Tensor2s<T>::Dot2Product(const Tensor2s<T>& o_) const
{
	T result = 0.;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result += Rows[i][j] * o_[j][i];
	return result;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor2s<T>::CrossProduct(const Tensor2s<T>& o_, Tensor3s<T>& result_) const
{
	const Tensor3s<T> levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						result_[i][j][k] += Rows[i][l] * o_[m][k] * levyChivita[l][m][j];
			}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor2s<T>::DirectProduct(const Tensor2s<T>& o_, Tensor4s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = Rows[i][j] * o_[k][l];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor2s<T>::Dot2Product(const Tensor3s<T>& o_, Tensor1s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i] += Rows[j][k] * o_[k][j][i];
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor2s<T>::CrossProduct(const Tensor3s<T>& o_, Tensor4s<T>& result_) const
{
	const Tensor3s<T> levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					result_[i][j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							result_[i][j][k][l] += Rows[i][m] * o_[n][k][l] * levyChivita[m][n][j];
				}
	return result_;
}
//------------------------------------------------------------------------------
// Tensor2s<T>& Tensor2s<T>::Dot2Multiply (const Tensor4& o_, small_t index1No_, small_t index2No_)
template<typename T>
Tensor2s<T>& Tensor2s<T>::Dot2Multiply(const SymmetricTensor4s<T>& o_, small_t index1No_, small_t index2No_)
{
#ifdef STRONGCHECK
	Assert(
		index1No_ != index2No_ && index1No_ > 0 && index1No_ <= 4 && index2No_ > 0 &&
			index2No_ <= 4,
		"invalid index # in Tensor2s<T>::Dot2Multiply(SymmetricTensor4s<T>&,small_t,small_t)");
#endif
	typename Tensor4<T>::tIndex ijkl;
	const Tensor2s<T> copyOfThis(*this);
	const small_t freeIndex1No =
					  index1No_ + index2No_ == 3 ? 3 : (index1No_ == 1 || index2No_ == 1 ? 2 : 1),
				  freeIndex2No = 10 - freeIndex1No - index1No_ - index2No_;
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		ijkl[freeIndex1No] = i;
		for (j = 0; j < 3; ++j)
		{
			ijkl[freeIndex2No] = j;
			Rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
			{
				ijkl[index1No_] = k;
				for (l = 0; l < 3; ++l)
				{
					ijkl[index2No_] = l;
					Rows[i][j] += copyOfThis[k][l] * o_[ijkl];
				}
			}
		}
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::Dot2Multiply(const Tensor4s<T>& o_)
{
	const Tensor2s<T> copyOfThis(*this);
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			Rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) Rows[i][j] += copyOfThis[k][l] * o_[l][k][i][j];
		}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor2s<T>::Dot2Multiply_left(const Tensor4s<T>& o_)
{
	const Tensor2s<T> copyOfThis(*this);
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			Rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) Rows[i][j] += o_[i][j][k][l] * copyOfThis[l][k];
		}
	return *this;
}
//------------------------------------------------------------------------------
/*Tensor2& Tensor2s<T>::ApplyFunction (const Tensor2s<T>::fClassicFunction& function_)
{//according to (Lurie, 1980), page 436.
 Tensor1s<T> eigVals, funOfEigVals;   EigenValues(eigVals);
 try{
	 funOfEigVals[0] = function_(eigVals[0]);
	 funOfEigVals[1] = function_(eigVals[1]);
	 funOfEigVals[2] = function_(eigVals[2]);
	}
 catch (...) {throw tMessage("Invalid arg in Tensor2s<T>::ApplyFunction (maybe divide by 0 etc.)");}
 Tensor2s<T> currentEigDiade, result(0.);
 small_t i=0,j,k,l;   T tmp;//   Index index;
 try{
	 for (; i<3; ++i)
		 {
		  j = i==0?1:0; k = 3-i-j;
		  (currentEigDiade.Assign1(-eigVals[j]-eigVals[k]) += *this) *= *this;
		  tmp = eigVals[j]*eigVals[k];
		  for (l=0; l<3; ++l) currentEigDiade[l][l] += (tmp);
		  tmp = funOfEigVals[i] / ((eigVals[i]-eigVals[j])*(eigVals[i]-eigVals[k]));
		  currentEigDiade *= tmp;//Direct product of i-th eigenvector to i-th eigenvector of the
mutual basis is multiplied by the function of i-th eigenvalue result += currentEigDiade;
		 }
	}
 catch (tMessage) {throw;}
 catch (...)//if 2- or 3-fold eigen value
	{//if 2-fold eigen value
	 if (!isSymmetric()) throw tMessage("Fold eigen values for non-simmetric tensor");
	 j=0;//i - single eigvalue's number, j - number of one of same eigvalues. Finding of i,j:
	 if (fabs(eigVals[0]-eigVals[1]) < fabs(eigVals[0]-eigVals[2]))   i=2;
	 else                                                           i=1;
	 if (fabs(eigVals[1]-eigVals[2]) < fabs(eigVals[0]-eigVals[3-i])){i=0; j=1;}
	 tmp = eigVals[j];//2-fold eigenvalue
	 (result.Assign1(-tmp-tmp) += *this) *= *this;
	 tmp *= tmp;
	 for (l=0; l<3; ++l) result[l][l] += tmp;
	 tmp = (eigVals[i]-eigVals[j]); tmp*= tmp;
	 if (NE0(tmp))
	   {
		result /= tmp;
		tmp = funOfEigVals[j];
		result *= (funOfEigVals[i] - tmp);
		for (l=0; l<3; ++l) result[l][l] += tmp;
	   }//result = E*f(fold_eigval) +
(f(not_fold_eigval)-f(fold_eigval))*eigvect_diade_of_not_fold_eigval else//if 3-fold eigen value
		 result.Assign1(funOfEigVals[1]);
	}
 return operator=(result);
} */
//------------------------------------------------------------------------------
//==========================Tensor2a<T>=============================================
//------------------------------------------------------------------------------
/*ostream& operator<< (ostream& out_, const Tensor2a<T>& o_)
{

 for (small_t i=1,j; i<=3; ++i)
   {
	out_ << "|";
	out_ << setw(26) << setprecision(19) << o_(i,1);
	for (j=2; j<=3;  ++j)
		   out_ << ' ' << setw(26) << setprecision(19) << o_(i,j);
	out_ << "|\n";
   }
 out_ << '\n';
 return out_.flush();
}*/
//------------------------------------------------------------------------------
template<typename T>
Tensor2a<T>& Tensor2a<T>::operator-=(const Tensor2a<T>& o_)
{
#ifdef STRONGCHECK
	Assert(this != &o_, "the same parameter in Tensor2a<T>::operator-=");
#endif
	// if (this == &o_) return Assign0();
	if (__AdjustAdd(o_))
		pData->operator-=(*o_.pData);
	else if (!is0())
		pData->Negate();
	return *this;
}
//------------------------------------------------------------------------------
//==========================Tensor3s<T>==========================================
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::operator-=(const Tensor3s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 27; ++pd, ++po) *pd -= *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::operator+=(const Tensor3s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 27; ++pd, ++po) *pd += *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::Assign1(T v_)	// Levy-Chivita multiplied by v_
{
	Assign0();
	(*Rows[0])[1][2] = (*Rows[1])[2][0] = (*Rows[2])[0][1] = v_;
	(*Rows[0])[2][1] = (*Rows[2])[1][0] = (*Rows[1])[0][2] = -v_;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor3s<T>::Contraction(small_t n1_, small_t n2_, Tensor1s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(
		n1_ > 0 && n2_ > 0 && n1_ <= 3 && n2_ <= 3 && n1_ != n2_,
		"invalid number of index in Tensor3s<T>::Contraction");
#endif
	const small_t k = 6 - n1_ - n2_;
	tIndex index;
	for (small_t i = 0, j; i < 3; ++i)
	{
		result_[i] = 0.;
		index[k] = i;
		for (j = 0; j < 3; ++j)
		{
			index[n1_] = index[n2_] = j;
			result_[i] += operator[](index);
		}
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor3s<T>::ScalarProductWithBasisVector(small_t k_, Tensor2s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(k_ > 0 && k_ <= 3, "invalid index in Tensor3s<T>::ScalarProductWithBasisVector");
#endif
	--k_;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = (*Rows[i])[j][k_];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor3s<T>::ScalarProductWithBasisVector_left(small_t k_, Tensor2s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(k_ > 0 && k_ <= 3, "invalid index in Tensor3s<T>::ScalarProductWithBasisVector_left");
#endif
	--k_;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = (*Rows[k_])[i][j];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor3s<T>::DotProduct(const Tensor1s<T>& o_, Tensor2s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k) result_[i][j] += (*Rows[i])[j][k] * o_[k];
		}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::CrossMultiply(const Tensor1s<T>& o_)
{
	const Tensor3s<T> copyOfThis(*this), levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*Rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						(*Rows[i])[j][k] += copyOfThis[i][j][l] * o_[m] * levyChivita[l][m][k];
			}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::CrossMultiply_left(const Tensor1s<T>& o_)
{
	const Tensor3s<T> copyOfThis(*this), levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*Rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						(*Rows[i])[j][k] += o_[l] * copyOfThis[m][j][k] * levyChivita[l][m][i];
			}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor3s<T>::DirectProduct(const Tensor1s<T>& o_, Tensor4s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = (*Rows[i])[j][k] * o_[l];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor3s<T>::Dot2Product(const Tensor2s<T>& o_, Tensor1s<T>& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i] += (*Rows[i])[j][k] * o_[k][j];
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::DotMultiply(const Tensor2s<T>& o_)
{
	T data_ij0, data_ij1;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			data_ij0 = (*Rows[i])[j][0];
			data_ij1 = (*Rows[i])[j][1];
			(*Rows[i])[j][0] =
				data_ij0 * o_[0][0] + data_ij1 * o_[1][0] + (*Rows[i])[j][2] * o_[2][0];
			(*Rows[i])[j][1] =
				data_ij0 * o_[0][1] + data_ij1 * o_[1][1] + (*Rows[i])[j][2] * o_[2][1];
			(*Rows[i])[j][2] =
				data_ij0 * o_[0][2] + data_ij1 * o_[1][2] + (*Rows[i])[j][2] * o_[2][2];
		}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::DotMultiply_left(const Tensor2s<T>& o_)
{
	T data_0jk, data_1jk;
	for (small_t j = 0, k; j < 3; ++j)
		for (k = 0; k < 3; ++k)
		{
			data_0jk = (*Rows[0])[j][k];
			data_1jk = (*Rows[1])[j][k];
			(*Rows[0])[j][k] =
				o_[0][0] * data_0jk + o_[0][1] * data_1jk + o_[0][2] * (*Rows[2])[j][k];
			(*Rows[1])[j][k] =
				o_[1][0] * data_0jk + o_[1][1] * data_1jk + o_[1][2] * (*Rows[2])[j][k];
			(*Rows[2])[j][k] =
				o_[2][0] * data_0jk + o_[2][1] * data_1jk + o_[2][2] * (*Rows[2])[j][k];
		}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
T Tensor3s<T>::Dot3Product(const Tensor3s<T>& o_) const
{
	T result = 0.;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result += (*Rows[i])[j][k] * o_[k][j][i];
	return result;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor3s<T>::Dot2Product(const Tensor3s<T>& o_, Tensor2s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j] += (*Rows[i])[k][l] * o_[l][k][j];
		}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor3s<T>::DotProduct(const Tensor3s<T>& o_, Tensor4s<T>& result_) const
{
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					result_[i][j][k][l] = 0.;
					for (m = 0; m < 3; ++m) result_[i][j][k][l] += (*Rows[i])[j][m] * o_[m][k][l];
				}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor3s<T>::Dot3Product(const Tensor4s<T>& o_, Tensor1s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i] += (*Rows[j])[k][l] * o_[l][k][j][i];
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::Dot2Multiply(const Tensor4s<T>& o_)
{
	Tensor2s<T> tmp;
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) tmp[j][k] = (*Rows[i])[j][k];
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*Rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) (*Rows[i])[j][k] += tmp[l][m] * o_[m][l][j][k];
			}
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor3s<T>::Dot2Multiply_left(const Tensor4s<T>& o_)
{
	Tensor2s<T> tmp;
	for (small_t i, j, k = 0, l, m; k < 3; ++k)
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j) tmp[i][j] = (*Rows[i])[j][k];
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				(*Rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) (*Rows[i])[j][k] += o_[i][j][l][m] * tmp[m][l];
			}
	}
	return *this;
}
//------------------------------------------------------------------------------
//==========================Tensor4s<T>==========================================
//------------------------------------------------------------------------------
/*Tensor4s<T>& Tensor4s<T>::operator= (const Tensor4& o_)
{
 for (small_t i=0,j,k,l; i<3; ++i)
  for   (j=0; j<3; ++j)
   for  (k=0; k<3; ++k)
	for (l=0; l<3; ++l)
	 (*Rows[i])[j][k][l] = o_[tIndex(i,j,k,l)];
//     (*Rows[i])[j][k][l] = o_(i+1,j+1,k+1,l+1);
 return *this;
} */
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::operator=(const SymmetricTensor4s<T>& o_)
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		for (k = 0; k < 3; ++k)
		{
			(*Rows[i])[i][k][k] = o_[tIndex(i, i, k, k)];
			for (l = 0; l < k; ++l)
				(*Rows[i])[i][k][l] = (*Rows[i])[i][l][k] = o_[tIndex(i, i, k, l)];
		}
		for (j = 0; j < i; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*Rows[i])[j][k][k] = (*Rows[j])[i][k][k] = o_[tIndex(i, j, k, k)];
				for (l = 0; l < k; ++l)
					(*Rows[i])[j][k][l] = (*Rows[i])[j][l][k] = (*Rows[j])[i][k][l] =
						(*Rows[j])[i][l][k] = o_[tIndex(i, j, k, l)];
			}
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::operator-=(const Tensor4s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 81; ++pd, ++po) *pd -= *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::operator+=(const Tensor4s<T>& o_)
{
	const T* po = o_.Data;
	for (T* pd = Data; pd < Data + 81; ++pd, ++po) *pd += *po;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::Assign1(small_t kind_, T v_)  // Isotropic multiplied by v_
{
#ifdef STRONGCHECK
	Assert(kind_ > 0 && kind_ <= 3, "invalid of isotropic Tensor4s<T> in Tensor4s<T>::Assign1");
#endif
	Assign0();
	small_t i = 0;
	small_t j = 0;
	const small_t &k = (kind_ == 1) ? i : j, &l = (kind_ == 2) ? i : j, &m = (kind_ == 3) ? i : j;
	for (; i < 3; ++i)
		for (j = 0; j < 3; ++j) (*Rows[i])[k][l][m] = v_;
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor4s<T>::Contraction(small_t n1_, small_t n2_, Tensor2s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(
		n1_ > 0 && n2_ > 0 && n1_ <= 4 && n2_ <= 4 && n1_ != n2_,
		"invalid number of index in Tensor4s<T>::Contraction");
#endif
	small_t k1 = 10 - n1_ - n2_, k2 = k1;
	tIndex index4;
	if (k1 > 5)
	{
		k2 = 4;
		k1 -= k2;
	}
	else if (k1 < 5)
	{
		k1 = 1;
		k2 -= k1;
	}
	else if (n1_ == 1 || n2_ == 1)
	{
		k1 = 2;
		k2 = 3;
	}
	else
	{
		k1 = 1;
		k2 = 4;
	}
	for (small_t i = 0, j, k; i < 3; ++i)
	{
		index4[k1] = i;
		for (j = 0; j < 3; ++j)
		{
			index4[k2] = j;
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k)
			{
				index4[n1_] = index4[n2_] = k;
				result_[i][j] += operator[](index4);
			}
		}
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor4s<T>::ScalarProductWithBasisVector(small_t l_, Tensor3s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(l_ > 0 && l_ <= 3, "invalid index in Tensor4s<T>::ScalarProductWithBasisVector");
#endif
	--l_;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = (*Rows[i])[j][k][l_];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor4s<T>::ScalarProductWithBasisVector_left(small_t l_, Tensor3s<T>& result_) const
{
#ifdef STRONGCHECK
	Assert(l_ > 0 && l_ <= 3, "invalid index in Tensor4s<T>::ScalarProductWithBasisVector_left");
#endif
	--l_;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = (*Rows[l_])[i][j][k];
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor3s<T>& Tensor4s<T>::DotProduct(const Tensor1s<T>& o_, Tensor3s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l) result_[i][j][k] += (*Rows[i])[j][k][l] * o_[l];
			}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::CrossMultiply(const Tensor1s<T>& o_)
{
	const Tensor4s<T> copyOfThis(*this);
	const Tensor3s<T> levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*Rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							(*Rows[i])[j][k][l] +=
								copyOfThis[i][j][k][m] * o_[n] * levyChivita[m][n][l];
				}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::CrossMultiply_left(const Tensor1s<T>& o_)
{
	const Tensor4s<T> copyOfThis(*this);
	const Tensor3s<T> levyChivita(TENSOR_KIND::UNIT);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*Rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							(*Rows[i])[j][k][l] +=
								o_[m] * copyOfThis[n][j][k][l] * levyChivita[m][n][i];
				}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::DotMultiply(const Tensor2s<T>& o_, small_t ownIndexNo_, small_t facIndexNo_)
{
#ifdef STRONGCHECK
	Assert(
		ownIndexNo_ > 0 && ownIndexNo_ <= 4 && facIndexNo_ > 0 && facIndexNo_ <= 2,
		"invalid index # in Tensor4s<T>::DotMultiply(Tensor4s<T>&,small_t,small_t)");
#endif
	T data_ijk0, data_ijk1;
	typename Tensor2s<T>::tIndex index2;
	typename Tensor2s<T>::tIndex index4;
	const small_t ownFreeNo1 = ownIndexNo_ == 1 ? 2 : 1,
				  ownFreeNo2 = ownIndexNo_ + ownFreeNo1 == 3 ? 3 : 2,
				  ownFreeNo3 = 10 - ownIndexNo_ - ownFreeNo1 - ownFreeNo2,
				  facFreeNo = 3 - facIndexNo_;
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		index4[ownFreeNo1] = i;
		for (j = 0; j < 3; ++j)
		{
			index4[ownFreeNo2] = j;
			for (k = 0; k < 3; ++k)
			{
				index4[ownFreeNo3] = k;
				data_ijk0 = operator[]((index4[ownIndexNo_] = 0, index4));
				data_ijk1 = operator[]((index4[ownIndexNo_] = 1, index4));
				for (l = 0; l < 3; ++l)
				{
					index2[facFreeNo] = l;
					operator[]((index4[ownIndexNo_] = l, index4)) =
						data_ijk0 * o_[(index2[facIndexNo_] = 0, index2)] +
						data_ijk1 * o_[(index2[facIndexNo_] = 1, index2)] +
						operator[]((index4[ownIndexNo_] = 2, index4)) *
							o_[(index2[facIndexNo_] = 2, index2)];
				}
			}
		}
	}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::DotMultiply(const Tensor2s<T>& o_)
{
	T data_ijk0, data_ijk1;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				data_ijk0 = (*Rows[i])[j][k][0];
				data_ijk1 = (*Rows[i])[j][k][1];
				(*Rows[i])[j][k][0] =
					data_ijk0 * o_[0][0] + data_ijk1 * o_[1][0] + (*Rows[i])[j][k][2] * o_[2][0];
				(*Rows[i])[j][k][1] =
					data_ijk0 * o_[0][1] + data_ijk1 * o_[1][1] + (*Rows[i])[j][k][2] * o_[2][1];
				(*Rows[i])[j][k][2] =
					data_ijk0 * o_[0][2] + data_ijk1 * o_[1][2] + (*Rows[i])[j][k][2] * o_[2][2];
			}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::DotMultiply_left(const Tensor2s<T>& o_)
{
	T data_0jkl, data_1jkl;
	for (small_t j = 0, k, l; j < 3; ++j)
		for (k = 0; k < 3; ++k)
			for (l = 0; l < 3; ++l)
			{
				data_0jkl = (*Rows[0])[j][k][l];
				data_1jkl = (*Rows[1])[j][k][l];
				(*Rows[0])[j][k][l] =
					o_[0][0] * data_0jkl + o_[0][1] * data_1jkl + o_[0][2] * (*Rows[2])[j][k][l];
				(*Rows[1])[j][k][l] =
					o_[1][0] * data_0jkl + o_[1][1] * data_1jkl + o_[1][2] * (*Rows[2])[j][k][l];
				(*Rows[2])[j][k][l] =
					o_[2][0] * data_0jkl + o_[2][1] * data_1jkl + o_[2][2] * (*Rows[2])[j][k][l];
			}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor1s<T>& Tensor4s<T>::Dot3Product(const Tensor3s<T>& o_, Tensor1s<T>& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i] += (*Rows[i])[j][k][l] * o_[l][k][j];
	}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
T Tensor4s<T>::Dot4Product(const Tensor4s<T>& o_) const
{
	T result = 0.;
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result += (*Rows[i])[j][k][l] * o_[l][k][j][i];
	return result;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor2s<T>& Tensor4s<T>::Dot3Product(const Tensor4s<T>& o_, Tensor2s<T>& result_) const
{
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) result_[i][j] += (*Rows[i])[k][l][m] * o_[m][l][k][j];
		}
	return result_;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::Dot2Multiply(const Tensor4s<T>& o_)
{
	Tensor2s<T> tmp;
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) tmp[k][l] = (*Rows[i])[j][k][l];
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*Rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n) (*Rows[i])[j][k][l] += tmp[m][n] * o_[n][m][k][l];
				}
		}
	return *this;
}
//------------------------------------------------------------------------------
template<typename T>
Tensor4s<T>& Tensor4s<T>::Dot2Multiply_left(const Tensor4s<T>& o_)
{
	Tensor2s<T> tmp;
	for (small_t i, j, k = 0, l, m, n; k < 3; ++k)
		for (l = 0; l < 3; ++l)
		{
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j) tmp[i][j] = (*Rows[i])[j][k][l];
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					(*Rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n) (*Rows[i])[j][k][l] += o_[i][j][m][n] * tmp[n][m];
				}
		}
	return *this;
}
//------------------------------------------------------------------------------
/*ostream& operator<< (ostream& out_, Tensor4& t_)
{
 for (small_t i=0,j,k; i<3; ++i)
 for (j=0; j<3; ++j)
   {
	out_<<"| ";
	for (k=0; k<3; ++k)
	   out_<<t_[i][j][k][0]<<' '<<t_[i][j][k][1]<<' '<<t_[i][j][k][2]<<"\t\t";
	out_<<"|\n";
   }
 return out_;
}*/
//------------------------------------------------------------------------------
