
#include "Tensors.hpp"
#include "ThrowMessage.hpp"

using namespace std;

//==========================Tensor1s============================================

Tensor2s& Tensor1s::DirectProduct(const Tensor1s& o_, Tensor2s& result_) const
{
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = m_data[i] * o_[j];
	return result_;
}

Tensor1s& Tensor1s::DotMultiply(const Tensor2s& o_)
{
	real_t data_0 = m_data[0], data_1 = m_data[1];
	for (small_t i = 0; i < 3; ++i)
		m_data[i] = data_0 * o_[0][i] + data_1 * o_[1][i] + m_data[2] * o_[2][i];
	return *this;
}

Tensor1s& Tensor1s::DotMultiply_left(const Tensor2s& o_)
{
	real_t data_0 = m_data[0], data_1 = m_data[1];
	for (small_t i = 0; i < 3; ++i)
		m_data[i] = o_[i][0] * data_0 + o_[i][1] * data_1 + o_[i][2] * m_data[2];
	return *this;
}

Tensor3s& Tensor1s::DirectProduct(const Tensor2s& o_, Tensor3s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = m_data[i] * o_[j][k];
	return result_;
}

Tensor2s& Tensor1s::DotProduct(const Tensor3s& o_, Tensor2s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k) result_[i][j] += m_data[k] * o_[k][i][j];
		}
	return result_;
}

Tensor4s& Tensor1s::DirectProduct(const Tensor3s& o_, Tensor4s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = m_data[i] * o_[j][k][l];
	return result_;
}

Tensor3s& Tensor1s::DotProduct(const Tensor4s& o_, Tensor3s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l) result_[i][j][k] += m_data[l] * o_[l][i][j][k];
			}
	return result_;
}

//==========================Tensor1a=============================================

Tensor1a& Tensor1a::operator-=(const Tensor1a& o_)
{
#ifdef STRONGCHECK
	Assert(this != &o_, "the same parameter in Tensor1a::operator-=");
#endif
	// if (&o_ == this) return Assign0();
	if (__AdjustAdd(o_))
		m_pData->operator-=(*o_.m_pData);
	else if (!is0())
		m_pData->Negate();
	return *this;
}

//==========================Tensor2s=======================================

Tensor2s::Tensor2s(const Tensor1s& col1_, const Tensor1s& col2_, const Tensor1s& col3_) : m_data()
{
	__SET_ROWS;
	for (small_t i = 0; i < 3; ++i)
	{
		m_rows[i][0] = col1_.m_data[i];
		m_rows[i][1] = col2_.m_data[i];
		m_rows[i][2] = col3_.m_data[i];
	}
}

Tensor2s::Tensor2s(const Tensor1s& col1_, const Tensor1s& col2_, const Tensor1s& col3_, bool)
	: m_data()
{
	__SET_ROWS;
	for (small_t i = 0; i < 3; ++i)
	{
		m_rows[0][i] = col1_.m_data[i];
		m_rows[1][i] = col2_.m_data[i];
		m_rows[2][i] = col3_.m_data[i];
	}
}

Tensor2s& Tensor2s::operator-=(const Tensor2s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 9; ++pd, ++po) *pd -= *po;
	return *this;
}

Tensor2s& Tensor2s::operator+=(const Tensor2s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 9; ++pd, ++po) *pd += *po;
	return *this;
}

Tensor2s& Tensor2s::AssignRows(const Tensor1s& row1_, const Tensor1s& row2_, const Tensor1s& row3_)
{
	real_t *p0 = *m_rows, *p1 = m_rows[1], *p2 = m_rows[2];
	const real_t *psrc0 = row1_.m_data, *psrc1 = row2_.m_data, *psrc2 = row3_.m_data;
	for (; p0 < m_rows[1]; ++p0, ++p1, ++p2, ++psrc0, ++psrc1, ++psrc2)
	{
		*p0 = *psrc0;
		*p1 = *psrc1;
		*p2 = *psrc2;
	}
	return *this;
}

Tensor2s& Tensor2s::AssignCols(const Tensor1s& col1_, const Tensor1s& col2_, const Tensor1s& col3_)
{
	real_t *p0 = *m_rows, *p1 = m_rows[1], *p2 = m_rows[2];
	const real_t *psrc0 = col1_.m_data, *psrc1 = col2_.m_data, *psrc2 = col3_.m_data;
	for (; p0 <= m_rows[2]; p0 += 3, p1 += 3, p2 += 3, ++psrc0, ++psrc1, ++psrc2)
	{
		*p0 = *psrc0;
		*p1 = *psrc1;
		*p2 = *psrc2;
	}
	return *this;
}

real_t Tensor2s::Minor(small_t i_, small_t j_) const
{
#ifdef STRONGCHECK
	Assert(i_ > 0 && j_ > 0 && i_ <= 3 && j_ <= 3, "invalid index in Tensor2s::Minor");
#endif
	small_t i1 = 4 - i_, j1 = i1 == 1 ? 1 : 2, i2 = 4 - j_, j2 = i2 == 1 ? 1 : 2;
	i1 -= j1;
	i2 -= j2;
	return m_rows[i1][i2] * m_rows[j1][j2] - m_rows[i1][j2] * m_rows[j1][i2];
}

Tensor2s& Tensor2s::Invert()
{
	real_t det1 = Det();
#ifdef STRONGCHECK
	Assert(NE0(det1), tMessage("singular tensor in Tensor2s::Invert\ndet = ") << det1);
#endif
	det1 = 1. / det1;
	real_t tmp[9] = {(m_data[4] * m_data[8] - m_data[5] * m_data[7]) * det1,
					 (m_data[2] * m_data[7] - m_data[1] * m_data[8]) * det1,
					 (m_data[1] * m_data[5] - m_data[2] * m_data[4]) * det1,
					 (m_data[5] * m_data[6] - m_data[3] * m_data[8]) * det1,
					 (m_data[0] * m_data[8] - m_data[2] * m_data[6]) * det1,
					 (m_data[2] * m_data[3] - m_data[0] * m_data[5]) * det1,
					 (m_data[3] * m_data[7] - m_data[4] * m_data[6]) * det1,
					 (m_data[1] * m_data[6] - m_data[0] * m_data[7]) * det1,
					 (m_data[0] * m_data[4] - m_data[1] * m_data[3]) * det1};
	copy(tmp, tmp + 9, m_data);
	return *this;
}

void Tensor2s::EigenValues(real_t eigVals_[3]) const
{
	// Characteristic equation eigVal^3 - I1*eigVal^2 + I2*eigVal - I3 =0
	// reduced to x^3 + p*x + q =0  by substitution {x = eigVal-I1/3}:
	const real_t i1 = I1(), i2 = I2(), i3 = I3(), i1_3 = i1 / 3., p = i2 - i1 * i1_3,
				 q = i1_3 * (i2 - 2. * i1_3 * i1_3) - i3;  // coefs of reduced eqn
	Assert(mathdef::is_lte(p, 0.), "Eigen values with Im!=0 in Tensor2s::EigenValues(real_t[])");
	const real_t ro = sqrt(-p / 3.);  // auxiliary param
	if (mathdef::is_not_zero(ro))
	{
		const real_t fi = acos(-q / (2. * ro * ro * ro)) / 3.,	// auxiliary param
			pi23 = acos(-0.5);									//=2Pi/3
		eigVals_[0] = 2. * ro * cos(fi) + i1_3;
		eigVals_[1] = 2. * ro * cos(fi + pi23) + i1_3;
		eigVals_[2] = 2. * ro * cos(fi + pi23 * 2.) + i1_3;
	}
	else
		eigVals_[0] = eigVals_[1] = eigVals_[2] = (sign(pow(fabs(q), real_t(1. / 3.)), -q) + i1_3);
}

struct less_in_abs
{
	bool operator()(const real_t& a, const real_t& b) const
	{
		return fabs(a) < fabs(b);
	}
};

void Tensor2s::MaxAbsIndex(small_t& row_num_, small_t& col_num_) const
{
	const small_t lin_num =
		static_cast<small_t>(max_element(m_data, m_data + 9, less_in_abs()) - m_data);
	row_num_ = lin_num / 3;
	col_num_ = lin_num % 3;
}

small_t Tensor2s::EigenVectors(real_t eigVal_, Tensor2s& eigVects_, small_t begPos_) const
{
#ifdef STRONGCHECK
	Assert(
		begPos_ > 0 && begPos_ <= 3, "invalid position of eigenvector in Tensor2s::EigenVectors");
#endif
	Tensor2s singular(*this);
	real_t tmp;
	small_t i = 0, j, i1, i2, j1, j2, multiplicity = 1;
	const small_t posOfCurrent = begPos_ - 1;
	for (; i < 3; ++i) singular[i][i] -= eigVal_;  // i.e. singular = *this - eigval*E
#ifdef STRONGCHECK
	Assert(
		EQ0(singular.Det(), PRECISION * 50.),
		tMessage(
			"value expected to be eigen is not eigen in Tensor2s::EigenVectors\ndet(A-lambdaE) = ")
			<< singular.Det() << " (abs must be <= " << PRECISION * 20. << ')');
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
	Assert(begPos_ <= 2, "invalid position of eigenvector in Tensor2s::EigenVectors");
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
		Assert(begPos_ == 1, "invalid position of eigenvector in Tensor2s::EigenVectors");
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

void Tensor2s::Eigen(real_t eigvals_[3], Tensor2s& eigvects_) const
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
		"eigen vectors are not correspond to eigen values in Tensor2s::Eigen");
#endif
	if (eig_vects_count < 3)
		EigenVectors(eigvals_[2], eigvects_, 3);
}

Tensor2s& Tensor2s::RotateToEigenSystem(Tensor2s& rotation_)
{  // assigns columns of arg to eigenvectors and assigns *this to diagonal form
	real_t eigenvals[3];
	Eigen(eigenvals, rotation_);
	Assign0();
	for (small_t i = 0; i < 3; ++i) m_data[3 * i] = eigenvals[i];
	return *this;
}

Tensor2s& Tensor2s::DotMultiply(const Tensor2s& o_)	 //(T*F)[ij]=T[ip]*F[pj]
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "mult by itself in Tensor2s::DotMultiply(const Tensor2s&)");
#endif
	// if (&o_ == this)
	//       return DotMultiply(Tensor2s(*this));
	real_t rows_i0, rows_i1;
	for (small_t i = 0; i < 3; ++i)
	{
		rows_i0 = m_rows[i][0];
		rows_i1 = m_rows[i][1];
		m_rows[i][0] = rows_i0 * o_[0][0] + rows_i1 * o_[1][0] + m_rows[i][2] * o_[2][0];
		m_rows[i][1] = rows_i0 * o_[0][1] + rows_i1 * o_[1][1] + m_rows[i][2] * o_[2][1];
		m_rows[i][2] = rows_i0 * o_[0][2] + rows_i1 * o_[1][2] + m_rows[i][2] * o_[2][2];
	}
	return *this;
}

Tensor2s& Tensor2s::DotMultiply_left(const Tensor2s& o_)
{
#ifdef STRONGCHECK
	Assert(&o_ != this, "mult by itself in Tensor2s::DotMultiply_left(const Tensor2s&)");
#endif
	// if (&o_ == this)
	//        return DotMultiply_left(Tensor2s(Data));
	real_t rows_0j, rows_1j;
	for (small_t j = 0; j < 3; ++j)
	{
		rows_0j = m_rows[0][j];
		rows_1j = m_rows[1][j];
		m_rows[0][j] = o_[0][0] * rows_0j + o_[0][1] * rows_1j + o_[0][2] * m_rows[2][j];
		m_rows[1][j] = o_[1][0] * rows_0j + o_[1][1] * rows_1j + o_[1][2] * m_rows[2][j];
		m_rows[2][j] = o_[2][0] * rows_0j + o_[2][1] * rows_1j + o_[2][2] * m_rows[2][j];
	}
	return *this;
}

Tensor2s& Tensor2s::CrossMultiply(const Tensor1s& o_)  //(T*f)[ij]=T[ip]*f[q]*LevyChivita[pqj]
{
	real_t rows_i0, rows_i1;
	for (small_t i = 0; i < 3; ++i)
	{
		rows_i0 = m_rows[i][0];
		rows_i1 = m_rows[i][1];
		m_rows[i][0] = rows_i1 * o_[2] - m_rows[i][2] * o_[1];
		m_rows[i][1] = m_rows[i][2] * o_[0] - rows_i0 * o_[2];
		m_rows[i][2] = rows_i0 * o_[1] - rows_i1 * o_[0];
	}
	return *this;
}

Tensor2s& Tensor2s::CrossMultiply_left(const Tensor1s& o_)	//(t*F)[ij]=t[p]*F[qj]*LevyChivita[pqi]
{
	real_t rows_0j, rows_1j;
	for (small_t j = 0; j < 3; ++j)
	{
		rows_0j = m_rows[0][j];
		rows_1j = m_rows[1][j];
		m_rows[0][j] = o_[1] * m_rows[2][j] - o_[2] * rows_1j;
		m_rows[1][j] = o_[2] * rows_0j - o_[0] * m_rows[2][j];
		m_rows[2][j] = o_[0] * rows_1j - o_[1] * rows_0j;
	}
	return *this;
}

Tensor3s& Tensor2s::DirectProduct(const Tensor1s& o_, Tensor3s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = m_rows[i][j] * o_[k];
	return result_;
}

real_t Tensor2s::Dot2Product(const Tensor2s& o_) const
{
	real_t result = 0.;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result += m_rows[i][j] * o_[j][i];
	return result;
}

Tensor3s& Tensor2s::CrossProduct(const Tensor2s& o_, Tensor3s& result_) const
{
	const Tensor3s levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						result_[i][j][k] += m_rows[i][l] * o_[m][k] * levyChivita[l][m][j];
			}
	return result_;
}

Tensor4s& Tensor2s::DirectProduct(const Tensor2s& o_, Tensor4s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = m_rows[i][j] * o_[k][l];
	return result_;
}

Tensor1s& Tensor2s::Dot2Product(const Tensor3s& o_, Tensor1s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i] += m_rows[j][k] * o_[k][j][i];
	}
	return result_;
}

Tensor4s& Tensor2s::CrossProduct(const Tensor3s& o_, Tensor4s& result_) const
{
	const Tensor3s levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					result_[i][j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							result_[i][j][k][l] +=
								m_rows[i][m] * o_[n][k][l] * levyChivita[m][n][j];
				}
	return result_;
}

// Tensor2s& Tensor2s::Dot2Multiply (const Tensor4& o_, small_t index1No_, small_t index2No_)
Tensor2s& Tensor2s::Dot2Multiply(const SymmetricTensor4s& o_, small_t index1No_, small_t index2No_)
{
#ifdef STRONGCHECK
	Assert(
		index1No_ != index2No_ && index1No_ > 0 && index1No_ <= 4 && index2No_ > 0 &&
			index2No_ <= 4,
		"invalid index # in Tensor2s::Dot2Multiply(SymmetricTensor4s&,small_t,small_t)");
#endif
	Tensor4::tIndex ijkl;
	const Tensor2s copyOfThis(*this);
	const small_t freeIndex1No =
					  index1No_ + index2No_ == 3 ? 3 : (index1No_ == 1 || index2No_ == 1 ? 2 : 1),
				  freeIndex2No = 10 - freeIndex1No - index1No_ - index2No_;
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		ijkl[freeIndex1No] = i;
		for (j = 0; j < 3; ++j)
		{
			ijkl[freeIndex2No] = j;
			m_rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
			{
				ijkl[index1No_] = k;
				for (l = 0; l < 3; ++l)
				{
					ijkl[index2No_] = l;
					m_rows[i][j] += copyOfThis[k][l] * o_[ijkl];
				}
			}
		}
	}
	return *this;
}

Tensor2s& Tensor2s::Dot2Multiply(const Tensor4s& o_)
{
	const Tensor2s copyOfThis(*this);
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			m_rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) m_rows[i][j] += copyOfThis[k][l] * o_[l][k][i][j];
		}
	return *this;
}

Tensor2s& Tensor2s::Dot2Multiply_left(const Tensor4s& o_)
{
	const Tensor2s copyOfThis(*this);
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			m_rows[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) m_rows[i][j] += o_[i][j][k][l] * copyOfThis[l][k];
		}
	return *this;
}

/*Tensor2& Tensor2s::ApplyFunction (const Tensor2s::fClassicFunction& function_)
{//according to (Lurie, 1980), page 436.
 Tensor1s eigVals, funOfEigVals;   EigenValues(eigVals);
 try{
	 funOfEigVals[0] = function_(eigVals[0]);
	 funOfEigVals[1] = function_(eigVals[1]);
	 funOfEigVals[2] = function_(eigVals[2]);
	}
 catch (...) {throw tMessage("Invalid arg in Tensor2s::ApplyFunction (maybe divide by 0 etc.)");}
 Tensor2s currentEigDiade, result(0.);
 small_t i=0,j,k,l;   real_t tmp;//   Index index;
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

//==========================Tensor2a=============================================

/*ostream& operator<< (ostream& out_, const Tensor2a& o_)
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

Tensor2a& Tensor2a::operator-=(const Tensor2a& o_)
{
#ifdef STRONGCHECK
	Assert(this != &o_, "the same parameter in Tensor2a::operator-=");
#endif
	// if (this == &o_) return Assign0();
	if (__AdjustAdd(o_))
		m_pData->operator-=(*o_.m_pData);
	else if (!is0())
		m_pData->Negate();
	return *this;
}

//==========================Tensor3s==========================================

Tensor3s& Tensor3s::operator-=(const Tensor3s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 27; ++pd, ++po) *pd -= *po;
	return *this;
}

Tensor3s& Tensor3s::operator+=(const Tensor3s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 27; ++pd, ++po) *pd += *po;
	return *this;
}

Tensor3s& Tensor3s::Assign1(real_t v_)	// Levy-Chivita multiplied by v_
{
	Assign0();
	(*m_rows[0])[1][2] = (*m_rows[1])[2][0] = (*m_rows[2])[0][1] = v_;
	(*m_rows[0])[2][1] = (*m_rows[2])[1][0] = (*m_rows[1])[0][2] = -v_;
	return *this;
}

Tensor1s& Tensor3s::Contraction(small_t n1_, small_t n2_, Tensor1s& result_) const
{
#ifdef STRONGCHECK
	Assert(
		n1_ > 0 && n2_ > 0 && n1_ <= 3 && n2_ <= 3 && n1_ != n2_,
		"invalid number of index in Tensor3s::Contraction");
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

Tensor2s& Tensor3s::ScalarProductWithBasisVector(small_t k_, Tensor2s& result_) const
{
#ifdef STRONGCHECK
	Assert(k_ > 0 && k_ <= 3, "invalid index in Tensor3s::ScalarProductWithBasisVector");
#endif
	--k_;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = (*m_rows[i])[j][k_];
	return result_;
}

Tensor2s& Tensor3s::ScalarProductWithBasisVector_left(small_t k_, Tensor2s& result_) const
{
#ifdef STRONGCHECK
	Assert(k_ > 0 && k_ <= 3, "invalid index in Tensor3s::ScalarProductWithBasisVector_left");
#endif
	--k_;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j) result_[i][j] = (*m_rows[k_])[i][j];
	return result_;
}

Tensor2s& Tensor3s::DotProduct(const Tensor1s& o_, Tensor2s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k) result_[i][j] += (*m_rows[i])[j][k] * o_[k];
		}
	return result_;
}

Tensor3s& Tensor3s::CrossMultiply(const Tensor1s& o_)
{
	const Tensor3s copyOfThis(*this), levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*m_rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						(*m_rows[i])[j][k] += copyOfThis[i][j][l] * o_[m] * levyChivita[l][m][k];
			}
	return *this;
}

Tensor3s& Tensor3s::CrossMultiply_left(const Tensor1s& o_)
{
	const Tensor3s copyOfThis(*this), levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*m_rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m)
						(*m_rows[i])[j][k] += o_[l] * copyOfThis[m][j][k] * levyChivita[l][m][i];
			}
	return *this;
}

Tensor4s& Tensor3s::DirectProduct(const Tensor1s& o_, Tensor4s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j][k][l] = (*m_rows[i])[j][k] * o_[l];
	return result_;
}

Tensor1s& Tensor3s::Dot2Product(const Tensor2s& o_, Tensor1s& result_) const
{
	for (small_t i = 0, j, k; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i] += (*m_rows[i])[j][k] * o_[k][j];
	}
	return result_;
}

Tensor3s& Tensor3s::DotMultiply(const Tensor2s& o_)
{
	real_t data_ij0, data_ij1;
	for (small_t i = 0, j; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			data_ij0 = (*m_rows[i])[j][0];
			data_ij1 = (*m_rows[i])[j][1];
			(*m_rows[i])[j][0] =
				data_ij0 * o_[0][0] + data_ij1 * o_[1][0] + (*m_rows[i])[j][2] * o_[2][0];
			(*m_rows[i])[j][1] =
				data_ij0 * o_[0][1] + data_ij1 * o_[1][1] + (*m_rows[i])[j][2] * o_[2][1];
			(*m_rows[i])[j][2] =
				data_ij0 * o_[0][2] + data_ij1 * o_[1][2] + (*m_rows[i])[j][2] * o_[2][2];
		}
	return *this;
}

Tensor3s& Tensor3s::DotMultiply_left(const Tensor2s& o_)
{
	real_t data_0jk, data_1jk;
	for (small_t j = 0, k; j < 3; ++j)
		for (k = 0; k < 3; ++k)
		{
			data_0jk = (*m_rows[0])[j][k];
			data_1jk = (*m_rows[1])[j][k];
			(*m_rows[0])[j][k] =
				o_[0][0] * data_0jk + o_[0][1] * data_1jk + o_[0][2] * (*m_rows[2])[j][k];
			(*m_rows[1])[j][k] =
				o_[1][0] * data_0jk + o_[1][1] * data_1jk + o_[1][2] * (*m_rows[2])[j][k];
			(*m_rows[2])[j][k] =
				o_[2][0] * data_0jk + o_[2][1] * data_1jk + o_[2][2] * (*m_rows[2])[j][k];
		}
	return *this;
}

real_t Tensor3s::Dot3Product(const Tensor3s& o_) const
{
	real_t result = 0.;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result += (*m_rows[i])[j][k] * o_[k][j][i];
	return result;
}

Tensor2s& Tensor3s::Dot2Product(const Tensor3s& o_, Tensor2s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i][j] += (*m_rows[i])[k][l] * o_[l][k][j];
		}
	return result_;
}

Tensor4s& Tensor3s::DotProduct(const Tensor3s& o_, Tensor4s& result_) const
{
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					result_[i][j][k][l] = 0.;
					for (m = 0; m < 3; ++m) result_[i][j][k][l] += (*m_rows[i])[j][m] * o_[m][k][l];
				}
	return result_;
}

Tensor1s& Tensor3s::Dot3Product(const Tensor4s& o_, Tensor1s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i] += (*m_rows[j])[k][l] * o_[l][k][j][i];
	}
	return result_;
}

Tensor3s& Tensor3s::Dot2Multiply(const Tensor4s& o_)
{
	Tensor2s tmp;
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
	{
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) tmp[j][k] = (*m_rows[i])[j][k];
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*m_rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) (*m_rows[i])[j][k] += tmp[l][m] * o_[m][l][j][k];
			}
	}
	return *this;
}

Tensor3s& Tensor3s::Dot2Multiply_left(const Tensor4s& o_)
{
	Tensor2s tmp;
	for (small_t i, j, k = 0, l, m; k < 3; ++k)
	{
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j) tmp[i][j] = (*m_rows[i])[j][k];
		for (i = 0; i < 3; ++i)
			for (j = 0; j < 3; ++j)
			{
				(*m_rows[i])[j][k] = 0.;
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) (*m_rows[i])[j][k] += o_[i][j][l][m] * tmp[m][l];
			}
	}
	return *this;
}

//==========================Tensor4s==========================================

/*Tensor4s& Tensor4s::operator= (const Tensor4& o_)
{
 for (small_t i=0,j,k,l; i<3; ++i)
  for   (j=0; j<3; ++j)
   for  (k=0; k<3; ++k)
	for (l=0; l<3; ++l)
	 (*Rows[i])[j][k][l] = o_[tIndex(i,j,k,l)];
//     (*Rows[i])[j][k][l] = o_(i+1,j+1,k+1,l+1);
 return *this;
} */

Tensor4s& Tensor4s::operator=(const SymmetricTensor4s& o_)
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		for (k = 0; k < 3; ++k)
		{
			(*m_rows[i])[i][k][k] = o_[tIndex(i, i, k, k)];
			for (l = 0; l < k; ++l)
				(*m_rows[i])[i][k][l] = (*m_rows[i])[i][l][k] = o_[tIndex(i, i, k, l)];
		}
		for (j = 0; j < i; ++j)
			for (k = 0; k < 3; ++k)
			{
				(*m_rows[i])[j][k][k] = (*m_rows[j])[i][k][k] = o_[tIndex(i, j, k, k)];
				for (l = 0; l < k; ++l)
					(*m_rows[i])[j][k][l] = (*m_rows[i])[j][l][k] = (*m_rows[j])[i][k][l] =
						(*m_rows[j])[i][l][k] = o_[tIndex(i, j, k, l)];
			}
	}
	return *this;
}

Tensor4s& Tensor4s::operator-=(const Tensor4s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 81; ++pd, ++po) *pd -= *po;
	return *this;
}

Tensor4s& Tensor4s::operator+=(const Tensor4s& o_)
{
	const real_t* po = o_.m_data;
	for (real_t* pd = m_data; pd < m_data + 81; ++pd, ++po) *pd += *po;
	return *this;
}

Tensor4s& Tensor4s::Assign1(small_t kind_, real_t v_)  // Isotropic multiplied by v_
{
#ifdef STRONGCHECK
	Assert(kind_ > 0 && kind_ <= 3, "invalid of isotropic Tensor4s in Tensor4s::Assign1");
#endif
	Assign0();
	small_t i = 0, j = 0;
	const small_t &k = (kind_ == 1) ? i : j, &l = (kind_ == 2) ? i : j, &m = (kind_ == 3) ? i : j;
	for (; i < 3; ++i)
		for (j = 0; j < 3; ++j) (*m_rows[i])[k][l][m] = v_;
	return *this;
}

Tensor2s& Tensor4s::Contraction(small_t n1_, small_t n2_, Tensor2s& result_) const
{
#ifdef STRONGCHECK
	Assert(
		n1_ > 0 && n2_ > 0 && n1_ <= 4 && n2_ <= 4 && n1_ != n2_,
		"invalid number of index in Tensor4s::Contraction");
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

Tensor3s& Tensor4s::ScalarProductWithBasisVector(small_t l_, Tensor3s& result_) const
{
#ifdef STRONGCHECK
	Assert(l_ > 0 && l_ <= 3, "invalid index in Tensor4s::ScalarProductWithBasisVector");
#endif
	--l_;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = (*m_rows[i])[j][k][l_];
	return result_;
}

Tensor3s& Tensor4s::ScalarProductWithBasisVector_left(small_t l_, Tensor3s& result_) const
{
#ifdef STRONGCHECK
	Assert(l_ > 0 && l_ <= 3, "invalid index in Tensor4s::ScalarProductWithBasisVector_left");
#endif
	--l_;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k) result_[i][j][k] = (*m_rows[l_])[i][j][k];
	return result_;
}

Tensor3s& Tensor4s::DotProduct(const Tensor1s& o_, Tensor3s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				result_[i][j][k] = 0.;
				for (l = 0; l < 3; ++l) result_[i][j][k] += (*m_rows[i])[j][k][l] * o_[l];
			}
	return result_;
}

Tensor4s& Tensor4s::CrossMultiply(const Tensor1s& o_)
{
	const Tensor4s copyOfThis(*this);
	const Tensor3s levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*m_rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							(*m_rows[i])[j][k][l] +=
								copyOfThis[i][j][k][m] * o_[n] * levyChivita[m][n][l];
				}
	return *this;
}

Tensor4s& Tensor4s::CrossMultiply_left(const Tensor1s& o_)
{
	const Tensor4s copyOfThis(*this);
	const Tensor3s levyChivita(eUnit);
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*m_rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n)
							(*m_rows[i])[j][k][l] +=
								o_[m] * copyOfThis[n][j][k][l] * levyChivita[m][n][i];
				}
	return *this;
}

Tensor4s& Tensor4s::DotMultiply(const Tensor2s& o_, small_t ownIndexNo_, small_t facIndexNo_)
{
#ifdef STRONGCHECK
	Assert(
		ownIndexNo_ > 0 && ownIndexNo_ <= 4 && facIndexNo_ > 0 && facIndexNo_ <= 2,
		"invalid index # in Tensor4s::DotMultiply(Tensor4s&,small_t,small_t)");
#endif
	real_t data_ijk0, data_ijk1;
	Tensor2s::tIndex index2;
	tIndex index4;
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

Tensor4s& Tensor4s::DotMultiply(const Tensor2s& o_)
{
	real_t data_ijk0, data_ijk1;
	for (small_t i = 0, j, k; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
			{
				data_ijk0 = (*m_rows[i])[j][k][0];
				data_ijk1 = (*m_rows[i])[j][k][1];
				(*m_rows[i])[j][k][0] =
					data_ijk0 * o_[0][0] + data_ijk1 * o_[1][0] + (*m_rows[i])[j][k][2] * o_[2][0];
				(*m_rows[i])[j][k][1] =
					data_ijk0 * o_[0][1] + data_ijk1 * o_[1][1] + (*m_rows[i])[j][k][2] * o_[2][1];
				(*m_rows[i])[j][k][2] =
					data_ijk0 * o_[0][2] + data_ijk1 * o_[1][2] + (*m_rows[i])[j][k][2] * o_[2][2];
			}
	return *this;
}

Tensor4s& Tensor4s::DotMultiply_left(const Tensor2s& o_)
{
	real_t data_0jkl, data_1jkl;
	for (small_t j = 0, k, l; j < 3; ++j)
		for (k = 0; k < 3; ++k)
			for (l = 0; l < 3; ++l)
			{
				data_0jkl = (*m_rows[0])[j][k][l];
				data_1jkl = (*m_rows[1])[j][k][l];
				(*m_rows[0])[j][k][l] =
					o_[0][0] * data_0jkl + o_[0][1] * data_1jkl + o_[0][2] * (*m_rows[2])[j][k][l];
				(*m_rows[1])[j][k][l] =
					o_[1][0] * data_0jkl + o_[1][1] * data_1jkl + o_[1][2] * (*m_rows[2])[j][k][l];
				(*m_rows[2])[j][k][l] =
					o_[2][0] * data_0jkl + o_[2][1] * data_1jkl + o_[2][2] * (*m_rows[2])[j][k][l];
			}
	return *this;
}

Tensor1s& Tensor4s::Dot3Product(const Tensor3s& o_, Tensor1s& result_) const
{
	for (small_t i = 0, j, k, l; i < 3; ++i)
	{
		result_[i] = 0.;
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result_[i] += (*m_rows[i])[j][k][l] * o_[l][k][j];
	}
	return result_;
}

real_t Tensor4s::Dot4Product(const Tensor4s& o_) const
{
	real_t result = 0.;
	for (small_t i = 0, j, k, l; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) result += (*m_rows[i])[j][k][l] * o_[l][k][j][i];
	return result;
}

Tensor2s& Tensor4s::Dot3Product(const Tensor4s& o_, Tensor2s& result_) const
{
	for (small_t i = 0, j, k, l, m; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			result_[i][j] = 0.;
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
					for (m = 0; m < 3; ++m) result_[i][j] += (*m_rows[i])[k][l][m] * o_[m][l][k][j];
		}
	return result_;
}

Tensor4s& Tensor4s::Dot2Multiply(const Tensor4s& o_)
{
	Tensor2s tmp;
	for (small_t i = 0, j, k, l, m, n; i < 3; ++i)
		for (j = 0; j < 3; ++j)
		{
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l) tmp[k][l] = (*m_rows[i])[j][k][l];
			for (k = 0; k < 3; ++k)
				for (l = 0; l < 3; ++l)
				{
					(*m_rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n) (*m_rows[i])[j][k][l] += tmp[m][n] * o_[n][m][k][l];
				}
		}
	return *this;
}

Tensor4s& Tensor4s::Dot2Multiply_left(const Tensor4s& o_)
{
	Tensor2s tmp;
	for (small_t i, j, k = 0, l, m, n; k < 3; ++k)
		for (l = 0; l < 3; ++l)
		{
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j) tmp[i][j] = (*m_rows[i])[j][k][l];
			for (i = 0; i < 3; ++i)
				for (j = 0; j < 3; ++j)
				{
					(*m_rows[i])[j][k][l] = 0.;
					for (m = 0; m < 3; ++m)
						for (n = 0; n < 3; ++n) (*m_rows[i])[j][k][l] += o_[i][j][m][n] * tmp[n][m];
				}
		}
	return *this;
}

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
