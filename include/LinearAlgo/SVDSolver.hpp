// This file is part of COPT, a light-weight C++ based optimization open source library
//
// Copyright (C) 2015 Ruimin Wang <ruimin.wang13@gmail.com>
// Copyright (C) 2015 MathU
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.


#ifndef SVD_SOLVER_HPP__
#define SVD_SOLVER_HPP__


namespace COPT
{

template<class Matrix>
class SVD
{

private:

	typedef typename Matrix::scalar 			scalar;
	typedef typename Matrix::DMatrix 			DMatrix;
	typedef typename Matrix::DVector 			DVector;

	DMatrix __u;
	DMatrix __vt;
	DMatrix __s;

	scalar *__a;


	int __info;
public:
	SVD();
	~SVD();
	SVD(const Matrix& mat);
	void compute(const Matrix& mat);
	const DMatrix& U() const;
	const DMatrix& VT() const;
	const DMatrix& S() const;
};

template<class Matrix>
SVD<Matrix>::SVD()
	:
	__a(nullptr)
{
}

template<class Matrix>
SVD<Matrix>::SVD(const Matrix& mat)
	:
	__a(NULL)
{
	this->compute(mat);
}

template<class Matrix>
SVD<Matrix>::~SVD()
{
	SAFE_DELETE_ARRAY(__a);
}

template<class Matrix>
void SVD<Matrix>::compute(const Matrix& mat)
{
	SAFE_DELETE_ARRAY(__a);
	__a = new scalar[mat.size()];
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),mat.interval(),__a,1);
	__u.resize(mat.rows(),mat.rows());
	__vt.resize(mat.cols(),mat.cols());
	DVector s(std::min(mat.rows(),mat.cols()));
	copt_lapack_gesvd('A','A',mat.rows(),mat.cols(),__a,mat.lda(),s.dataPtr(),__u.dataPtr(),__u.lda(),__vt.dataPtr(),__vt.lda(),&__info);
	__s = DMatrix::diag(mat.rows(),mat.cols(),s);
}

template<class Matrix>
const typename SVD<Matrix>::DMatrix& SVD<Matrix>::U() const
{
	return __u;
}

template<class Matrix>
const typename SVD<Matrix>::DMatrix& SVD<Matrix>::VT() const
{
	return __vt;
}

template<class Matrix>
const typename SVD<Matrix>::DMatrix& SVD<Matrix>::S() const
{
	return __s;
}

}// End of namespace COPT

#endif