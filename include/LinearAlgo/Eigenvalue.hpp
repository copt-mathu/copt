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


#ifndef EIGEN_VALUE_HPP__
#define EIGEN_VALUE_HPP__

namespace COPT
{

/** 		Reduce symmetric matrix to tridiagonal matrix.
  *			As introduced in 'Lapack', in order to solve the eigenproblem
  *			of a symmetric matrix the matrix is first reducd to a tridiagonal
  *			matrix consisting of diagonal and sub-diagonal elements.
  *			And the orthogonal matrix Q is also stored in '__a' and '__tau'.
  *			All routines are derived from 'Lapack'.
  */
template<class Matrix>
class SymmetricToTridiagonal
{
private:
	typedef typename Matrix::Kernel 				kernel;
	typedef typename kernel::scalar 				scalar;
	typedef typename kernel::podscalar 				podscalar;
	typedef typename kernel::index 					index;

	/** uplo */
	char 			__uplo;
	/** the dimension of symmetric matrix */
	index 			__n;
	/** the A */
	scalar 			*__a;
	/** the diagonal elements */
	podscalar 		*__d;
	/** the sub-diagonal elements */
	podscalar		*__e;
	/** tau */
	scalar 			*__tau;
	/** information */
	index 			__info;
public:
	/** constructor and deconstructor */
	//%{
	SymmetricToTridiagonal();
	SymmetricToTridiagonal( const Matrix& mat );
	~SymmetricToTridiagonal();
	//%}

	/** kernel function */
	void compute(const Matrix& mat);

	/** getter */
	//%{
	index n() const;
	podscalar* a();
	podscalar* d();
	podscalar* e();
	scalar* tau();
	//%}
};

/**			EigenSolver solves eigenproblem of a general symmetric
  * 		matrix. 'Lapack' routine is used.
  *
  */
template<class Matrix>
class EigenSolver
{

private:

	typedef typename Matrix::Kernel 			kernel;
	typedef typename Matrix::index 				index;
	typedef typename Matrix::scalar 			scalar;
	typedef typename kernel::podscalar 			podscalar;
	typedef SymmetricToTridiagonal<Matrix>		tri_generator;
	typedef typename Matrix::DVector 			DVector;
	typedef typename Matrix::DMatrix 			DMatrix;

	/** the triagonal */
	tri_generator 		__tri;
	/** the dimension of the matrix */
	index 				__dim;
	/** eigenvalue */
	DVector 			__eigenvalue;
	/** eigenvector */
	DMatrix 			__eigenvector;
	/** information */
	index 				__info;
	/** whether is solved */
	bool 				__is_solved;
	/** isuppz */
	index 				*__isuppz;

public:

	EigenSolver();
	EigenSolver( const Matrix& mat );
	~EigenSolver( );

	/** compute the symmetric matrix */
	void compute(const Matrix& mat);
	/** get the eigen vector */
	const DMatrix& eigenVector() const;
	/** get the eigen value */
	const DVector& eigenValue() const;
};

/*			'ParitialEigenSolver' can partially compute the eigenvalues of a dense
 *			matrix. The algorithms takes a tridiagonal as input which can be achieved
 *			by previous class 'SymmetricToTridiagonal'. If no SymmetricToTridiagonal is
 *			used as input, the class itself computes it at first. The solver can compute 
 *			the eigenvalues in a range of indices like il<=i<=iu as well as a range of 
 *			values like vi<=v<=vu.
 *
 */
template<class Matrix>
class PartialEigenSolver
{
private:
	typedef typename Matrix::Kernel 			kernel;
	typedef typename kernel::index 				index;
	typedef typename kernel::scalar 			scalar;
	typedef typename kernel::podscalar 			podscalar;
	typedef SymmetricToTridiagonal<Matrix>		tri_generator;

	/** the triagonal */
	tri_generator 		__tri;
	/** the dimension of the matrix */
	index 				__dim;
	/** the actual number of eigenvalues fount */
	index 				__m;
	/** storing the eigenvalues */
	podscalar 			*__w;
	/** block infomation */
	index 				*__iblock;
	/**	nsplit */
	index 				__nsplit;
	/** the split array */
	index 				*__isplit;
	/** the information */
	index 				__info;
	/** whether is solved */
	bool 				__is_solved;

public:
	PartialEigenSolver( );
	PartialEigenSolver( const Matrix& mat );
	~PartialEigenSolver( );

	void compute(const Matrix& mat);

	/** solve partial eigenvalue problem in range [il,iu] */
	void solveIndex( const index il, const index iu );
	/** solve partial eigenvalue problem in value range [vl,vu] */
	void solveValueRange( const podscalar vl, const podscalar vu);

	podscalar computeLargestEigenvalue();
};

/********************Implementation of 'SymmetricToTridiagonal'****************/
template<class Matrix>
SymmetricToTridiagonal<Matrix>::SymmetricToTridiagonal()
	:
	__uplo('U'),
	__a(nullptr),
	__d(NULL),
	__e(NULL),
	__tau(NULL)
{
}

template<class Matrix>
SymmetricToTridiagonal<Matrix>::SymmetricToTridiagonal( const Matrix& mat )
	:
	__uplo('U'),
	__a(nullptr),
	__d(NULL),
	__e(NULL),
	__tau(NULL)
{
	compute(mat);
}

template<class Matrix>
void SymmetricToTridiagonal<Matrix>::compute(const Matrix& mat)
{
	if (!mat.isSymmetric())
		std::cerr<<" SymmetricToTridiagonal Warning: Please make sure that whether the input matrix is symmetric or not!"<<std::endl;
	__n = mat.cols();
	SAFE_DELETE_ARRAY(__d);
	SAFE_DELETE_ARRAY(__e);
	SAFE_DELETE_ARRAY(__tau);
	SAFE_DELETE_ARRAY(__a);
	__d = new podscalar[__n];
	__e = new podscalar[__n-1];
	__tau = new scalar[__n-1];
	__a = new scalar[mat.size()];
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	if ( is_real<scalar>::value )
		copt_lapack_sytrd(__uplo,__n,__a,mat.lda(),__d,__e,__tau,&__info);
	else if( is_complex<scalar>::value )
		copt_lapack_hetrd(__uplo,__n,__a,mat.lda(),__d,__e,__tau,&__info);
	else
		throw COException("Unknown scalar type for symmetric matrix calculation!");
}

template<class Matrix>
SymmetricToTridiagonal<Matrix>::~SymmetricToTridiagonal()
{
	SAFE_DELETE_ARRAY(__d);
	SAFE_DELETE_ARRAY(__e);
	SAFE_DELETE_ARRAY(__tau);
	SAFE_DELETE_ARRAY(__a);
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::index SymmetricToTridiagonal<Matrix>::n() const
{
	return __n;
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::podscalar* SymmetricToTridiagonal<Matrix>::a()
{
	return __a;
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::podscalar* SymmetricToTridiagonal<Matrix>::d()
{
	return __d;
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::podscalar* SymmetricToTridiagonal<Matrix>::e()
{
	return __e;
}

template<class Matrix>
typename SymmetricToTridiagonal<Matrix>::scalar* SymmetricToTridiagonal<Matrix>::tau()
{
	return __tau;
}
//////////////End of implementation of 'SymmetricToTridiagonal'

/*********************Implementation of 'EigenSolver'***********************/
template<class Matrix>
EigenSolver<Matrix>::EigenSolver()
	:
	__isuppz(nullptr)
{}

template<class Matrix>
EigenSolver<Matrix>::EigenSolver(const Matrix& mat)
	:
	__isuppz(nullptr)
{
	this->compute(mat);
}

template<class Matrix>
EigenSolver<Matrix>::~EigenSolver()
{
	SAFE_DELETE_ARRAY(__isuppz);
}

template<class Matrix>
void EigenSolver<Matrix>::compute(const Matrix& mat)
{
	if(!mat.isSymmetric())
		std::cerr<<" SymmetricToTridiagonal Warning: Please make sure that whether the input matrix is symmetric or not!"<<std::endl;
	__dim = mat.rows();
	__tri.compute(mat);
	std::cout<<"here"<<std::endl;
	SAFE_DELETE_ARRAY(__isuppz);
	__isuppz = new index[2*__dim];
	__eigenvalue.setArray(__dim,__tri.d());
	__eigenvector.resize(__dim,__dim);
	blas::copt_blas_copy(__eigenvector.size(),__tri.a(),1,__eigenvector.dataPtr(),1);
	std::cout<<"here"<<std::endl;
	copt_lapack_orgtr('U',__dim,__eigenvector.dataPtr(),__eigenvector.lda(),const_cast<scalar*>(__tri.tau()),&__info);
	copt_lapack_steqr('V',__dim,__eigenvalue.dataPtr(),const_cast<scalar*>(__tri.e()),__eigenvector.dataPtr(),__eigenvector.lda(),&__info);
	// copt_lapack_stegr('V','A',__dim,const_cast<scalar*>(__tri.d()),const_cast<scalar*>(__tri.e()),0.0,0.0,0,0,0.0,__dim,__eigenvalue.dataPtr(),__eigenvector.dataPtr(),__eigenvector.lda(),__isuppz,&__info);
}

template<class Matrix>
const typename EigenSolver<Matrix>::DVector& EigenSolver<Matrix>::eigenValue() const
{
	return __eigenvalue;
}

template<class Matrix>
const typename EigenSolver<Matrix>::DMatrix& EigenSolver<Matrix>::eigenVector() const
{
	return __eigenvector;
}

/*********************Implementation of 'PartialEigenSolver'****************/ 
template<class Matrix>
PartialEigenSolver<Matrix>::PartialEigenSolver()
	:
	__w(NULL),
	__iblock(NULL),
	__isplit(NULL),
	__is_solved(false)
{
}

template<class Matrix>
PartialEigenSolver<Matrix>::PartialEigenSolver(const Matrix& mat)
	:
	__w(NULL),
	__iblock(NULL),
	__isplit(NULL),
	__is_solved(false)
{
	compute(mat);
}

template<class Matrix>
PartialEigenSolver<Matrix>::~PartialEigenSolver()
{
	SAFE_DELETE_ARRAY(__w);
	SAFE_DELETE_ARRAY(__iblock);
	SAFE_DELETE_ARRAY(__isplit);
}

template<class Matrix>
void PartialEigenSolver<Matrix>::compute( const Matrix& mat )
{
	if (!mat.isSymmetric())
		std::cerr<<" SymmetricToTridiagonal Warning: Please make sure that whether the input matrix is symmetric or not!"<<std::endl;
	__tri.compute(mat);
	__dim = mat.cols();
	SAFE_DELETE_ARRAY(__w);
	SAFE_DELETE_ARRAY(__iblock);
	SAFE_DELETE_ARRAY(__isplit);
	__w = new podscalar[__dim];
	__iblock = new index[__dim];
	__isplit = new index[__dim];
	__is_solved = false;
}

template<class Matrix>
void PartialEigenSolver<Matrix>::solveIndex( index il , index iu)
{
	podscalar *work = new podscalar[4*__dim];
	index *iwork = new index[3*__dim];
	copt_lapack_stebz('I','B',__dim,0.0,0.0,il,iu,0.0,const_cast<podscalar*>(__tri.d()),const_cast<podscalar*>(__tri.e()) ,&__m,&__nsplit,__w,__iblock,__isplit,work,iwork,&__info);
	delete[]work;
	delete[]iwork;
	__is_solved = true;
}

template<class Matrix>
typename PartialEigenSolver<Matrix>::podscalar PartialEigenSolver<Matrix>::computeLargestEigenvalue()
{
	solveIndex(__dim,__dim);
	return __w[0];
}


}
#endif