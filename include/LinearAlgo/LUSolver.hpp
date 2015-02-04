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


#ifndef LU_SOLVER_HPP__
#define LU_SOLVER_HPP__

namespace COPT
{

/**			LU solver of a general linear system. This is actually a wrapper of famous library 
  *			lapack. The solver aims to factorize any input matrix A with A=LU. 
  *
  *
  */
template<class Matrix>
class LU
	:
	public LinearSolver<Matrix>
{
private:
	typedef typename Matrix::scalar 			scalar;
	typedef typename Matrix::index 				index;

	/** the corresponding dynamic vector type for Matrix */
	typedef typename Matrix::DVector 			DVector;
	/** the corresponding dynamic matrix type for Matrix */
	typedef typename Matrix::DMatrix 			DMatrix;

	typedef typename Matrix::AbstractMatrix		AbstractMatrix;
	typedef typename Matrix::AbstractVector		AbstractVector;

	/** the array */
	scalar 						*__a;
	/** the size of array */
	index 						__size;
	/** piv */
	index 						*__piv;
	/** the factorization information */
	index 						__info;

	void doCompute (const Matrix& mat);

	virtual DVector doSolve( const AbstractVector& b );
	virtual DMatrix doSolve( const AbstractMatrix& b );

public:
	/** constructor and deconstructor */
	//%{
	LU ( );
	LU ( const Matrix& mat );
	~LU();
	//%}

	/** inverse matrix */
	DMatrix inverse ( );

	/** clear */
	void clear( );

};

/*************Implementation of class 'LU'*******************/

template<class Matrix>
LU<Matrix>::LU()
	:
	__a(NULL),
	__piv(NULL)
{
}

template<class Matrix>
LU<Matrix>::LU(const Matrix& mat)
	:
	__a(NULL),
	__piv(NULL)
{
	this->compute(mat);
}

template<class Matrix>
LU<Matrix>::~LU()
{
	clear();	
}

template<class Matrix>
void LU<Matrix>::doCompute( const Matrix& mat )
{
	clear();
	__a = new scalar[mat.size()];
	__piv = new index[std::min(mat.rows(),mat.cols())];
	this->setLDA(mat.lda());
	this->setRowNum(mat.rows());
	this->setColNum(mat.cols());
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	copt_lapack_getrf(this->rowNum(),this->colNum(),__a,this->lda(),__piv,&__info);
}

template<class Matrix>
typename LU<Matrix>::DVector LU<Matrix>::doSolve( const AbstractVector& b )
{
	this->squareValidation();
	if ( this->rowNum() != b.dimension() )
	{
		std::cerr<<"the order of matrix is "<<this->rowNum()<<" and the size of vector is "<<b.dimension()<<std::endl;
		throw COException("Linear system solving error: the size is not consistent!");
	}
	DVector result(b);
	copt_lapack_getrs('N',this->rowNum(),1,__a,this->lda(),__piv,result.dataPtr(),result.dimension(),&__info);
	return result;
}

template<class Matrix>
typename LU<Matrix>::DMatrix LU<Matrix>::doSolve( const AbstractMatrix& b )
{
	this->squareValidation();
	if( this->rowNum() != b.rows() )
	{
		std::cerr<<"the order of matrix is "<<this->rowNum()<<" and the size of right hand vectors are "<<b.rows()<<std::endl;
		throw COException("Linear system solving error: the size is not consistent!");
	}
	DMatrix result(b);
	copt_lapack_getrs('N',this->rowNum(),b.cols(),__a,this->lda(),__piv,result.dataPtr(),result.lda(),&__info);
	return result;
}

template<class Matrix>
typename LU<Matrix>::DMatrix LU<Matrix>::inverse( )
{
	this->squareValidation();
	DMatrix result(this->rowNum(),this->rowNum(),__a);
	copt_lapack_getri(this->rowNum(),result.dataPtr(),this->rowNum(),__piv,&__info);
	return result;
}

template<class Matrix>
void LU<Matrix>::clear()
{
	SAFE_DELETE_ARRAY(__a);
	SAFE_DELETE_ARRAY(__piv);
}

}

#endif