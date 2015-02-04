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


#ifndef QR_SOLVER_HPP__
#define QR_SOLVER_HPP__


namespace COPT
{
/*			QR factorization of a dense matrix. 
 *			The solver is kind of  Lapack wrapper which uses traditional 
 *			lapack solver to factorize a dense matrix. Unlike directly solving a linear system, this QR solver first factorize the matrix and stores necessary information. Then the coming right hand vectors are easily solved in a short time. What's more, the solver derives from 'LinearSolver'.
 *			
 *			Usage introduction:
 *				The QR solver takes the type of matrix as template. The first thing you have to do is to make sure what matrix is used. (real and complex matrices are supported )
 *			Kernel functions:
 *				void QR::compute( const Matrix& mat ): the function computes the input matrix and stores everthing need
 *				Vector QR::compute( const Vector& b ): solves a linear system Ax=b where A is the input matrix in compute.
 *				Matrix QR::compute( const Matrix& b ): similar to the former one but takes Matrix as input and output type meanning several right hand vectors are computed.
 *
 */
template<class Matrix>
class QR
	:
	public LinearSolver<Matrix>
{
private:
	typedef typename Matrix::scalar 			scalar;
	typedef typename Matrix::index 				index;

	typedef typename Matrix::DMatrix 			DMatrix;
	typedef typename Matrix::DVector 			DVector;

	typedef typename Matrix::AbstractVector 	AbstractVector;
	typedef typename Matrix::AbstractMatrix 	AbstractMatrix;

	/** private variables */
	//%{
	scalar 					*__a;
	/** Q part */
	scalar 					*__tau;
	/** information */
	index 					__info;
	//%}


	void doCompute ( const Matrix& mat );
	DVector doSolve( const AbstractVector& b );
	DMatrix doSolve( const AbstractMatrix& b );

public:

	/** constructor and deconstructor */
	//%{
	QR();
	QR(const Matrix& mat );
	~QR();
	//%}

	void clear();
};

template<class Matrix>
QR<Matrix>::QR()
	:
	__a(NULL),
	__tau(NULL)
{
}

template<class Matrix>
QR<Matrix>::QR( const Matrix& mat )
	:
	__a(NULL),
	__tau(NULL)
{
	this->compute(mat);
}

template<class Matrix>
QR<Matrix>::~QR()
{
	clear();
}

template<class Matrix>
void QR<Matrix>::doCompute( const Matrix& mat )
{
	clear();
	this->setLDA( mat.lda() );
	this->setRowNum( mat.rows() );
	this->setColNum( mat.cols() );
	__a = new scalar[mat.size()];
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	__tau = new scalar[std::min(mat.rows(),mat.cols())];
	copt_lapack_geqrf(this->rowNum(),this->colNum(),__a,this->lda(),__tau,&__info);
	if (__info != 0 )
		std::cerr<<"QR factorization warning: there is something wrong when doing QR factoriation in COPT!"<<std::endl;
}

template<class Matrix>
typename QR<Matrix>::DVector QR<Matrix>::doSolve( const AbstractVector& b )
{
	if( b.size() != this->colNum())
		throw COException("QR solving error: the size of matrix and vector are not consistent!");
	DVector result(b);
	copt_lapack_geqrs(this->rowNum(),this->colNum(),1,__a,this->lda(),__tau,result.dataPtr(),result.size(),&__info);
	return result;
}

template<class Matrix>
typename QR<Matrix>::DMatrix QR<Matrix>::doSolve( const AbstractMatrix& b )
{
	if ( b.rows() != this->colNum() )
		throw COException("QR solving error: the size of matrix and right hand vectors are not consistent!");
	DMatrix result(b);
	copt_lapack_geqrs(this->rowNum(),this->colNum(),result.cols(),__a,this->lda(),__tau,result.dataPtr(),
		result.lda(),&__info);
	return result;
}

template<class Matrix>
void QR<Matrix>::clear()
{
	SAFE_DELETE_ARRAY(__a);
	SAFE_DELETE_ARRAY(__tau);
}
}
#endif