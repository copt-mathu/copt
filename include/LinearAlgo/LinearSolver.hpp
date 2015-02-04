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


#ifndef LINEAR_SOLVER_HPP__
#define LINEAR_SOLVER_HPP__


namespace COPT
{
/*			A general design for solving linear system 'LinearSolver'.
 *			LinearSolver takes the type of matrix as template. The class
 *			stores the size and lda of input matrix.
 *			A derived class from LinearSolver has to implement several
 *			special functions as the computation of the linear system and
 *			the solving of the system.
 */
template<class Matrix>
class LinearSolver
	:
	COPTObject,
	noncopyable
{

	typedef typename Matrix::scalar 		scalar;
	typedef typename Matrix::index 			index;

	/** the corresponding dynamic vector type for Matrix */
	typedef typename Matrix::DVector 			DVector;
	/** the corresponding dynamic matrix type for Matrix */
	typedef typename Matrix::DMatrix 			DMatrix;

	typedef typename Matrix::AbstractMatrix		AbstractMatrix;
	typedef typename Matrix::AbstractVector		AbstractVector;

	/** private varibles */
	//%{
	/** the row number of matrix */
	index 					__m;
	/** the column number of matrix */
	index 					__n;
	/** lda of matrix */
	index 					__lda;
	//%}

	/** virtual function for implementation */
	//%{
	/** computation of the solver */
	virtual void doCompute( const Matrix& mat ) = 0;
	/** solve a linear system with a right hand vectors */
	virtual DVector doSolve( const AbstractVector& b ) = 0;
	virtual DMatrix doSolve( const AbstractMatrix& b ) = 0;
	//%}

public:

	virtual ~LinearSolver() {}

	/** compute the solver */
	virtual void compute ( const Matrix& mat );

	/** solving one right hand vector */
	virtual DVector solve ( const AbstractVector& b );

	/** solving several right hand vectors */
	virtual DMatrix solve ( const AbstractMatrix& b );

	/** clear the solver */
	virtual void clear( ) = 0;

	/** square validation */
	virtual void squareValidation() const;

	/** setter and getter */
	//%{
	/** lda of the matrix */
	void setLDA( const index lda );
	index lda( ) const;
	/** row number */
	void setRowNum( const index m );
	index rowNum() const;
	/** column number */
	void setColNum( const index n );
	index colNum() const;
	//%}
};

/***********Implementation of class 'LinearSolver'****************/
template<class Matrix>
void LinearSolver<Matrix>::compute( const Matrix& mat )
{
	this->doCompute(mat);
}

template<class Matrix>
typename LinearSolver<Matrix>::DVector LinearSolver<Matrix>::solve ( const AbstractVector& b )
{
	return this->doSolve( b );
}

template<class Matrix>
typename LinearSolver<Matrix>::DMatrix LinearSolver<Matrix>::solve ( const AbstractMatrix& b )
{
	return this->doSolve(b);
}

template<class Matrix>
void LinearSolver<Matrix>::squareValidation() const
{
	if ( __n != __m )
		throw COException("Linear system solving error: the matrix is not square! ");
}

template<class Matrix>
void LinearSolver<Matrix>::setLDA( const index lda )
{
	__lda = lda;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::lda() const
{
	return __lda;
}

template<class Matrix>
void LinearSolver<Matrix>::setRowNum( const index m )
{
	__m = m;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::rowNum() const
{
	return __m;
}

template<class Matrix>
void LinearSolver<Matrix>::setColNum( const index n )
{
	__n = n;
}

template<class Matrix>
typename LinearSolver<Matrix>::index LinearSolver<Matrix>::colNum() const
{
	return __n;
}
/////////////End of implementation of class 'LinearSolver'

/* 
 *		some construction of COPT linear solvers 
 */
template<class Matrix> class LU;
template<class Matrix> class QR;
template<class Matrix> class Cholesky;
template<class Matrix> class EigenSolver;



}
#endif