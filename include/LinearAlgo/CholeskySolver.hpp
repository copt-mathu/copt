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


#ifndef CHOLESKY_SOLVER_HPP__
#define CHOLESKY_SOLVER_HPP__

namespace COPT
{

/*			Cholesky factorization for symmetric semidefine matrix.
 *
 *
 */
template<class Matrix>
class Cholesky
	:
	public LinearSolver<Matrix>
{
private:
	typedef typename Matrix::index 				index;
	typedef typename Matrix::scalar 			scalar;
	typedef VectorBase<scalar,index>			Vector;

	typedef typename Matrix::DVector 			DVector;
	typedef typename Matrix::DMatrix 			DMatrix;

	typedef typename Matrix::AbstractVector 	AbstractVector;
	typedef typename Matrix::AbstractMatrix 	AbstractMatrix;

	/** the data */
	scalar 							*__a;

	/** info */
	index 							__info;

	/** since the input matrix has to be symmetric semidefinite */
	bool validation( const Matrix& mat ) const;

	/** implementation of virtual functions */
	//%{
	void doCompute( const Matrix& mat );
	DVector doSolve( const AbstractVector& b );
	DMatrix doSolve( const AbstractMatrix& b );
	//%}

public:

	/** constructor and deconstructor */
	//%{
	/** default constructor */
	Cholesky();
	Cholesky(const Matrix& mat);
	~Cholesky();
	//%}

	/** inver matrix */
	Matrix inverse ();

	/** clear */
	void clear();
};

template<class Matrix>
Cholesky<Matrix>::Cholesky()
	:
	__a(NULL)
{
}

template<class Matrix>
Cholesky<Matrix>::Cholesky( const Matrix& mat )
	:
	__a(NULL)
{
	this->compute(mat);
}

template<class Matrix>
Cholesky<Matrix>::~Cholesky()
{
	clear();
}

template<class Matrix>
void Cholesky<Matrix>::doCompute( const Matrix& mat )
{
	if (!validation(mat))
	{
		std::cerr<<"Warning in cholesky solver: the input matrix is not symmetric. Please confirm that!"<<std::endl;
		return;
	}
	clear();
	__a = new scalar[mat.size()];
	this->setRowNum(mat.rows());
	this->setColNum(mat.cols());
	this->setLDA(mat.lda());
	blas::copt_blas_copy(mat.size(),mat.dataPtr(),1,__a,1);
	copt_lapack_potrf('U',this->rowNum(),__a,this->lda(),&__info);
	if( __info != 0 )
	{
		std::cout<<__info<<std::endl;
		std::cerr<<"Warning in Cholesky solver: something computation is wrong!"<<std::endl;
	}
}

template<class Matrix>
typename Cholesky<Matrix>::DVector Cholesky<Matrix>::doSolve( const AbstractVector& b )
{
	if ( this->rowNum() != b.size() )
	{
		std::cerr<<"The order of matrix is "<<this->rowNum()<<" and the dimension of vector is "<<b.size()<<std::endl;
		throw COException("Cholesky solving error: the size if not consistent!");
	}
	DVector result(b);
	copt_lapack_potrs('U',this->rowNum(),1,__a,this->lda(),result.dataPtr(),result.size(),&__info);
	if ( __info != 0 )
		std::cerr<<"Warning in Cholesky solver: solving is wrong!"<<std::endl;
	return result;
}

template<class Matrix>
typename Cholesky<Matrix>::DMatrix Cholesky<Matrix>::doSolve( const AbstractMatrix& b )
{
	if ( this->rowNum() != b.rows() )
	{
		std::cerr<<"The order of matrix is "<<this->rowNum()<<" and the dimension of vector is "<<b.rows()<<std::endl;
		throw COException("Cholesky solving error: the size if not consistent!");
	}
	DMatrix result(b);
	copt_lapack_potrs('U',this->rowNum(),b.cols(),__a,this->lda(),result.dataPtr(),result.lda(),&__info);
	if ( __info != 0 )
		std::cerr<<"Warning in Cholesky solver: solving is wrong!"<<std::endl;
	return result;
}

template<class Matrix>
bool Cholesky<Matrix>::validation( const Matrix& mat )const
{
	return mat.isSymmetric();
}

template<class Matrix>
Matrix Cholesky<Matrix>::inverse( )
{
	this->squareValidation();
	Matrix result(this->rowNum(),this->rowNum(),__a);
	copt_lapack_potri('U',this->rowNum(),result.dataPtr(),this->lda(),&__info);
	result.setSymmetricFlag(true);
	return result;
}

template<class Matrix>
void Cholesky<Matrix>::clear()
{
	SAFE_DELETE_ARRAY(__a);
}

}

#endif