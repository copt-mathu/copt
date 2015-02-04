#include "LeastSquares.h"

// for dev
#include <iostream>

namespace COpt
{

// Least Square Solver
LeastSquareSolver::LeastSquareSolver()
	:
	__lstype(Tra),
	__dim(1)
{

}
LeastSquareSolver::LeastSquareSolver(int dim)
	:
	__lstype(Tra),
	__dim(dim)
{
	__v.resize(dim);
	std::cout<<"success !"<<std::endl;
}

void LeastSquareSolver::setDimension(int dim)
{
	__dim = dim;
	__v.resize(dim);
}

void LeastSquareSolver::solve(const Matrix& A,const Vector& b)
{
	// Assert the combination
	if (A.rows()!=b.size())
	{
		std::cerr<<"Warning: The shape of input matrix 'A' and vector 'b' is not consistent! Least square solver ends"<<std::endl;
		return;
	}
	else if (A.cols() != __dim)
	{
		std::cerr<<"Warning: The shape of input matrix 'A' is not consistent with the dimension of the problem. The dimension is changed automatically"<<std::endl;
		setDimension(A.cols());
	}
	
	switch (__lstype)
	{
	case Tra:
		solveTraditional(A,b);
		break;
	case LMS:
		solveLMS(A,b);
		break;
	case RLS:
		solveRLS(A,b);
		break;
	default:
		break;
	}

}

void LeastSquareSolver::solveTraditional(const Matrix& A,const Vector& b)
{
	// accurate algorithm that solves least squares
	// solving (A^TA)^-1A^Tb
	Matrix M = A.transpose()*A;
	Vector mb = A.transpose()*b;

	// using ldlt with high speed to solve the problem
	__v = M.ldlt().solve(mb);

	// test
	// std::cout<<"result is "<<__v<<std::endl;
}

void LeastSquareSolver::solveLMS(const Matrix& A,const Vector& b,double length)
{
	//least mean squares method
	// simple iteration
}

/*
	initialization of LMS
*/
void LeastSquareSolver::initializeLMS( ) 
{
	__v.resize(__dim);
	__v.setZero();
}
void LeastSquareSolver::initializeLMS( int dim )
{
	__dim = dim;
	__v.resize(__dim);
	__v.setZero();
}
void LeastSquareSolver::initializeLMS( int dim , const Vector& v)
{
	// give an initial weight vector v;
	__dim = dim;
	__v = v;
}
/*
	update of LMS
*/
void LeastSquareSolver::updateLMS(const Vector& a,double e)
{
}

void LeastSquareSolver::solveRLS(const Matrix& A,const Vector& b)
{
}

};