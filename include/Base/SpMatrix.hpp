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


#ifndef SP_MATRIX_HPP__
#define SP_MATRIX_HPP__

namespace COPT
{


/*
 *		class Triplet for sparse matrix assignment
 */
template<class Scalar,class I>
struct TripletBase
{

public:
	typedef Scalar 								scalar;
	typedef I 									index;
	typedef KernelTrait<Scalar,I>				kernel;
	typedef typename kernel::size 				size;
private:
	/** private variables */
	//%{

	/** row index */
	index				__r;

	/** column index */
	index 				__c;

	/** value */
	Scalar 				__v;

	//%}

	/** private constructors */
	//%{
	TripletBase();
	//%}

public:

	/** constructor and deconstructor */
	//%{

	/** the only constructor */
	TripletBase(
		const index r,
		const index c,
		const Scalar v);

	~TripletBase();
	//%}

	/** getter (no other setter is allowed) */
	//%{

	/** row index */
	const index& rowIndex() const;

	/** column index */
	const index& columnIndex() const;

	/** value */
	const scalar& value() const;
	//%}
};

/*		compare two triplets according to row index
 *
 */
template<class Triplet>
struct rowComparison
{
	bool operator()(const Triplet& t1,const Triplet& t2);
};

/*		compare two triplets accordSize to column index
 *
 */
template<class Triplet>
struct columnComparison
{
	bool operator()(const Triplet& t1,const Triplet& t2);
};

/*		Sparse matrix class
 *		the sparse matrix is designed for solve sparse linear systems
 */

template<class SpMatrix>
class UMFLinearSolver;

template<class T,class I>
class SpMatrixBase
{
public:
	typedef 	T 						scalar;
	typedef  	I 	 					index;
	typedef 	TripletBase<T,I>		Triplet;
	typedef 	sp_matrix_object 		ObjectCategory;
	typedef 	KernelTrait<T,I>		kernel;
private:

	typedef 	VectorBase<T,I>			Vector;
	/** private variables */
	//%{

	/** the number of rows */
	index 				__rows;

	/** the number of columns */
	index 				__cols;

	/** the number of elements */
	index 				__elesize;

	/** the col poSizeers */
	index*				__colptr;

	/** the indices of the rows */
	index*				__rowind;

	/** the values */
	scalar*		 		__vals;

	/** static zero */
	static const scalar __zero;

	//%}

	/** private functions */
	//%{

	/** judge the rationality */
	void judgeRationality();

	//%}

public:

	/** constructor and deconstructor */
	//%{

	/** default constructor */
	SpMatrixBase();

	SpMatrixBase(
		const index 				rows,
		const index 				cols,
		const index 				elesize,
		const index*			 	colptr,
		const index*				rowind,
		const scalar*				vals);

	SpMatrixBase(
		const SpMatrixBase&);

	/** deconstructor */
	~SpMatrixBase();
	//*}

	/** getter and setter */
	//%{

	/** traditional setter of sparse matrix*/
	void setSparseMatrix(
		const index 					rows,
		const index 					cols,
		const index 					elesize,
		const index*			 		colptr,
		const index*			 		rowind,
		const scalar*			 		vals);

	/** overload of operator = */
	SpMatrixBase& operator = (const SpMatrixBase& );

	/** set from triplets */
	void setFromTriplets(
		const index rows,
		const index cols,
		std::vector<Triplet>& triplets);

	/** set from triplet iterator */
	template<class InputIterator>
	void setFromTriplets(
		const index rows,
		const index cols,
		const InputIterator& begin,
		const InputIterator& end);

	/** fast setting from triplets iterator (like mtx file) */
	template<class InputIterator>
	void fastSetFromTriplets(
		const index rows,
		const index cols,
		const InputIterator& begin,
		const InputIterator& end);

	/** clear the data */
	void clear();

	/** get the row number */
	const index& rows() const;

	/** get the column number */
	const index& cols() const;

	/** get the element size */
	const index& elementSize() const;

	/** get the column poSizeer */
	const index* columnPointer() const; 

	/** get the row indices */
	const index* rowIndex() const;

	/** get the values */
	const scalar* values() const;

	/** scale with a scalar s */
	void scale(const scalar s);

	/** negative sign of the matrix */
	void neg();

	//%}

	/** element access */
	//%{

	/** only const access is allowed */
	const scalar& operator()(
		const index i,
		const index j) const;

	//%}

	/** operations */
	//%{

	/** summation */
	SpMatrixBase operator+(const SpMatrixBase& mat) const;

	/** subtraction */
	SpMatrixBase operator-(const SpMatrixBase& mat) const;

	/** negative sign */
	SpMatrixBase operator-() const;

	/** multiplication with another sparse matrix */
	SpMatrixBase operator*(const SpMatrixBase& mat) const;

	/** multiplication with vector */
	Vector operator*(const Vector& vec) const;

	/** multiplication with scalar */
	SpMatrixBase operator*(const scalar s) const;

	/** transform to a dense matrix */
	MatrixBase<scalar,index> toDenseMatrix() const;

	/** solve a linear system */
	VectorBase<scalar,index> solve(const VectorBase<scalar,index>& vec);

	//%}

}; // end of class SpMatrixBase


/** Sparse matrix related operator */
template<class scalar,class index>
SpMatrixBase<scalar,index> operator* (const scalar s, const SpMatrixBase<scalar,index>& mat);
template<class scalar,class index,class T>
SpMatrixBase<scalar,index> operator* (const T s, const SpMatrixBase<scalar,index>& mat);

}

#endif