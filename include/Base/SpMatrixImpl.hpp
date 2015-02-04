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


#ifndef SP_MATRIX_IMPL_HPP__
#define SP_MATRIX_IMPL_HPP__

namespace COPT
{

/*******************Implementation of Triplet******************/
template<class scalar,class index>
TripletBase<scalar,index>::TripletBase(
	const index r,
	const index c,
	const scalar v)
	:
	__r(r),
	__c(c),
	__v(v)
{
}

template<class scalar,class index>
TripletBase<scalar,index>::~TripletBase()
{
}

template<class scalar,class index>
const index& TripletBase<scalar,index>::rowIndex() const
{
	return __r;
}

template<class scalar,class index>
const index& TripletBase<scalar,index>::columnIndex() const
{
	return __c;
}

template<class scalar,class index>
const scalar& TripletBase<scalar,index>::value() const
{
	return __v;
}

/** 		Implementation of rowComparison 	*/
template<class Triplet>
bool rowComparison<Triplet>::operator() (const Triplet& t1,const Triplet& t2)
{
	return t1.rowIndex()<t2.rowIndex();
}

/**			Implementation of columnComparison	*/
template<class Triplet>
bool columnComparison<Triplet>::operator() (const Triplet& t1,const Triplet& t2)
{
	return t1.columnIndex()<t2.columnIndex();
}

/**			Implementation of SpMatrixBase 		*/

template<class scalar,class index>
void SpMatrixBase<scalar,index>::judgeRationality()
{
	for ( index i = 0 ; i < __elesize ; ++ i )
	{
		if ( __rowind[i] >= __rows )
			throw COException("Sparse matrix not rational: row index out of range!");
	}
}

template<class scalar,class index>
const typename SpMatrixBase<scalar,index>::scalar SpMatrixBase<scalar,index>::__zero = static_cast<scalar>(0.0);

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase()
	:
	__rows(0),
	__cols(0),
	__elesize(0),
	__colptr(NULL),
	__rowind(NULL),
	__vals(NULL)
{
}

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase(
	const index 						rows,
	const index 						cols,
	const index 						elesize,
	const index*						colptr,
	const index*						rowind,
	const scalar*						vals)
	:
	__rows(rows),
	__cols(cols),
	__elesize(elesize),
	__colptr(new index[cols+1]),
	__rowind(new index[elesize]),
	__vals(new scalar[elesize])
{

	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);
	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);

	judgeRationality();
}

template<class scalar,class index>
SpMatrixBase<scalar,index>::SpMatrixBase(
	const SpMatrixBase& mat)
	:
	__colptr(NULL),
	__rowind(NULL),
	__vals(NULL)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.columnPointer(),
		mat.rowIndex(),
		mat.values());
}

template<class scalar,class index>
SpMatrixBase<scalar,index>::~SpMatrixBase()
{
	if(__rowind)
		SAFE_DELETE_ARRAY(__rowind);
	if(__vals)
		SAFE_DELETE_ARRAY(__vals);
	if(__colptr)
		SAFE_DELETE_ARRAY(__colptr);
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::setSparseMatrix(
	const index 					rows,
	const index 					cols,
	const index 					elesize,
	const index*			 		colptr,
	const index*					rowind,
	const scalar*			 		vals)
{
	clear();
	__rows = rows;
	__cols = cols;
	__elesize = elesize;
	
	__rowind = new index[__elesize];
	__vals = new scalar[__elesize];
	__colptr = new index[__cols+1];

	blas::copt_blas_copy(__elesize,rowind,1,__rowind,1);
	blas::copt_blas_copy(__elesize,vals,1,__vals,1);
	blas::copt_blas_copy(__cols+1,colptr,1,__colptr,1);

	judgeRationality();
}

template<class scalar,class index>
SpMatrixBase<scalar,index>& SpMatrixBase<scalar,index>::operator=(const SpMatrixBase& mat)
{
	setSparseMatrix(
		mat.rows(),
		mat.cols(),
		mat.elementSize(),
		mat.columnPointer(),
		mat.rowIndex(),
		mat.values());
	return *this;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::setFromTriplets(
	const index rows,
	const index cols,
	std::vector<Triplet>& triplets)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(triplets.begin(),triplets.end(),columnComparison<Triplet>());
	// compute how many elements there are in one column
	index colind = 0 , ip = 0;
	for ( index i = 0 ; i < triplets.size() ; ++ i )
	{
		if(triplets[i].columnIndex()>colind)
		{
			colind = triplets[i].columnIndex();
			std::sort(triplets.begin()+ip,triplets.begin()+i,rowComparison<Triplet>());
			ip = i;
		}
	}
	// last sort
	std::sort(triplets.begin()+ip,triplets.end(),rowComparison<Triplet>());
	std::vector<index> 	colcounts(__cols,0);
	std::list<index> 		rowinds;
	std::list<scalar>	vals;
	colind = 0 , ip = 0;
	index count = 0;
	for ( index i = 0 ; i < triplets.size() ; ++ i )
	{
		if(triplets[i].columnIndex()>colind)
		{
			colcounts[colind] = count;
			colind = triplets[i].columnIndex();
			count = 0;
		}
		
		if(i==0)
		{
			rowinds.push_back(triplets[i].rowIndex());
			vals.push_back(triplets[i].value());
			++count;
		}
		else if(triplets[i].rowIndex()==triplets[i-1].rowIndex())
		{
			vals.back() += triplets[i].value();
		}
		else
		{
			rowinds.push_back(triplets[i].rowIndex());
			vals.push_back(triplets[i].value());
			++count;
		}
	}
	// last column
	colcounts[colind] = count;

	__colptr = new index[__cols+1];
	index columncount = 0;
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new index[rowinds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new scalar[vals.size()];
	i = 0;
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class scalar,class index> template<class InputIterator>
void SpMatrixBase<scalar,index>::setFromTriplets(
	const index rows,
	const index cols,
	const InputIterator& begin,
	const InputIterator& end)
{
	clear();
	__rows = rows;
	__cols = cols;
	// sort the triplets according to the column index at first
	std::sort(begin,end,columnComparison<Triplet>());
	// compute how many elements there are in one column
	index colind = 0 , ip = 0;
	InputIterator previter = begin;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
	{
		if(iter->columnIndex()>colind)
		{
			colind = iter->columnIndex();
			std::sort(previter,iter,rowComparison<Triplet>());
			previter = iter;
		}
	}
	// last sort
	std::sort(previter,end,rowComparison<Triplet>());
	std::vector<index> 		colcounts(__cols,0);
	std::list<index> 		rowinds;
	std::list<scalar>	vals;
	colind = 0 , ip = 0;
	previter = begin;
	index count = 0;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
	{
		if(iter->columnIndex()>colind)
		{
			colcounts[colind] = count;
			colind = iter->columnIndex();
			count = 0;
		}
		if(iter == begin)
		{
			rowinds.push_back(iter->rowIndex());
			vals.push_back(iter->value());
			++count;
		}
		else if(iter->rowIndex()==previter->rowIndex()&&iter->columnIndex()==previter->columnIndex())
		{
			vals.back() += iter->value();
		}
		else
		{
			rowinds.push_back(iter->rowIndex());
			vals.push_back(iter->value());
			++count;
		}
		previter = iter;
	}
	// last column
	colcounts[colind] = count;

	__colptr = new index[__cols+1];
	index columncount = 0;
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		__colptr[i] = columncount;
		columncount += colcounts[i];
	}
	__colptr[__cols] = columncount;

	__rowind = new index[rowinds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = rowinds.begin() ; iter != rowinds.end() ; ++ iter , ++ i )
		__rowind[i] = *iter;

	__vals = new scalar[vals.size()];
	i = 0;
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		__vals[i] = *iter;

	__elesize = rowinds.size();
}

template<class scalar,class index> template<class InputIterator>
void SpMatrixBase<scalar,index>::fastSetFromTriplets(
	const index rows,
	const index cols,
	const InputIterator& begin,
	const InputIterator& end)
{
	clear();
	__rows = rows;
	__cols = cols;
	__colptr = new index[__cols+1];
	// compute nnz
	index nnz = 0;
	for ( InputIterator iter = begin ; iter != end ; ++ iter )
		++ nnz;
	__elesize = nnz;
	__rowind = new index[nnz];
	__vals = new scalar[nnz];
	InputIterator previter = begin;
	index count = 0 , i = 0;
	std::vector<index> counts (__cols,0);
	for ( InputIterator iter = begin ; iter != end ; ++ iter , ++ i )
	{
		if(iter->columnIndex()!=previter->columnIndex())
		{
			counts[previter->columnIndex()] = count;
			count = 0;
			previter = iter;
		}
		__rowind[i] = iter->rowIndex();
		__vals[i] = iter->value();
		++ count;
	}
	counts[previter->columnIndex()] = count;
	count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		__colptr[c] = count;
		count += counts[c];
	}
	__colptr[__cols] = count;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::clear()
{
	__rows = 0;
	__cols = 0;
	__elesize = 0;

	if(__rowind)
	{
		SAFE_DELETE_ARRAY(__rowind);
		__rowind = NULL;
	}
	if(__vals)
	{
		SAFE_DELETE_ARRAY(__vals);
		__vals = NULL;
	}
	if(__colptr)
	{
		SAFE_DELETE_ARRAY(__colptr);
		__colptr = NULL;
	}
}

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::rows() const
{
	return __rows;
}

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::cols() const
{
	return __cols;
}

template<class scalar,class index>
const index& SpMatrixBase<scalar,index>::elementSize() const
{
	return __elesize;
}

template<class scalar,class index>
const index* SpMatrixBase<scalar,index>::columnPointer() const
{
	return __colptr;
}

template<class scalar,class index>
const index* SpMatrixBase<scalar,index>::rowIndex() const
{
	return __rowind;
}

template<class scalar,class index>
const scalar* SpMatrixBase<scalar,index>::values() const
{
	return __vals;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::scale(const scalar s)
{
	for ( index i = 0 ; i < __elesize ; ++ i )
		__vals[i] *= s;
}

template<class scalar,class index>
void SpMatrixBase<scalar,index>::neg()
{
	for ( index i = 0 ; i < __elesize ; ++ i )
		__vals[i] = -__vals[i];
}

template<class scalar,class index>
const scalar& SpMatrixBase<scalar,index>::operator()(
	const index i,
	const index j) const
{
	if(i>=__rows||j>=__cols)
		throw COException("Sparse Matrix error, index out of range!");
	index ip = __colptr[j],in=__colptr[j+1];
	for ( index ind = ip ; ind < in ; ++ ind )
	{
		if( i == __rowind[ind] )
			return __vals[ind];
	}
	return __zero;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator*(const SpMatrixBase& mat) const
{
	if(__cols != mat.rows() )
		throw COException("Multiplication error: matrix index does not fit!");
	std::list<Triplet> tris;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		index ci = __colptr[c] , cn = __colptr[c+1];
		for ( index r = ci ; r < cn ; ++ r )
		{
			index rind = __rowind[r];
			// traverse mat
			for ( index mc = 0 ; mc < mat.cols() ; ++ mc )
			{
				index mci = mat.columnPointer()[mc], mcn = mat.columnPointer()[mc+1];
				for ( index mr = mci ; mr < mcn ; ++ mr )
				{
					if ( c == mat.rowIndex()[mr] )
						tris.push_back(Triplet(rind,mc,__vals[r]*mat.values()[mr]));
				}
			}
		}
	}
	std::vector<Triplet> triplets;
	triplets.reserve(tris.size());
	for ( typename std::list<Triplet>::iterator iter = tris.begin() ; iter != tris.end() ; ++ iter )
		triplets.push_back(*iter);
	SpMatrixBase result;
	result.setFromTriplets(__rows,mat.cols(),triplets);
	return result;
}


template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator+ ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: index does not fit!");
	std::list<index> inds;
	std::list<scalar> vals;
	index *colptr = new index[__cols+1];
	index count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		index ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		index i1 = ci1,i2=ci2;
		while( i1<cn1 || i2<cn2 )
		{
			if( i1 < cn1 && i2 < cn2 )
			{
				if(__rowind[i1]==mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]+mat.values()[i2]);
					++ i1;
					++ i2;
				}
				else if(__rowind[i1]<mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]);
					++ i1;
				}
				else
				{
					inds.push_back(mat.rowIndex()[i2]);
					vals.push_back(mat.values()[i2]);
					++ i2;
				}
			}
			else if ( i1 < cn1 )
			{
				inds.push_back(__rowind[i1]);
				vals.push_back(__vals[i1]);
				++ i1;
			}
			else{
				inds.push_back(mat.rowIndex()[i2]);
				vals.push_back(mat.values()[i2]);
				++ i2;
			}
			++ count;
		}
	}
	colptr[__cols] = count;
	index *rowind = new index[inds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	scalar *vs = new scalar[vals.size()];
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);

	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator- ( const SpMatrixBase& mat ) const
{
	if( __cols != mat.cols() || __rows != mat.rows() )
		throw COException("Sparse matrix summation error: index does not fit!");
	std::list<index> inds;
	std::list<scalar> vals;
	index *colptr = new index[__cols+1];
	index count = 0;
	for ( index c = 0 ; c < __cols ; ++ c )
	{
		colptr[c] = count;
		index ci1 = __colptr[c],cn1 = __colptr[c+1],ci2 = mat.columnPointer()[c],cn2=mat.columnPointer()[c+1];
		index i1 = ci1,i2=ci2;
		while( i1<cn1 || i2<cn2 )
		{
			if( i1 < cn1 && i2 < cn2 )
			{
				if(__rowind[i1]==mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]-mat.values()[i2]);
					++ i1;
					++ i2;
				}
				else if(__rowind[i1]<mat.rowIndex()[i2])
				{
					inds.push_back(__rowind[i1]);
					vals.push_back(__vals[i1]);
					++ i1;
				}
				else
				{
					inds.push_back(mat.rowIndex()[i2]);
					vals.push_back(-mat.values()[i2]);
					++ i2;
				}
			}
			else if ( i1 < cn1 )
			{
				inds.push_back(__rowind[i1]);
				vals.push_back(__vals[i1]);
				++ i1;
			}
			else{
				inds.push_back(mat.rowIndex()[i2]);
				vals.push_back(-mat.values()[i2]);
				++ i2;
			}
			++ count;
		}
	}
	colptr[__cols] = count;
	index *rowind = new index[inds.size()];
	index i = 0;
	for ( typename std::list<index>::iterator iter = inds.begin() ; iter != inds.end() ; ++ iter , ++i )
		rowind[i] = *iter;
	i = 0;
	scalar *vs = new scalar[vals.size()];
	for ( typename std::list<scalar>::iterator iter = vals.begin() ; iter != vals.end() ; ++ iter , ++ i )
		vs[i] = *iter;
	SpMatrixBase result;
	result.setSparseMatrix(__rows,__cols,inds.size(),colptr,rowind,vs);
	delete[]rowind;
	delete[]colptr;
	delete[]vs;
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator- () const
{
	SpMatrixBase result(*this);
	result.neg();
	return result;
}

template<class scalar,class index>
typename SpMatrixBase<scalar,index>::Vector SpMatrixBase<scalar,index>::operator*(const Vector& vec) const
{
	if(__cols!=vec.size() )
		throw COException("Multiplication error: index does not fit!");
	Vector result(__rows);
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		index ip = __colptr[i] , in = __colptr[i+1];
		for ( index r = ip ; r < in ; ++ r ){
			result[__rowind[r]] += __vals[r]*vec[i];
		}
	}
	return result;
}


template<class scalar,class index, class T>
SpMatrixBase<scalar,index> operator* (const T s,const SpMatrixBase<scalar,index>& mat)
{
	return mat.operator*(static_cast<scalar>(s));
}

template<class scalar,class index>
MatrixBase<scalar,index> SpMatrixBase<scalar,index>::toDenseMatrix() const
{
	MatrixBase<scalar,index> result(__rows,__cols);
	for ( index i = 0 ; i < __cols ; ++ i )
	{
		index ip = __colptr[i] , in = __colptr[i+1];
		for ( index r = ip ; r < in ; ++ r )
		{
			result(__rowind[r],i) = __vals[r];
		}
	}
	return result;
}

template<class scalar,class index>
VectorBase<scalar,index> SpMatrixBase<scalar,index>::solve(const VectorBase<scalar,index>& vec)
{
	return UMFLinearSolver<SpMatrixBase>(*this).solve(vec);
}

template<class scalar,class index>
SpMatrixBase<scalar,index> SpMatrixBase<scalar,index>::operator* ( const scalar s ) const
{
	SpMatrixBase result(*this);
	result.scale(s);
	return result;
}

template<class scalar,class index>
SpMatrixBase<scalar,index> operator*(const scalar s,const SpMatrixBase<scalar,index>& mat)
{
	return mat.operator*(s);
}

}

#endif