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


#ifndef MATRIX_IMPL_HPP__
#define MATRIX_IMPL_HPP__

namespace COPT
{

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase()
	:
	Array(),
	__rows((RowAtCompileTime==Dynamic)?0:RowAtCompileTime),
	__cols((ColAtCompileTime==Dynamic)?0:ColAtCompileTime),
	__sym(false),
	__trans(false),
	__is_row_dynamic((RowAtCompileTime==Dynamic)?true:false),
	__is_col_dynamic((ColAtCompileTime==Dynamic)?true:false)
{
	if(RowAtCompileTime!=Dynamic)
		__lda = RowAtCompileTime+1;
	else
		__lda = 0;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	index m,
	index n,
	const scalar* data,
	bool trans)
	:
	Array(),
	__rows(m),
	__cols(n),
	__sym(false),
	__trans(trans),
	__lda(trans?(n+1):(m+1)),
	__is_row_dynamic(false),
	__is_col_dynamic(false)
{
	assert(RowAtCompileTime==Dynamic);
	assert(ColAtCompileTime==Dynamic);
	this->setArray(trans?(m*(n+1)):((m+1)*n),data);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	const index m,
	const scalar *data)
	:
	Array(),
	__sym(false),
	__trans(false)
{
	if(RowAtCompileTime!=Dynamic&&ColAtCompileTime!=Dynamic)
		throw COException("Matrix constructing error: One dimension constructor for Matrix requires Row number of column number is fixed at compile time");
	if(RowAtCompileTime==Dynamic)
	{
		__rows = m;
		__cols = ColAtCompileTime;
		this->setArray((m+1)*ColAtCompileTime,data);
		__is_row_dynamic = false;
		__is_col_dynamic = true;
		__lda = m+1;
	}
	else
	{
		__rows = RowAtCompileTime;
		__cols = m;
		this->setArray((RowAtCompileTime+1)*m,data);
		__is_row_dynamic = true;
		__is_col_dynamic = false;
		__lda = RowAtCompileTime+1;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	const MatrixBase& mat)
	:
	Array(mat.size(),mat.dataPtr()),
	__rows(mat.rows()),
	__cols(mat.cols()),
	__sym(mat.isSymmetric()),
	__trans(mat.isTranspose()),
	__lda(mat.lda())
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::MatrixBase(
	const AbstractMatrix& mat)
{
	assert(!(RowAtCompileTime!=Dynamic&&RowAtCompileTime!=mat.rowAtCompileTime()&&RowAtCompileTime!=mat.rows()));
	assert(!(ColAtCompileTime!=Dynamic&&ColAtCompileTime!=mat.colAtCompileTime()&&ColAtCompileTime!=mat.cols()));
	this->setArray(mat.size(),mat.dataPtr());
	__rows = mat.rows();
	__cols = mat.cols();
	__sym = mat.isSymmetric();
	__trans = mat.isTranspose();
	__lda = mat.lda();
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::~MatrixBase()
{
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
index MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rows() const
{
	return __rows;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isRowDynamic() const
{
	return __is_row_dynamic;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
index MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::cols() const
{
	return __cols;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isColumnDynamic() const
{
	return __is_col_dynamic;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator() (const index i,const index j)
{
	if(i<0||j<0)
		throw COException("MatrixBase error: index is less than zero!");
	else if (i>=__rows||j>=__cols)
		throw COException("MatrixBase error: index is out of range!");
	else if (__sym)
	{
		if ( i <= j )
			return this->operator[](j*__lda+i);
		else
			return this->operator[](i*__lda+j);
	}
	else if(__trans)
	{
		return this->operator[](i*__lda+j);
	}
	else{
		return this->operator[](j*__lda+i);
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator() ( const index i , const index j ) const
{
	return const_cast<MatrixBase&>(*this).operator()(i,j);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::set(const index i , const scalar value)
{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= this->size() )
		throw COException("MatrixBase error: index is out of range!");
	else
		this->operator[](i) = value;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::scalar& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::data ( const index i ) const
{
	if ( i < 0 )
		throw COException("MatrixBase error: index is less that zero!");
	else if ( i >= this->size() )
		throw COException("MatrixBase error: index is out of range!");
	else
		return this->operator[](i);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::row(const index num){
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range when getting a row of the matrix!");
	else if(__trans)
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->lda(),1);
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->lda());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::row(const index num )const {
	if ( num >= __rows || num < 0 )
		throw COException("MatrixBase error: row index out of range when getting a row of the matrix!");
	else if(__trans)
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->lda(),1);
	else
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->lda());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::col(const index num){
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range when getting a column of the matrix!");
	else if(__trans)
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->lda());
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->lda(),1);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const VectorBase<scalar,index> MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::col(const index num) const{
	if ( num >= __cols || num < 0 )
		throw COException("MatrixBase error: col index out of range when getting a column of the matrix!");
	else if(__trans)
		return Vector(this->cols(),referred_array(),this->dataPtr()+num,this->rows());
	else
		return Vector(this->rows(),referred_array(),this->dataPtr()+num*this->rows(),1);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
index MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::lda() const
{
	return __lda;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
scalar* MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::dataPtr()
{
	return Array::dataPtr();
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
const scalar* MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::dataPtr() const
{
	return Array::dataPtr();
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
index MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::size() const
{
	return Array::size();
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
int MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowAtCompileTime() const
{
	return RowAtCompileTime;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
int MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::colAtCompileTime() const
{
	return ColAtCompileTime;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setSymmetricFlag( bool sym )
{
	__sym = sym;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isSymmetric() const
{
	return __sym;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setTransposeFlag( bool trans )
{
	__trans = trans;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
bool MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::isTranspose() const
{
	return __trans;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::resize(index m,index n)
{
	if(__is_row_dynamic&&__is_col_dynamic)
	{
		__rows = m;
		__cols = n;
		if( __trans )
		{
			__lda = __cols+1;
			this->reset(__rows*__lda);
		}
		else
		{
			__lda = __rows+1;
			this->reset(__lda*__cols);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::resize(const MSize<index>& s)
{
	this->resize(s.m(),s.n());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::binaryCheck(const Mat&mat) const
{
	if(!is_same<typename MatrixBase::scalar,typename Mat::scalar>::value||!is_same<typename MatrixBase::index,typename Mat::index>::value)
		throw COException("Matrix binary check failed: the scalar type or the index type of two matrices is not the same! Please check it");
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>& MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator=(const Mat& mat)
{
	binaryCheck(mat);
	// if this matrix is a pure dynamic matrix
	if(__is_row_dynamic&&__is_col_dynamic)
	{
		__rows = mat.rows();
		__cols = mat.cols();
		this->setArray(mat.size(),mat.dataPtr());
	}
	// the column number is dynamic, row number is fixed
	else if(__is_col_dynamic)
	{
		if(__rows!=mat.rows())
			throw COException("Copy operator error: The matrix's row number is fixed and the row number of the input matrix is not the same!");
		else{
			__cols = mat.cols();
			this->setArray(mat.size(),mat.dataPtr());
		}
	}
	// the row number is dynamic, column number is fixed
	else if(__is_row_dynamic)
	{
		if(__cols!=mat.cols())
			throw COException("Copy operator error: The matrix's column number is fixed at compile time and the column number of the input matrix is not the same!");
		else
		{
			__rows = mat.rows();
			this->setArray(mat.size(),mat.dataPtr());
		}
	}
	// both row and column number are fixed
	else
	{
		if(__rows!=mat.rows()||__cols!=mat.cols())
			throw COException("Copy operator error: The matrix's size is fixed at compile time and the row number or column number of the input matrix is not the same!");
		else
		{
			this->setArray(mat.size(),mat.dataPtr());
		}
	}
	__lda = mat.lda();
	__trans = mat.isTranspose();
	__sym = mat.isSymmetric();
	return *this;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator+(const Mat& mat ) const
{
	binaryCheck(mat);
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	DMatrix result(*this);
	blas::copt_blas_axpy(this->size(),1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator-(const Mat& mat ) const
{
	binaryCheck(mat);
	if ( __rows != mat.rows() || __cols != mat.cols() )
	{
		throw COException("MatrixBase summation error: the size of two matrices are not consistent!");
	}
	DMatrix result(*this);
	blas::copt_blas_axpy(this->size(),-1.0,mat.dataPtr(),1,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Vec>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DVector MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::multi(const Vec &vec, const vector_object&)const
{
	binaryCheck(vec);
	if ( __cols != vec.size() )
		throw COException("MatrixBase multiply error: the size of MatrixBase and vector are not consistent!");
	DVector result(__rows);
	if (__sym)
		blas::copt_blas_symv(CblasColMajor,CblasUpper,__rows,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	else if (__trans)
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemv(CblasColMajor,CblasTrans,__cols,__rows,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
		else if ( is_complex<scalar>::value )
			blas::copt_blas_gemv(CblasColMajor,CblasConjTrans,__cols,__rows,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	}
	else
	{
		blas::copt_blas_gemv(CblasColMajor,CblasNoTrans,__rows,__cols,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	}
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class T>
typename T::DType MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator*(const T &t )const
{
	return multi(t,typename T::ObjectCategory());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::multi(const Mat& mat,const matrix_object&)const
{
	binaryCheck(mat);
	if ( __cols != mat.rows() )
		throw COException("MatrixBase multiply error: the size of two matrices are not consistent!");
	DMatrix result(__rows,mat.cols());
	if(__trans&&mat.isTranspose())
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
		else if( is_complex<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasConjTrans,CblasConjTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	}
	else if( __trans )
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
		else if ( is_complex<scalar>::value ){
			blas::copt_blas_gemm(CblasColMajor,CblasConjTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
		}
	}
	else if (mat.isTranspose())
	{
		if( is_real<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
		else if( is_complex<scalar>::value )
			blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasConjTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	}
	else
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,__rows,mat.cols(),__cols,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operator*(const scalar s)const
{
	DMatrix result(rows(),cols());
	for ( index i = 0 ; i < rows() ; ++ i )
		for ( index j = 0 ; j < cols() ; ++ j )
			result(i,j) = this->operator()(i,j)*s;
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transpose() const
{
	if(!__trans)
		return DMatrix(__cols,__rows,this->dataPtr(),true);
	else
		return DMatrix(__cols,__rows,this->dataPtr());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class T>
typename T::DType MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transMulti(const T &t)const
{
	return transMulti(t,typename T::ObjectCategory());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Vec>
typename  MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DVector MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transMulti( const Vec& vec, const vector_object& ) const
{
	if(__rows != vec.size() )
		throw COException("transpose multiplication error: the size of vector and matrix is not consistent!");
	DVector result(__cols);
	if(__trans)
		blas::copt_blas_gemv(CblasColMajor,CblasNoTrans,__cols,__rows,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	else
		blas::copt_blas_gemv(CblasColMajor,CblasTrans,__rows,__cols,1.0,this->dataPtr(),lda(),vec.dataPtr(),vec.interval(),0.0,result.dataPtr(),1);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::transMulti(const Mat &mat,const matrix_object&) const
{
	if(__rows != mat.rows() )
		throw COException("transpose multiplication error: the size of two matrices are not consistent!");
	DMatrix result(__cols,mat.cols());
	if (__trans&&mat.isTranspose())
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	else if(__trans)
		blas::copt_blas_gemm(CblasColMajor,CblasNoTrans,CblasNoTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	else if(mat.isTranspose())
		blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasTrans,__cols,mat.cols(),__rows,1.0,this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	else
		blas::copt_blas_gemm(CblasColMajor,CblasTrans,CblasNoTrans,__cols,mat.cols(),__rows,1.0,
			this->dataPtr(),lda(),mat.dataPtr(),mat.lda(),0.0,result.dataPtr(),result.lda());
	return result;
}

template<class Matrix>
class LU;

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Vec>
typename Vec::DType MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::solve(const Vec& vec)
{
	LU<MatrixBase> lu(*this);
	return lu.solve(vec);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::identity(
	index m,
	index n,
	const scalar s)
{
	DMatrix result(m,n);
	index min = std::min(m,n);
	for ( index i = 0 ; i < min ; ++ i )
		result(i,i) = s;
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Vec>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::diag(index m, index n, const Vec& vec)
{
	DMatrix result(m,n);
	index mi = std::min(m,n);
	if(vec.size()!=mi)
	{
		throw COException("Matrix diag error: please check the size of input vector!");
	}
	for (index i=0; i < mi; ++ i)
		result(i,i) = vec(i);
	return result;
}


template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums,const std::set<index>& colnums)
{
	if (*rownums.rbegin()>=mat.rows()||*colnums.rbegin()>=mat.cols()){
		std::cerr<<*rownums.rbegin()<<' '<<*colnums.rbegin()<<std::endl;
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),colnums.size());
	index r = 0,c = 0;
	for( typename std::set<index>::iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter)
	{
		c = 0;
		for ( typename std::set<index>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
		{
			this->operator()(r,c) = mat(*riter,*citer);
			++c;
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(const MatrixBase& mat,const std::set<index>& colnums)
{
	if(*colnums.rbegin()>=mat.cols()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(mat.rows(),colnums.size());
	index c = 0;
	for ( typename std::set<index>::const_iterator citer = colnums.begin() ; citer != colnums.end() ; ++ citer )
	{
		for (index r = 0 ; r < mat.rows() ; ++ r )
		{
			this->operator()(r,c) = mat(r,*citer);
		}
		++ c;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(const MatrixBase& mat,const std::set<index>& rownums)
{
	if(*rownums.rbegin()>=mat.rows()){
		throw COException("Index out of range in matrix blocking!");
		return;
	}

	this->resize(rownums.size(),mat.cols());
	index r = 0;
	for ( typename std::set<index>::const_iterator riter = rownums.begin() ; riter != rownums.end() ; ++ riter )
	{
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*riter,c);
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums,const std::vector<index>& colnums)
{
	this->resize(rownums.size(),colnums.size());
	for ( index r = 0 ; r < rownums.size() ; ++ r )
	{
		if(rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(rownums[r],colnums[c]);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime> template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::blockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns and rows at first
	index rownum = 0, colnum = 0;
	for (InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	for (InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(rownum,colnum);
	index r = 0, c = 0;
	for ( InputIterator ri = rowbegin ; ri != rowend ; ++ ri )
	{
		c = 0;
		if ( *ri >= mat.rows () )
			throw COException("Index out of range in matrix blocking!");
		for ( InputIterator ci = colbegin ; ci != colend ; ++ ci )
		{
			if ( *ci >= mat.cols() )
				throw COException("Column index out of range in matrix blocking!");
			this->operator()(r,c) = mat(*ri,*ci);
			++ c;
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& colnums)
{
	this->resize(mat.rows(),colnums.size());
	
	for ( index r = 0 ; r < mat.rows() ; ++ r )
	{
		for ( index c = 0 ; c < colnums.size() ; ++ c )
		{
			if(colnums[c]>=mat.cols())
				throw COException("Index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,colnums[c]);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime> template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::columnBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& colbegin,
	const InputIterator& colend)
{
	// count the number of columns
	index colnum = 0;
	for ( InputIterator iter = colbegin ; iter != colend ; ++ iter )
		++ colnum;
	this->resize(mat.rows(),colnum);
	for ( index r = 0 ; r < mat.rows() ; ++ r )
	{
		index c = 0;
		for ( InputIterator ci = colbegin ; ci != colend ; ++ ci )
		{
			if( *ci >= mat.cols())
				throw COException("Column index out of range in matrix blocking!");
			this->operator()(r,c) = mat(r,*ci);
			++ c;
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(const MatrixBase& mat,const std::vector<index>& rownums)
{
	this->resize(rownums.size(),mat.cols());
	for ( index r = 0 ; r < rownums.size() ; ++ r )
	{
		if (rownums[r]>=mat.rows())
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c)=mat(rownums[r],c);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime> 
template<class InputIterator>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::rowBlockFromMatrix(
	const MatrixBase& mat,
	const InputIterator& rowbegin,
	const InputIterator& rowend)
{
	// count the number of rows
	index rownum = 0;
	for ( InputIterator iter = rowbegin ; iter != rowend ; ++ iter )
		++ rownum;
	this->resize(rownum,mat.cols());
	index r = 0;
	for ( InputIterator ri = rowbegin ; ri != rowend ; ++ ri )
	{
		if ( *ri >= mat.rows() || *ri < 0 )
			throw COException("Index out of range in matrix blocking!");
		for ( index c = 0 ; c < mat.cols() ; ++ c )
		{
			this->operator()(r,c) = mat(*ri,c);
		}
		++ r;
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::combineAlongRow(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongRow(m1,m2,*this);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::combineAlongColumn(const MatrixBase& m1,const MatrixBase& m2)
{
	MatrixBase::stCombineAlongColumn(m1,m2,*this);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::stCombineAlongRow(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.cols()!=m2.cols())
		throw COException("Please make sure the column number of two matrices is the same before combination!");
	m.resize(m1.rows()+m2.rows(),m1.cols());
	for ( index i = 0 ; i < m1.cols() ; ++ i ){
		for ( index j = 0 ; j < m1.rows() ; ++ j ){
			m(j,i) = m1(j,i);
		}
		index n = m1.rows();
		for ( index j = 0 ; j < m2.rows() ; ++ j ){
			m(j+n,i) = m2(j,i);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::stCombineAlongColumn(const MatrixBase& m1,const MatrixBase& m2,MatrixBase& m)
{
	if(m1.rows()!=m2.rows())
		throw COException("Please make sure the row number of two matrices is the same before combination!");
	m.resize(m1.rows(),m1.cols()+m2.cols());
	for ( index i = 0 ; i < m1.rows() ; ++ i ){
		for (index j = 0 ; j < m1.cols() ; ++ j ){
			m(i,j) = m1(i,j);
		}
		index n = m1.cols();
		for (index j = 0 ; j < m2.cols() ; ++ j ){
			m(i,j+n)=m2(i,j);
		}
	}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::setRandom(const index rows,const index cols)
{
	if(rows<0||cols<0)
		throw COException("Please make sure that the number of row and column is bigger than zero!");
	std::uniform_real_distribution<typename get_pod_type<scalar>::type> unif(0.0,1.0);
	this->resize(rows,cols);
	for ( int i = 0 ; i < rows ; ++ i )
		for ( int j = 0 ; j < cols ; ++ j )
		{
			if(is_real<scalar>::value)
				ForceAssignment(unif(copt_rand_eng),this->operator()(i,j));
			else
				ForceAssignment(std::complex<typename get_pod_type<scalar>::type>(unif(copt_rand_eng),unif(copt_rand_eng)),this->operator()(i,j));
		}
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::DMatrix MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::random(const index rows,const index cols)
{
	DMatrix result(rows,cols);
	std::uniform_real_distribution<typename get_pod_type<scalar>::type> unif(0.0,1.0);
	for ( int i = 0 ; i < rows; ++ i )
		for ( int j = 0 ; j < cols ; ++ j )
		{
			if(is_real<scalar>::value)
				ForceAssignment(unif(copt_rand_eng),result(i,j));
			else
				ForceAssignment(std::complex<podscalar>(unif(copt_rand_eng),unif(copt_rand_eng)),result(i,j));
		}
	// result.setRandom(rows,cols);
	return result;
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
template<class Mat>
void MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::mtm(Mat &mat ) const
{
	int m = this->rows();
	int n = this->cols();
	if(mat.rows()!=n||mat.cols()!=n)
		mat.resize(n,n);
	if ( is_real<scalar>::value )
		blas::copt_blas_syrk(CblasColMajor,CblasUpper,CblasTrans,n,m,1.0,this->dataPtr(),lda(),0.0,mat.dataPtr(),mat.lda());
	else
		blas::copt_blas_herk(CblasColMajor,CblasUpper,CblasConjTrans,n,m,1.0,this->dataPtr(),lda(),0.0,mat.dataPtr(),mat.lda());
	for ( int i = 0 ; i < mat.rows() ; ++ i )
		for ( int j = 0 ; j < i ; ++ j )
			mat(i,j) = mat(j,i);
	mat.setSymmetricFlag(true);
}

template<class Matrix>
class PartialEigenSolver;

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::operationNorm() const
{
	DMatrix mtm;
	this->mtm(mtm);
	PartialEigenSolver<MatrixBase> solver(mtm);
	podscalar e = solver.computeLargestEigenvalue();
	return std::sqrt(e);
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::frobeniusNorm() const
{
	return copt_lapack_lange('f',__rows,__cols,const_cast<scalar *>(this->dataPtr()),this->lda());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::maxNorm() const
{
	return copt_lapack_lange('m',__rows,__cols,const_cast<scalar *>(this->dataPtr()),this->lda());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::oneNorm() const
{
	return copt_lapack_lange('o',__rows,__cols,const_cast<scalar *>(this->dataPtr()),this->lda());
}

template<class scalar,class index,int RowAtCompileTime,int ColAtCompileTime>
typename MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::podscalar MatrixBase<scalar,index,RowAtCompileTime,ColAtCompileTime>::infinityNorm() const
{
	return copt_lapack_lange('i',__rows,__cols,const_cast<scalar *>(this->dataPtr()),this->lda());
}
///////////////End of implementation of 'MatrixBase'



}// End of namespace COPT

#endif