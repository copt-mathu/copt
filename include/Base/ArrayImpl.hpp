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


#ifndef ARRAY_IMPL_HPP__
#define ARRAY_IMPL_HPP__

namespace COPT
{

/*********************Implementation of 'Array'************************/

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array()
	:
	__data_ptr(nullptr),
	__referred(false)
{
	// for dynamic array a zero-dimensional, non-referred array is constructed
	if(SizeAtCompileTime == Dynamic)
	{
		__size = 0;
		__inter = 1;
		__data_ptr = nullptr;
	}
	// for non-dynamic array, a zero-constant array is created
	else
	{
		__size = SizeAtCompileTime;
		__inter = 1;
		__data_ptr = new scalar[__size*__inter];
		for(index i=0; i<__size; ++i)
			__data_ptr[i]=0.0;
	}
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array(const index size, const scalar* data, const index inter)
	:
	__size(size),
	__inter(inter),
	__data_ptr(new scalar[size*inter]),
	__referred(false)
{
	if(SizeAtCompileTime!=Dynamic)
		throw COException("Size specified constructor only works for dynamic array");
	if ( data ){
		blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
	}
	else{
		for ( index i = 0 ; i < __size ; ++ i )
			__data_ptr[i] = static_cast<scalar>(0.0);
	}
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array(const scalar* data, const index inter)
	:
	__size(SizeAtCompileTime),
	__inter(inter),
	__data_ptr(nullptr),
	__referred(false)
{
	if(SizeAtCompileTime==Dynamic)
		throw COException("Wrong constructor for dynamic array. You should try to use size specified constructor or give the dimension of the array");
	__data_ptr = new scalar[__inter*__size];
	if(data)
		blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
	else{
		for (index i = 0; i < __size; ++ i)
			__data_ptr[i] = static_cast<scalar>(0.0);
	}
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array(const index size, const referred_array&, scalar *data, const index inter)
	:
	__size(size),
	__inter(inter),
	__data_ptr(data),
	__referred(true)
{
	if(SizeAtCompileTime!=Dynamic)
		throw COException("Size specified constructor only works for dynamic array");
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array(const referred_array&, scalar *data, const index inter)
	:
	__size(SizeAtCompileTime),
	__inter(inter),
	__data_ptr(data),
	__referred(true)
{
	if(SizeAtCompileTime==Dynamic)
		throw COException("Wrong constructor for dynamic array. You should try to use size specified constructor or give the dimension of the array");
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::Array(const Array& arr)
	:
	__data_ptr(nullptr)
{
	if(arr.isReferred())
		setReferredArray(arr.size(),arr.dataPtr(),arr.interval());
	else
		setArray(arr);
	__is_dynamic = (SizeAtCompileTime==Dynamic)?true:false;
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>::~Array()
{
	if (__referred){
		
		__data_ptr = nullptr;
	}
	else{
		SAFE_DELETE_ARRAY(__data_ptr);
	}
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::clear()
{
	if (__referred)
	{
		__data_ptr = NULL;
		__size = 0;
		__referred = false;
	}
	else
	{
		SAFE_DELETE_ARRAY(__data_ptr);
		__size = 0;
	}
}

template<class scalar,class index,int SizeAtCompileTime>
scalar* Array<scalar,index,SizeAtCompileTime>::dataPtr()
{
	return __data_ptr;
}

template<class scalar,class index,int SizeAtCompileTime>
const scalar* Array<scalar,index,SizeAtCompileTime>::dataPtr() const
{
	return __data_ptr;
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::copy(const Array& arr)
{
	if(this->isReferred())
		__referred = false;
	reset(arr.size(),arr.interval());
	blas::copt_blas_copy(__size,arr.dataPtr(),1,__data_ptr,arr.interval());
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::swap(Array& arr)
{
	if(arr.size()!=size())
		throw COException("the size of two arrays must be the same if anyone wants to swap them!");
	else if(arr.isReferred()||this->isReferred())
		throw COException("Swap requires that both arrays are not referred array");
	else
		blas::copt_blas_swap(__size,const_cast<scalar*>(arr.dataPtr()),1,__data_ptr,1);
}

template<class scalar,class index,int SizeAtCompileTime>
const index& Array<scalar,index,SizeAtCompileTime>::size() const
{
	return __size;
}

template<class scalar,class index,int SizeAtCompileTime>
bool Array<scalar,index,SizeAtCompileTime>::isReferred() const
{
	return __referred;
}

template<class scalar,class index,int SizeAtCompileTime>
bool Array<scalar,index,SizeAtCompileTime>::isDynamic() const
{
	return __is_dynamic;
}

template<class scalar,class index,int SizeAtCompileTime>
index Array<scalar,index,SizeAtCompileTime>::interval() const
{
	return __inter;
}

template<class scalar,class index,int SizeAtCompileTime>
typename Array<scalar,index,SizeAtCompileTime>::iterator Array<scalar,index,SizeAtCompileTime>::begin()
{
	return &__data_ptr[0];
}

template<class scalar,class index,int SizeAtCompileTime>
typename Array<scalar,index,SizeAtCompileTime>::const_iterator Array<scalar,index,SizeAtCompileTime>::begin() const
{
	return &__data_ptr[0];
}

template<class scalar,class index,int SizeAtCompileTime>
typename Array<scalar,index,SizeAtCompileTime>::iterator Array<scalar,index,SizeAtCompileTime>::end()
{
	return &__data_ptr[__size];
}

template<class scalar,class index,int SizeAtCompileTime>
typename Array<scalar,index,SizeAtCompileTime>::const_iterator Array<scalar,index,SizeAtCompileTime>::end() const
{
	return &__data_ptr[__size];
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::reset(const index size,const index inter)
{
	if(SizeAtCompileTime!=Dynamic)
		throw COException("Size specified array is not able to be resized");
	if(this->isReferred())
	{
		std::cerr<<"COPT Warning: it is dangerous to resize a reffered array. The previous pointer will not be deleted here"<<std::endl;
		clear();
		__referred = false;
		__data_ptr = new scalar[__size*__inter];
		for ( index i = 0 ; i < __size ; ++ i )
			__data_ptr[i*__inter] = static_cast<scalar>(0.0);
	}
	else
	{
		if(__size != size || __inter != inter ){
			clear();
			__size = size;
			__inter = inter;
			__data_ptr = new scalar[__size*__inter];
		}
		this->setZeros();
	}
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::setZeros()
{
	std::for_each(this->begin(),this->end(),[](scalar& s){s=0.0;});
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::setArray(const index size, const scalar *data, const index inter)
{
	if(SizeAtCompileTime!=Dynamic&&SizeAtCompileTime!=size)
		throw COException("Size specified array is not able to be resized");
	if (__referred)
		throw COException("referred array is not allowed to be reset ");
	else{
		if(SizeAtCompileTime==Dynamic) reset(size,inter);
		if(data)
			blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
	}
}

template<class scalar,class index,int SizeAtCompileTime>
template<int Size>
void Array<scalar,index,SizeAtCompileTime>::setFromArray(const Array<scalar,index,Size> &arr)
{
	if(SizeAtCompileTime!=Dynamic&&SizeAtCompileTime!=arr.size())
		throw COException("Size specified array is not able to be resized");
	if ( __referred )
		throw COException("referred array is not allowed to be reset ");
	else{
		reset(arr.size(),arr.interval());
		blas::copt_blas_copy(__size,arr.dataPtr(),arr.interval(),__data_ptr,__inter);
	}
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::setArray(const scalar *data)
{
	blas::copt_blas_copy(__size,data,1,__data_ptr,__inter);
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::setReferredArray(const index size, scalar* data, const index inter)
{
	if(SizeAtCompileTime!=Dynamic)
		throw COException("Size specified array is not able to be resized");
	__referred = true;
	__size = size;
	if(!__referred)
		SAFE_DELETE_ARRAY(__data_ptr);
	__data_ptr = data;
	__inter = inter;
}

template<class scalar,class index,int SizeAtCompileTime>
void Array<scalar,index,SizeAtCompileTime>::setReferredArray(scalar *data, const index inter)
{
	if(!__referred)
		SAFE_DELETE_ARRAY(__data_ptr);
	__referred = true;
	__data_ptr = data;
	__inter = inter;
}

template<class scalar,class index,int SizeAtCompileTime>
bool Array<scalar,index,SizeAtCompileTime>::isValid() const
{
	return is_scalar<scalar>::value;
}

template<class scalar,class index,int SizeAtCompileTime>
scalar& Array<scalar,index,SizeAtCompileTime>::operator[](index i)
{
	if ( i < 0 ){
		// index less than zero
		throw COException("Vector error, index less than zero.");
	}
	else if ( i >= __size ){
		// out of range
		throw COException("Vector error, index larger than the length.");
	}
	else
		return __data_ptr[i*__inter];
}

template<class scalar,class index,int SizeAtCompileTime>
const scalar& Array<scalar,index,SizeAtCompileTime>::operator[](index i)const
{
	return const_cast<Array&>(*this).operator[](i);
}

template<class scalar,class index,int SizeAtCompileTime>
Array<scalar,index,SizeAtCompileTime>& Array<scalar,index,SizeAtCompileTime>::operator=(const Array& arr)
{
	if (arr.isReferred())
	{
		__referred = true;
		__data_ptr = const_cast<scalar *>(arr.dataPtr());
		__size = arr.size();
		__inter = arr.interval();
	}
	else{
		__referred = false;
		copy(arr);
	}
	return *this;
}

}

#endif