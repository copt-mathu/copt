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


#ifndef ARITH_OPERATION_HPP__
#define ARITH_OPERATION_HPP__

namespace COPT
{

/** 		This header file introduces basic arithmetic operations like distance computation,
  *  		norm computation, computing mean vector and so on  for mathematical types
  * 		in open source library COPT.
  */

/** 		Check the dimension of two mathematical types */
template<class T1,class T2>
void checkDimension(const T1& t1, const T2& t2)
{
	checkDimension(t1,t2,typename T1::ObjectCategory(),typename T2::ObjectCategory());
}

/** 		Assert if the dimension of two vectors are not the same
  */
template<class V1,class V2>
void checkDimension(const V1& v1, const V2& v2, const vector_object&, const vector_object& )
{
	if(v1.dimension()!=v2.dimension())
		throw COException("COPT error: The size of two vectors are not consistent!");
}

template<class M1,class M2>
void checkDimension(const M1& m1, const M2& m2, const matrix_object&, const matrix_object& )
{
	if(m1.rows()!=m2.rows()||m1.cols()!=m2.cols())
		throw COException("COPT error: The size of two matrices are not consistent!");
}

/** 		Computation of the norm of a mathematical types including vector and matrix.
  */
template<class T>
typename T::podscalar norm(const T& t, char c)
{
	return norm(t,c,typename T::ObjectCategory());
}

/** 		Norm computation of a vector.
  */
template<class Vector>
typename Vector::podscalar norm(const Vector& vec, char c, const vector_object&)
{
	typedef typename Vector::index 		ind;
	typedef typename Vector::scalar 	scalar;
	typedef typename Vector::podscalar 	podscalar;
	switch (c)
	{
		// l0-norm
		case '0':
		{
			podscalar norm = 0.0;
			std::for_each(vec.begin(),vec.end(),[&norm](const scalar& s){norm+=(IS_ZERO(s)?0.0:1.0);});
			return norm;
		}
		break;
		// l1-norm
		case '1':
		{
			return vec.absNorm();
		}
		break;
		// l2-norm
		case '2':
		{
			return vec.norm();
		}
		break;
		// max-norm
		case 'm':
		case 'M':
		{
			podscalar norm = -1.0;
			std::for_each(vec.begin(),vec.end(),[&norm](const scalar& s){podscalar as=std::abs(s); if(norm<as) norm=as;});
			return norm;
		}
		break;
		default:
		{
			throw COException("Unknown type of norm!");
			return 0.0;
		}
		break;
	}
}

/** 		Norm computation of a matrix 
  */
template<class Matrix>
typename Matrix::podscalar norm(const Matrix& mat, char c, const matrix_object& )
{
	typedef typename Matrix::scalar 	scalar;
	typedef typename Matrix::podscalar 	podscalar;
	switch(c)
	{
		// l1-norm, max column sum
		case '1':
		{
			return mat.oneNorm();
		}
		break;
		// l2-norm, max singular value
		case '2':
		{
			return mat.operationNorm();
		}
		break;
		// max norm, max abs element
		case 'm':
		case 'M':
		{
			return mat.maxNorm();
		}
		break;
		// frobenium norm
		case 'f':
		case 'F':
		{
			return mat.frobeniusNorm();
		}
		break;
		// infinity norm, maximum row sum
		case 'i':
		case 'I':
		{
			return mat.infinityNorm();
		}
		break;
		default:
		{
			throw COException("Unknown type of matrix norm");
		}
		break;
	}
}

/** 		Compute the distance between two vectors 
  */
template<class T1,class T2>
typename T1::podscalar distance(const T1& t1, const T2& t2)
{
	return distance(t1,t2,typename T1::ObjectCategory(),typename T2::ObjectCategory());
}

/** 		Distance between two vectors.
  */
template<class Vector>
typename Vector::podscalar distance(const Vector &v1, const Vector &v2)
{
	checkDimension(v1,v2);
	return (v1-v2).norm();
}

/** 		Compute the mean vector of a set of vectors 
  */
template<class VectorIterator>
typename VectorIterator::value_type mean(VectorIterator begin,VectorIterator end)
{
	typename VectorIterator::value_type vec((*begin).dimension());
	int n=0;
	std::for_each(begin,end,[&vec,&n](typename VectorIterator::value_type& v){vec=vec+v;++n;});
	vec.scale(1.0/n);
	return vec;
}

/** 		Compute the sgn of a given type 
  */
template<class T>
T sgn(const T& t)
{
	typedef typename T::scalar scalar;
	T result(t);
	for_each(t.begin(),t.end(),[](scalar& s){s=(s<0)?-1:(s>0?1:0);});
}

/** 		Add sparse noise 
  */
template<class T,class S,class I>
void addSparseNoise(T& t,const I sp,const S n)
{
	std::vector<int> tt(t.size());
	for ( int i = 0 ; i < tt.size() ; ++i ) tt[i]=i;
	std::random_shuffle(tt.begin(),tt.end());
	std::uniform_real_distribution<typename T::podscalar> unif(-1.0,1.0);
	for ( int i = 0 ; i < sp ; ++ i ) t[i]+=n*unif(copt_rand_eng);
}


}

#endif