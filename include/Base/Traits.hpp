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


#ifndef COPT_TRAITS_HPP__
#define COPT_TRAITS_HPP__

/*		In this file, we introduce several traits that are used in COPT.
 *		
 *
 */
namespace COPT
{

/** type judgetments and definition */
//%{
template<class T1>
struct get_pod_type
{ typedef T1 type;};

template<class T2>
struct get_pod_type<std::complex<T2> >
{ typedef T2 type;};

template<class T1,class T2>
struct is_same
{ static const bool value = false;};

template<class T>
struct is_same<T,T>
{ static const bool value = true;};

template<class T1>
struct is_float
{ static const bool value = false;};

template<>
struct is_float<float>
{ static const bool value = true;};

template<class T1>
struct is_double
{ static const bool value = false;};

template<>
struct is_double<double>
{ static const bool value = true;};

template<class T2>
struct is_real
{ static const bool value = false;};

template<>
struct is_real<float>
{ static const bool value = true;};

template<>
struct is_real<double>
{ static const bool value = true;};

template<class T1>
struct is_complex
{ static const bool value = false;};

template<class eT>
struct is_complex<std::complex<eT> >
{ static const bool value = false;};

template<>
struct is_complex<std::complex<float> >
{ static const bool value = true;};

template<>
struct is_complex<std::complex<double> >
{ static const bool value = true;};

template<class eT>
struct is_complex_float
{ static const bool value = false;};

template<>
struct is_complex_float<std::complex<float> >
{ static const bool value = true;};

template<class eT>
struct is_complex_double
{ static const bool value = false;};

template<>
struct is_complex_double<std::complex<double> >
{ static const bool value = true;};

template<class T>
struct is_size
{ static const bool value = false;};

template<class T>
struct is_index
{ 
	typedef T size;
	static const bool value = false;
};

template<>
struct is_index<int>
{ 
	typedef unsigned int size;
	static const bool value = true;
};

template<>
struct is_index<long>
{ 
	typedef unsigned long size;
	static const bool value = true;
};

template<class T>
struct is_unsigned_size
{ 
	typedef T index;
	static const bool value = false;
};

template<>
struct is_unsigned_size<unsigned int>
{ 
	typedef 	int 			index;
	static const bool value = true;
};

template<>
struct is_unsigned_size<unsigned long>
{ 
	typedef 	long 			index;
	static const bool value = true;
};

template<class T>
struct is_scalar
{ static const bool value = 
					is_scalar<T>::value||
					is_complex<T>::value;
};
//%}

/** kernel of COPT */
//%{
template<class T,class I,int SizeAtCompileTime>
class Array;
template<class T,class I,int SizeAtCompileTime>
class VectorBase;
template<class T,class I,int RowAtCompileTime,int ColAtCompileTime>
class MatrixBase;
template<class T,class I>
class SpMatrixBase;
/*		A trait class describing basic types that might be
 *		used in a numerical solver. A solver should take trait
 *		as template for flexibility.
 */
template<class T,class I = int >
class KernelTrait
{
public:
	typedef T 											scalar;
	typedef typename get_pod_type<scalar>::type 		podscalar;
	typedef I 											index;
	typedef typename is_index<I>::size					size;
	// whether the kernel is valid:
	static const bool valid  = is_scalar<T>::value&&is_index<I>::value;

	typedef COPT::Array<T,I,Dynamic> 					Array;
	typedef VectorBase<T,I,Dynamic>						Vector;
	typedef MatrixBase<T,I,Dynamic,Dynamic>				Matrix;
	typedef SpMatrixBase<T,I>							SpMatrix;
};
//%}

/** tags */
//%{
/** base tag for copt object */
struct copt_object{};

/** objects storing data */
struct data_object:copt_object{}; 				// the data object
struct referred_array:data_object{};			// refered array
struct array_object:data_object{}; 				// array
struct matrix_object:data_object{};				// normal matrix	
struct vector_object:data_object{};				// normal vector
struct sp_matrix_object:data_object{};			// sparse matrix

/** function objects */
struct scalar_func_object:copt_object{};		// scalar function
struct vector_func_object:copt_object{};		// vector function

/** problem object */
struct problem_object:copt_object{};			// problem
struct linear_programming:problem_object{};
struct lasso_problem:problem_object{};			// lasso problem

/** solver for problems */
struct solver_object:copt_object{};				// solver
struct linear_solver:solver_object{};			// linear solver
struct admm_solver:solver_object{};				// admm solver
struct proximal_solver:solver_object{}; 		// proximal solver
struct fista_solver:solver_object{};			// FISTA solver

/** constraint objects */
struct constraint_object:copt_object{};			// constraint

/** no time statistics tag */
struct time_stat_object:copt_object{};			// time statistics object
struct no_time_stat_tag:time_stat_object{};
struct solver_time_stat_tag:time_stat_object{};

/** traits of constraints and functions*/
struct linear_constraint:constraint_object{};
struct linear_leq_constraint:linear_constraint{};
struct linear_neq_constraint:linear_constraint{};
struct linear_eq_constraint:linear_constraint{};
struct quadratic_constraint_tag:constraint_object{};
struct non_linear_constraint_tag:constraint_object{};
//%}

template<class Constraint>
struct constraint_trait{
	typedef typename Constraint::constraint_category	constraint_category;
};



/** trais of functions */
struct scalar_function_tag{};
struct log_scalar_function_tag:public scalar_function_tag{};
struct abs_scalar_function_tag:public scalar_function_tag{};



/** force assignment between complex to real numbers */
template<class T>
void ForceAssignment(const T& t1,std::complex<T>& t2)
{
	// real to complex
	if( is_real<T>::value )
		t2 = std::complex<T>(t1,0);
}

template<class T>
void ForceAssignment(const std::complex<T>&t1 , T& t2 )
{
	if ( is_real<T>::value )
		t2 = t1.real();
}

template<class T>
void ForceAssignment(const T& t1, T& t2)
{
	t2 = t1;
}

/** get the dynamic type of mathematical type */
template<class T>
struct DynamicType{
	typedef T Type;
};
template<class scalar,class index,int Row>
struct DynamicType<MatrixBase<scalar,index,Row,Dynamic> >{
	typedef MatrixBase<scalar,index,Dynamic,Dynamic> Type;
};
template<class scalar,class index,int Col>
struct DynamicType<MatrixBase<scalar,index,Dynamic,Col> >{
	typedef MatrixBase<scalar,index,Dynamic,Dynamic> Type;
};
template<class scalar,class index,int Row,int Col>
struct DynamicType<MatrixBase<scalar,index,Row,Col> >{
	typedef MatrixBase<scalar,index,Dynamic,Dynamic> Type;
};

template<class FT,class I>
class AbstractMatrix;
template<class FT,class I>
class AbstractVector;

/** copt traits of numerical types */
template<class T>
struct copt_traits{
	typedef float 	scalar;
	typedef long 	index;
};

// template<class T>
// struct copt_traits<AbstractMatrix<T> >{
// 	typedef typename copt_traits<T>::scalar scalar;
// 	typedef typename copt_traits<T>::index index;
// };

// template<class T>
// struct copt_traits<AbstractVector<T> >{
// 	typedef typename copt_traits<T>::scalar scalar;
// 	typedef typename copt_traits<T>::index index;
// };

// template<class s,class i,int R,int C>
// struct copt_traits<MatrixBase<s,i,R,C> >{
// 	typedef s scalar;
// 	typedef i index;
// 	const static int RowAtCompileTime = R;
// 	const static int ColAtCompileTime = C;
// 	const static int SizeAtCompileTime = (RowAtCompileTime==Dynamic||ColAtCompileTime==Dynamic)?Dynamic:(RowAtCompileTime+1)*(ColAtCompileTime);
// };

// template<class s,class i,int Si>
// struct copt_traits<VectorBase<s,i,Si> >{
// 	typedef s scalar;
// 	typedef i index;
// 	const static int SizeAtCompileTime = Si;
// };

// template<class s,class i,int Si>
// struct copt_traits<Array<s,i,Si> >{
// 	typedef s scalar;
// 	typedef i index;
// };

}// End of namespace COPT


#endif