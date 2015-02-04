/*****************************************************************************
 *
 *		This file introduces the basic usage of COPT's arithmetic type 
 *		'VectorBase'. 'VectorBase' is included in Header file 'Vector.hpp'
 * 		and feel free to read the header file for a detailed description.
 *		
 *		'VectorBase' describes widely-used mathematical type vector. 
 *		User can use simple syntax to initialize, assign and operate vectors
 * 		in COPT. Matlab-like assignment is also implemented for people who
 *		are familiar with Matlab. Basic operations like summation, subtraction
 *		and dot operation are easy to be called.
 *
 ****************************************************************************/

#include "Core"
#include <iostream>

/** basic declaration for COPT */
typedef double 							scalar;
typedef int 							ind;
typedef COPT::VectorBase<scalar,ind> 	Vector;

/** the size of vector */
constexpr int n = 5;

int main(int argc,char *argv[])
{
	scalar s[n]{1.0,4.0,5.0,2.1,3.2};

	/** initialize vector v1 based on a given array s */
	Vector v1(5,s);
	std::cout<<"v1 is "<<v1<<std::endl; // [1.0,4.0,5.0,2.1,3.2]

	/** any initialized vector is assumed to be zero */
	Vector v2(5);
	std::cout<<"v2 is "<<v2<<std::endl; // [0,0,0,0,0]

	/** Matlab-like assignment */
	v2(0) = 0.1;
	v2(3) = 5.1;
	std::cout<<"v2 becomes "<<v2<<std::endl; // [0.1,0,0,5.1,0]

	/** basic operations */
	std::cout<<"v1+v2="<<v1+v2<<std::endl; 		// [1.1,4.0,5.0,7.1,3.2]
	std::cout<<"(v1,v2)="<<v1.dot(v2)<<std::endl; // 10.81

	/** computation of l2 norm of v1 */
	std::cout<<"l2 norm of v1 is "<<v1.norm()<<std::endl; // 7.52662
	std::cout<<"l0 norm of v1 is "<<COPT::norm(v1,'0')<<std::endl; // 5

	/** vector blocking of v1 with indices 0 and 2 */
	ind inds[]{0,2};
	std::cout<<"blocking of v1 is "<<v1.block(2,inds,inds+2)<<std::endl; // [1.0,5.0]
	v2.blockFromVector(v1,2,inds,inds+2);
	std::cout<<"v2 becomes "<<v2<<std::endl; // [1.0,5.0]
}