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


#ifndef NUM_DIFFERENTIAL_HPP__
#define NUM_DIFFERENTIAL_HPP__

/*
	compute the differential of functions
*/
namespace COPT{

/*		Fast but inaccurate method to compute differential
 *		The main thought is to use a small step length to compute differential quantities
 *		The method is described in Chapter 8 in 'Numerical Optimization' 
 *		/param func:		the input function
 *		/param x:			current point
 *		/param episilon: 	the step length
 *		return value: the approximating difference of the scalar function at x.
 */
template<class SFunc>
typename SFunc::ScalarType fastDifference(const SFunc& func,const typename SFunc::ScalarType x,const typename SFunc::ScalarType epsilon)
{
	return (func(x+epsilon)-func(x-epsilon)/(2*epsilon));
}

/*		Fast computation of second order difference
 *		/param func:		the input function
 *		/param x:			the input point
 *		/param epsilon:		the step length
 *			f''(x) = (f(x+h)+f(x-h)-f(x))/h^2
 */
template<class SFunc>
typename SFunc::ScalarType fastSecondDifference( const SFunc& func, const typename SFunc::ScalarType x,const typename SFunc::ScalarType epsilon)
{
	return (func(x+epsilon)+func(x-epsilon)-2*func(x))/(epsilon*epsilon);
}

/*		Fast computaton of partial difference of a vector function
 *		/param func:		the input function
 *		/param i:			the index of the partial variable
 *		/param x:			the input point
 *		/param episilon:	the step length
 *		return value: the approximating partial difference of the vector function at x
 */
template<class VFunc>
typename VFunc::ScalarType fastPartialDifference(const VFunc& func,const int i,const typename VFunc::Vector x,const typename VFunc::ScalarType epsilon)
{
	typedef typename VFunc::Vector 			Vector;
	Vector e = Vector::vecE(x.size(),i,epsilon);
	return (func(x+e)-func(x-e))/(2*epsilon);
}

/*		Fast computation of gradient of a vector function at certian point
 *		/param func:		the input function
 *		/param x:			the input point
 *		/param episilon:	the step length
 *		return value: the approximating gradient of a vector function at x
 */
template<class VFunc>
typename VFunc::Vector fastGradient(const VFunc& func,const typename VFunc::Vector x,const typename VFunc::ScalarType epsilon)
{
	typedef typename VFunc::Vector 			Vector;
	Vector g(x.size());
	for ( int i = 0 ; i < x.size() ; ++ i )
	{
		g[i] = fastPartialDifference(func,i,x,epsilon);
	}
	return g;
}
/*		Fast computation of second order partial difference
 *
 *
 */
template<class VFunc>
typename VFunc::ScalarType fastSecondPartialDifference(const VFunc& func,const int i,const int j,const typename VFunc::Vector& x,const typename VFunc::ScalarType epsilon)
{
	typedef typename VFunc::Vector 			Vector;
	Vector e1 = Vector::vecE(x.size(),i,epsilon);
	Vector e2 = Vector::vecE(x.size(),j,epsilon);
	return (func(x+e1+e2)-func(x+e1)-func(x+e2)+func(x))/(epsilon*epsilon);
}

template<class VFunc>
typename VFunc::Matrix fastHessianMatrix(
	const VFunc& func,
	const typename VFunc::Vector& x,
	const typename VFunc::ScalarType epsilon)
{
	typedef typename VFunc::Vector 			Vector;
	typedef typename VFunc::Matrix 			Matrix;
	Matrix result(x.size(),x.size());
	for ( int i = 0 ; i < x.size() ; ++ i )
		for ( int j = 0 ; j < x.size() ; ++ j )
			result(i,j) = fastSecondPartialDifference(func,i,j,x,epsilon);
	return result;
}


/*		class 'ScalarDifferential' taking scalar function as its template
 *
 *
 */
template<class SFunc>
class ScalarDifferential{
private:
	// the type of float number
	typedef typename SFunc::FT 			FT;

	// the reference to the scalar function
	const SFunc&						__func;

	// the error threshold
	FT 									__epsilon;
	// the iteration number that is used
	// this parameter is only used in dev version
	int 								__iterused;

	/*
		private Functions
	*/

	// central differentce
	FT centerDifference (FT x,FT h){
		return (__func(x+h)-__func(x-h))/(2*h);
	}

	// second order differential
	FT secondDifference (FT x,FT h){
		return (__func(x+h)+__func(x-h)-2*__func(x))/(h*h);
	}

public:
	/*
	 *	Constructor:
	 *		no default constructor is allowed
	 *		one must specify the function that is used
	 */
	ScalarDifferential(const SFunc& func,FT epsilon = 1e-5)
		:
		__func(func),
		__epsilon(epsilon),
		__iterused(0)
	{}

	// compute the differential
	FT diff(FT x, FT h = 0.01){
		FT error 		= __epsilon+1.0;
		FT forwarddiff 	= centerDifference(x,h);
		FT currdiff		= forwarddiff;
		__iterused		= 0;
		while ( error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= centerDifference(x,h);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

	// compute second order differential
	FT sDiff(FT x,FT h = 0.01){
		FT error 			= __epsilon + 1.0;
		FT forwarddiff 		= secondDifference(x,h);
		FT currdiff 		= forwarddiff;
		__iterused 			= 0;
		while (error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= secondDifference(x,h);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

};


/*			class 'VectorDifferential' is desigend to compute difference of a vector function
 *
 */
template<class VFunc>
class VectorDifferential{
public:
	enum Type{
		FAST,		// fast computation
		ACC 		// more accurate computation
	};
private:
	typedef typename VFunc::Vector 					Vector;
	typedef typename VFunc::Matrix 					Matrix;
	typedef typename Vector::ScalarType 					FT;
	/*
	 *		const reference to the function
	 */
	const VFunc& 									__vfunc;
	// 
	FT 												__epsilon;
	//
	int 											__iterused;
	//			the type of algorithm
	Type 											__type;

	/*
	 *				Private Function
	 */
	FT centerDifference(const Vector& vec,FT h,int i){
		Vector vp(vec),vm(vec);
		vp[i] += h;
		vm[i] -= h;
		return (__vfunc(vp)-__vfunc(vm))/(2*h);
	}

	/*
	 *
	 *
	 *
	 */
	FT computeDiff(const Vector& vec,int i,FT h = 0.01){
		FT error 		= __epsilon+1.0;
		FT forwarddiff 	= centerDifference(vec,h,i);
		FT currdiff		= forwarddiff;
		__iterused		= 0;
		while ( error > __epsilon ){
			h /= 2.0;
			forwarddiff = currdiff;
			currdiff 	= centerDifference(vec,h,i);
			error 		= fabs(currdiff - forwarddiff);
			if ( ++ __iterused > 20 ){
				std::cerr<<"Reach the max iteration number, differential result might be inaccurate"<<std::endl;
				break;
			}
		}
		return currdiff;
	}

	/*
	 *				compute differential of second order
	 */
	// FT secondDifference(const Vector& vec,FT h,int i,int j){
	// 	Vector vp(vec),vm(vec);
	// 	vp[i] +
	// }
public:

	

	VectorDifferential(const VFunc& func,FT epsilon = 1e-5,Type t=FAST) 
		: 
		__vfunc(func),
		__epsilon(epsilon),
		__type(t)
		{}
	/*
	 *			Main function, compute the gradient of '__func'
	 */
	Vector gradient(const Vector& vec,FT h = 1e-4){
		switch(__type){
		case FAST:
		{
			return fastGradient(__vfunc,vec,h);
		}
			break;
		case ACC:
		{
			Vector result(vec.size());
	 		for ( int i = 0 ; i < vec.size() ; ++ i ){
	 			result[i] = computeDiff(vec,i,h);
	 		}
	 		return result;
	 	}
	 		break;
	 	default:
	 	{
	 		throw COException("Unknown type in difference computation!");
	 	}
	 		break;
		}
	 }

	 Matrix hessian(const Vector& vec,FT h = 1e-5){
	 	switch(__type){
	 	case FAST:
	 	{
	 		return fastHessianMatrix(__vfunc,vec,h);
	 	}
	 		break;
	 	case ACC:
	 	{
	 		throw COException("No algorithm for accurate compuation of hessian matrix yet");
	 	}
	 		break;
	 	default:
	 	{
	 		throw COException("Unknown type in difference computation!");
	 	}
	 		break;
	 	}
	 }
};
};

#endif