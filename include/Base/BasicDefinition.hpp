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


#ifndef BASIC_DEFINITION_HPP__
#define BASIC_DEFINITION_HPP__

/*			This file introduces basic definitions of open source library
 *			COPT. 
 *
 */
namespace COPT
{

#ifdef _WIN64
typedef __int64 		COPTlong;
#else
typedef long 			COPTlong;
#endif

typedef unsigned long 	longsize;

const int Dynamic = -1;

/** the string that is used */
typedef std::string 			ostring;


const double ZERO_THRESH = 1e-10;				// the threshold to judge whether a scalar is zero
const int    MAX_SEARCH = 10000;				// default maximum number of search

const double	DEFAULT_CONVERGE_ERROR = 1e-5; 	// default converge error
const double 	DEFAULT_STEP_FOR_DIFFERENTIAL = 1e-5;

template<class T>
struct Infty{
	static inline T maximal(){
		return std::numeric_limits<T>::has_infinity()?std::numeric_limits<T>::infinity(): std::numeric_limits<T>::max();
	}
};

/*
 *				Judge that whether a scalar is zero
 */
template<class T>
inline bool IS_ZERO( T data )
{
	return fabs(data) < ZERO_THRESH ? true : false;
}

template<class T>
inline void SAFE_DELETE(T* value)
{
	if ( value ) { delete value; }
}

template<class T>
inline void SAFE_DELETE_ARRAY(T* array)
{
	if ( array ) {delete[] array; }
}

/**		base class who is not copyable */
class noncopyable
{
	noncopyable(const noncopyable& );
	const noncopyable& operator=(const noncopyable&);
protected:
	noncopyable(){}
	~noncopyable(){}
};


/** random engine */
static std::mt19937 copt_rand_eng(time(NULL));

/** the boolean operations on scalar types */
template<class T>
bool StrictLessThan(const T t1, const T t2)
{
	return t1<t2;
}

template<class T>
bool StrictLessThan(const std::complex<T>& t1, const std::complex<T>& t2 )
{
	std::cerr<<"Warning: no less than comparison on complex numbers is possible!"<<std::endl;
	return 0;
}

template<class T>
bool LessThan(const T t1, const T t2)
{
	return t1<=t2;
}

template<class T>
bool LessThan(const std::complex<T>& t1, const std::complex<T>& t2)
{
	std::cerr<<"Warning: no less than comparison on complex numbers is possible!"<<std::endl;
	return 0;
}

template<class T>
bool StrictLargerThan(const T t1, const T t2)
{
	return t1>t2;
}

template<class T>
bool StrictLargerThan(const std::complex<T>& t1, const std::complex<T>& t2 )
{
	std::cerr<<"Warning: the boolean opeations on complex numbers is not valid!"<<std::endl;
	return 0;
}

template<class T>
bool LargerThan(const T t1, const T t2)
{
	return t1>=t2;
}

template<class T>
bool LargerThan(const std::complex<T>& t1, const std::complex<T>& t2)
{
	std::cerr<<"Warning: the boolean opeations on complex numbers is not valid!"<<std::endl;
	return 0;
}


/** the types of linear solvers that COPT contains now */
enum LinearSolverType
{
LUSolver,				// LU solver for matrices (lapack)
QRSolver,				// QR solver for matrices (lapack)
CholeskySolver,			// Cholesky solver for symmetric semidefine matrix (lapack)
};

} // End of namespace COPT
#endif