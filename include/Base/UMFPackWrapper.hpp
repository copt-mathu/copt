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


#ifndef UMFPACK_WRAPPER_HPP__
#define UMFPACK_WRAPPER_HPP__

#ifdef USE_UMFPACK
namespace COPT
{

/** wrappers for umf functions */
//%{
inline void umfpack_free_numeric(void **Numeric,double,int)
{
	umfpack_di_free_numeric(Numeric);
}

inline void umfpack_free_numeric(void **Numeric,double,COPTlong)
{ 
	umfpack_dl_free_numeric(Numeric); 
	*Numeric = NULL;
}

// inline void umfpack_free_numeric(void **Numeric,std::complex<double>)
// {
// 	umfpack_zl_free_numeric(Numeric);
// 	*Numeric = NULL;
// }
inline void umfpack_free_symbolic(void **Symbolic,double,int)
{
	umfpack_di_free_symbolic(Symbolic);
}

inline void umfpack_free_symbolic(void **Symbolic,double,COPTlong)
{
	umfpack_dl_free_symbolic(Symbolic);
	*Symbolic = NULL;
}

// inline void umfpack_free_symbolic(void **Symbolic,std::complex<double>)
// {
// 	umfpack_zl_free_symbolic(Symbolic);
// 	*Symbolic = NULL;
// }

inline int umfpack_symbolic(
	int rows,int cols,
	const int* colptr,const int* rowind, const double*vals,void **Symbolic,
	const double* Control,double* Info)
{
	return umfpack_di_symbolic(rows,cols,colptr,rowind,vals,Symbolic,Control,Info);
}

inline COPTlong umfpack_symbolic(
	COPTlong rows,COPTlong cols,
	const COPTlong* colptr,const COPTlong* rowind,const double* vals,void **Symbolic,
	const double* Control,double* Info)
{
	return umfpack_dl_symbolic(rows,cols,colptr,rowind,vals,Symbolic,Control,Info);
}

inline COPTlong umfpack_symbolic(
	longsize rows,longsize cols,
	const longsize* colptr,const longsize* rowind,const double*vals,void **Symbolic,
	const double* Control , double *Info)
{
	return umfpack_dl_symbolic((COPTlong)rows,(COPTlong)cols,(const COPTlong*)colptr,(const COPTlong*)rowind,vals,Symbolic,Control,Info);
}

// inline int umfpack_symbolic(
// 	int rows,int cols,
// 	const COPTlong* colptr,const COPTlong* rowind,const std::complex<double>* vals,void **Symbolic,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_symbolic(rows,cols,colptr,rowind,vals,0,Symbolic,Control,Info);
// }

inline int umfpack_numeric(
	const int* colptr, const int* rowind, const double*vals,
	void *Symbolic, void **Numeric,
	const double* Control, double* Info)
{
	return umfpack_di_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
}

inline COPTlong umfpack_numeric(
	const COPTlong *colptr, const COPTlong* rowind, const double* vals,
	void *Symbolic, void **Numeric,
	const double* Control,double *Info)
{
	return umfpack_dl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
}

inline COPTlong umfpack_numeric(
	const longsize* colptr , const longsize* rowind , const double* vals,
	void *Symbolic, void **Numeric,
	const double* Control,double *Info)
{
	return umfpack_dl_numeric((const COPTlong*)colptr,(const COPTlong*)rowind,vals,Symbolic,Numeric,Control,Info);
}

// inline int umfpack_numeric(
// 	const COPTlong* colptr , const COPTlong* rowind ,const std::complex<double>* vals,
// 	void *Symbolic, void **Numeric,
// 	const double* Control,double *Info)
// {
// 	return umfpack_zl_numeric(colptr,rowind,vals,Symbolic,Numeric,Control,Info);
// }

inline int umfpack_solve(
	int sys,const int* colptr,const int* rowind, const double* vals,
	double* X, const double* B,void *Numeric,
	const double* Control, double* Info)
{
	return umfpack_di_solve(sys,colptr,rowind,vals,X,B,Numeric,Control,Info);
}
inline COPTlong umfpack_solve(
	int sys,const COPTlong *colptr,const COPTlong* rowind, const double* vals,
	double* X, const double* B, void *Numeric,
	const double *Control,double *Info)
{
	return umfpack_dl_solve(sys,colptr,rowind,vals,X,B,Numeric,Control,Info);
}

inline COPTlong umfpack_solve(
	int sys,const longsize* colptr, const longsize* rowind, const double* vals,
	double* X, const double* B, void *Numeric,
	const double* Control, double* Info)
{
	return umfpack_dl_solve(sys,(const COPTlong*)colptr,(const COPTlong*)rowind,vals,X,B,Numeric,Control,Info);
}

// inline int umfpack_solve(
// 	int sys,const COPTlong *colptr,const COPTlong* rowind, const double*vals,
// 	std::complex<double>* X,const std::complex<double>* B,void *Numeric,
// 	const double*Control,double *Info)
// {
// 	return umfpack_zl_solve(sys,colptr,rowind,vals,0,X,0,B,0,Numeric,Control,Info);
// }

/** default umfpack control */
inline void umfpack_defaults(
	double Control[UMFPACK_CONTROL], double, int)
{
	umfpack_di_defaults(Control);
}

inline void umfpack_defaults(
	double Control[UMFPACK_CONTROL],double, COPTlong)
{
	umfpack_dl_defaults(Control);
}
	

//%}

/*			A wrapper for UMFPack solving sparse linear system.
 *			The input of the method is a COPT sparse matrix type.
 *			The wrapper first analyze the input matrix and factorize it.
 *			Warning and error information is provided
 */
template<class SpMatrix>
class UMFLinearSolver
	:
	noncopyable
{
private:
	typedef		typename SpMatrix::scalar			scalar;
	typedef 	typename SpMatrix::index 			index;
	typedef VectorBase<scalar,index,Dynamic>	 	DVector;

	/**	private variables */
	//%{

	/** the sparse matrix */
	const SpMatrix& 	__mat;

	/** umfpack numeric */
	void*				__symbolic;

	/** umfpack numeric */
	void*				__numeric;

	/** UMFPack Control */
	double* 			__control;

	/** UMFPack Info */
	double*				__info;

	/** the last computation time */
	double 				__compute_time;

	/** the last solving time */
	double 				__solve_time;

	/** whether warning happens */
	bool				__haswarning;

	/** whether error happens */
	bool				__haserror;

	/** whether the matrix is square */
	bool 				__issquare;

	/** status information */
	std::string 		__info_str;

	//%}

	UMFLinearSolver();

	void analyzeStatus(const int status);
public:

	UMFLinearSolver( const SpMatrix& mat);
	~UMFLinearSolver();

	/** initialization of umfpack solver */
	void init();

	/** given the matrix first analyze the symbolic */
	void analyzeSymbolic();

	/** analyze numeric */
	void analyzeNumeric();

	/** solve a linear system */
	template<class Vec>
	DVector solve(const Vec& vec);

	/** print information */
	void printInfo();

	/** compute time in seconds */
	double computeTime();

	/** solve time in seconds */
	double solveTime();


};

/** 		Implementation			*/

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeStatus( const int status )
{
	if ( status == 0 )
	{
		// everythign is ok
		return;
	}
	else if ( status > 0 )
	{
		// there is some warning
		__haswarning = true;
		switch (status)
		{
			case 1:
			{
				__info_str.append("UMFPack Warning: matrix is singular. There are exact zeros on the diagonal of U.\n");
			}
			break;
			case 2:
			{
				__info_str.append("UMFPack Warning: the determinant is nonzero, but smaller in magnitude than the smallest positive floating-point number.\n");
			}
			break;
			case 3:
			{
				__info_str.append("UMFPack Warning: the determinant is larger in magnitude than the largest positive floating-point number (IEEE Inf).\n");
			}
			break;
			default:
			{
				__info_str.append("UMFPack Warning: unknown warning happes!\n");
			}
			break;
		}
	}
	else
	{
		// error happens
		__haserror = true;
		switch (status)
		{
			case -1:
			{
				__info_str.append("UMFPack Error: out of memory!\n");
			}
			break;
			case -3:
			{
				__info_str.append("UMFPack Error: Numeric object is invalid. You should check umfpack user guide for further infomation!\n");
			}
			break;
			case -4:
			{
				__info_str.append("UMFPack Error: Symbolic object is invalid. You should check umfpack user guide for further infomation!\n");
			}
			break;
			case -5:
			{
				__info_str.append("UMFPack Error: A NULL pointer is passed while it need to be present. Check that whether Symbolic or Numeric is NULL or not!\n");
			}
			break;
			case -6:
			{
				__info_str.append("UMFPack Error: the number of rows and columns should be greater than zero!\n");
			}
			break;
			case -8:
			{
				__info_str.append("UMFPack Error: the matrix is invalid!\n");
			}
			break;
			case -11:
			{
				__info_str.append("UMFPack Error: different pattern now. Check that whether that you change the matrix between symbolic and numeric factorization.\n");
			}
			break;
			case -13:
			{
				__info_str.append("UMFPack Error: sys system argument is invalid!\n");
			}
			break;
			case -15:
			{
				__info_str.append("UMFPack Error: the provided permutation vector is invalid!\n");
			}
			break;
			case -17:
			{
				__info_str.append("UMFPack Error: file I/O error happens when saving or loading Numeric or Symbolic object!\n");
			}
			break;
			case -18:
			{
				__info_str.append("UMFPack Error: the ordering method failed!\n");
			}
			break;
			case -911:
			{
				__info_str.append("UMFPack Error: an internal error has occured!\n");
			}
			break;
			default:
			{
				__info_str.append("UMFPack Error: unknown error occured!\n");
			}
			break;
		}
	}
}

template<class SpMatrix>
UMFLinearSolver<SpMatrix>::UMFLinearSolver(
	const SpMatrix& mat)
	:
	__mat(mat),
	__symbolic(NULL),
	__numeric(NULL),
	__control(NULL),
	__info(NULL),
	__compute_time(0.0),
	__solve_time(0.0),
	__haswarning(false),
	__haserror(false)
{
	init();

	clock_t compute_start = clock() , compute_end;
	analyzeSymbolic();
	analyzeNumeric();
	compute_end = clock();

	__compute_time = static_cast<double>(compute_end-compute_start)/CLOCKS_PER_SEC;

	if(__haswarning||__haserror)
		printInfo();
	else
		__info_str.append("factorization success!");
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::init()
{
	if (__mat.rows()==__mat.cols())
		__issquare = true;
	else
		__issquare = false;
	__control = new double[UMFPACK_CONTROL];
	__info = new double[UMFPACK_INFO];
	umfpack_defaults(__control,scalar(),index());
}

template<class SpMatrix>
UMFLinearSolver<SpMatrix>::~UMFLinearSolver()
{
	umfpack_free_symbolic(&__symbolic,scalar(),index());
	umfpack_free_numeric(&__numeric,scalar(),index());
	SAFE_DELETE_ARRAY(__info);
	SAFE_DELETE_ARRAY(__control);
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeSymbolic()
{
	int status = umfpack_symbolic(__mat.rows(),__mat.cols(),__mat.columnPointer(),__mat.rowIndex(),__mat.values(),&__symbolic,__control,__info);
	analyzeStatus(status);
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::analyzeNumeric()
{
	int status = umfpack_numeric(__mat.columnPointer(),__mat.rowIndex(),__mat.values(),__symbolic,&__numeric,__control,__info);
	analyzeStatus(status);
}

template<class SpMatrix>
template<class Vec>
typename UMFLinearSolver<SpMatrix>::DVector UMFLinearSolver<SpMatrix>::solve(const Vec &vec)
{
	DVector result(vec.size());
	clock_t solve_start = clock() , solve_end;
	umfpack_solve(UMFPACK_A,__mat.columnPointer(),__mat.rowIndex(),__mat.values(),result.dataPtr(),vec.dataPtr(),__numeric,__control,__info);
	solve_end = clock();
	__solve_time = static_cast<double>(solve_end-solve_start)/CLOCKS_PER_SEC;
	return result;
}

template<class SpMatrix>
void UMFLinearSolver<SpMatrix>::printInfo()
{
	std::cerr<<__info_str<<std::endl;
	for ( int i = 0 ; i < UMFPACK_INFO ; ++ i )
	{
		std::cout<<i<<" "<<__info[i]<<std::endl;
	}
}

template<class SpMatrix>
double UMFLinearSolver<SpMatrix>::computeTime()
{
	return __compute_time;
}

template<class SpMatrix>
double UMFLinearSolver<SpMatrix>::solveTime()
{
	return __solve_time;
}

}// End of namespace COPT
#endif

#endif