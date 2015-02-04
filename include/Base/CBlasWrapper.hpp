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


#ifndef CBLAS_WRAPPER_HPP__
#define CBLAS_WRAPPER_HPP__

// 		This file defines C++ wrapper for library cblas

#ifdef CBLAS
namespace blas{
//		level 1
template<class eT>
inline void copt_blas_copy( const int N,const eT* X,const int incX, eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		Y[i*incY]=X[i*incX];
	}
}

inline void copt_blas_copy( const int N , const double* X, const int incX , double* Y , const int incY)
{
	cblas_dcopy(N,X,incX,Y,incY);
}

inline void copt_blas_copy( const int N, const float* X , const int incX , float* Y , const int incY )
{
	cblas_scopy(N,X,incX,Y,incY);
}

inline void copt_blas_copy(const int N, const std::complex<float>* X,const int incX, std::complex<float>* Y,const int incY)
{
	cblas_ccopy(N,X,incX,Y,incY);
}

inline void copt_blas_copy(const int N, const std::complex<double>* X,const int incX,std::complex<double>* Y,const int incY)
{
	cblas_zcopy(N,X,incX,Y,incY);
}

/*			overload of swap operation
 *
 *
 */
template<class eT>
inline void copt_blas_swap (const int N,eT* X,const int incX,eT* Y,const int incY)
{
}
inline void copt_blas_swap (const int N,float* X,const int incX,float* Y,const int incY)
{
	cblas_sswap(N,X,incX,Y,incY);
}
inline void copt_blas_swap (const int N,double* X,const int incX,double* Y,const int incY)
{
	cblas_dswap(N,X,incX,Y,incY);
}
inline void copt_blas_swap (const int N,std::complex<float>* X,const int incX,std::complex<float>* Y,const int incY)
{
	cblas_cswap(N,X,incX,Y,incY);
}
inline void copt_blas_swap (const int N,std::complex<double>* X,const int incX,std::complex<double>* Y,const int incY)
{
	cblas_zswap(N,X,incX,Y,incY);
}

/** 			summation 			*/

inline void copt_blas_axpy( const int N , const float alpha , const float* X , const int incX , float *Y , const int incY)
{
	cblas_saxpy(N,alpha,X,incX,Y,incY);
}
inline void copt_blas_axpy( const int N , const double alpha , const double* X , const int incX , double *Y , const int incY )
{
	cblas_daxpy(N,alpha,X,incX,Y,incY);
}
inline void copt_blas_axpy( const int N , const std::complex<float>& alpha , const std::complex<float> *X , const int incX , std::complex<float> *Y , const int incY)
{
	cblas_caxpy(N,&alpha,X,incX,Y,incY);
}
inline void copt_blas_axpy( const int N , const std::complex<double>& alpha , const std::complex<double> *X , const int incX , std::complex<double> *Y , const int incY )
{
	cblas_zaxpy(N,&alpha,X,incX,Y,incY);
}

/**				dot operation 			*/
template<class eT>
inline eT copt_blas_dot(const int N,const eT* X,const int incX,const eT* Y,const int incY)
{
	return static_cast<eT>(0.0);
}
inline double copt_blas_dot(const int N,const double* X,const int incX,const double* Y,const int incY)
{
	return cblas_ddot(N,X,incX,Y,incY);
}
inline float copt_blas_dot(const int N,const float* X,const int incX,const float* Y,const int incY)
{
	return cblas_sdot(N,X,incX,Y,incY);
}
inline std::complex<float> copt_blas_dot(const int N,const std::complex<float>* X,const int incX,const std::complex<float>* Y,int incY)
{
	std::complex<float> re;
	cblas_cdotu_sub(N,X,incX,Y,incY,&re);
	return re;
}
inline std::complex<double> copt_blas_dot(const int N,const std::complex<double>* X,const int incX,const std::complex<double>* Y,int incY)
{
	std::complex<double> re;
	cblas_zdotu_sub(N,X,incX,Y,incY,&re);
	return re;
}

/*				scale operation
 */
template<class eT>
inline void copt_blas_scal(const int N,const eT alpha,eT* X,const int incX)
{
}
inline void copt_blas_scal(const int N,const float alpha,float* X,const int incX)
{
	cblas_sscal(N,alpha,X,incX);
}
inline void copt_blas_scal(const int N,const double alpha,double* X,const int incX)
{
	cblas_dscal(N,alpha,X,incX);
}
inline void copt_blas_scal(const int N,const std::complex<float>& alpha,std::complex<float>* X,const int incX)
{
	cblas_cscal(N,&alpha,X,incX);
}
inline void copt_blas_scal(const int N,const std::complex<double>& alpha,std::complex<double>* X,const int incX)
{
	cblas_zscal(N,&alpha,X,incX);
}
inline void copt_blas_scal(const int N,const float alpha,std::complex<float>* X,const int incX)
{
	cblas_csscal(N,alpha,X,incX);
}
inline void copt_blas_scal(const int N,const double alpha,std::complex<double>* X,const int incX)
{
	cblas_zdscal(N,alpha,X,incX);
}

template<class T>
inline T copt_blas_nrm2(const int N , const T* X , const int incX )
{
	throw COPT::COException("Unsupported type for nrm2 computation!");
}
inline float copt_blas_nrm2(const int N , const float *X , const int incX )
{
	return cblas_snrm2(N,X,incX);
}
inline double copt_blas_nrm2(const int N , const double *X , const int incX )
{
	return cblas_dnrm2(N,X,incX);
}
inline float copt_blas_nrm2(const int N , const std::complex<float> *X , const int incX )
{
	return cblas_scnrm2(N,X,incX);
}
inline double copt_blas_nrm2(const int N , const std::complex<double> *X , const int incX )
{
	return cblas_dznrm2(N,X,incX);
}

/** level 2 operations */

inline void copt_blas_gemv(const enum CBLAS_ORDER order,const enum CBLAS_TRANSPOSE TransA , const int M , const int N, const float alpha, const float *A , const int lda , const float*X , const int incX , const float beta , float *Y , const int incY)
{
	cblas_sgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}

inline void copt_blas_gemv(const enum CBLAS_ORDER order, const enum CBLAS_TRANSPOSE TransA , const int M , const int N,const double alpha , const double* A , const int lda , const double*X , const int incX , const double beta , double *Y , const int incY )
{
	cblas_dgemv(order,TransA,M,N,alpha,A,lda,X,incX,beta,Y,incY);
}

inline void copt_blas_gemv(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA , const int M , const int N , const std::complex<float>& alpha , const std::complex<float>* A , const int lda , const std::complex<float>* X , const int incX , const std::complex<float>& beta , std::complex<float>* Y , const int incY )
{
	cblas_cgemv(order,TransA,M,N,&alpha,A,lda,X,incX,&beta,Y,incY);
}

inline void copt_blas_gemv(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA , const int M , const int N , const std::complex<double>& alpha , const std::complex<double>* A, const int lda , const std::complex<double>* X , const int incX , const std::complex<double>& beta , std::complex<double>* Y , const int incY )
{
	cblas_zgemv(order,TransA,M,N,&alpha,A,lda,X,incX,&beta,Y,incY);
}

/** symmetric multiplication */

template<class index,class scalar>
inline void copt_blas_symv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
					const index N, const scalar alpha, const scalar *A,
					const index lda, const scalar *X , const index incX,
					const scalar beta, scalar *Y , const index incY )
{
	throw COPT::COException("Unknown type for blas wrapper!");
}

template<>
inline void copt_blas_symv(const enum CBLAS_ORDER order , const enum CBLAS_UPLO Uplo,
					const int N , const float alpha , const float *A ,
					const int lda , const float *X , const int incX,
					const float beta , float *Y , const int incY)
{
	cblas_ssymv(order,Uplo,N,alpha,A,lda,X,incX,beta,Y,incY);
}

template<>
inline void copt_blas_symv(const enum CBLAS_ORDER order , const enum CBLAS_UPLO Uplo,
					const int N , const double alpha , const double *A , 
					const int lda , const double *X , const int incX,
					const double beta , double *Y , const int incY)
{
	cblas_dsymv(order,Uplo,N,alpha,A,lda,X,incX,beta,Y,incY);
}

inline void copt_blas_symv(const enum CBLAS_ORDER order , const enum CBLAS_UPLO Uplo,
					const int N , const std::complex<float>&alpha , const std::complex<float>*A,
					const int lda , const std::complex<float> *X , const int incX,
					const std::complex<float>& beta , std::complex<float> *Y , const int incY)
{
	return;
	// cblas_csymv(order,Uplo,N,&alpha,A,lda,X,incX,&beta,Y,incY);
}

inline void copt_blas_symv(const enum CBLAS_ORDER order , const enum CBLAS_UPLO Uplo,
					const int N , const std::complex<double>& alpha , const std::complex<double> *A,
					const int lda , const std::complex<double>* X , const int incX , 
					const std::complex<double>& beta , std::complex<double> *Y , const int incY)
{
	return;
	// cblas_zsymv(order,Uplo,N,&alpha,A,lda,X,incX,&beta,Y,incY);
}

/** level three operations */

inline void copt_blas_gemm(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA,
	const enum CBLAS_TRANSPOSE TransB , const int M , const int N , const int K , const float alpha , const float*A ,
	const int lda , const float* B , const int ldb,
	const float beta , float *C , const int ldc )
{
	cblas_sgemm(order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
}

inline void copt_blas_gemm(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA,
	const enum CBLAS_TRANSPOSE TransB , const int M , const int N , const int K , const double alpha , const double* A,
	const int lda , const double *B , const int ldb,
	const double beta , double *C , const int ldc )
{
	cblas_dgemm(order,TransA,TransB,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);
}

inline void copt_blas_gemm(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA,
					const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
					const int K, const std::complex<float>& alpha, const std::complex<float>* A,
					const int lda, const std::complex<float>* B, const int ldb,
					const std::complex<float>& beta, std::complex<float> *C , const int ldc )
{
	cblas_cgemm(order,TransA,TransB,M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc);
}

inline void copt_blas_gemm(const enum CBLAS_ORDER order , const enum CBLAS_TRANSPOSE TransA,
					const enum CBLAS_TRANSPOSE TransB, const int M, const int N,
					const int K, const std::complex<double>& alpha, const std::complex<double>* A,
					const int lda, const std::complex<double>* B, const int ldb,
					const std::complex<double>& beta, std::complex<double> *C , const int ldc )
{
	cblas_zgemm(order,TransA,TransB,M,N,K,&alpha,A,lda,B,ldb,&beta,C,ldc);
}

/** symm suit for s,d,c,z */
inline void copt_blas_symm(const CBLAS_ORDER Order , const enum CBLAS_SIDE Side,
					const CBLAS_UPLO Uplo , const int M , const int N,
					const float alpha , const float* A ,const int lda , 
					const float *B , const int ldb , const float beta,
					float *C , const int ldc )
{
	cblas_ssymm(Order,Side,Uplo,M,N,alpha,A,lda,B,ldb,beta,C,ldc);
}

inline void copt_blas_syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int n, const int k, const float alpha,const float*A , const int lda , const float beta , float* C, const int ldc)
{
	cblas_ssyrk(Order,Uplo,Trans,n,k,alpha,A,lda,beta,C,ldc);
}

inline void copt_blas_syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int n, const int k, const double alpha,const double* A, const int lda, const double beta , double* C, const int ldc )
{
	cblas_dsyrk(Order,Uplo,Trans,n,k,alpha,A,lda,beta,C,ldc);
}

inline void copt_blas_syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans , const int n, const int k , const std::complex<float>& alpha , const std::complex<float>* A , const int lda , const std::complex<float>& beta , std::complex<float>* C , const int ldc )
{
	cblas_csyrk(Order,Uplo,Trans,n,k,&alpha,A,lda,&beta,C,ldc);
}

inline void copt_blas_syrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const int n , const int k , const std::complex<double>& alpha , const std::complex<double>* A , const int lda , const std::complex<double>& beta , std::complex<double>* C , const int ldc )
{
	cblas_zsyrk(Order,Uplo,Trans,n,k,&alpha,A,lda,&beta,C,ldc);
}

inline void copt_blas_herk(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
				const enum CBLAS_TRANSPOSE Trans, const int n, const int k,
				const float alpha , const float *A , const int lda,
				const float beta, float *C , const int ldc )
{
	throw COPT::COException("Unknown type for blas wrapper!");
}

inline void copt_blas_herk(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
				const enum CBLAS_TRANSPOSE Trans, const int n, const int k,
				const double alpha , const double *A , const int lda,
				const double beta, double *C , const int ldc )
{
	throw COPT::COException("Unknown type for blas wrapper!");
}

inline void copt_blas_herk(const enum CBLAS_ORDER Order,const enum CBLAS_UPLO Uplo,
				const enum CBLAS_TRANSPOSE Trans, const int n, const int k,
				const float alpha , const std::complex<float> *A , const int lda,
				const float beta, std::complex<float> *C , const int ldc )
{
	cblas_cherk( Order, Uplo, Trans, n, k, alpha, A, lda, beta, C, ldc);
}

inline void copt_blas_herk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
				const enum CBLAS_TRANSPOSE Trans, const int n, const int k,
				const double alpha, const std::complex<double>* A, const int lda,
				const double beta, std::complex<double> *C , const int ldc)
{
	cblas_zherk(Order,Uplo,Trans,n,k,alpha,A,lda,beta,C,ldc);
}

/*				summation
*/

//			level 2
/*				multipliation
 */
// inline void copt_blas_gemm()


}// End of namespace blas

#else
/*					simple blas operations
 */
namespace blas{
template<class eT>
inline void copt_blas_copy( const int N,const eT* X,const int incX, eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		Y[i*incY] = X[i*incX];
	}
}
template<class eT>
inline void copt_blas_swap (const int N,eT* X,const int incX,eT* Y,const int incY)
{
	for ( int i = 0 ; i < N ; ++ i ){
		eT t = X[i*incX];
		X[i*incX] = Y[i*incY];
		Y[i*incY] = t;
	}
}
template<class eT>
eT copt_blas_dot(const int N,const eT* X,const int incX,const eT* Y,const int incY)
{
	eT t(0.0);
	for ( int i = 0 ; i < N ; ++ i ){
		t += X[i*incX]*Y[i*incY];
	}
	return t;
}
template<class eT>
inline void copt_blas_scal(const int N,const eT alpha,eT* X,const int incX)
{
	for (int i = 0 ; i < N ; ++ i ){
		X[i*incX] = X[i*incX]*alpha;
	}
}
}

// End of ifdef CBLAS
#endif


#endif