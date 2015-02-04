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


#ifndef MATRIX_BASE_HPP__
#define MATRIX_BASE_HPP__

/*
 *				class 'MatrixBase'
 *					dervied from base class 'Array'
 */
namespace COPT
{

/** 	Abstract matrix class for special consideration. the AbstractMatrix has no
 *  	actual functionality and any matrix that derives from AbstractMatrix is able
 *  	to be used in the COPT.
 */
template<class FT,class I>
class AbstractMatrix
	:
	DataObject<FT,I>
{
public:
	/** the category of matrix */
	typedef 		matrix_object 						ObjectCategory;
	/** Dynamic type */
	typedef 		AbstractMatrix 						DType;
	typedef FT 											scalar;
	typedef I 											index;

private:
	
public:
	virtual ~AbstractMatrix(){}
	/** the row and column number */
	virtual index rows() const = 0;
	virtual index cols() const = 0;
	/** the lda of matrix */
	virtual index lda() const = 0;
	virtual scalar* dataPtr() = 0;
	virtual const scalar *dataPtr() const = 0;

	/** the size */
	virtual index size() const = 0;
	virtual bool isSymmetric() const = 0;
	virtual bool isTranspose() const = 0;

	virtual int rowAtCompileTime() const = 0;
	virtual int colAtCompileTime() const = 0;

};
/*
 *	Class of 'MatrixBase'
 *		the data is stored column by column
 *		the matrix is assumed not to be transpose or symmetric
 *		but symmetric and transpose flag is designed
 */
template<class FT,class I = int,int RowAtCompileTime=Dynamic,int ColAtCompileTime=Dynamic>
class MatrixBase
	:
	public Array<FT,I,(RowAtCompileTime==Dynamic||ColAtCompileTime==Dynamic)?Dynamic:(RowAtCompileTime+1)*(ColAtCompileTime)>,
	public AbstractMatrix<FT,I> 
{
public:

	/** the scalar type */
	typedef 			FT 				 				scalar;
	/** pod scalar type */
	typedef typename get_pod_type<scalar>::type 		podscalar;
	/** the size type used */
	typedef 			I 								index;						
	/** define fthe category */
	typedef 			matrix_object					ObjectCategory;
	/** define the trait */
	typedef 			KernelTrait<FT,index>			Kernel;
	/** the dynamic object */
	typedef MatrixBase<scalar,index,Dynamic,Dynamic> 	DType;
	/** the dynamic matrix */
	typedef MatrixBase<FT,I,Dynamic,Dynamic> 			DMatrix;
	/** the dynamic vector */
	typedef VectorBase<scalar,index,Dynamic> 			DVector;
	/** the corresponding Abstract Matrix */


	typedef COPT::AbstractMatrix<FT,I> 						AbstractMatrix;

	typedef COPT::AbstractVector<FT,I>						AbstractVector;

private:
	
	/** definition used in implementation */
	typedef 			VectorBase<FT,index>			Vector;
	typedef 			COPT::Array<FT,I,(RowAtCompileTime==Dynamic||ColAtCompileTime==Dynamic)?Dynamic:(RowAtCompileTime+1)*(ColAtCompileTime)>			
														Array;


	/**			private variables			*/
	//%{
	/** the size of rows */
	index					__rows;
	/** the size of columns */
	index 					__cols;
	/** whether the matrix is symmetric */
	bool					__sym;
	/** whether the matrix is transpose */
	bool 					__trans;
	/** lda. in current version lda = __rows+1 */
	index 					__lda;
	/** is row number dynamic */
	bool 					__is_row_dynamic;
	/** is column number dynamic */
	bool 					__is_col_dynamic;
	//%}


	/** private functions */
	template<class Mat>
	DMatrix multi(const Mat &mat, const matrix_object&) const;

	template<class Vec>
	DVector multi(const Vec &vec,const vector_object&) const;

public:
	/** struct for vector iterator */
	class RowVectorIterator
		:
		public std::iterator<std::input_iterator_tag,DVector>
	{
		
		MatrixBase& __m;
		index __i;
		DVector* __v;

	public:
		RowVectorIterator(const index i,MatrixBase& mat)
			:
			__m(mat),
			__i(i),
			__v(nullptr)
		{
			if(__i<0||__i>=__m.rows())
				__v = nullptr;
			else
				__v = new DVector(__m.row(i));
		}
		RowVectorIterator(const RowVectorIterator& iter)
			:
			__m(iter.__m),
			__i(iter.__i),
			__v(nullptr)
		{
			if(__i<0||__i>=__m.rows())
				__v = nullptr;
			else
				__v = new DVector(__m.row(__i));
		}
		~RowVectorIterator()
		{
			SAFE_DELETE(__v);
			__v = nullptr;
		}
		RowVectorIterator& operator++() 
		{
			if(__i>=0&&__i<__m.rows())
			{
				SAFE_DELETE(__v);
				++__i;
				if(__i<__m.rows())
					__v = new DVector(__m.row(__i));
				else
					__v = nullptr;
			}
			else
			{
				SAFE_DELETE(__v);
				__v = nullptr;
			}
			return (*this);
		}
		RowVectorIterator operator++(int)
		{
			RowVectorIterator tmp (*this);
			if(__i>=0&&__i<__m.rows())
			{
				SAFE_DELETE(__v);
				++__i;
				if(__i<__m.rows())
					__v = new DVector(__m.row(__i));
				else
					__v = nullptr;
			}
			else
			{
				SAFE_DELETE(__v);
				__v = nullptr;
			}
			return tmp;
		}
		bool operator==(const RowVectorIterator& rhs){
			return (__i==rhs.__i&& &__m==&rhs.__m);
		}
		bool operator!=(const RowVectorIterator& rhs){
			return (!operator==(rhs));
		}
		DVector& operator*(){return *__v;}
	};

	class ColVectorIterator
		:
		public std::iterator<std::input_iterator_tag,DVector>
	{
		
		MatrixBase& __m;
		index __i;
		DVector* __v;

	public:
		ColVectorIterator(const index i,MatrixBase& mat)
			:
			__m(mat),
			__i(i),
			__v(nullptr)
		{
			if(__i<0||__i>=__m.cols())
				__v = nullptr;
			else
				__v = new DVector(__m.col(i));
		}
		ColVectorIterator(const ColVectorIterator& iter)
			:
			__m(iter.__m),
			__i(iter.__i),
			__v(nullptr)
		{
			if(__i<0||__i>=__m.cols())
				__v = nullptr;
			else
				__v = new DVector(__m.col(__i));
		}
		~ColVectorIterator()
		{

			SAFE_DELETE(__v);
			__v = nullptr;
		}
		ColVectorIterator& operator++() 
		{
			if(__i>=0&&__i<__m.cols())
			{
				SAFE_DELETE(__v);
				++__i;
				if(__i<__m.cols())
					__v = new DVector(__m.col(__i));
				else
					__v = nullptr;
			}
			else
			{
				SAFE_DELETE(__v);
				__v = nullptr;
			}
			return (*this);
		}
		ColVectorIterator operator++(int)
		{
			ColVectorIterator tmp (*this);
			if(__i>=0&&__i<__m.cols())
			{
				SAFE_DELETE(__v);
				++__i;
				if(__i<__m.cols())
					__v = new DVector(__m.col(__i));
				else
					__v = nullptr;
			}
			else
			{
				SAFE_DELETE(__v);
				__v = nullptr;
			}
			return tmp;
		}
		bool operator==(const ColVectorIterator& rhs){
			return (__i==rhs.__i&& &__m==&rhs.__m);
		}
		bool operator!=(const ColVectorIterator& rhs){
			return (!operator==(rhs));
		}
		DVector& operator*(){return *__v;}
	};
	/** constructor and deconstructor */
	//%{
	/** default constructor */
	MatrixBase();
	/** construct a matrix with input dimension m*n */
	MatrixBase(const index m, const index n, const scalar *data=nullptr,bool trans = false);
	/** constructor for Matrix with row number of column number specified */
	MatrixBase(const index m, const scalar *data=nullptr);
	/** Copy assignment */
	MatrixBase(const MatrixBase& mat);
	/** Copy assignment */
	MatrixBase(const AbstractMatrix& mat);
	/** deconstructor */
	~MatrixBase();
	//%} end of constructor and deconstructor

	/**	getters and setters*/
	//%{
	/** get the number of rows */
	index		rows() const;
	/** is the row number dynamic */
	bool isRowDynamic() const;
	/** get the number of columns */
	index 		cols() const;
	bool isColumnDynamic() const;
	/** lda of the matrix */
	index lda() const;
	/** overload of dataPtr() */
	scalar* dataPtr();
	const scalar* dataPtr() const;
	/** overload of size */
	index size() const;
	/** overload of rowAtCompileTime */
	int rowAtCompileTime() const;
	/** overload of colAtCompileTime */
	int colAtCompileTime() const;

	/**	matlab-like (i,j) element getter */
	scalar& operator() (const index i, const index j);
	const scalar& operator() (const index i, const index j) const;

	/** get the element using Arr */
	const scalar& data(const index i) const;

	/**	obtain the i-th row */
	Vector row(const index num);
	const Vector row(const index num) const;
	/**	obtain the i-th column */
	Vector col(const index num);
	const Vector col(const index num) const;

	/** the first iterator */
	RowVectorIterator rowBegin(){return RowVectorIterator(0,*this);}
	RowVectorIterator rowEnd(){return RowVectorIterator(__rows,*this);}
	ColVectorIterator colBegin(){return ColVectorIterator(0,*this);}
	ColVectorIterator colEnd(){return ColVectorIterator(__cols,*this);}


	/** set element using Arr */
	void set (const index i, const scalar value);

	/** resize the matrix */
	void resize (index m, index n);
	/** resize for one dimension fixed matrix */
	void resize (index m);
	/** resize the matrix */
	void resize (const MSize<index>& s);

	/** set the matrix to be symmetric */
	void setSymmetricFlag(bool sym);
	bool isSymmetric() const;

	/** set the matrix to be transpose mode */
	void setTransposeFlag(bool sym);
	bool isTranspose() const;

	/** basic check for two matrices */
	template<class Mat>
	void binaryCheck(const Mat& mat) const;

	/** Copy operation */
	template<class Mat>
	MatrixBase& operator= (const Mat& mat);

	/** Mathematical operations */
	/** summation */
	template<class Mat>
	DMatrix operator+ (const Mat& mat) const;

	/** subtraction */
	template<class Mat>
	DMatrix operator- (const Mat& mat) const;

	/** matrix multiplications */
	template<class T>
	typename T::DType operator* (const T &vec) const;

	/** multiplication with a scalar */
	DMatrix operator* (const scalar s) const;

	/** multiplication between a scalar and a matrix */
	friend DMatrix operator* (const scalar s, const MatrixBase& mat)
	{
		DMatrix result(mat.rows(),mat.cols());
		for ( index i = 0 ; i < mat.rows() ; ++ i )
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				result(i,j) = mat(i,j)*s;
		return result;
	}

	/** transpose */
	DMatrix transpose() const;

	/** transpose muliplication */
	template<class T>
	typename T::DType transMulti(const T &t) const;

	template<class Vec>
	DVector transMulti(const Vec& vec, const vector_object&) const;

	template<class Mat>
	DMatrix transMulti(const Mat &mat, const matrix_object&) const;

	/*				overload of ostream
	 */
	friend std::ostream& operator<<(std::ostream& os, const MatrixBase& mat)
	{
		for ( index i = 0 ; i < mat.rows() ; ++ i ){
			for ( index j = 0 ; j < mat.cols() ; ++ j )
				os<<mat(i,j)<<' ';
			os<<std::endl;
		}
		return os;
	}
	
	template<class Vec>
	typename Vec::DType solve(const Vec& vec);

	Vector leastSquareSolve(const Vector& vec){
		if( __cols != vec.size() )
			throw COException("least square error: the size is not consistent!");
		scalar *a = new scalar[__rows*__cols];
		scalar *b = new scalar[vec.size()];
		scalar *s = new scalar[std::min(__rows,__cols)];
		index info = -1;
		index rank;
		blas::copt_blas_copy(__rows*__cols,this->dataPtr(),1,a,1);
		blas::copt_blas_copy(vec.size(),vec.dataPtr(),1,b,1);
		copt_lapack_gelss(__rows,__cols,1,a,__rows,b,vec.size(),s,-1.0,&rank,&info);
		Vector result(__cols,b);
		delete[]a;
		delete[]b;
		return result;
	}
	


	/*			Special MatrixBase
	 *
	 */
	static DMatrix identity(index m, index n){
		DMatrix result(m,n);
		index min = std::min(m,n);
		for ( index i = 0 ; i < min ; ++ i )
			result(i,i) = static_cast<scalar>(1.0);
		return result;
	}

	static DMatrix identity(index m, index n, const scalar s);

	template<class Vec>
	static DMatrix diag(index m, index n, const Vec& vec);

	/** Blocking methods */
	//%{
	
	/** Blocking matrix from a given matrix with specific row numbers and column numbers*/
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just columns */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& colnums);
	
	/** Blocking matrix from a given matrix with just rows */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::set<index>& rownums);

	/** Blocking matrix from a given matrix, order is not considered */
	void blockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums,
		const std::vector<index>& colnums);

	/** blocking matrix using iterator */
	template<class InputIterator>
	void blockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend,
		const InputIterator& colbegin,
		const InputIterator& colend );

	/** Blocking matrix from a given matrix, order is not considered */
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& colnums);

	/** column blocking from a given matrix */
	template<class InputIterator>
	void columnBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& colbegin,
		const InputIterator& colend);

	/** Blocking matrix from a given matrix, order is not considered */
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const std::vector<index>& rownums);

	template<class InputIterator>
	void rowBlockFromMatrix(
		const MatrixBase& mat,
		const InputIterator& rowbegin,
		const InputIterator& rowend);
	//%}

	/** combining methods */
	//%{

	/** combining along row direction */
	void combineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2);
	/** combining along column direction */
	void combineAlongColumn(
		const MatrixBase& m1,
		const MatrixBase& m2);

	/** combine along row direction taking matrix as parameter */
	static inline void stCombineAlongRow(
		const MatrixBase& m1,
		const MatrixBase& m2,
		MatrixBase& m);
	/** combine along column direction taking matrix as parameter */
	static inline void stCombineAlongColumn(
		const MatrixBase&m1,
		const MatrixBase& m2,
		MatrixBase& m);
	//%}

	/** set a random matrix */
	void setRandom(const index rows=RowAtCompileTime, const index cols=ColAtCompileTime);
	static inline DMatrix random(const index rows, const index cols);
	/** compute A^TA of a given matrix */
	template<class Mat>
	void mtm(Mat &mat) const;

	/** norms */
	//%{
	/** compute the operation norm */
	podscalar operationNorm() const;
	/** frobenius norm */
	podscalar frobeniusNorm() const;
	/** max element, actually not a norm */
	podscalar maxNorm() const;
	/** l_1 norm, maximum column sum */
	podscalar oneNorm() const;
	/** infinity norm, maximum row sum */
	podscalar infinityNorm() const;
	//%}

};// End of class MatrixBase




}// End of namespace COPT
#endif