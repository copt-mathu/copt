#ifndef SIZE_HPP__
#define SIZE_HPP__

namespace COPT
{
template<class T>
class MSize
{
private:
	T __m;
	T __n;
public:
	MSize(){throw COException("Error: unknown size type! Check your template carefully!");}
	MSize(T m,T n){throw COException("Error: unknown size type! Check your template carefully!");}
};

template<>
class MSize<int>
{
private:
	int __m;
	int __n;
public:
	MSize():__m(0),__n(0){}
	MSize(int m,int n):__m(m),__n(n){}
	template<class Matrix>
	MSize(const Matrix& mat):__m(mat.rows()),__n(mat.cols()){}

	int m() const {return __m;}
	int n() const {return __n;}

	MSize& operator=(const MSize& s){__m=s.m();__n=s.n();return *this;}
	MSize& operator=(MSize&& s){__m=s.m();__n=s.n();return *this;}
};

template<>
class MSize<long>
{
private:
	long __m;
	long __n;
public:
	MSize():__m(0),__n(0){}
	MSize(long m,long n):__m(m),__n(n){}
	template<class Matrix>
	MSize(const Matrix& mat):__m(mat.rows()),__n(mat.cols()){}

	long m() const {return __m;}
	long n() const {return __n;}

	MSize& operator=(const MSize& s){__m=s.m();__n=s.n();return *this;}
	MSize& operator=(MSize&& s){__m=s.m();__n=s.n();return *this;}
};
}

#endif