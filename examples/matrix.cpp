#include "Core"


using namespace COPT;

//Dynamic Matrix
typedef MatrixBase<double,int,Dynamic,Dynamic>		DMatrix;
typedef MatrixBase<double,int,Dynamic,3> 			MatrixD3;
typedef MatrixBase<double,int,5,1> 					Matrix51;

int main(int argc, char *argv[])
{
	DMatrix mat;
	// Matrix51 m(3,2);
	MatrixD3 md(2);

	std::cout<<"md is "<<std::endl<<md<<std::endl;
	std::cout<<"lda of md is "<<md.lda()<<std::endl;
	std::cout<<"the actual size is "<<md.size()<<std::endl;

	MatrixBase<double,int,2,Dynamic> m3(4);
	std::cout<<"lda of m3 is "<<m3.lda()<<std::endl;
	std::cout<<"size of m3 is "<<m3.rows()<<' '<<m3.cols()<<std::endl;
	m3(1,2) = 1.0;
	m3(1,3) = 3.0;
	std::cout<<m3<<std::endl;
	mat = m3;
	std::cout<<mat<<std::endl;

	std::cout<<5.0*mat<<std::endl;

	DMatrix iden = DMatrix::identity(4,4);
	DMatrix iden2 = DMatrix::identity(2,2);
	std::cout<<m3*iden<<std::endl;

	std::cout<<m3.lda()<<std::endl;
	std::cout<<m3.transpose()*iden2<<std::endl;
	std::cout<<m3.transMulti(iden2)<<std::endl;
	std::cout<<m3.transpose().transMulti(iden)<<std::endl;
	std::cout<<m3.transpose().transpose()*iden<<std::endl;

	MatrixBase<double,int,4,4> mtm;
	m3.mtm(mtm);
	std::cout<<mtm<<std::endl;

	std::cout<<COPT::mean(mtm.rowBegin(),mtm.rowEnd())<<std::endl;
	std::cout<<COPT::mean(mtm.colBegin(),mtm.colEnd())<<std::endl;
	std::cout<<"here"<<std::endl;
}