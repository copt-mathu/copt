/** 		In this example file, we use conjugate gradient method used in 
  * 		solving linear system to introduce how our optimization solver
  *			is used. 
  *			The objective function is:
  * 				min_x \phi(x) = 1/2*x^TAx - b^Tx
  * 
  */


#include "Core"

typedef COPT::KernelTrait<double> 				kernel;
typedef kernel::Vector 							Vector;
typedef kernel::Matrix 							Matrix;
typedef COPT::BasicParameter<double> 			BasicParameter;
typedef COPT::BasicOption<double> 				BasicOption;

/** 	Parameter of the solver 
  * 		Since we have a specific parameter 'A' and 'b' in the problem, we inherit
  * 		basic parameter introduced in COPT and put 'A' and 'b' in it.
  */
struct Parameter:public BasicParameter
{
 	/** A matrix used in the problem */
 	Matrix 		A;
 	/** b vector */
 	Vector 		b;
 	/** current residual */
 	Vector 		r;
 	/** the current direction p */
 	Vector 		p;
};

/** 	the objective function: 1/2*x^TAx - b^Tx */
double phi(const Vector& x,const Parameter& p)
{
	return 0.5*x.dot(p.A*x)-p.b.dot(x);
}

/** 	initialization of the solver */
void initialization(Vector& x, Parameter& para)
{
	x.resize(para.A.rows());
	para.r = para.A*x-para.b;
	para.p = -para.r;
}

/** 	one iteration of conjugate gradient method */
double iteration(Vector& x,Parameter& para)
{
	double alpha = para.r.dot(para.r)/(para.p.dot(para.A*para.p));
	x = x+alpha*para.p;
	Vector r_n = para.r + alpha*para.A*para.p;
	double beta = (r_n.dot(r_n))/(para.r.dot(para.r));
	para.r = r_n;
	para.p = -para.r+beta*para.p;
	return (para.r.norm());
}

/** 	determination of the terminal condition */
bool terminal(const Vector& x,const BasicOption& o,const Parameter& para)
{
	if(para.r.norm()<=1e-10)
		return true;
	else
		return false;
}

typedef COPT::Solver<double,Vector,Vector,BasicOption,Parameter> 			Solver;
int main(int argc, char *argv[])
{
	Parameter para;
	Matrix A = Matrix::identity(5,5);
	A(0,1) = 1.0;A(0,2)=1.0;A(0,3)=1.0;A(0,4)=1.0;
	A(1,0) = -1.0;
	para.A = A.transpose()*A;
	para.b = Vector(5);
	para.b[0]=1.0;
	para.b[3]=3.0;
	Solver sol(para,&phi,&initialization,&iteration,&terminal);
	sol.solve();
	std::cout<<sol.result()<<std::endl;
	auto f = sol.objectiveFunction();
	std::cout<<f(sol.result(),sol.parameter())<<std::endl;
}