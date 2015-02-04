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


#ifndef SOLVER_HPP__
#define SOLVER_HPP__

/**		"""
  *		A procedural-based framework for optimization solver.
  * 	Since Optimization is currently a foundamental and widely-used mathematical tool 
  * 	in many academic and industrial applications, it is an important problem to develop
  * 	and validate optimization solvers in a short time. Thus it is one of our core interests to 
  * 	propose an easy-to-validate framework for researchers to test their own algorithms.
  * 	In this header file, we introduce an easy-to-use procedural-based solver system.
  * 	Users can write their own functions that compute the objective function, how the solver
  * 	iterates, when the solver terminates and so on. Then then can easily separate an
  * 	optimization problem into several module which is a nice property of optimization problem.
  * 	This design helps researchers focus on the problem not how to code and debug. 
  * 	"""
  */
namespace COPT{

/** 		Description for objective function. 
  *
  * 		Input:
  * 			const ArgType& 		x:		Input argument of the objective function.
  * 			const Parameter& 	para: 	Input parameter which contains necessary
  *			parameters in an optimization solver.
  *
  * 		Simple use case:
  * 			consider a simple objective function that computes the Frobenius
  * 		norm of a given matrix A and the input argument x.
  *
  * 		struct Parameter:public BasicParameter
  * 		{
  *				Matrix A;
  *			};
  *			double func(const Matrix& x,const Parameter&p){
  *				return (p.A*mat).frobeniusNorm();
  * 		}
  *			ObjectiveFunction f = func;
  *			Matrix m = Matrix::identity(4,4);
  * 		Parameter p;
  * 		p.A = Matrix::identity(4.4);
  *			std::cout<<f(m,para); // 2.0;
  */
template<class Scalar,class ArgType,class Para>
using ObjectiveFunction = Scalar(*)(const ArgType&,const Para&);

/** 		Iteration function for one iteration.
  *			Input and Output:
  * 			ArgType& 	x: 		Main argument that is used and updated in the 
  * 		iteration.
  * 			Parameter& 	para:	Parameter in an optimization solver.
  */
template<class Scalar,class ArgType,class Para>
using OneIteration = Scalar(*)(ArgType& x,Para& para);

/** 		Initialization function for a solver. You can set the initial parameters,
  * 		prepare for solving the problem etc.
  * 		Input and Output:
  *				ArgType& 	x: 		Main argument that is initialized and will be used
  * 		in the whole procedure.
  * 			Para& 		para: 	Parameter that is initialized.
  */
template<class ArgType,class Parameter>
using SolverInitialization = void(*)(ArgType& x, Parameter& para);

/** 		Check whether the solver is terminated.
  * 		Input:
  * 			const ArgType& 		x: 		The argument that is used in iterations
  * 			const Option& 		opt:	The option of the problem like maximum iteration
  *			number, error threshold and so on.
  * 			const Parameter& 	para:	The parameters of the solver
  */
template<class ArgType,class Option,class Parameter>
using CheckTermination = bool(*)(const ArgType& x, const Option& opt, const Parameter& para);


/** 		Obtain the result from argument 
  * 		Input:
  * 			const ArgType& 	x:		The solved x.
  * 		Output:
  * 			OutputType& 	result: The result.
  *
  * 		Generally, you do not need to overload this function. But sometimes, one might add
  * 		redundent information in 'ArgType'. Then you will need to extract the result from x. 
  */
template<class ArgType,class OutputType>
using GetResult = void(*)(const ArgType& x, OutputType& result);

/** 		Normal print function for solver which can be modified.
  * 		Input: 
  * 			int 			level: 	The print level.
  * 		Input and Output
  * 			std::ostream 	os:		The output stream
  */
using NormalPrintFunction = void(*)(int level,std::ostream& os);

/** 		Print function for the solver, can be modified
  * 		Input:
  * 			const ArgType& 		x: 		The argument
  * 			const Option& 		opt: 	The option
  * 			const Parameter&	para:	The parameter
  * 			int 				level:	Print level
  *			Input and Output:
  * 			std::ostream& 		os: 	Output stream
  */
template<class ArgType,class Option,class Parameter>
using PrintFunction = void(*)(const ArgType& x, const Option& opt, const Parameter& para, int level, std::ostream& os);


/** 		Types of termination 
  */
enum TerminalType
{
	C, 		// converged
	U,		// Unbound
	N,		// not feasible
	M 		// maximum iteration is reached
};

/** 		Basic option structure.
  * 		One might define their own option structure by inheriting from this structure.
  * 		Note that the default maximum iteration number is 1000 and the default
  * 		threshold is set as 1e-7.
  */
template<class Scalar>
struct BasicOption
{
	/** the maximum iteration number */
	int 					MaxIter;
	/** the error threshold */
	Scalar 					Threshold;

	BasicOption():MaxIter(1000),Threshold(1e-7){}
};

/**			Basic parameter structure.
  * 		One might define their own parameter structure by inheriting from the following
  * 		structure. At least there must be iteration number, current estimated error, current
  * 		objective value, computation time and termination type.
  */
template<class Scalar>
struct BasicParameter
{
	/** the iteration number */
	int 					IterNum;
	/** the iterative error */
	Scalar 					Error;
	/** current objective value */
	Scalar 					Object;
	/** the computation time */
	double 					ComputationTime;
	/** the terminal type */
	TerminalType 			Termination;

	BasicParameter()
		:
		IterNum(0),
		Error(static_cast<Scalar>(0.0)),
		Object(static_cast<Scalar>(0.0)),
		ComputationTime(0.0),
		Termination(N)
	{}
};

class Timer
{
private:
	clock_t 	__begin;
	clock_t	 	__end;

	bool 		__is_begined;
	bool		__is_ended;
	double 		__time;

public:

	/** default constructor */
	Timer():__is_begined(false),__is_ended(false),__time(0.0){}

	void tic(){ 
		__begin = clock();
		__is_begined = true; 
		__is_ended = false;
	}

	void toc(){ 
		__end = clock(); 
		__time = static_cast<double>(__end-__begin)/CLOCKS_PER_SEC;
		if(!__is_begined) std::cerr<<"COPT Timer Warning: the timer ends without beginning!"<<std::endl;
	}

	void end(){
		__is_begined = false;
		__is_ended = true;
	}

	friend std::ostream& operator<<(std::ostream& os,const Timer& timer){
		os<<"Time elapsed "<<timer.__time<<" seconds"<<std::endl;
		return os;
	}

	double time() const{
		return __time;
	}
};

/** default print function at the beginning of a solver */
void solverBeginPrint(int level,std::ostream& os)
{
	if(level>0){
		os<<std::fixed<<std::setw(8)<<"Iters"<<"\t"<<std::setw(12)<<"Objective"<<"\t"<<std::setw(12)<<"Error"<<std::endl;
	}
}

/** default print function for iteration */
template<class ArgType,class Option,class Parameter>
void iterationPrint(const ArgType& x, const Option& o, const Parameter& para, int level, std::ostream& os)
{
	if(level==0) // nothing will be printed
		return;
	else if(level==1)
	{
		os<<std::fixed<<std::setw(8)<<para.IterNum<<"\t";
		os<<std::setw(12)<<std::setprecision(4)<<std::scientific<<para.Object<<"\t";
		os<<std::setw(12)<<para.Error<<std::endl;
		os<<std::defaultfloat;
	}
}

/** default print function at the end of solving the problem */
template<class ArgType,class Option,class Parameter>
void solverEndPrint(const ArgType& x, const Option& o, const Parameter& para, int level, std::ostream& os)
{
	if(level==0)
		return;
	else if(level==1)
	{
		switch(para.Termination)
		{
			case C:
			{
				os<<"Solver successfully converges."<<std::endl;
			}
			break;
			case U:
			{
				os<<"The problem is unbounded."<<std::endl;
			}
			break;
			case N:
			{
				os<<"The problem is not feasible."<<std::endl;
			}
			break;
			case M:
			{
				os<<"Maximum itreation number is reached."<<std::endl;
			}
			break;
			default:
			{
				os<<"Unknown terminal type for the solver."<<std::endl;
			}
			break;
		}
		os<<"The whole solver costs "<<para.ComputationTime<<" s."<<std::endl;
		os<<"The final objective value is "<<para.Object<<"."<<std::endl;
		os<<para.IterNum<<" iterations are used."<<std::endl;
	}
}

template<class ArgType,class OutputType>
void argEqualToResult(const ArgType& x, OutputType& result)
{
	result = x;
}

/** 		The main class for this header file. A procedural-based framework for 
  *			optimization problems. The solver is used by setting different functions
  * 		that is used. In another word, we separate an optimization problem into
  * 		simple parts. We design a framework like this is because that when we are
  * 		doing optimization research and writing optimization codes we find that there
  * 		are a lot of repeative works and we want to avoid it.
  * 		
  */
template<class Scalar,class ArgType,class OutputType=ArgType,class Option=BasicOption<Scalar>,class Parameter = BasicParameter<Scalar> >
class Solver
{
private:

	typedef COPT::ObjectiveFunction<Scalar,ArgType,Parameter> 			ObjectiveFunction;
	typedef COPT::OneIteration<Scalar,ArgType,Parameter> 				IterationFunction;
	typedef COPT::CheckTermination<ArgType,Option,Parameter> 			TerminationFunction;
	typedef COPT::SolverInitialization<ArgType,Parameter> 				InitializationFunction;
	typedef COPT::GetResult<ArgType,OutputType> 						ArgToResultFunction;
	typedef COPT::PrintFunction<ArgType,Option,Parameter> 				PrintFunction;
	
	/** the argument x */
	ArgType 				__x;
	/** the result */
	OutputType 				__result;
	/** the option of the solver */
	Option 					__op;
	/** the parameter */
	Parameter 				__para;
	/** print level of the solver */
	int 					__print_level;
	/** the call for objective function */
	ObjectiveFunction 		__ob_func;
	/** initialization function */
	InitializationFunction	__init_func;
	/** the iteration for the solver */
	IterationFunction 		__iter_func;
	/** the determination of the terminal condition */
	TerminationFunction 	__ter_func;
	/** function describing how argument is transformed into result */
	ArgToResultFunction 	__arg_to_re_func;
	/** the output stream */
	std::ostream& 			__ostream;
	/** print functions when the solver begins */
	NormalPrintFunction		__begin_print_func;
	/** print function for iteration */
	PrintFunction 			__iter_print_func;
	/** print function at the end */
	PrintFunction 			__end_print_func;

	virtual void doSolve();

public:

	/** constructor */
	Solver(
		const Parameter& para,									// the input parameter has be been given
		ObjectiveFunction obfunc = nullptr, 					
		InitializationFunction initfunc=nullptr,
		IterationFunction iterfunc = nullptr,
		TerminationFunction terfunc = nullptr,
		int printlevel = 1);

	/** solve the problem */
	void solve();

	/** return the result of the solver */
	OutputType result() const;

	/** validate the solver */
	bool validation() const;

	/** setter and getter */
	//%{
	/** 	Set the print level of the solver.
	  * 	The print level of the solver determines what information will be
	  * 	printed. 
	  */
	void setPrintLevel(int printlevel);
	/** set the objective function */
	void setObjectiveFunction(ObjectiveFunction func);
	/** return the current objective value */
	Scalar objective() const;
	/** how to compute objective function */
	ObjectiveFunction objectiveFunction() const;
	/** set the initialization function */
	void setInitializationFunction(InitializationFunction func);
	/** obtain the initialization function */
	InitializationFunction initializationFunction() const;
	/** set the iteration function */
	void setIterationFunction(IterationFunction func);
	/** obtain the iteration function */
	IterationFunction iterationFunction() const;
	/** set the terminal function */
	void setTerminationFunction(TerminationFunction func);
	/** obtain the terminal function */
	TerminationFunction terminationFunction() const;
	/** set the argument to result function */
	void setArgToResultFunction(ArgToResultFunction func);
	/** set the print function for the begining */
	void setBeginningPrintFunction(NormalPrintFunction func);
	/** set the iterative print function */
	void setIterationPrintFunction(PrintFunction func);
	//%}

	/** get the option */
	const Option& option() const;

	/** get the parameter */
	const Parameter& parameter() const;
};

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
Solver<Scalar,ArgType,OutputType,Option,Parameter>::Solver(
		const Parameter& para,
		ObjectiveFunction obfunc, 					
		InitializationFunction initfunc,
		IterationFunction iterfunc,
		TerminationFunction terfunc,
		int printlevel)
	:
	__para(para),
	__print_level(printlevel),
	__ob_func(obfunc),
	__init_func(initfunc),
	__iter_func(iterfunc),
	__ter_func(terfunc),
	__arg_to_re_func(argEqualToResult<ArgType,OutputType>),
	__ostream(std::cout),
	__begin_print_func(&solverBeginPrint),
	__iter_print_func(&iterationPrint<ArgType,Option,Parameter>),
	__end_print_func(&solverEndPrint<ArgType,Option,Parameter>)
{
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::solve()
{
	this->doSolve();
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
OutputType Solver<Scalar,ArgType,OutputType,Option,Parameter>::result() const
{
	return __result;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
bool Solver<Scalar,ArgType,OutputType,Option,Parameter>::validation() const
{
	return __ob_func&&__iter_func&&__init_func&&__iter_func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setPrintLevel(int printlevel)
{
	__print_level = printlevel;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setObjectiveFunction(ObjectiveFunction func)
{
	__ob_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
Scalar Solver<Scalar,ArgType,OutputType,Option,Parameter>::objective() const
{
	if (__ob_func)
		return __ob_func(__x,__para);
	else
		throw COException("Solver error: Objective Function is not defined yet!");
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
typename Solver<Scalar,ArgType,OutputType,Option,Parameter>::ObjectiveFunction Solver<Scalar,ArgType,OutputType,Option,Parameter>::objectiveFunction() const
{
	return __ob_func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setInitializationFunction(InitializationFunction func)
{
	__init_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
typename Solver<Scalar,ArgType,OutputType,Option,Parameter>::InitializationFunction Solver<Scalar,ArgType,OutputType,Option,Parameter>::initializationFunction() const
{
	return __init_func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setIterationFunction(IterationFunction func)
{
	__iter_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
typename Solver<Scalar,ArgType,OutputType,Option,Parameter>::IterationFunction Solver<Scalar,ArgType,OutputType,Option,Parameter>::iterationFunction() const
{
	return __iter_func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setTerminationFunction(TerminationFunction func)
{
	__ter_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
typename Solver<Scalar,ArgType,OutputType,Option,Parameter>::TerminationFunction Solver<Scalar,ArgType,OutputType,Option,Parameter>::terminationFunction() const
{
	return __ter_func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setArgToResultFunction(ArgToResultFunction func)
{
	__arg_to_re_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setBeginningPrintFunction(NormalPrintFunction func)
{
	__begin_print_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::setIterationPrintFunction(PrintFunction func)
{
	__iter_print_func = func;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
void Solver<Scalar,ArgType,OutputType,Option,Parameter>::doSolve()
{
	if(!this->validation())
	{
		std::cerr<<"The solver is not valid in current state. Please make sure that every function needed is set."<<std::endl;
		return;
	}
	// the beginning of copt solver
	__ostream<<"*******************COPT solver begins*******************"<<std::endl;
	__begin_print_func(__print_level,__ostream);
	Timer timer;
	timer.tic();
	__init_func(__x,__para);
	__para.IterNum = 0;
	int i;
	for(i=0; i<__op.MaxIter; ++i)
	{
		__para.Error = __iter_func(__x,__para);
		__para.IterNum ++;
		if(__print_level>0)
			__para.Object = __ob_func(__x,__para);
		__iter_print_func(__x,__op,__para,__print_level,__ostream);
		if(__ter_func(__x,__op,__para))
		{
			__para.Termination = C;
			break;
		}
	}
	// operations when solver ends
	__arg_to_re_func(__x,__result);
	timer.toc();
	if(i==__op.MaxIter) 
		__para.Termination = M; // max iteration number
	__para.ComputationTime = timer.time();
	// print the end of the solver
	__end_print_func(__x,__op,__para,__print_level,__ostream);
	__ostream<<"********************COPT solver ends********************"<<std::endl;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
const Option& Solver<Scalar,ArgType,OutputType,Option,Parameter>::option() const
{
	return __op;
}

template<class Scalar,class ArgType,class OutputType,class Option,class Parameter>
const Parameter& Solver<Scalar,ArgType,OutputType,Option,Parameter>::parameter() const
{
	return __para;
}

};

#endif