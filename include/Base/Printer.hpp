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


#ifndef PRINTER_HPP__
#define PRINTER_HPP__

namespace COPT
{


/** 	The printer of COPT library
 *
 *
 */
class Printer
	:
	public COPTObject
{
	int 		__print_level;

	/** data objects */
	static inline std::ostream& printObjectTag( const data_object& 			, std::ostream& os);
	static inline std::ostream& printObjectTag( const matrix_object& 		, std::ostream& os);
	static inline std::ostream& printObjectTag( const sp_matrix_object& 	, std::ostream& os );
	static inline std::ostream& printObjectTag( const vector_object& 		, std::ostream& os );
	/** functions */
	static inline std::ostream& printObjectTag( const scalar_func_object& 	, std::ostream& os );
	static inline std::ostream& printObjectTag( const vector_func_object& 	, std::ostream& os );
	/** problems */
	static inline std::ostream& printObjectTag( const lasso_problem& 		, std::ostream& os );
	/** solvers */
	static inline std::ostream& printObjectTag( const solver_object& 		, std::ostream& os );
	static inline std::ostream& printObjectTag( const proximal_solver& 		, std::ostream& os );
	static inline std::ostream& printObjectTag( const admm_solver& 			, std::ostream& os );
	static inline std::ostream& printObjectTag( const fista_solver& 		, std::ostream& os );

public:
	Printer()
		:
		__print_level(0)
	{
	}

	std::ostream& print(const COPTObject& obj , std::ostream& os );


	/*		Static function for printing the information of an object from COPT.
	 *		Note that if no specific tag is set, the information of base class
	 *		will be printed.
	 *
	 *		An example:
	 *
	 *			class Array is a derived class from COPTObject.
	 *			Array a;
	 *			Printer::printType(a,std::cout);
	 */
	template<class Object>
	static inline std::ostream& printType( const Object& , std::ostream& os );
	// taking std::cout as output stream
	template<class Object>
	static inline std::ostream& printType( const Object& );
};

template<class Object>
std::ostream& Printer::printType( const Object& obj , std::ostream& os )
{
	return printObjectTag(typename Object::ObjectCategory(),os);
}

template<class Object>
std::ostream& Printer::printType( const Object& obj )
{
	return printType(obj,std::cout);
}

std::ostream& Printer::printObjectTag( const data_object& , std::ostream& os )
{
	os<<"-----------------This is a data object in COPT------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const matrix_object& , std::ostream& os )
{
	os<<"--------------------This is a matrix in COPT--------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const sp_matrix_object& , std::ostream& os )
{
	os<<"------------------This is a sparse matrix in COPT---------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const vector_object& , std::ostream& os )
{
	os<<"---------------------This is a vector in COPT-------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const scalar_func_object& , std::ostream& os )
{
	os<<"-------------------This is a scalar function in COPT------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const vector_func_object& , std::ostream& os )
{
	os<<"-----------------This is a vector function in COPT--------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const lasso_problem& , std::ostream& os )
{
	os<<"------------------This is a lasso problem in COPT---------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const solver_object& , std::ostream& os )
{
	os<<"---------------------This is a solver int COPT------------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const proximal_solver& , std::ostream& os )
{
	os<<"-----------------This is a proximal problem in COPT-------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const admm_solver& , std::ostream& os )
{
	os<<"------------------This is an ADMM solver in COPT----------------"<<std::endl;
	return os;
}

std::ostream& Printer::printObjectTag( const fista_solver& , std::ostream& os )
{
	os<<"------------------This is an FISTA solver in COPT---------------"<<std::endl;
	return os;
}

}

#endif