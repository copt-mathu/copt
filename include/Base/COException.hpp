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


#ifndef COEXCEPTION_HPP__
#define COEXCEPTION_HPP__

#include <exception>
#include <string>

// the base class of the exception that is used in the library
namespace COPT
{
class COException : public std::exception
{
public:
	COException(const char* str):__str("COpt Exception: "){ __str.append(str); }
	virtual const char* what() const throw()
	{
		return __str.c_str();
	}
	~COException() throw(){}
private:
	std::string			__str;
};
};

#endif