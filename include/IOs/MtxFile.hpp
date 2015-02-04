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


#ifndef MTX_FILE_HPP__
#define MTX_FILE_HPP__


namespace COPT
{

/************************Mtx File reader*********************/

template<class T>
inline void readMtxFile(const std::string&filename , T& t)
{
	readMtxFile(filename,t,typename T::ObjectCategory());
}

template<class Matrix>
inline void readMtxFile( const std::string& filename , Matrix& mat , const matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		typedef typename Matrix::scalar 			scalar;
		typedef typename get_pod_type<scalar>::type podtype;		
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename Matrix::index rows,cols,i=0;
		bool first = true;
		while (fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%' )
				continue;
			else if(temp.empty())
				continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if(first)
				{
					if(tokens.size()!=2)
					{
						std::cerr<<"the format of matrix market file is not right !"<<std::endl;
						return;
					}
					if(sizeof(typename Matrix::index)==4)
					{
						rows = atoi(tokens[0].c_str());
						cols = atoi(tokens[1].c_str());
					}
					else
					{
						rows = atol(tokens[0].c_str());
						cols = atol(tokens[1].c_str());
					}
					mat.resize(rows,cols);
					first = false;
				}
				else
				{
					// real matrix
					if (tokens.size() == 1 )
					{
						if(is_real<scalar>::value)
							mat.dataPtr()[i++]=static_cast<podtype>(atof(tokens[0].c_str()));
						else
						{
							ForceAssignment(std::complex<podtype>(static_cast<podtype>(atof(tokens[0].c_str())),0.0),mat.dataPtr()[i++]);
						}
					}
					else if(tokens.size() == 2 )
					{
						if(is_complex<scalar>::value)
							ForceAssignment(std::complex<podtype>(static_cast<podtype>(atof(tokens[0].c_str())),static_cast<podtype>(atof(tokens[1].c_str()))),mat.dataPtr()[i++]);
						else
							throw COException("Mtx reading error: try to assign complex value to a real matrix!");
					}
				}
			}
		}
	}
	else
	{	
		std::cerr<<"Mtx file reader warning: File extension is not mtx!"<<std::endl;
		return;
	}
}

template<class SpMatrix>
inline void readMtxFile( const std::string& filename , SpMatrix& mat , const sp_matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		typedef typename SpMatrix::scalar 				scalar;
		typedef typename get_pod_type<scalar>::type		podtype;
		// deal the matrix
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename SpMatrix::index rows,cols,nnz;
		bool first = true;
		std::vector<typename SpMatrix::Triplet> tris;
		while(fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%')
				continue;
			else if(temp.empty())
				continue;
			else if(temp.c_str()[0] == ' ')
				continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if( first )
				{
					if(tokens.size() == 3)
					{
						// matrix
						if(sizeof(typename SpMatrix::index)==4)
						{
							rows = atoi(tokens[0].c_str());
							cols = atoi(tokens[1].c_str());
							nnz = atoi(tokens[2].c_str());
						}
						else
						{
							
							rows = atol(tokens[0].c_str());
							cols = atol(tokens[1].c_str());
							nnz = atol(tokens[2].c_str());
						}
						tris.reserve(nnz);
					}
					else
					{
						throw COException("Error, the input should be a matrix not a vector!");
					}
					first = false;
				}
				else
				{
					if(sizeof(typename SpMatrix::index)==4)
					{
						tris.push_back(typename SpMatrix::Triplet(atoi(tokens[0].c_str())-1,atoi(tokens[1].c_str())-1,atof(tokens[2].c_str())));
					}
					else
					{
						tris.push_back(typename SpMatrix::Triplet(atol(tokens[0].c_str())-1,atol(tokens[1].c_str())-1,atof(tokens[2].c_str())));
					}
				}
			}
		}
		mat.fastSetFromTriplets(rows,cols,tris.begin(),tris.end());
		fin.close();
	}
	else
		return;
}

template<class Vector>
inline void readMtxFile( const std::string& filename , Vector& vec , const vector_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)=="mtx")
	{
		typedef typename Vector::scalar 			scalar;
		typedef typename get_pod_type<scalar>::type podtype;
		// deal the matrix
		std::ifstream fin(filename);
		if(!fin)
		{
			std::cerr<<"file does not exist!"<<std::endl;
			return;
		}
		std::string temp;
		typename Vector::index size,i = 0;
		bool first = true;
		while(fin)
		{
			std::getline(fin,temp);
			if(temp.c_str()[0] == '%')
				continue;
			else if(temp.empty())
				continue;
			// else if(temp.c_str()[0] == ' ')
			// 	continue;
			else
			{
				std::istringstream iss(temp);
				std::vector<std::string> tokens;
				std::copy(std::istream_iterator<std::string>(iss),
					std::istream_iterator<std::string>(),
					std::back_inserter(tokens));
				if( first )
				{
					if(tokens.size() == 2)
					{
						// matrix
						if(sizeof(typename Vector::index)==4)
						{
							if(atoi(tokens[1].c_str()) != 1 )
								throw COException("Error, the input is a matrix! Please check it out!");
							size = atoi(tokens[0].c_str());
						}
						else
						{
							if(atol(tokens[1].c_str()) != 1 )
								throw COException("Error, the input is a matrix! Please check it out!");
							size = atol(tokens[0].c_str());
						}
						vec.resize(size);
					}
					else
					{
						throw COException("Error, the input should be a vector not a sparse matrix!");
					}
					first = false;
				}
				else
				{
					if( tokens.size() == 1 )
					{
						// real vector
						if(is_real<scalar>::value)
							vec(i++) = static_cast<scalar>(atof(tokens[0].c_str()));
						else
							ForceAssignment(std::complex<podtype>(static_cast<podtype>(atof(tokens[0].c_str())),0.0),vec(i++));
					}
					else
					{
						// complex vector
						if(is_complex<scalar>::value)
						{
							ForceAssignment(std::complex<podtype>(static_cast<podtype>(atof(tokens[0].c_str())),static_cast<podtype>(atof(tokens[1].c_str()))),vec(i++));
						}
						else
							throw COException("Mtx file reading error: try to assign complex value to a real vector!");
					}
				}
			}
		}
		fin.close();
	}
	else
		return;
}

/****************************Mtx file writer*********************************/
template<class T>
void writeMtxFile( const std::string& filename , const T& t )
{
	writeMtxFile(filename,t,typename T::ObjectCategory());
}

template<class Vector>
void writeMtxFile( const std::string& filename , const Vector& vec , const vector_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket vector data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<vec.size()<<' '<<1<<std::endl;
	for ( int i = 0 ; i < vec.size() ; ++ i )
	{
		fout<<vec(i)<<std::endl;
	}
}

template<class Matrix>
void writeMtxFile( const std::string& filename , const Matrix& mat , const matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket matrix data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<mat.rows()<<' '<<mat.cols()<<std::endl;
	for ( int i = 0 ; i < mat.rows()*mat.cols() ; ++ i )
	{
		fout<<mat.dataPtr()[i]<<std::endl;
	}
}

template<class SpMatrix>
void writeMtxFile( const std::string& filename , const SpMatrix& spmat , const sp_matrix_object& )
{
	if(filename.substr(filename.find_last_of(".")+1)!="mtx")
	{
		std::cerr<<"Mtx file writing warning: File extension is not mtx!"<<std::endl;
	}
	std::ofstream fout(filename);
	fout<<"%% MatrixMarket sparse matrix data"<<std::endl;
	fout<<"%% Generated by open source library COPT"<<std::endl;
	fout<<spmat.rows()<<' '<<spmat.cols()<<' '<<spmat.elementSize()<<std::endl;

}

}

#endif