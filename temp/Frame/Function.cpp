#include "Function.h"
#include <iostream>

// // Implementation of 'Function.h'


// // Default constructor for 'Function'
// //	a function's differential is not cacluable as default
// Function::Function()
// 	:__exact_diff(false)
// {
// }
// /*
// 	Cosine Function
// */
// ScalarCosineFunc::ScalarCosineFunc( double l )
// {
// 	__fp.resize(1);
// 	__fp[0] = l;
// }

// double ScalarCosineFunc::operator()(double x)
// {
// 	return cos(__fp[0]*x);
// }

// double ScalarCosineFunc::diff( double x )
// {
// 	return -__fp[0]*sin(__fp[0]*x);
// }

// /*
// 	'Vector Cosine Function'
// */
// VectorCosineFunc::VectorCosineFunc()
// 	:__dim(0)
// {

// }
// VectorCosineFunc::VectorCosineFunc(int dim, double *vec)
// 	:__dim(dim)
// {
// 	__fp.setVector(dim,vec);
// }

// Vector VectorCosineFunc::operator() ( const Vector& vec )
// {
// 	return Vector(10);
// }
