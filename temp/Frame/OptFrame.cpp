#include "OptFrame.h"

namespace COpt
{
// "Parameter"
VectorOptimization::Parameter::Parameter()
	:__iterthresh(100),__errorthresh(1e-7)
{
}

int VectorOptimization::Parameter::iterationThreshold()
{
	return __iterthresh;
}

void VectorOptimization::Parameter::setIterationThreshold( int iter )
{
	__iterthresh = iter;
}

double VectorOptimization::Parameter::errorThreshold()
{
	return __errorthresh;
}

void VectorOptimization::Parameter::setErrorThreshold(double e)
{
	__errorthresh = e;
}

// "OptOutput"

VectorOptimization::OptOutput::OptOutput()
	:__iterationnumber(0)
{
}

void VectorOptimization::OptOutput::append(double e)
{
	++ __iterationnumber;
	__errorvec.push_back(e);
}

// "VectorOptimization"
VectorOptimization::VectorOptimization()
	:__v(1)
{
}

VectorOptimization::~VectorOptimization()
{

}

};