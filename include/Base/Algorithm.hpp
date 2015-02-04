//		Copyright (C) Ruimin Wang, ruimin.wang13@gmail.com
//		Copyright (C) MathU

#ifndef ALGORITHM_HPP__
#define ALGORITHM_HPP__

namespace COPT
{
/*			General algorithm in COPT. 
 *
 *
 */
template<class kernel>
class Algorithm
	:
	public COPTObject,
	noncopyable
{
private:
	typedef typename kernel::scalar 			scalar;
	typedef typename kernel::podscalar 			podscalar;
	typedef typename kernel::index 				index;
	typedef typename kernel::Vector 			Vector;
	typedef typename kernel::Matrix 			Matrix;

	
public:

	Algorithm();
};
}

#endif