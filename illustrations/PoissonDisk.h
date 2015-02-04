#ifndef POISSON_DISK_H__
#define POISSON_DISK_H__

#include "../include/Core"

typedef double 							scalar;
typedef int 							ind;
typedef COPT::KernelTrait<scalar,ind>	kernel;
typedef kernel::Vector 					Vector;


class PoissonDisk
{
private:
	
	struct Pt{
		Vector 		coord;
		bool 		active;
	};

	/** dimension of Vector */
public:
};

#endif