#ifndef POISSON_DISK_H__
#define POISSON_DISK_H__

#include "../../include/Core"
#include <vector>

typedef double 							scalar;
typedef int 							ind;
typedef COPT::KernelTrait<scalar,ind>	kernel;
typedef kernel::Vector 					Vector;

const double pi = 3.1415926535897;

void poissonRandomSampling(ind dim, scalar r, Vector &sample);

/**		Poisson Disk produces a poisson sampling
 *
 *
 *
 */
class PoissonDisk
{
	
private:

	struct Pt{
		Vector 		coord;
		bool 		active;
	};

	/** dimension of Vector */
	ind 						__dim;
	/** the list of points */
	std::list<Pt> 				__pts;
	/** radius */
	scalar 						__r;
	/** k samples per iteration */
	ind 						__k;
	/** radius/sqrt(dim) */
	scalar 						__sr;
	/** the rectangle region */
	Vector 						__min_vec;
	Vector 						__max_vec;
	/** the mark of list */
	std::vector<ind>			__m;
	/** the current number of points */
	ind 						__num;
	/** active set */
	std::list<ind>				__active_set;
	/** the grid information */
	Vector 						__grid_num;

	PoissonDisk();

public:

	
	PoissonDisk(ind dim, scalar r, ind k);

	void setRange(scalar min, scalar max);
	void setRange(const Vector &minvec, const Vector &maxvec);

	void constructGrid();

	ind getIndex(const Vector &c);

	void initialize();

	bool generateNewPoint();

	bool checkPoint(const Vector &c);

	void checkTraverse(ind level, const std::vector<ind>& lc, const Vector& c, std::vector<ind> &v, bool& found);

	void generate();

	void output(const std::string& str);

	const std::list<Pt>& points()const;

};


#endif