#include "PoissonDisk.h"
#include <iostream>
#include <fstream>

void poissonRandomSampling(ind dim, scalar r, Vector &sample)
{
	scalar radius = std::uniform_real_distribution<>(r,r*2)(COPT::copt_rand_eng);
	sample.resize(dim);
	for (ind i = 0; i < dim; ++ i)
	{
		sample[i]= std::uniform_real_distribution<>(-1,1)(COPT::copt_rand_eng);
	}
	sample.normalize();
	sample = sample*radius;
}

PoissonDisk::PoissonDisk(ind dim, scalar r, ind k)
	:
	__dim(dim),
	__r(r),
	__k(k),
	__min_vec(dim),
	__max_vec(dim),
	__grid_num(dim)
{
	__sr = __r/(std::sqrt(__dim));
}

void PoissonDisk::setRange(scalar min, scalar max)
{
	for(int i = 0; i < __dim ; ++i)
	{
		__min_vec = min;
		__max_vec = max;
	}
}

void PoissonDisk::setRange(const Vector &minvec, const Vector &maxvec)
{
	if(minvec.size()!=__dim||maxvec.size()!=__dim)
		throw COPT::COException("make sure that the dimensions of input vector is equal to the dimension of space!");
	__min_vec = minvec;
	__max_vec = maxvec;
}

void PoissonDisk::constructGrid()
{
	// get the size
	ind gridsize = 1;
	for (ind i = 0; i < __dim; ++ i)
	{
		__grid_num[i] = (static_cast<ind>((__max_vec[i]-__min_vec[i])/__sr)+1);
		gridsize *= __grid_num[i];
		
	}
	__m.resize(gridsize,-1);
	for ( int i = 0 ;i < gridsize; ++i)
	{
		__m[i]=-1;
	}
}

ind PoissonDisk::getIndex(const Vector &c)
{
	ind index = 0, shift = 1;
	for (ind i = 0; i < __dim; ++ i)
	{
		if(c[i]<__min_vec[i]||c[i]>__max_vec[i])
			return -1;
		ind j = static_cast<ind>((c[i]-__min_vec[i])/__sr);
		if(i==0)
			index += j;
		else
			index += j*shift;

		if(i<__dim-1)
			shift *= (static_cast<ind>((__max_vec[i]-__min_vec[i])/__sr)+1);
	}
	return index;
}

void PoissonDisk::initialize()
{
	// initialize a poisson disk
	Vector vec(__dim);
	for (ind i = 0; i < __dim; ++ i)
	{
		vec[i] = std::uniform_real_distribution<>(__min_vec[i],__max_vec[i])(COPT::copt_rand_eng);
	}
	__pts.push_back({vec,true});
	__active_set.push_back(0);
	ind index = getIndex(vec);
	__m[index] = 0;
}

bool PoissonDisk::generateNewPoint()
{
	auto dis = __active_set.begin();
	std::advance(dis, std::uniform_int_distribution<>(0,__active_set.size()-1)(COPT::copt_rand_eng));
	auto centeriter = __pts.begin();
	std::advance(centeriter, *dis);
	auto center = centeriter->coord;

	for(ind i = 0; i < __k; ++i)
	{
		Vector direction;
		poissonRandomSampling(__dim,__r,direction);
		ind index = getIndex(center+direction);
		if(index==-1)
			continue;
		if(__m[index]==-1)
		{
			if(!checkPoint(center+direction))
			{
				continue;
			}
			__pts.push_back({center+direction,true});
			__active_set.push_back(__num);
			__m[index]=__num;
			++ __num;
			return true;
		}
	}
	__active_set.erase(dis);
	return false;
}

bool PoissonDisk::checkPoint(const Vector &c)
{
	std::vector<ind> lc(__dim);
	for (ind i = 0; i < __dim; ++ i)
	{
		lc[i] = static_cast<ind>((c[i]-__min_vec[i])/__sr);
	}
	ind currind;
	ind dim = __dim;
	const std::list<Pt>& pts = __pts;
	std::vector<ind> v(__dim);
	bool found=true;;
	checkTraverse(0,lc,c,v,found);
	return found;
}

void PoissonDisk::checkTraverse(ind level, const std::vector<ind>& lc, const Vector& c, std::vector<ind> &vec, bool &found)
{
	if(level==__dim)
	{
		ind index = 0;
		ind mul = 1;
		for (ind j=0;j<__dim;++j){
			index += mul*vec[j];
			mul *= __grid_num[j];
		}
		if(__m[index]==-1)
			return;
		else{
			auto pt = __pts.begin();
			std::advance(pt, __m[index]);
			if(COPT::distance(pt->coord,c)<__r)
			{
				found = false;
				return;
			}
		}
	}
	else
	{
		for (ind i=-1;i<=1;++i){
			if(lc[level]==0&&i==-1)
			{
				continue;
			}
			else if(lc[level]==__grid_num[level]-1&&i==1)
			{
				continue;
			}
			else
			{
				vec[level]=lc[level]+i;
				checkTraverse(level+1,lc,c,vec,found);
			}
		}
	}
}

void PoissonDisk::generate()
{
	__active_set.clear();
	__pts.clear();
	std::cout<<"1"<<std::endl;
	constructGrid();
	std::cout<<"2"<<std::endl;
	// initialize a point
	initialize();
	std::cout<<"3"<<std::endl;
	__num = 1;

	while(!__active_set.empty())
	{
		generateNewPoint();
	}
	std::cout<<"there are "<<__pts.size()<<" points"<<std::endl;
}

void PoissonDisk::output(const std::string& str)
{
	std::ofstream fout(str);
	for ( auto iter = __pts.begin(); iter != __pts.end(); ++ iter)
	{
		for (ind i = 0; i < iter->coord.size(); ++ i)
		{
			fout<<iter->coord(i)<<"\t";
		}
		fout<<std::endl;
	}
}

const std::list<PoissonDisk::Pt>& PoissonDisk::points()const
{
	return __pts;
}