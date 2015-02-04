#include "Core"

using namespace COPT;
typedef Array<double,int,Dynamic> DArray;

int main(int argc, char *argv[])
{
	/** a size-specified array */
	Array<double,int,3> a;	// an array with size 3
	a[0]=1.0;
	a[1]=2.0;
	a[2]=3.0;
	std::cout<<"a is "<<a<<std::endl;	// [1,2,3]

	/** set the array */
	double s[3]{3.0,2.0,1.0};
	a.setArray(s);
	std::cout<<"a becomes "<<a<<std::endl; // [3,2,1]

	/** Dynamic array */
	DArray da; 		// an array with dynamic size
	da.setFromArray(a); // set from 'a'
	std::cout<<"da is"<<da<<std::endl; //[3,2,1]

	try
	{
		a.setFromArray(da); //not allowed even with the same size
	}
	catch(COException &e)
	{
		std::cout<<e.what()<<std::endl;
	}

	try
	{
		a.reset(5); //not allowed to be resized
	}
	catch(COException &e)
	{
		std::cout<<e.what()<<std::endl;
	}

	/** dynamic array can be changed at any moment */
	double ds[5]{1.0,1.1,1.2,1.3,1.4};
	da.setArray(5,ds);
	std::cout<<"da becomes "<<da<<std::endl; //[1,1.1,1.2,1.3,1.4]

	/** normal array-like element access */
	std::cout<<"the third element of a is "<<a[2]<<std::endl; //1.0
	std::cout<<"the fourth element of da is "<<da[3]<<std::endl; //1.3

	/** test */
	double &&e=0;
	std::for_each(a.begin(),a.end(),[&](double x){e+=x;});
	std::cout<<e<<std::endl;
}