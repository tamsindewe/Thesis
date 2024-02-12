// include guard
#ifndef __MGE_ADAPTED_H_INCLUDED__
#define __MGE_ADAPTED_H_INCLUDED__

//inclusions
#include <iostream>
#include <random>
#include <boost/ptr_container/clone_allocator.hpp> //for clones
#include <boost/ptr_container/ptr_vector.hpp>
#include "Templates_adapted.h"

using namespace std;

//////////////
class MGEType : public Counter<MGEType*>{ //using MGE as 'type + suffix' can be reserved keywords
public:
	string MGETypeName;
	double MGETypeBeta;
	double MGETypeDefence; //probability BEATS host defences
	double MGETypeChangeFitness; 
	double MGETypeChangeGamma;
	double MGETypeChangeBeta;
	int MGETypeSite;
	int MGETypeIndex;
	int MGETypeIntro;
	bool MGETypeResistant;
	MGEType(string, double, double, double,  double, double, int,  int, bool); //declare constructor #should take all parms
};

MGEType::MGEType(string a, double b, double c,double d, double g, double gg, int f,  int h, bool e) //constructor
	:MGETypeName(a), MGETypeBeta(b), MGETypeDefence(c), MGETypeChangeFitness(d), MGETypeChangeGamma(g), MGETypeChangeBeta(gg), MGETypeSite(f), MGETypeIntro(h), MGETypeResistant(e)
{
}

//////////////////////
class MGE : public Bug, public Counter<MGE>{ //: boost::noncopyable

public:
	shared_ptr<MGEType> OfType;

	using Counter<MGE>::getCount;

	//clear transfers that shouldn't have happened
	//######using this?!?!!?!?
	unsigned int bactdefence(MGE a, int thread) {
		if (ranf_mt(thread) < a.OfType->MGETypeDefence){ //MGE does not transfer due to 
			return 0; //
		}
		else {
			return 1;
		}
	}

	virtual Bug* do_clone() const
	{
		return new MGE(*this);
	}

	MGE(shared_ptr<MGEType>); //this is the constructor - same name, no output
	~MGE(){ } 
};

MGE::MGE(shared_ptr<MGEType> a) //this is constructor
	: Bug(), OfType(a) //copies existing shared pointer
{

}

#endif