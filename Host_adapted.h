// include guard
#ifndef __HOST_ADAPTED_H_INCLUDED__
#define __HOST_ADAPTED_H_INCLUDED__

//inclusions
#include <iostream>
#include <random>
#include <boost/ptr_container/ptr_vector.hpp>
#include "Bacterium_adapted.h"  //because includes doing stuff to Bacteria

using namespace std;

//specify class
class Host{ //# : boost::noncopyable

public:
	boost::ptr_vector<Bacterium> WithBact; 
	boost::ptr_vector<Bug> CandBact;
	//#check http://stackoverflow.com/questions/9469968/stl-container-with-stdunique-ptrs-vs-boostptr-container

	void add_bact(Strain& const a){ //create new Bacteria and add
		WithBact.push_back(new Bacterium(a));
	}

	size_t getMGE(void); //gives number of MGEs in host
	bool hasMGE(void); //says whether there is at least one MGE in host

	Host(); //constructor
};

Host::Host() //this is constructor: (void) means takes no argument
{
	WithBact.reserve(4); CandBact.reserve(2);

}

inline size_t Host::getMGE(void){  //#see if any inbuilt function to count pointers
	if (WithBact.size() == 0){
		return 0;
	}
	else {
		size_t mges = 0; //counter
		for (size_t i = 0; i < WithBact.size(); ++i){ //for every bacterium in host 
			mges += WithBact[i].WithMGE.size(); //adds number of MGE in bacterium to counter
			}
		return mges;
	}
}

inline bool Host::hasMGE(void){
	if (WithBact.size() == 0){
		return 0;
	}
	else {
		for (size_t i = 0; i < WithBact.size(); ++i){ //for every bacterium in host 
			if (WithBact[i].WithMGE.size() > 0){ //if there are any MGEs in bacterium
				goto hasMGE;
			}
		}
		return 0; //no MGE if hasn't jumped to hasMGE
		
	hasMGE: 
		return 1;
	}
}

#endif