// include guard
#ifndef __TEMPLATES_ADAPTED_H_INCLUDED__
#define __TEMPLATES_ADAPTED_H_INCLUDED__

//inclusions
#include <iostream>
#include <random>
#include <boost/ptr_container/clone_allocator.hpp> //for clones

using namespace std;

////////////////
template <class T> //this is a counter for anything
class Counter
{
private:
	static int count;
public:
	Counter()
	{
		count++;
	}
	Counter(const Counter &c)
	{
		count++;
	}
	~Counter() //counter goes down when destructor is called
	{
		count--;
	}
	static int getCount() {

		return count;
	}
};
template <class T> int Counter<T>::count = 0;


template <typename T>
void printvector(const vector<T> &a){ //#could write this as template function
	cout << "print vector: ";
	size_t n = a.size();
	for (size_t i = 0; i < n; ++i){
		cout << a[i] << " ";
	}
	cout << endl;
}
///////////////
// abstract base class for Bacterium and MGE
class Bug : public Counter<Bug>{//: boost::noncopyable { 
	/*http://www.boost.org/doc/libs/1_55_0/libs/ptr_container/doc/examples.html */

private:
	virtual Bug* do_clone() const = 0;
	/*this makes a pure virtual function - means that base class only acts as base class, 
	and can't create object purely of this type */


public:
	Bug(){}
	virtual ~Bug() {} //{ cout << "Bug destructor called\n"; } //destructor MUST be virtual

	Bug* clone() const
	{
		return do_clone();
	}


};

Bug* new_clone(const Bug& a) //this is the command that goes in main()
{
	return a.clone();
}


/////////////
//functions

void addbool(const vector<bool> &a, vector<int> &b){
	size_t n = a.size();
	for (size_t i = 0; i < n; ++i) {
		//b[i] = int(a[i]) + b[i]; //change bool to int
		b[i] += int(a[i]);
	}
}

//sample from uniform distribution
//http://stackoverflow.com/questions/2254498/c-how-i-can-get-random-value-from-01-to-12/2254535#2254535

//#possible problems?!?!?!?!?!?!
int sampleint(int limit) {
	//return a random number between 0 and limit inclusive.
	cout << "limit is " << limit;
	cout << "RAND_MAX is " << RAND_MAX;
	int divisor = RAND_MAX / (limit + 1);
	cout << " divisor is " << divisor;
	int retval;

	do {
		retval = rand() / divisor;
	} while (retval > limit);

	return retval;
}

#endif