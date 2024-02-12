// include guard
#ifndef __BACTERIUM_ADAPTED_H_INCLUDED__
#define __BACTERIUM_ADAPTED_H_INCLUDED__

//inclusions
#include <iostream>
#include <random>
#include <boost/ptr_container/clone_allocator.hpp> //for clones
#include <boost/ptr_container/ptr_vector.hpp>
#include "MGE_adapted.h"


using namespace std;

//////////////
class Strain { //: public Counter<Strain>{ //#shouldn't need counter as have use_count for shared pointer
public:
	vector<double> StrainComp;
	string StrainName;
	double StrainPrev;
	double StrainBeta;
	double StrainGamma;
	double StrainFitness;
	bool StrainResistant;
	int StrainIntro;
	int StrainIndex;

	Strain(string,vector<double>, double, double, double, double, bool, int); //declare constructor #should take all parms
	//using Counter<Strain>::getCount;

	void printstrain() {
		cerr<< "Strain Name " << StrainName << "\n"
			<< "StrainIntro " << StrainIntro << "\n"
			<< "StrainPrev " << StrainPrev << "\n"
			<< "StrainBeta " << StrainBeta << "\n"
			<< "StrainGamma " << StrainGamma << "\n"
			<< "StrainFitness " << StrainFitness << "\n"
			<< "StrainResistant " << StrainResistant << "\n";
	}

};

Strain::Strain(string a, vector<double> g, double b, double c, double d, double e, bool f, int intro) //constructor
	:StrainName(a), StrainComp(g), StrainPrev(b), StrainBeta(c), StrainGamma(d), StrainFitness(e), StrainResistant(f), StrainIntro(intro)
{
//	printstrain();
}

///////////////
class Bacterium : public Bug, public Counter<Bacterium>{ //: boost::noncopyable
public: 
	vector<MGEType*> WithMGE; //this is vector of all MGEs inside bacterium
	Strain* OfStrain;
	double BactGamma;
	double BactBeta;
	double BactFitness;
	int BactCarriage;
	bool BactResistant; 
	
	void add_MGE(MGEType& const a){
			//cout << "This MGE has site " << a.MGETypeSite << "\n";

		if (WithMGE[a.MGETypeSite] == nullptr) {
			WithMGE[a.MGETypeSite] = & a;

			//update bacterial resistance if MGE is resistant
			if (a.MGETypeResistant == 1) {
				//if new MGE is resistant, make bacteria resistant
				BactResistant = 1;
			}

			//update bacterial fitness
			BactFitness *= a.MGETypeChangeFitness;
			//BactGamma += (1 - a.MGETypeChangeGamma)*(1 - BactGamma);
			BactGamma = 1 - pow((1 - BactGamma),a.MGETypeChangeGamma); //1-(1-BactGamma)^a.MGETypeChangeGamma
			BactBeta = 1 - pow((1 - BactBeta), 1/a.MGETypeChangeBeta);

		} 
	}

	using Counter<Bacterium>::getCount;

	//competition function
	unsigned int competition_cont(Bacterium& a, int b2index, int thread) { //competition loop
		
		//see if competition occurs
		if (ranf_mt(thread) > OfStrain->StrainComp[b2index]){ //b is index giving strain of second bact
			return 0; //no competition
		}
			double fit = BactFitness; //fitness score of first bact (invader) is halved
			double denom = fit + a.BactFitness; //sum of fitness scores of first and second bact
			if (denom == 0) { //avoids error of denominator 0 producing NaN --> bias to 2nd bacteria
				return floor(ranf_mt(1)*2)+1;//returns 1 or 2 at random 
			}
				if (ranf_mt(thread) < (fit / denom)) { //if random number less than proportion
					return 1; //first bact wins
				}
					return 2; //second bact wins
				}

	unsigned int competition_acq(Bacterium& a, int b2index, int thread, double penalty) { //competition loop
		//see if competition occurs
		if (ranf_mt(thread) > OfStrain->StrainComp[b2index]) { //b is index giving strain of second bact
			return 0; //no competition
		}
		double fit = BactFitness; 
		double denom = fit + a.BactFitness + penalty; //sum of fitness scores of first and second bact
		if (denom == 0) { //avoids error of denominator 0 producing NaN --> bias to 2nd bacteria
			return floor(ranf_mt(1) * 2) + 1;//returns 1 or 2 at random 
		}
		if (ranf_mt(thread) < (fit / denom)) { //if random number less than proportion
											   //if (ranf_mt(thread) < 0.5){ //if random number less than proportion
			return 1; //first bact wins
		}
		return 2; //second bact wins
	}

	//competition function - decreasing fitness if acquired later
	unsigned int competition_lastin(Bacterium& a, int b2index, int thread, int b, int b2) { //b is index of the more recent addition, b2 is index of the other one
		//see if competition occurs
		if (ranf_mt(thread) > OfStrain->StrainComp[b2index]) { //b2index is index giving strain of second bact
			return 0; //no competition
		}

		//fitness of each bacterium is multipled by PDF of gamma-distribution of position in vector
		//distribution is gamma(x, alpha, beta) where alpha=1 and beta =2
		//because alpha =1, PDF is same as PDF of exponential distribution = lambda*exp(-lambda*x), where lambda = 1/beta
		//so PDF is (1/beta)*exp(-x/beta)

		double fit = BactFitness * (1/2)*exp(-b/2); //fitness score of first bact, by gamma distribution
		double denom = fit + a.BactFitness* (1 / 2)*exp(-b2 / 2); //sum of fitness scores of first and second bact
		if (denom == 0) { //avoids error of denominator 0 producing NaN --> bias to 2nd bacteria
			return floor(ranf_mt(1) * 2) + 1;//returns 1 or 2 at random 
		}
		if (ranf_mt(thread) < (fit / denom)) { //if random number less than proportion
			return 1; //first bact wins
		}
		return 2; //second bact wins
	}
		
	//competition function - no denominator, just whether one has higher competition fitness than other
	unsigned int competition_lessfit(Bacterium& a, int b2, int thread, double penalty) { //competition loop
		//see if competition occurs
		if (ranf_mt(thread) > OfStrain->StrainComp[b2]) { //b is index giving strain of second bact
			return 0; //no competition
		}
		if (BactFitness > a.BactFitness+penalty) {	//if first bact (invader's) competition fitness is more than second bact (other invader, or resident)'s fitness
			return 1; // first bact wins
		}
		else {
			return 2; //first bact doesn't win
		}
	}

	virtual Bug* do_clone() const
	{
		return new Bacterium(*this);
	}
	Bacterium(Strain& const); //this is the destructor - same name, no output
	~Bacterium() {
	} 
};

Bacterium::Bacterium(Strain& const a) //this is constructor - not called when cloning?
	: Bug()
{
	OfStrain = &a;
	BactFitness = a.StrainFitness; //gets value from Strain at time of construction
	BactResistant = a.StrainResistant; //takes Strain default when constructed
	BactGamma = a.StrainGamma;
	BactBeta = a.StrainBeta;
}





#endif
