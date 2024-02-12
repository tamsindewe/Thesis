#include <iostream>
#include <sstream>		// for stringstream
#include <vector>
#include <list>
#include <random>
#include <algorithm>	// for count_if, shuffle, boost, min
#include <ctime>		// for time()
#include <cstdlib>		// for srand() and rand() #check this: rand not so good?!
#include <memory>		// for smart pointers
#include <numeric>		// for iota
#include <fstream>		// for writing files
#include <string>		// for reading CSV
#include <iterator>		// to iterate through vectors
#include <chrono>		// for profiling
#include <omp.h>		// parallel processing
#include <windows.h>	// get file name in Windows
#include <limits>
#include <cstddef>

#include <boost/lambda/lambda.hpp> 
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/ptr_container/clone_allocator.hpp> //for clones
#include <boost/variant.hpp> //for complex data structures
#include <boost/algorithm/string/trim.hpp> //for trimming strings

#include "randlib_par.h" //Wes's random sampling file

//include classes
#include "Templates_adapted.h" //#do I need all of these if included within header files?
#include "MGE_adapted.h"
#include "Bacterium_adapted.h" //have these after parameters if need them in .h files
#include "Host_adapted.h"

//set 'death' rate = rate of aging-out
const double HostDeathYr = 0.2;  //based on 'lifespan' of 5 years
const double HostBirthYr = 0.2;	//make same as 'deaths'

//stuff for ELFI fit - generating summary statistics
const double SumStatNoStrains = 10;
double SumStat = 999999; //will be obvious if wrong

//function to get path in Windows
string ExePath() { 
	char buffer[MAX_PATH];
	GetModuleFileName(NULL, buffer, MAX_PATH);
	string::size_type pos = string(buffer).find_last_of("\\/");
	return string(buffer).substr(0, pos);
}

//start main
int main(int argc, char* argv[]) {
	//clock
	auto clock_start = chrono::steady_clock::now();

	//---------------------
	//declare variables
	//---------------------

	//baseline parameters
	int Runs;
	int Years;
	int BurnIn;
	int StepsYr; //number of timesteps/year - get from InputOptions.txt
	int OutputSteps; //how frequently need to print outputs - get from InputOptions.txt

	//modifications from batch file 
	double ParmsMGECompetitionCost;
	double ParmsMGEGammaCost;
	double ParmsMGEBetaCost;
	double ParmsIntraPcomp;
	double ParmsInterPcomp;
	double ParmsGammaFactor;
	double ParmsR0;
	int YearMGEIntro;
	string Page;
	int JobNumber;
	string StrainFile;
	string OptionsFile;

	//ELFI targets - only relevant for ELFI fits
	double TargetPrev;
	double TargetUnique;
	double TargetRes;
	vector<double> vStrainTargetDuration; //this is duration in years
	vector<double> vStrainTargetProportion; //proportion of isolates
	vector<double> vMGETargetPrev; //prevalence of MGE
	double MGETargetChi; //chi squared statistic
	string GetGamma; //now deprecated, allowed for different options for obtaining gamma

	//host parameters
	int HostStartNum;
	double HostABcYr; //this is consumption rate in courses/year

	//---------------------
	//initialise variables
	//---------------------
	if (argc == 13) { //if correct amount of arguments in batchfile

		//use %d for integers, as uses base 10
		//use %lf for double
		sscanf_s(argv[1], "%d", &Runs);
		Page = argv[2];
		sscanf_s(argv[3], "%d", &JobNumber);
		StrainFile = argv[4];
		OptionsFile = argv[5];
		sscanf_s(argv[6], "%lf", &HostABcYr);
		sscanf_s(argv[7], "%lf", &ParmsMGECompetitionCost);
		sscanf_s(argv[8], "%lf", &ParmsMGEGammaCost);
		sscanf_s(argv[9], "%lf", &ParmsMGEBetaCost);
		sscanf_s(argv[10], "%lf", &ParmsR0);
		sscanf_s(argv[11], "%lf", &ParmsIntraPcomp);
		sscanf_s(argv[12], "%lf", &ParmsInterPcomp);

		if ((ParmsIntraPcomp < ParmsInterPcomp)) {
			//cerr <<"Intra: " << ParmsIntraPcomp << "; inter: " << ParmsInterPcomp << " - FAIL\n";
			printf("%.6f", SumStatNoStrains); //this is maximum error - occurs when no strains surviving - allow if resistance
			return 0;
		}
		/*else {
			cerr << "Intra: " << ParmsIntraPcomp << "; inter: " << ParmsInterPcomp << "\n";
		}*/
	}
	else {
		cerr << "Arguments are missing, using default" << '\n';
		Runs = 1;
		Page = "blah";
		JobNumber = 1;
		StrainFile = "blah_InputStrains.txt";
		OptionsFile = "InputOptions.txt";
		HostABcYr = 0.2;
		ParmsMGECompetitionCost = 0.0;
		ParmsMGEGammaCost = 0;
		ParmsMGEBetaCost = 0;
		ParmsR0 = 0;
		ParmsIntraPcomp = 0.15;
		ParmsInterPcomp = 0.03;
	}

	//---------------------
	//set up Wes randlib sampling
	//---------------------
	initSeeds(JobNumber, JobNumber + 1);
	int thread = 0;

	//---------------------
	//get run options
	//---------------------
	int in; //input for loop
	int jn; //input for loop

	//output options
	bool OutputBactAllSites; //count how many bacteria have all MGE sites filled
	bool OutputBactClearances; //how bacteria are cleared
	bool OutputBactCarriage; //output carriage duration
	bool OutputHostBactDist; //distribution of number of isolates per host
	bool OutputResistance; // distribution of number of resistant isolates per host

	//competition options
	char CompetitionCont; //"N" means normal: fitness doesn't change depending on element in vector; "L" means last - reduced fitness for later acquisitions
	char CompetitionBlock; //"F" means fitness of a competing bacterium is fraction of the sum of the competition fitnesses of the two bacteria; "D" means it is binary, either more or less
	double InvadeFactor; //factor for multiplying fitness of invading bacteria (F and D Block options only)

	//migration rate
	double BactIntro;  //rate of migration (rate/person/year)

	//---------------------
	//read in strain and mge data
	//---------------------
	//strain details
	int numstrains;
	vector <string> vStrainName;
	vector <double> vStrainPrev;
	vector <double> vStrainBeta;
	vector <double> vStrainGamma;
	vector <double> vStrainR0;
	vector <double> vStrainDelta;
	vector <double> vStrainFitness;
	vector <int> vStrainIntro; //time in years at which strains introduced
	vector <bool> vStrainResistant;

	//mge details
	int nummges;
	vector <string> vMGETypeName;
	vector <double> vMGETypeBeta;
	vector <double> vMGETypeGamma;
	vector <double> vMGETypeDefence; //#this should be lookup depending on strain
	vector <bool> vMGETypeResistant;
	vector <int> vMGETypeSite;
	vector <int> vMGETypeIntro;
	vector <double> vMGETypeCompCost;
	vector <double> vMGETypeBetaCost;
	vector <double> vMGETypeGammaCost;

	//strain-mge interactions
	vector<vector<double> > vStrainMGEPrev;

	//get constant inputs
	ifstream optionsfile(OptionsFile); //this has been created manually

	if (optionsfile) { //error checking: make sure file exists
		string line;

		while (getline(optionsfile, line)) {

			istringstream iss(line);
			string s1, s2;
			if (iss >> s1 >> s2 && s1 != "#") {

				if (s1 == "HostStartNum") {
					HostStartNum = ::atoi(s2.c_str());
				}
				else if (s1 == "Years") {
					Years = ::atoi(s2.c_str());
				}
				else if (s1 == "BurnIn") {
					BurnIn = ::atoi(s2.c_str());
					if (BurnIn >= Years) {
						cerr << "Burn in " << BurnIn;
						cerr << "Years " << Years;
						cerr << "Burn-in period exceeds total run time - please correct\n";
						return 0; //terminate program
					}
				}
				else if (s1 == "StepsYr") {
					StepsYr = ::atoi(s2.c_str());
				}
				else if (s1 == "OutputSteps") {
					OutputSteps = ::atoi(s2.c_str());
				}
				else if (s1 == "GetGamma") {
					GetGamma = s2;
				}
				else if (s1 == "TargetPrev") {
					TargetPrev = ::atof(s2.c_str());
				}
				else if (s1 == "TargetUnique") {
					TargetUnique = ::atof(s2.c_str());
				}
				else if (s1 == "TargetRes") {
					TargetRes = ::atof(s2.c_str());
				}
				else if (s1 == "TargetDuration") {
					double tempdouble = ::atof(s2.c_str());
					vStrainTargetDuration.push_back(tempdouble);
					while (iss >> tempdouble) {
						vStrainTargetDuration.push_back(tempdouble);
					}
				}
				else if (s1 == "TargetProportion") {
					double tempdouble = ::atof(s2.c_str());
					vStrainTargetProportion.push_back(tempdouble);
					while (iss >> tempdouble) {
						vStrainTargetProportion.push_back(tempdouble);
					}
				}
				else if (s1 == "TargetMGEPrev") {
					double tempdouble = ::atof(s2.c_str());
					vMGETargetPrev.push_back(tempdouble);
					while (iss >> tempdouble) {
						vMGETargetPrev.push_back(tempdouble);
					}
				}
				else if (s1 == "TargetMGEChi") {
					MGETargetChi = ::atof(s2.c_str());
				}
				else if (s1 == "BactResistant") {
					OutputResistance = ::atoi(s2.c_str()); //convert to bool
				}
				else if (s1 == "BactAllSites") {
					OutputBactAllSites = ::atoi(s2.c_str()); //convert to bool
				}
				else if (s1 == "BactClearances") {
					OutputBactClearances = ::atoi(s2.c_str()); //convert to bool
				}
				else if (s1 == "BactCarriage") {
					OutputBactCarriage = ::atoi(s2.c_str()); //convert to bool
				}
				else if (s1 == "HostBactDist") {
					OutputHostBactDist = ::atoi(s2.c_str()); //convert to bool
				}
				else if (s1 == "CompetitionContinuous") {
					CompetitionCont = s2[0]; //convert to char? - only takes first character
				}
				else if (s1 == "CompetitionBlock") {
					CompetitionBlock = s2[0]; //convert to char? - only takes first character
				}
				else if (s1 == "Migration") {
					double tempdouble = ::atof(s2.c_str()); //changed to rate input
					BactIntro = 1 - exp(-tempdouble / StepsYr); //now convert to probability of intro into host in single timestep
				}
				else if (s1 == "InvasionFactor") {
					InvadeFactor = ::atof(s2.c_str());
				}
			}
		}
	}
	else {
		cerr << "Options file missing\n";
		return 0; //terminate program
	}

	ifstream strainfile(StrainFile); //generate by Python script for ELFI

	if (strainfile) { //error checking: make sure file exists

		vector<double> tempvec1;
		vector<double> tempvec2;
		string line;

		while (getline(strainfile, line)) {

			istringstream iss(line);
			string s1, s2;
			if (iss >> s1 >> s2 && s1 != "#") {
				if (s1 == "numstrains") {
					numstrains = ::atof(s2.c_str()); //atof turns string into double

					vStrainName.resize(numstrains);
					vStrainIntro.resize(numstrains);
					vStrainPrev.resize(numstrains);
					vStrainBeta.resize(numstrains);
					vStrainGamma.resize(numstrains);
					vStrainR0.resize(numstrains);
					vStrainDelta.resize(numstrains);
					vStrainFitness.resize(numstrains);
					vStrainResistant.resize(numstrains);

					//vStrainComp.resize(numstrains);
					//vStrainMGEPrev.resize(numstrains);
					tempvec1.resize(numstrains);

				}
				else if (s1 == "StrainName") {
					vStrainName[0] = s2; //no conversion, because this is string
					for (in = 1; in < numstrains; ++in) {
						iss >> vStrainName[in];
					}
				}
				else if (s1 == "StrainIntro") {
					int tempint = ::atoi(s2.c_str()); //convert to integer
					vStrainIntro[0] = tempint * StepsYr;
					for (in = 1; in < numstrains; ++in) {
						int tempint;
						iss >> tempint;
						vStrainIntro[in] = tempint * StepsYr;
					}
				}
				else if (s1 == "StrainPrev") {
					vStrainPrev[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < numstrains; ++in) {
						iss >> vStrainPrev[in]; 
					}
				}
				else if (s1 == "StrainBeta") {
					double tempdouble = ::atof(s2.c_str()); //this is annual rate, independent of population size
					vStrainBeta[0] = tempdouble; //NOT converting to probability per timestep yet
					for (in = 1; in < numstrains; ++in) {
						double tempdouble;
						iss >> tempdouble;
						vStrainBeta[in] = tempdouble;
					}
				}
				else if (s1 == "StrainGamma") {
					double tempdouble = ::atof(s2.c_str()); //this is annual rate
					vStrainGamma[0] = tempdouble; //NOT converting to probability per timestep yet
					for (in = 1; in < numstrains; ++in) {
						double tempdouble;
						iss >> tempdouble;
						vStrainGamma[in] = tempdouble;
					}
				}
				//else if (s1 == "StrainR0") {
				   // vStrainR0[0] = ::atof(s2.c_str()); //convert to double
				   // for (in = 1; in < numstrains; ++in) {
				   //	 iss >> vStrainR0[in]; //#is this converting to double???
				   // }
				//}
				else if (s1 == "StrainDelta") {
					vStrainDelta[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < numstrains; ++in) {
						iss >> vStrainDelta[in];
					}
				}
				else if (s1 == "StrainFitness") {
					vStrainFitness[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < numstrains; ++in) {
						iss >> vStrainFitness[in]; 
					}
				}
				else if (s1 == "StrainResistant") {
					vStrainResistant[0] = ::atoi(s2.c_str()); //convert to bool
					for (in = 1; in < numstrains; ++in) {
						bool tempbool;
						iss >> tempbool; //doesn't work if do it directly as for other types
						vStrainResistant[in] = tempbool;
					}
				}
				else if (s1 == "nummges") {
					nummges = ::atof(s2.c_str());

					//resize MGE type vectors - necessary to add values by index
					vMGETypeName.resize(nummges);
					vMGETypeCompCost.resize(nummges);
					vMGETypeBetaCost.resize(nummges);
					vMGETypeGammaCost.resize(nummges);
					vMGETypeIntro.resize(nummges);
					vMGETypeBeta.resize(nummges);
					vMGETypeGamma.resize(nummges);
					vMGETypeDefence.resize(nummges); //#this should be lookup depending on strain
													 /*vMGETypeChangeFitness.resize(nummges);
													 fill(vMGETypeChangeFitness.begin(), vMGETypeChangeFitness.end(), 1);
													 vMGETypeChangeGamma.resize(nummges);
													 fill(vMGETypeChangeGamma.begin(), vMGETypeChangeGamma.end(), 1);*/
					vMGETypeResistant.resize(nummges);
					vMGETypeSite.resize(nummges);
					vStrainMGEPrev.resize(numstrains, vector<double>(nummges));

					tempvec2.resize(nummges);
				}
				else if (s1 == "MGETypeName") {
					vMGETypeName[0] = s2; //no conversion, because this is string
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeName[in];
					}
				}
				else if (s1 == "MGETypeCompetitionCost") {
					vMGETypeCompCost[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeCompCost[in];
					}
				}
				else if (s1 == "MGETypeBetaCost") {
					vMGETypeBetaCost[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeBetaCost[in];
					}
				}
				else if (s1 == "MGETypeGammaCost") {
					vMGETypeGammaCost[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeGammaCost[in];
					}
				}
				else if (s1 == "MGETypeIntro") {
					vMGETypeIntro[0] = ::atoi(s2.c_str()) * StepsYr; //convert to integer
					for (in = 1; in < nummges; ++in) {
						int tempint;
						iss >> tempint;
						vMGETypeIntro[in] = tempint * StepsYr;
					}
				}
				else if (s1 == "MGETypeBeta") {
					double tempdouble = ::atof(s2.c_str()); //annual rate
					vMGETypeBeta[0] = 1 - exp(-tempdouble / StepsYr); //converted to probability per timestep
					for (in = 1; in < nummges; ++in) {
						double tempdouble;
						iss >> tempdouble;
						vMGETypeBeta[in] = 1 - exp(-tempdouble / StepsYr);
					}
				}
				else if (s1 == "MGETypeDefence") {
					vMGETypeDefence[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeDefence[in]; 
					}
				}
				else if (s1 == "MGETypeGamma") { 
					double tempdouble = ::atof(s2.c_str()); //this is annual rate
					vMGETypeGamma[0] = 1 - exp(-tempdouble / StepsYr); //converted to probability per timestep
					for (in = 1; in < nummges; ++in) {
						double tempdouble;
						iss >> tempdouble;
						vMGETypeGamma[in] = 1 - exp(-tempdouble / StepsYr);
					}
				}
				else if (s1 == "MGETypeResistant") {
					vMGETypeResistant[0] = ::atoi(s2.c_str()); //convert to bool
					for (in = 1; in < nummges; ++in) {
						bool tempbool;
						iss >> tempbool; //doesn't work if do it directly as for other types
						vMGETypeResistant[in] = tempbool;
					}
				}
				else if (s1 == "MGETypeSite") {
					vMGETypeSite[0] = ::atof(s2.c_str()); //convert to double
					for (in = 1; in < nummges; ++in) {
						iss >> vMGETypeSite[in]; 
					}
				}
				else if (s1 == "StrainMGEPrev") {
					for (jn = 0; jn < nummges; ++jn) {
						if (s2 == vMGETypeName[jn]) {
							for (in = 0; in < numstrains; ++in) {
								iss >> vStrainMGEPrev[in][jn]; 
							}
						}
					}
				}
			}
		}

		if (numstrains != vStrainTargetDuration.size()) {
			cerr << "Number of strains in " << StrainFile << " does not match number of strain durations given in " << OptionsFile << " - please correct\n";
			return 0; //terminate program
		}
	}
	else {
		cerr << "Strain file missing: " << StrainFile << "\n";
		return 0; //terminate program
	}

	if (!((CompetitionCont == 'N') | (CompetitionCont == 'L') | (CompetitionBlock == 'D') | (CompetitionBlock == 'F'))) { //if no competition specified
		return 0; //need error message
	}

	//convert timescales
	const double HostDeath = 1 - exp(-HostDeathYr / StepsYr); //individual's probability of death in a single time step
	const double HostBirth = 1 - exp(-HostBirthYr / StepsYr); //individual's probability of birth in a single time step
	const double HostABc = 1 - exp(-HostABcYr / StepsYr); //convert abc from weekly treatments/yr per host to probability of treatments/timestep
	int Timesteps = round(Years * StepsYr);
	int Burnsteps = round(BurnIn * StepsYr);

	//--------------------------
	//loop through Strains 
	//#can't concatenate strings to dynamically name variables, so use vectors...

	//create strains
	vector<Strain> Strains; //#note is not boost::ptr_vector because can't push_back(make_shared)
	Strains.reserve(numstrains);

	if (GetGamma == "proportional") {// did have other options here, now deprecated
		for (in = 0; in < numstrains; ++in) {
			vStrainBeta[in] = ParmsR0 * vStrainDelta[in]; //vStrainR0[in] * vStrainDelta[in]; 
		}
	}

	for (in = 0; in < numstrains; ++in) {
		vector<double> tempvec(numstrains, ParmsInterPcomp);
		tempvec[in] = ParmsIntraPcomp;

		Strain newstrain(vStrainName[in], tempvec, vStrainPrev[in], 1 - exp(-vStrainBeta[in] / (HostStartNum * StepsYr)), 1 - exp(-vStrainGamma[in] / StepsYr), vStrainFitness[in], vStrainResistant[in], vStrainIntro[in]);
		newstrain.StrainIndex = Strains.size();
		Strains.push_back(newstrain);
	}

	if (GetGamma == "proportional") {
		//use this for error-checking
		string summaryname = Page + "_InputSummary.csv";
		ifstream summaryfile(summaryname);
		if (!summaryfile) //if this file doesn't exist
		{

			ofstream summary;
			summary.open(summaryname);
			summary << "StrainName,StrainBeta,StrainGamma\n";
			summary.close();
		}

		ofstream summary;
		summary.open(summaryname, ios::app); //appends to existing file
		for (in = 0; in < numstrains; ++in) {
			summary << Strains[in].StrainName << "," << -log(1 - Strains[in].StrainBeta) * StepsYr * HostStartNum << "," << -log(1 - Strains[in].StrainGamma) * StepsYr << "\n";
		}
		summary.close();
	}

	//--------------------------
	//loop through MGETypes 

	//create MGE types
	vector<MGEType> MGETypes;
	MGETypes.reserve(nummges);
	vector<int> MGESites;

	for (in = 0; in < nummges; ++in) {

		double MGEChangeFitness = 1;
		double MGEChangeGamma = 1;
		double MGEChangeBeta = 1;

		//add competition cost
		MGEChangeFitness = 1 - ParmsMGECompetitionCost * vMGETypeCompCost[in]; //competition cost from R wrapper is multipled by factor in InputStrains.txt

		//add beta cost
		MGEChangeBeta = 1 + ParmsMGEBetaCost * vMGETypeBetaCost[in]; //changed equation for way MGE changes bacterial gamma

		//add gamma cost
		MGEChangeGamma = 1 + ParmsMGEGammaCost * vMGETypeGammaCost[in]; //changed equation for way MGE changes bacterial gamma

		MGEType newmge(vMGETypeName[in], vMGETypeBeta[in], vMGETypeDefence[in], MGEChangeFitness, MGEChangeGamma, MGEChangeBeta, vMGETypeSite[in], vMGETypeIntro[in], vMGETypeResistant[in]);
		newmge.MGETypeIndex = MGETypes.size();
		MGETypes.push_back(newmge);

		//see if MGEs already exist that use this site
		auto result1 = find(begin(MGESites), end(MGESites), newmge.MGETypeSite);
		if (result1 == end(MGESites)) {

			//if they don't, add site to list of sites
			MGESites.push_back(newmge.MGETypeSite);
		}
	}

	//make global variable so don't have to keep looking up 
	int numMGESites = MGESites.size();

	//--------------------------
	//Delete read-in vectors
	vector<string>().swap(vStrainName);

	vector<double>().swap(vStrainPrev);

	vector<double>().swap(vStrainGamma);
	vector<double>().swap(vStrainFitness);
	vector<int>().swap(vStrainIntro);
	vector<bool>().swap(vStrainResistant);

	vector<string>().swap(vMGETypeName);
	vector<double>().swap(vMGETypeCompCost);
	vector<double>().swap(vMGETypeBetaCost);
	vector<double>().swap(vMGETypeGammaCost);
	vector<int>().swap(vMGETypeIntro);
	vector<double>().swap(vMGETypeBeta);
	vector<double>().swap(vMGETypeGamma);
	vector<double>().swap(vMGETypeDefence); 
	vector<bool>().swap(vMGETypeResistant);
	vector<int>().swap(vMGETypeSite);

	//---------------------
	//set up monitoring

		//cout << chrono::high_resolution_clock::period::den << endl;
		/*##profiling
		auto start_time = chrono::high_resolution_clock::now();

		long setup_time = 0;
		long run_time = 0;
		long write_time = 0;
		long cum_time = 0;

		long discard_time = 0;
		long rand_time = 0; //rand_max on every host
		long kill_time = 0;
		long getinf_time = 0; //getinf on every live host
		long clear_time = 0;
		long abc_time = 0;
		long immune_time = 0;
		long comp_time = 0;
		long mge_time = 0;
		long infect_time = 0;
		long count_time = 0;
		long chain_time = 0; //long index look-up of MGEs to add to counter
		long point_time = 0; //iterator look-up of MGEs to add to counter
		long binom_time = 0;
		long clone_time = 0;
		long infecting_time = 0;
		long demographic_time = 0;
		*/

		//threaded code
		//int cores = omp_get_max_threads(); //gets number of cores
	//	int cores = 4;
	//	omp_set_num_threads(cores); //states how many threads to use
	//	//cout << "max threads: " << cores << '\n';
	//	
	//	//int cores = set_NCores();
	//
	//	div_t divresult = div(Runs, cores); //divresult.quot is answer, divresult.rem is remainder
	//	vector<int> coreruns(cores + 1, 0); // create vector of number of runs per core
	//	for (in = 0; in < cores + 1; ++in){
	//		coreruns[in + 1] = coreruns[in] + divresult.quot;
	//		if (in < divresult.rem){
	//			++coreruns[in + 1];
	//		}
	//	}
	//
	//#pragma omp parallel for private(thread) schedule(static,1)
	//for (thread = 0; thread<cores; ++thread){
	//

	//things for loops
	unsigned int i; //for for loops
	unsigned int j; //for for loops
	int h; //must be UNSIGNED so loop through hosts works correctly
	int b; //must be unsigned for backwards loop
	int b2; //must be unsigned for backwards loop
	int m; //must be unsigned for backwards loop
	int m2; //must be unsigned for backwards loop
	unsigned int k;
	double r; //for random numbers
	double r2;

	//check if output file exists, otherwise create
	string outputname = Page + string("_") + to_string(JobNumber) + string("_OutputCounts.csv");

	ifstream file(outputname);
	if (!file) //if this file doesn't exist
	{
		//create output file and headers
		ofstream OutputCounts;
		//	OutputCounts.open(Folder + "//" + to_string(Page) + "_OutputCounts.csv");
		OutputCounts.open(outputname);
		OutputCounts << "count,years,host_startnum,intrapcomp, interpcomp, ICECompetitionCost, PhageGammaCost, PlasmidBetaCost, host_abc_yr,run,R0,";

		//now changing 0,1,2, iterator to time elapsed
		for (unsigned int t = 0; t < (Timesteps - 1); t += OutputSteps) { //#MAKE SURE THIS MATCHES OUTPUTS
			OutputCounts << t << ","; //adding 1 for time elapsed
		}
		OutputCounts << Timesteps;  //last item doesn't get comma. Timesteps -1 +1.

									//now close output file
		OutputCounts.close();
	}

	//this is for ELFI
	//check if summary statistic file exists, otherwise create
	string statname = "OutputStats.csv";
	ifstream file2(statname);
	if (!file2) //if this file doesn't exist
	{
		//create output file and headers
		ofstream OutputStats;
		OutputStats.open(statname);
		OutputStats << "file, R0, intrapcomp, interpcomp, ICECompetitionCost, PhageGammaCost, PlasmidBetaCost,strain_survival, SS_duration,prop_inf, SS_infection, prop_unique, SS_unique, prop_res, SS_res,SS_MGE_prev, SS_MGE_chi, summary_statistic, seconds";
		OutputStats.close();
	}

	for (int run = 0; run < Runs; ++run) {

		double Timestop = Timesteps; //this is when model stops running

		//summary statistics
		double strainsurvive = numstrains;
		double strainstat1 = 0.0;
		double strainstat2 = 0.0;
		double tempstat1;
		double tempstat2;
		double resstat;

		//rebuild counters
		//counter vectors for each strain over time 
		vector<vector<int> > vCountStrain(numstrains, vector<int>(Timesteps, 0));
		vector<vector<long long unsigned int> > vStrainClears(numstrains, vector<long long unsigned int>(Timesteps, 0));
		vector<vector<long long unsigned int> > vStrainCarriage(numstrains, vector<long long unsigned int>(Timesteps, 0));
		vector<vector<int> > vCountStrainRes(numstrains, vector<int>(Timesteps, 0)); //counter for resistant strains 

		//3D array with strains, mges, timesteps
		vector<vector<vector<int> > > vvStrainMGE(numstrains, vector<vector<int> >(nummges, vector<int>(Timesteps, 0)));
		vector<vector<vector<int> > > vvStrainMGERes(numstrains, vector<vector<int> >(nummges, vector<int>(Timesteps, 0)));

		//number of hosts infected with each strain
		vector<vector<int> > vCountHostStrain(numstrains, vector<int>(Timesteps, 0));

		//MGEs
		//counter for no. MGEs of each type in existence at each timestep
		vector<vector<int> > vCountMGE(nummges, vector<int>(Timesteps, 0));
		vector<vector<int> > vCountMGERes(nummges, vector<int>(Timesteps, 0)); //count number resistant MGEs ##error-checking

		vector<vector<vector<vector<int> > > > vvvSuccessHGT(nummges, vector<vector<vector<int> > >(numstrains, vector<vector<int> >(numstrains, vector<int>(Timesteps, 0))));
		vector<vector<int> > vFailedHGT(nummges, vector<int>(Timesteps, 0));

		//--------------------------
		//create hosts
		vector<Host> Hosts(HostStartNum); // vector of Host called Hosts with HostStartNum elements 
		Hosts.reserve(round(HostStartNum * 1.1)); //note that push_back on full vector doubles allocated memory

//--------------------------
// Get ready for timestep loop

		//vectors to store outputs 
		vector<unsigned int> CountHostsAlive(Timesteps, 0); //this will store number of hosts alive at the end of each day
		vector<unsigned int> CountHostsInf(Timesteps, 0); //this will store number of hosts infected with anything at the end of each day

		vector<int> CountBact(Timesteps, 0);
		vector<int> CountBactRes(Timesteps, 0);
		vector<int> CountMGE(Timesteps, 0);
		vector<int> CountMGERes(Timesteps, 0);
		vector<vector<int> > CountResSens(3, vector<int>(Timesteps, 0)); //each timestep has vector with counts of hosts infected with resistant only, sensitive only, or both 

		vector<vector<int> > CountGetInf(6, vector<int>(Timesteps, 0)); //count number of hosts infected with 0,1,2,... bacteria
		vector<vector<int> > CountGetInfRes(6, vector<int>(Timesteps, 0)); //count number of hosts infected with 0,1,2,... resistant bacteria
		vector<vector<int> > CountGetInfUnique(5, vector<int>(Timesteps, 0)); //count number of hosts infected with 1,2,..,5+ unique strains
		vector<vector<int> > vCountSingleStrain(numstrains, vector<int>(Timesteps, 0));

		//vector of integers of host position for later killing/infecting
		vector<int> HostNewInfect; HostNewInfect.reserve(HostStartNum);
		vector<int> HostNewDead; HostNewDead.reserve(HostDeath * HostStartNum * 2);
		vector<int> BactNewInfectMGE; BactNewInfectMGE.reserve(10);

		//ptr_vector to put bacterial clones into 
		vector<MGEType*> MGENewInfect;

		//counters 
		vector<int> clones(Timesteps, 0);
		vector<int> infecteds(numstrains, 0); //vector for host infected with each strain 
		vector<int> cleardead(Timesteps, 0);
		vector<int> clearabc(Timesteps, 0);
		vector<int> cleargamma(Timesteps, 0);
		vector<int> clearcomp(Timesteps, 0);
		vector<int> clearcompW(Timesteps, 0);
		vector<int> clearcompB(Timesteps, 0);
		vector<int> acceptcompRR(Timesteps, 0);
		vector<int> acceptcompRS(Timesteps, 0);
		vector<int> acceptcompSR(Timesteps, 0);
		vector<int> acceptcompSS(Timesteps, 0);
		vector<int> allsites(Timesteps, 0);

		//new measurement of res vs sens clearance
		//4 vectors in each of these, in reverse alphabetical order
		//first index is resistance status of winner, 2nd index is resistance status of loser
		vector<vector<vector<int> > > vvclearWithinStatus(2, vector<vector<int> >(2, vector<int>(Timesteps, 0)));
		vector<vector<vector<int> > > vvclearBetweenStatus(2, vector<vector<int> >(2, vector<int>(Timesteps, 0)));

		//counting the other clearances
		vector<vector<int> > vcleargammastrain(numstrains, vector<int>(Timesteps, 0));
		vector<vector<int> > vclearcompstrain(numstrains, vector<int>(Timesteps, 0));
		vector<vector<int> > vcleardeathstrain(numstrains, vector<int>(Timesteps, 0));
		vector<vector<int> > vclonestrain(numstrains, vector<int>(Timesteps, 0));

		vector<int> noseed(Timesteps, 0); //seedings that don't happen, due to antibiotics
		vector<int> noinfect(Timesteps, 0); //infections that don't happen, due to antibiotics
		//#I'm actually going to remove this - it is complicating the scheduling

		vector<int> trycomp(Timesteps, 0); //number of calls to competition function
		vector<int> nocomp(Timesteps, 0); //number of competition events that do not result in clearance of bacteria
		vector<int> residentlosescomp(Timesteps, 0); //competition events that resident bacteria loses (i.e. do not result in clearance event)

		//#problems with the below?
		vector<int> tryinfectcarrier(Timesteps, 0);
		vector<int> succeedinfectcarrier(Timesteps, 0);

		int failedinfection = 0;

		//initialise vector of antibiotic consumption
		//generate vector with booleans for hosts receiving antibiotics or not
		vector <bool> HostConsumeNext(HostStartNum); //zeroes by default
		for (i = 0; i < HostStartNum; ++i) {
			if (ranf_mt(thread) < HostABc) {
				HostConsumeNext[i] = 1; //change appropriate proportion to 1's
			}
		}

		///////////////
		//Loop through timesteps
		unsigned int t; //more efficient to declare before loop
		for (t = 0; t < Timesteps; ++t) {

			int numHosts = Hosts.size();

			//check we still have hosts
			if (numHosts == 0) {
				//	cout << "no more hosts left" << '\n';
				Timestop = t;
				break; //finish running model and go to outputs
			}

			//seed strains
			for (i = 0; i < numstrains; ++i) {
				if (Strains[i].StrainIntro == t) { //if strains introduced at this timestep
					//create vector of ints to encode Hosts
					vector<int> HostMarker(numHosts);
					iota(HostMarker.begin(), HostMarker.end(), 0); //fills with increasing numbers from 0 to numHosts

					//now get how many hosts we are infecting
					int howmany = ignpoi_mt(Strains[i].StrainPrev * numHosts, thread); //Poisson approximation to binomial 

					if (howmany > numHosts) {
						howmany = numHosts; //error-catching to stop 'infecting' more hosts than exist
					}

					//infect each of these hosts
					if (howmany > 0) { //avoid shuffle for no purpose
						random_shuffle(HostMarker.begin(), HostMarker.end()); //shuffle marker vector
						for (j = 0; j < howmany; ++j) {

							//choose host and add bacterium
							Hosts[HostMarker[j]].add_bact(Strains[i]);
							for (k = 0; k < numMGESites; ++k) {
								Hosts[HostMarker[j]].WithBact.back().WithMGE.push_back(nullptr);
							}

							//add MGEs if they are introduced at this time/before this time
							for (k = 0; k < nummges; ++k) {
								if ((MGETypes[k].MGETypeIntro <= t) && (ranf_mt(thread) < vStrainMGEPrev[i][k])) {
									Hosts[HostMarker[j]].WithBact.back().add_MGE(MGETypes[k]);
								}
							}

							//now delete sensitive bacteria if host on antibiotics (have to do it like this so don't erase resistant-MGE-bearing introductions before they get the MGE!
							//#may need to remove this step?
							if ((HostConsumeNext[HostMarker[j]] == 1) && (Hosts[HostMarker[j]].WithBact.back().BactResistant == 0)) {
								Hosts[HostMarker[j]].WithBact.pop_back();
								++noseed[t];
							}
						} //end seeding these hosts
					} //end howmany>0
				} //end strain introduction at this timestep
			} //end of strain introduction

			//COMPETITION HERE FOR ACQUISITION

			 //seed MGEs	- what about sensitive bacteria and hosts consuming antibiotics? Kill first? YES.
			//Actually, should go through hosts one by one, add mges if appropriate, otherwise have to cycle through them for each MGE seeded at this timestep
			for (k = 0; k < nummges; ++k) {
				if (MGETypes[k].MGETypeIntro == t) {
					//cerr << "Introducing MGE " << MGETypes[k].MGETypeName << " at time = " << t << "\n";
					for (h = 0; h < numHosts; ++h) {
						for (b = 0; b < Hosts[h].WithBact.size(); ++b) {
							int s = Hosts[h].WithBact[b].OfStrain->StrainIndex;

							//don't add MGEs to strains that have only been added in this timestep - already been added
							if ((t != Strains[s].StrainIntro) &&
								(ranf_mt(thread) < vStrainMGEPrev[s][k])) {
								Hosts[h].WithBact[b].add_MGE(MGETypes[k]);
							}
						}
					}
				}
			}

			//check we still have bacteria 
			if (t > 0 && CountBact[t - 1] == 0) {
				cerr << "no more bacteria left" << '\n';

				Timestop = t;

				//add in bit that gets missed with jump to printoutputs
				//note that host_unique_1 only refers to hosts with more than one isolate, but one unique strain (i.e. dual/multi infection same strain) 
				//add together first!
				for (int t2 = 0; t2 < Timestop; ++t2) {
					CountGetInfUnique[0][t2] += CountGetInf[1][t2]; //now unique_1 includes hosts with single bacterium
				}

				goto printoutputs; //finish running model and go to outputs
			}

			//this is migration 2023
			//pick hosts
			if ((BactIntro > 0) & (t >= Burnsteps)) {
				for (i = 0; i < numstrains; ++i) {

					if ((t > Strains[i].StrainIntro) ||(Strains[i].StrainIntro==999*StepsYr)) { //only migrate after strain introduced - or if fudge '999' - this enables new strains to migrate in
						//THIS IS COPIED FROM STRAIN SEEDING
						//create vector of ints to encode Hosts
						vector<int> HostMarker(numHosts);
						iota(HostMarker.begin(), HostMarker.end(), 0); //fills with increasing numbers from 0 to numHosts

						//now get how many hosts we are infecting with this strain
						int howmany = ignpoi_mt(Strains[i].StrainPrev * BactIntro * numHosts, thread); //Poisson approximation to binomial 

						if (howmany > numHosts) {
							howmany = numHosts; //error-catching to stop 'infecting' more hosts than exist
						}

						//infect each of these hosts
						if (howmany > 0) { //avoid shuffle for no purpose
							random_shuffle(HostMarker.begin(), HostMarker.end()); //shuffle marker vector
							for (j = 0; j < howmany; ++j) {

								//choose host and add bacterium
								Hosts[HostMarker[j]].add_bact(Strains[i]);
								for (k = 0; k < numMGESites; ++k) {
									Hosts[HostMarker[j]].WithBact.back().WithMGE.push_back(nullptr);
								}

								//add MGEs if they are introduced at this time/before this time
								for (k = 0; k < nummges; ++k) {
									if ((MGETypes[k].MGETypeIntro <= t) && (ranf_mt(thread) < vStrainMGEPrev[i][k])) {
										Hosts[HostMarker[j]].WithBact.back().add_MGE(MGETypes[k]);
									}
								}

								//now delete sensitive bacteria if host on antibiotics (have to do it like this so don't erase resistant-MGE-bearing introductions before they get the MGE!
								if ((HostConsumeNext[HostMarker[j]] == 1) && (Hosts[HostMarker[j]].WithBact.back().BactResistant == 0)) {
									Hosts[HostMarker[j]].WithBact.pop_back();
									++noseed[t];
								}
							} //end seeding these hosts
						} //end howmany>0
					} //end of if statement 
				} //end this strain's migration
			} //end migration

			HostNewInfect.clear();
			HostNewDead.clear();

			fill(infecteds.begin(), infecteds.end(), 0);

			//ORIGINAL
			vector <bool> HostConsumeNow(numHosts);
			for (i = 0; i < numHosts; ++i) {
				if (ranf_mt(thread) < HostABc) {
					HostConsumeNow[i] = 1;
				}
			}
			//now swap them
			HostConsumeNext.swap(HostConsumeNow); //HostConsumeNext should be permanently changed

			for (h = 0; h < numHosts; ++h) {

				if (ranf_mt(thread) < HostDeath) { //kill host
					HostNewDead.push_back(h); //add to list of integers
				}
				else //not dead!
				{

					if (Hosts[h].WithBact.size() == 0) { //host has no bacteria 
						//update no. hosts with no bact
						++CountGetInf[0][t];

						//update no. hosts with no resistant bact
						++CountGetInfRes[0][t];
					}
					else
					{ //host has bacteria
					// THIS IS CLEANING UP END OF LAST TIMESTEP

						//-------------------------

						//loop through bacteria and clear by immune response only
						for (b = Hosts[h].WithBact.size() - 1; b > -1; b--) {
							//clear by immune response 
							if (ranf_mt(thread) < Hosts[h].WithBact[b].BactGamma) {
								++cleargamma[t];
								++vcleargammastrain[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];

								++vStrainClears[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
								vStrainCarriage[Hosts[h].WithBact[b].OfStrain->StrainIndex][t] += (t - Hosts[h].WithBact[b].BactCarriage);
								Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b);
							}
						}

						//antibiotic clearance
						//loop through bacteria in host and clear by antibiotics
						if (HostConsumeNow[h] == 1) {
							for (b = Hosts[h].WithBact.size() - 1; b > -1; b--) {

								//clear by antibiotics 
								if (Hosts[h].WithBact[b].BactResistant == 0) {

									//update clear and carriage
									++vStrainClears[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
									vStrainCarriage[Hosts[h].WithBact[b].OfStrain->StrainIndex][t] += (t - (Hosts[h].WithBact[b].BactCarriage));
									Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b);
									++clearabc[t];
								}
							}
						}

						//check how many bacteria left after clearance
						if (Hosts[h].WithBact.size() == 0) {
							//update no. hosts with no bact
							++CountGetInf[0][t];

							//update no. hosts with no resistant bact
							++CountGetInfRes[0][t];

							//	cout << "updated zero counters\n";
						}
						else {

							//do competition loop
							if (CompetitionCont == 'N') {
								int tempcount = Hosts[h].WithBact.size();

								//loop through bacteria from most recently added to least recent
								for (b = tempcount - 1; b > 0; b--) { //not 0 as then nothing to compare

									//get next bacteria in line
									for (b2 = b - 1; b2 > -1; b2--) {
										//cout << "b=" << b << ", b2=" << b2 << '\n';

										++trycomp[t]; //count pairwise interactions

										//get index of strain of second bact
										int index = Hosts[h].WithBact[b2].OfStrain->StrainIndex;
										//cout << "index = " << index; //have to do this here rather than inside function because of hierarchy

										//see if competition happens
										switch (Hosts[h].WithBact[b].competition_cont(Hosts[h].WithBact[b2], index, thread)) {
										case 0:
											++nocomp[t];
											break; //goes to next b2
										case 1:
											//element b2 is deleted

											//test if within- or between-strain competition
											if (Hosts[h].WithBact[b].OfStrain->StrainIndex == index) {
												++vvclearWithinStatus[Hosts[h].WithBact[b].BactResistant][Hosts[h].WithBact[b2].BactResistant][t];
											}
											else {
												++vvclearBetweenStatus[Hosts[h].WithBact[b].BactResistant][Hosts[h].WithBact[b2].BactResistant][t];
											}

											++vclearcompstrain[index][t];
											++vStrainClears[Hosts[h].WithBact[b2].OfStrain->StrainIndex][t];
											vStrainCarriage[Hosts[h].WithBact[b2].OfStrain->StrainIndex][t] += (t - Hosts[h].WithBact[b].BactCarriage);
											Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b2); //erase element b2. Now b points to nothing.

											++clearcomp[t];
											--b; //reduce b by one
											break; //goes to next b2
										case 2: //kill b

											//test if within- or between-strain competition
											if (Hosts[h].WithBact[b].OfStrain->StrainIndex == index) {
												++vvclearWithinStatus[Hosts[h].WithBact[b2].BactResistant][Hosts[h].WithBact[b].BactResistant][t];
											}
											else {
												++vvclearBetweenStatus[Hosts[h].WithBact[b2].BactResistant][Hosts[h].WithBact[b].BactResistant][t];
											}

											++vclearcompstrain[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
											++vStrainClears[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
											vStrainCarriage[Hosts[h].WithBact[b].OfStrain->StrainIndex][t] += (t - Hosts[h].WithBact[b].BactCarriage);
											Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b); //erase element b

											++clearcomp[t];
											goto newb; //go to next b - possible because limits on b hard-coded

										default:
											cout << '\n' << "error in competition function";
										}

									} //end of b2
								newb:; //goes to next b

								} //end of b and competition loop
							} //end of optoin for contnuous, normal competition
							else if (CompetitionCont == 'L') {
								//now do competition loop
								int tempcount = Hosts[h].WithBact.size();

								//loop through bacteria from most recently added to least recent
								for (b = tempcount - 1; b > 0; b--) { //not 0 as then nothing to compare

																	  //get next bacteria in line
									for (b2 = b - 1; b2 > -1; b2--) {
										//cout << "b=" << b << ", b2=" << b2 << '\n';

										++trycomp[t];

										//get index of strain of second bact
										int index = Hosts[h].WithBact[b2].OfStrain->StrainIndex;
										//cout << "index = " << index; //have to do this here rather than inside function because of hierarchy

										//see if competition happens
										switch (Hosts[h].WithBact[b].competition_lastin(Hosts[h].WithBact[b2], index, thread, b, b2)) {
										case 0:
											++nocomp[t];
											break; //goes to next b2
										case 1:

											//test if within- or between-strain competition
											if (Hosts[h].WithBact[b].OfStrain->StrainIndex == index) {
												++vvclearWithinStatus[Hosts[h].WithBact[b].BactResistant][Hosts[h].WithBact[b2].BactResistant][t];
											}
											else {
												++vvclearBetweenStatus[Hosts[h].WithBact[b].BactResistant][Hosts[h].WithBact[b2].BactResistant][t];
											}

											++vclearcompstrain[index][t];
											++vStrainClears[Hosts[h].WithBact[b2].OfStrain->StrainIndex][t];
											vStrainCarriage[Hosts[h].WithBact[b2].OfStrain->StrainIndex][t] += (t - Hosts[h].WithBact[b].BactCarriage);
											Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b2); //erase element b2. Now b points to nothing.

											++clearcomp[t];
											--b; //reduce b by one
											break; //goes to next b2
										case 2: //kill b
											//test if within- or between-strain competition
											if (Hosts[h].WithBact[b].OfStrain->StrainIndex == index) {
												++vvclearWithinStatus[Hosts[h].WithBact[b2].BactResistant][Hosts[h].WithBact[b].BactResistant][t];
											}
											else {
												++vvclearBetweenStatus[Hosts[h].WithBact[b2].BactResistant][Hosts[h].WithBact[b].BactResistant][t];
											}

											++vclearcompstrain[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
											++vStrainClears[Hosts[h].WithBact[b].OfStrain->StrainIndex][t];
											vStrainCarriage[Hosts[h].WithBact[b].OfStrain->StrainIndex][t] += (t - Hosts[h].WithBact[b].BactCarriage);
											Hosts[h].WithBact.erase(Hosts[h].WithBact.begin() + b); //erase element b

											++clearcomp[t];
											goto newbL; //go to next b - possible because limits on b hard-coded

										default:
											cout << '\n' << "error in competition function";
										}

									} //end of b2
								newbL:; //goes to next b

								} //end of b and competition loop
							} //else{} //end of option for continuous competition = "L"

						//create temporary vectors
							vector<bool> tempinf(numstrains, 0); //temporary vector to get infection status of host for each strain 
							vector<int> tempressens(2, 0); //temporary vector to get whether host infected with resistant or sensitive

							//increase counter for no. hosts infected with anything
							++CountHostsInf[t];

							//THIS IS EFFECTIVELY THE END OF THE LAST TIMESTEP

							//-------------------------

							//if single bacteria, count resistance and MGE status
							if (Hosts[h].WithBact.size() == 1) {

								//count resistance status here before MGE loop
								int tempindex = Hosts[h].WithBact[0].OfStrain->StrainIndex;

								//add to no coinfection
								++vCountSingleStrain[tempindex][t];

								if (Hosts[h].WithBact[0].BactResistant == 1) {
									++vCountStrainRes[tempindex][t];
									++CountBactRes[t];

									++tempressens[0]; //update counter for host being infected with resistant bact
								}
								else {
									++tempressens[1]; //update bool for host being infected with sensitive bact
								}

								//temporary counter for dual sites
								int tempsites = 0;

								//add any MGEs to counter, but don't do MGE transfer
								for (m = 0; m < numMGESites; ++m) {
									if (Hosts[h].WithBact[0].WithMGE[m]) {
										++tempsites;

										//update counts
										++CountMGE[t];
										++vCountMGE[Hosts[h].WithBact[0].WithMGE[m]->MGETypeIndex][t];
										++vvStrainMGE[Hosts[h].WithBact[0].OfStrain->StrainIndex][Hosts[h].WithBact[0].WithMGE[m]->MGETypeIndex][t];
										if (Hosts[h].WithBact[0].WithMGE[m]->MGETypeResistant == 1) { //### why not MGETypes[i] ###CHANGED TO 1
											++CountMGERes[t];
											++vvStrainMGERes[Hosts[h].WithBact[0].OfStrain->StrainIndex][Hosts[h].WithBact[0].WithMGE[m]->MGETypeIndex][t];
										}
									}
								}

								//update counter for all MGE sites being filled
								if (tempsites == numMGESites) {
									++allsites[t];
								}
							} //end of hosts infected with single bacteria after competition

							//if more than one bacteria, count status, allow MGE transfer
							if (Hosts[h].WithBact.size() > 1) { //there is more than 1 bact: allow MGE transfer
								//shuffle bacteria to ensure random for MGE transfer and future competition loop
								//Only do this if don't want order!
								if (CompetitionCont == 'N') {
									random_shuffle(Hosts[h].WithBact.begin(), Hosts[h].WithBact.end());
								}

								//loop through bacteria
								int tempsize = Hosts[h].WithBact.size();
								vector<bool> tempunique(numstrains, 0); //temporary vector with places for each strain

								//cout << "Host " << h << ": ";
								for (b = 0; b < tempsize; ++b) {

									int tempindex = Hosts[h].WithBact[b].OfStrain->StrainIndex;
									tempunique[tempindex] = 1;

									if (Hosts[h].WithBact[b].BactResistant == 1) {
										++vCountStrainRes[tempindex][t];
										++CountBactRes[t];
										++tempressens[0]; //update counter for host being infected with resistant bact
									}
									else {
										++tempressens[1]; //update bool for host being infected with sensitive bact
									}

									int tempsites = 0;

									//iterate through MGEs in bacteria
									for (m = 0; m < numMGESites; ++m) {
										if (Hosts[h].WithBact[b].WithMGE[m]) {
											++tempsites;

											//update counts 
											++CountMGE[t];

											++vCountMGE[Hosts[h].WithBact[b].WithMGE[m]->MGETypeIndex][t];
											++vvStrainMGE[Hosts[h].WithBact[b].OfStrain->StrainIndex][Hosts[h].WithBact[b].WithMGE[m]->MGETypeIndex][t];
											if (Hosts[h].WithBact[b].WithMGE[m]->MGETypeResistant == 1) { //### why not MGETypes[i] ###CHANGED TO 1
												++CountMGERes[t];
												++vvStrainMGERes[Hosts[h].WithBact[b].OfStrain->StrainIndex][Hosts[h].WithBact[b].WithMGE[m]->MGETypeIndex][t];
											}

											//choose candidate bacterium for MGE infection - not doing as binomial because of small numbers
											for (b2 = 0; b2 < tempsize; ++b2) { //not choosing recipient in random order, but will infect in random order

												if (b2 != b) { //if b2 is different from b: stops self-infection
													if (ranf_mt(thread) < Hosts[h].WithBact[b].WithMGE[m]->MGETypeBeta) //and MGE transfer is attempted

												//	(ranf_mt(thread) < Hosts[h].WithBact[b].WithMGE[m]->MGETypeDefence)) //and MGE is not beaten by bacterial defences...
													{
														// transfer occurs!
														BactNewInfectMGE.push_back(b2); //add bacterial index to vector
														MGENewInfect.push_back(Hosts[h].WithBact[b].WithMGE[m]); //clone MGE and put in vector of clones
														++vvvSuccessHGT[Hosts[h].WithBact[b].WithMGE[m]->MGETypeIndex][Hosts[h].WithBact[b].OfStrain->StrainIndex][Hosts[h].WithBact[b2].OfStrain->StrainIndex][t];
													}

												}
											} //finished with this b2
										} //finished with non-null m
									} //finished with this m

									  //update counter for all MGE sites being filled
									if (tempsites == numMGESites) {
										++allsites[t];
									}

								} //finished with this b
								//this can be loop for Sonja option  - no shuffle

								//now implement MGE transfer 
								//choose element at random, swap then pop_back
								while (BactNewInfectMGE.size() > 0) {
									int r = floor(ranf_mt(thread) * (BactNewInfectMGE.size()));//Wes' code: will give uniform integer distribution between 0 and x-1 inclusive 
									Hosts[h].WithBact[BactNewInfectMGE[r]].add_MGE(*MGENewInfect[r]); //giving MGEType pointed to from MGENewInfect
									//#####need something here to show failed transfer - make add_MGE not function?
									swap(BactNewInfectMGE[r], BactNewInfectMGE.back());
									BactNewInfectMGE.pop_back();
									swap(MGENewInfect[r], MGENewInfect.back());
									MGENewInfect.pop_back();
								}

								//clear vectors for transfer 
								vector<int>().swap(BactNewInfectMGE);
								vector<MGEType*>().swap(MGENewInfect);

								//now collate unique infections
								//set markers
								int tempsum = 0;
								int tempindex = 0;

								//loop through tempunique vector and add value to sum
								//#could just change this to accumulate vector?!
								while ((tempsum < 5) && (tempindex < numstrains)) {
									tempsum += tempunique[tempindex];
									++tempindex;
								}
								++CountGetInfUnique[tempsum - 1][t];
							} //end of more than one bacteria after competition

						//-------------------------
						//now get bacterial infections
							boost::ptr_vector<Bacterium>::iterator bit;
							auto start = Hosts[h].WithBact.begin();
							auto end = Hosts[h].WithBact.end();
							for (bit = start; bit < end; ++bit) {
								++CountBact[t];

								int tempindex = bit->OfStrain->StrainIndex;
								tempinf[tempindex] = 1;
								++vCountStrain[tempindex][t];

								//adding in denominator here so Beta is proportional to number of bacteria in host
								//i.e. total Beta determined by host not by bacteria
								int howmany = ignpoi_mt(bit->BactBeta * (numHosts - 1) / Hosts[h].WithBact.size(), thread); //Poisson approximation to binomial 
								if (howmany > (numHosts - 1)) {
									howmany = numHosts - 1; //error-catching to stop 'infecting' more hosts than exist
								}

								for (i = 0; i < howmany; ++i) { //choose hosts to infect 
									//with replacement because infection is independent of host status

									int pick = floor(ranf_mt(thread) * (numHosts - 1));
									//Wes' code: will give uniform integer distribution between 0 and x-1 inclusive (here, x=Hosts.size()-1)
									if (pick == h) { //if number chosen is current host
										pick = numHosts - 1; //swap to last host
									}

									HostNewInfect.push_back(pick); //put host in vector of new infections
									Hosts[pick].CandBact.push_back(new_clone(*bit)); //put into list of candidate infectors
									++vclonestrain[bit->OfStrain->StrainIndex][t];
									++clones[t];
								}
							} //finish with bit
							//add to count of number of hosts infected
							addbool(tempinf, infecteds);

							//and if hosts infected with res, sens, or both
							if (tempressens[0] > 0) { //if at least one resistant bact in this host...
								if (tempressens[1] > 0) { //...and at least one sensitive bact as well
									++CountResSens[2][t]; //increase counter for resistant/sensitive coinfection
								}
								else { //resistant bact only
									++CountResSens[0][t]; //increase counter for resistant only
								}
							}
							else { //no resistant bact, can only be sensitive
								if (tempressens[1] > 0) { //if have sensitive bact
									++CountResSens[1][t]; //increase counter for sensitive only
								}
							}

							int tempint = tempressens[0] + tempressens[1];
							int resint = tempressens[0];
			
							//now update counts of multiple infections
							if (tempint < 5) {
								++CountGetInf[tempint][t]; //increase counter for number of bacteria in host
							}
							else {
								++CountGetInf[5][t]; //increase counter for >8 bacteria in host
							}

							//now update counts of multiple infections
							if (resint < 5) {
								++CountGetInfRes[resint][t]; //increase counter for number of bacteria in host
							}
							else {
								++CountGetInfRes[5][t]; //increase counter for >5 bacteria in host
							}
						} //host has bacteria after gamma/antibiotic clearance
					} //host has bacteria before clearance
				} //this finishes with alive host
			} //this finishes loop through hosts

			for (h = HostNewInfect.size() - 1; h > -1; h--) {
				//check if more than 1 to be infected
				int tempcount = Hosts[HostNewInfect[h]].CandBact.size();
				int tempindexresident = Hosts[HostNewInfect[h]].WithBact.size() - 1; //index for last element of resident bacteria
				int tempindexall = tempindexresident + tempcount;
				
				//DO THIS LIKE MGES: random choice and pop_back - checked - does not make any faster!!!
				//then could do competition against residents one-by-one without all the other crap!
				//transfer all candidate bacteria (have to do first before shuffling, due to abstract base class)
				Hosts[HostNewInfect[h]].WithBact.transfer(Hosts[HostNewInfect[h]].WithBact.end(), //put at end of WithBact....
					Hosts[HostNewInfect[h]].CandBact.begin(), Hosts[HostNewInfect[h]].CandBact.end(), //...whole CandBact vector
					Hosts[HostNewInfect[h]].CandBact);

				//set Carriage to current t for all new bacteria
				for (b = tempindexresident + 1; b < tempindexall + 1; ++b) {
					Hosts[HostNewInfect[h]].WithBact[b].BactCarriage = t;
					//make sure clock started on all new bacteria
				}


				//if statement for tempindexall>0
				if (tempindexall > 0) { //only do competition if there is more than one candidate or at least one resident

				//if more than 1 addition, shuffle additions
					if (tempcount > 1) {
						random_shuffle(Hosts[HostNewInfect[h]].WithBact.end() - tempcount, Hosts[HostNewInfect[h]].WithBact.end());
					}

					if (CompetitionBlock == 'F') { //competition fitness determined as fraction of sum of two competition fitnesses

						for (b = tempindexall; b > tempindexresident; b--) { //only new additions compete
							//count attempted infection of host who already has resident bacteria
							if (tempindexresident > -1) { //if host has at least one bacteria
								++tryinfectcarrier[t]; //increase count of attempted infection for every candidate bacteria
							}

							//get next bacteria in line
							for (b2 = tempindexresident; b2 > -1; b2--) { //new candidate invaders must compete against residents first - opportunity to get established earlier in timestep than other candidate
								++trycomp[t];

								int index = Hosts[HostNewInfect[h]].WithBact[b2].OfStrain->StrainIndex;
								//have to do this here rather than inside function because of hierarchy

								//see if competition happens
								switch (Hosts[HostNewInfect[h]].WithBact[b].competition_acq(Hosts[HostNewInfect[h]].WithBact[b2], index, thread, InvadeFactor)) {

								case 0:
									//no competition occurs
									++nocomp[t];
									break; //goes to next b2
								case 1:
									++residentlosescomp[t];
									break; //goes to next b2
								case 2:
								++vclearcompstrain[Hosts[HostNewInfect[h]].WithBact[b].OfStrain->StrainIndex][t];
									Hosts[HostNewInfect[h]].WithBact.erase(Hosts[HostNewInfect[h]].WithBact.begin() + b); //erase element b
									++clearcomp[t];
									goto newchallengerF;
								default:
									cout << '\n' << "error in competition function";
								} //end of switch statement (competition function)
							} // end of b2
							//if b is still present now, needs to become resident

							//if more than one invader, swap positions
							if (tempcount > 1) {
								//swap() function makes copy of first object, which is then cleared when goes out of scope					
								swap(Hosts[HostNewInfect[h]].WithBact[b], Hosts[HostNewInfect[h]].WithBact[tempindexresident + 1]); //swap with first invader
							}

							if (tempindexresident > -1) { //if host had at least one bacteria
								++succeedinfectcarrier[t]; //increase count of attempted infection for every candidate bacteria
							}

							++tempindexresident; //increase index of last resident
							//b needs to stay the same, as has been replaced by another invader
							++b; //then will de-increment when gets to end of loop, subtract 1, and start again with same b
						newchallengerF:;
						} //end of b and competition loop

					} //end of CompetitionBlock
					else if (CompetitionBlock == 'D') { //competition determined simply by whether competition fitness of one is greater than the other

						for (b = tempindexall; b > tempindexresident; b--) { //only new additions compete
							//count attempted infection of host who already has resident bacteria
							if (tempindexresident > -1) { //if host has at least one bacteria
								++tryinfectcarrier[t]; //increase count of attempted infection for every candidate bacteria
							}

							//get next bacteria in line
							for (b2 = tempindexresident; b2 > -1; b2--) { //new candidate invaders must compete against residents first - opportunity to get established earlier in timestep than other candidate
								++trycomp[t];

								int index = Hosts[HostNewInfect[h]].WithBact[b2].OfStrain->StrainIndex;
								//have to do this here rather than inside function because of hierarchy

								//see if competition happens
								switch (Hosts[HostNewInfect[h]].WithBact[b].competition_lessfit(Hosts[HostNewInfect[h]].WithBact[b2], index, thread, InvadeFactor)) {
								case 0:
									++nocomp[t];
									break; //goes to next b2
								case 1:
									//b wins
									++residentlosescomp[t];
									break; //goes to next b2
								case 2:
									//b2 wins
									//don't update strain clears or carriage, as these are new infections
									++vclearcompstrain[Hosts[HostNewInfect[h]].WithBact[b].OfStrain->StrainIndex][t];
									Hosts[HostNewInfect[h]].WithBact.erase(Hosts[HostNewInfect[h]].WithBact.begin() + b); //erase element b
									++clearcomp[t];

									goto newchallengerD;
								default:
									cout << '\n' << "error in competition function";
								} //end of switch statement (competition function)
							} // end of b2
							  //if b is still present now, needs to become resident
							//swap() function makes copy of first object, which is then cleared when goes out of scope
							swap(Hosts[HostNewInfect[h]].WithBact[b], Hosts[HostNewInfect[h]].WithBact[tempindexresident + 1]); //swap with first invader

							if (tempindexresident > -1) { //if host had at least one bacteria
								++succeedinfectcarrier[t]; //increase count of attempted infection for every candidate bacteria
							}

							++tempindexresident; //increase index of last resident
												 //b needs to stay the same, as has been replaced by another invader
							++b; //then will de-increment when gets to end of loop, subtract 1, and start again with same b
								 //cout << "b beat resident! New tempindexresident = " << tempindexresident << " and new b = " << b << "\n";
						newchallengerD:;
						} //end of b and competition loop
					} //end option for competition block
				} //end if statement for competition if more than one bacteria
			} //end of infections

			for (h = HostNewDead.size() - 1; h > -1; --h) {//iterate backwards over HostNewDead so do largest first
				//swap and pop_back much faster, but creates copy of first host
				//#need to do for HostConsumeNext as well

				for (b = 0; b < Hosts[HostNewDead[h]].WithBact.size(); ++b) {
					++cleardead[t];
					++vcleardeathstrain[Hosts[HostNewDead[h]].WithBact[b].OfStrain->StrainIndex][t];
					++vStrainClears[Hosts[HostNewDead[h]].WithBact[b].OfStrain->StrainIndex][t];
					vStrainCarriage[Hosts[HostNewDead[h]].WithBact[b].OfStrain->StrainIndex][t] += (t - Hosts[HostNewDead[h]].WithBact[b].BactCarriage);
				}

				swap(Hosts[HostNewDead[h]], Hosts.back());
				Hosts.pop_back();

				swap(HostConsumeNext[HostNewDead[h]], HostConsumeNext.back());
				HostConsumeNext.pop_back();

			}
			HostNewDead.clear(); //clear contents of HostNewDead 

			///MOVED HOSTS ALIVE
			CountHostsAlive[t] = Hosts.size(); //NOT numHosts as we have now made deaths

			//using Wes' function, Poisson approximation to binomial: lambda=N*p
			int births = ignpoi_mt(HostBirth * HostStartNum, thread);  //birth rate is constant, based on starting population size			

			//cout << births << " new births" << '\n';
			for (h = 0; h < births; ++h) {

				Hosts.push_back(Host());
				if (ranf_mt(thread) < HostABc) {
					HostConsumeNext.push_back(1); //change appropriate proportion to 1's
				}
				else {
					HostConsumeNext.push_back(0);
				}

			}

			//--------------------------


			//loop through strains and get counts

			for (i = 0; i < numstrains; ++i) {
				vCountHostStrain[i][t] = infecteds[i];
			}
		} //this is end of timestep


		//-------------------------
		//note that host_unique_1 only refers to hosts with more than one isolate, but one unique strain (i.e. dual/multi infection same strain) 
		//add together first!
		for (t = 0; t < Timestop; ++t) {
			CountGetInfUnique[0][t] += CountGetInf[1][t]; //now unique_1 includes hosts with single bacterium
		}

		//print outputfiles
		printoutputs:;
		ofstream OutputCounts;
		OutputCounts.open(outputname, ios::app); //appends to existing file
		

		//--------------------------------
		// Default outputs


		OutputCounts << '\n' << "hosts_alive," << Years << "," << HostStartNum << ","
			<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
			<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
		for (t = 0; t < (Timestop - 1); t += OutputSteps) { //don't need count for extra timestep
			OutputCounts << CountHostsAlive[t] << ",";
		}
		OutputCounts << CountHostsAlive[(Timestop - 1)];  //last item doesn't get comma

		//hosts infected
		OutputCounts << '\n' << "hosts_infected," << Years << "," << HostStartNum << ","
			<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
			<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
		for (t = 0; t < (Timestop - 1); t += OutputSteps) {
			OutputCounts << CountHostsInf[t] << ",";
		}
		OutputCounts << CountHostsInf[(Timestop - 1)]; //last item doesn't get comma

		OutputCounts << '\n' << "bact_total," << Years << "," << HostStartNum << ","
			<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
			<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
		for (t = 0; t < (Timestop - 1); t += OutputSteps) {
			OutputCounts << CountBact[t] << ",";
		}
		OutputCounts << CountBact[(Timestop - 1)];

		OutputCounts << '\n' << "clones," << Years << "," << HostStartNum << ","
			<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
			<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
		for (t = 0; t < (Timestop - 1); t += OutputSteps) {
			OutputCounts << clones[t] << ",";
		}
		OutputCounts << clones[(Timestop - 1)];

		for (i = 0; i < numstrains; ++i) {
			OutputCounts << '\n' << "hosts_inf_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << vCountHostStrain[i][t] << ",";
			}
			OutputCounts << vCountHostStrain[i][(Timestop - 1)]; //last item doesn't get comma

		}

		for (i = 0; i < numstrains; ++i) {
			OutputCounts << '\n' << "total_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << vCountStrain[i][t] << ",";
			}
			OutputCounts << vCountStrain[i][(Timestop - 1)]; //last item doesn't get comma

			OutputCounts << '\n' << "single_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << vCountSingleStrain[i][t] << ",";
			}
			OutputCounts << vCountSingleStrain[i][(Timestop - 1)]; //last item doesn't get comma
		}

		for (i = 0; i < numstrains; ++i) {
			for (j = 0; j < nummges; ++j) {
				if ((vStrainMGEPrev[i][j] > 0) || (MGETypes[j].MGETypeBeta > 0) || (BactIntro > 0)) { //only print MGEs for each strain if it's possible they could be in it (HGT>0 and/or starting prevalence >0, or migration)
					OutputCounts << '\n' << "strain_mge_" << Strains[i].StrainName << "_" << MGETypes[j].MGETypeName
						<< "," << Years << "," << HostStartNum << ","
						<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
						<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
					for (t = 0; t < (Timestop - 1); t += OutputSteps) {
						OutputCounts << vvStrainMGE[i][j][t] << ",";
					}
					OutputCounts << vvStrainMGE[i][j][(Timestop - 1)]; //last item doesn't get comma
				}
			}
		}


		for (j = 0; j < nummges; ++j) {
			for (i = 0; i < numstrains; ++i) {
				for (k = 0; k < numstrains; ++k) {
						if (MGETypes[j].MGETypeBeta > 0) { //only print if HGT>0
							OutputCounts << '\n' << "mge_from_to_" << MGETypes[j].MGETypeName << "_" << Strains[i].StrainName << "_" << Strains[k].StrainName
								<< "," << Years << "," << HostStartNum << ","
								<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
								<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
							for (t = 0; t < (Timestop - 1); t += OutputSteps) {
								OutputCounts << vvvSuccessHGT[j][i][k][t] << ",";
							}
							OutputCounts << vvvSuccessHGT[j][i][k][(Timestop - 1)]; //last item doesn't get comma
						}
				}
			}
		}

		//---------------
		// Optional outputs

		if (OutputBactAllSites == 1) {
			OutputCounts << '\n' << "bact_allsites," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << allsites[t] << ","; //printing number of bacteria with all MGE sites occupied
			}
			OutputCounts << allsites[(Timestop - 1)];
		}

		if (OutputHostBactDist == 1) {
			for (i = 1; i < 6; ++i) { //not printing zero counters
				OutputCounts << '\n' << "host_bact_" << i;
				if (i == 5) {
					OutputCounts << "+";
				}
				OutputCounts << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << CountGetInf[i][t] << ",";
				}
				OutputCounts << CountGetInf[i][(Timestop - 1)]; //last item doesn't get comma
			}

			for (i = 0; i < 5; ++i) { //not printing zero counters

				OutputCounts << '\n' << "host_unique_" << i + 1;
				OutputCounts << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";

				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << CountGetInfUnique[i][t] << ",";
				}
				OutputCounts << CountGetInfUnique[i][(Timestop - 1)]; //last item doesn't get comma
			}


		}

		if (OutputResistance == 1) {
			OutputCounts << '\n' << "bact_resistant," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << CountBactRes[t] << ",";
			}
			OutputCounts << CountBactRes[(Timestop - 1)];

			OutputCounts << '\n' << "host_resonly," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << CountResSens[0][t] << ","; //printing number of hosts infected with resistant only
			}
			OutputCounts << CountResSens[0][(Timestop - 1)];

			OutputCounts << '\n' << "host_sensonly," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << CountResSens[1][t] << ","; //printing number of resistant/sensitive coinfections
			}
			OutputCounts << CountResSens[1][(Timestop - 1)];

			OutputCounts << '\n' << "host_ressens," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << CountResSens[2][t] << ","; //printing number of resistant/sensitive coinfections
			}
			OutputCounts << CountResSens[2][(Timestop - 1)];

			for (i = 1; i < 6; ++i) { //not printing zero counters
				OutputCounts << '\n' << "host_resbact_" << i;
				if (i == 5) {
					OutputCounts << "+";
				}
				OutputCounts << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << CountGetInfRes[i][t] << ",";
				}
				OutputCounts << CountGetInfRes[i][(Timestop - 1)]; //no comma
			}
		}

		if (OutputBactClearances == 1) {
			OutputCounts << '\n' << "clear_dead," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << cleardead[t] << ",";
			}
			OutputCounts << cleardead[(Timestop - 1)]; //last item doesn't get comma

			OutputCounts << '\n' << "clear_abc," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << clearabc[t] << ",";
			}
			OutputCounts << clearabc[(Timestop - 1)]; //last item doesn't get comma

			OutputCounts << '\n' << "clear_gamma," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << cleargamma[t] << ",";
			}
			OutputCounts << cleargamma[(Timestop - 1)]; //last item doesn't get comma

			OutputCounts << '\n' << "clear_comp," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << clearcomp[t] << ",";
			}
			OutputCounts << clearcomp[(Timestop - 1)]; //last item doesn't get comma

			if (OutputResistance == 1) {
				OutputCounts << '\n' << "clear_comp_within_SkillS," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearWithinStatus[0][0][t] << ",";
				}
				OutputCounts << vvclearWithinStatus[0][0][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_within_SkillR," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearWithinStatus[0][1][t] << ",";
				}
				OutputCounts << vvclearWithinStatus[0][1][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_within_RkillS," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearWithinStatus[1][0][t] << ",";
				}
				OutputCounts << vvclearWithinStatus[1][0][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_within_RkillR," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearWithinStatus[1][1][t] << ",";
				}
				OutputCounts << vvclearWithinStatus[1][1][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_between_SS," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearBetweenStatus[0][0][t] << ",";
				}
				OutputCounts << vvclearBetweenStatus[0][0][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_between_SR," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearBetweenStatus[0][1][t] << ",";
				}
				OutputCounts << vvclearBetweenStatus[0][1][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_between_RS," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearBetweenStatus[1][0][t] << ",";
				}
				OutputCounts << vvclearBetweenStatus[1][0][(Timestop - 1)]; //last item doesn't get comma

				OutputCounts << '\n' << "clear_comp_between_RR," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vvclearBetweenStatus[1][1][t] << ",";
				}
				OutputCounts << vvclearBetweenStatus[1][1][(Timestop - 1)]; //last item doesn't get comma
			}

			for (i = 0; i < numstrains; ++i) {
				OutputCounts << '\n' << "clear_gamma_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vcleargammastrain[i][t] << ",";
				}
				OutputCounts << vcleargammastrain[i][(Timestop - 1)]; //last item doesn't get comma

			}

			for (i = 0; i < numstrains; ++i) {
				OutputCounts << '\n' << "clear_comp_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vclearcompstrain[i][t] << ",";
				}
				OutputCounts << vclearcompstrain[i][(Timestop - 1)]; //last item doesn't get comma

			}

			for (i = 0; i < numstrains; ++i) {
				OutputCounts << '\n' << "clear_dead_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
				for (t = 0; t < (Timestop - 1); t += OutputSteps) {
					OutputCounts << vcleardeathstrain[i][t] << ",";
				}
				OutputCounts << vcleardeathstrain[i][(Timestop - 1)]; //last item doesn't get comma

			}

		}

		if ((CompetitionBlock == 'F') || (CompetitionBlock == 'D')) {
			OutputCounts << '\n' << "attempt_inf_carrier," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << tryinfectcarrier[t] << ",";
			}
			OutputCounts << tryinfectcarrier[(Timestop - 1)]; //last item doesn't get comma

			OutputCounts << '\n' << "succeed_inf_carrier," << Years << "," << HostStartNum << ","
				<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
				<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";
			for (t = 0; t < (Timestop - 1); t += OutputSteps) {
				OutputCounts << succeedinfectcarrier[t] << ",";
			}
			OutputCounts << succeedinfectcarrier[(Timestop - 1)]; //last item doesn't get comma
		}

		//ELFI outputs - legacy
		string chiname = "Chitest.csv";
		double rawchi = 999.9;
		double samplechi = 999.9;
		double relvchi = 999.9;
		double censchi = 999.9;

		//return 'no strain' summary statistic if no bacteria
		//will go straight to 'nostrains' when i=0
		//THIS IS NOT CALIBRATED FOR MGE SS
		if (CountBact[(Timestop - 1)] == 0) {
			SumStat = SumStatNoStrains;
			goto nostrains;
		}

		//now calculate average duration of carriage for each strain
		for (i = 0; i < numstrains; ++i) {
			if (OutputBactCarriage == 1) {

				//carriage_mean is rolling mean
				OutputCounts << '\n' << "carriage_mean_strain_" << Strains[i].StrainName << "," << Years << "," << HostStartNum << ","
					<< ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << ","
					<< ParmsMGEGammaCost << "," << ParmsMGEBetaCost << "," << HostABcYr << "," << run << "," << ParmsR0 << ",";

				//burn-in period

				for (t = 0; t < Burnsteps; t += OutputSteps) {
					OutputCounts << ",";
				}

				//after burn-in period: reset counters
				for (t = Burnsteps; t < (Timestop - 1); t += OutputSteps) {
					if (vCountStrain[i][t] == 0) { //if no bacteria left, carriage is 0
						OutputCounts << 0 << ",";
					}
					else {
						int tempclears = accumulate(vStrainClears[i].begin() + Burnsteps, vStrainClears[i].begin() + t + 1, 0);
						int tempcarriage = accumulate(vStrainCarriage[i].begin() + Burnsteps, vStrainCarriage[i].begin() + t + 1, 0);
						if (tempclears == 0) {
							OutputCounts << tempcarriage / (double)vCountStrain[i][t] << ",";
							//if there have been no clearances since burnin, average duration is total length of carriage/no. of bacteria in existence
						}
						else {
							OutputCounts << tempcarriage / ((double)tempclears * StepsYr) << ",";
						}
					}
				}

				double tempstrainstat1;


				if (vCountStrain[i][(Timestop - 1)] < 100) { //if not bacteria of this strain left
					--strainsurvive;
					OutputCounts << 0; //last item doesn't get comma
					tempstrainstat1 = 1; //difference between observed and expected is expected --> relative difference is 1
					//tempstrainstat2 = 1;
				}
				else {
					//get average duration
					double mean = accumulate(vStrainCarriage[i].begin() + Burnsteps, vStrainCarriage[i].end(), 0) /
						((double)accumulate(vStrainClears[i].begin() + Burnsteps, vStrainClears[i].end(), 0) * StepsYr);
					OutputCounts << mean; //last item doesn't get comma
					tempstrainstat1 = abs(mean - vStrainTargetDuration[i]) / vStrainTargetDuration[i]; //relative difference
				}
				//update strainstat
				strainstat1 += tempstrainstat1;

			}
		}
		OutputCounts.close();

		//relative difference
		double inf_denom = (double)CountHostsInf[Timestop - 1] / CountHostsAlive[Timestop - 1];
		tempstat1 = abs(TargetPrev - inf_denom) / TargetPrev;
		double uniq_denom = (double)CountGetInfUnique[0][Timestop - 1] / CountHostsInf[Timestop - 1];
		tempstat2 = abs(TargetUnique - uniq_denom) / TargetUnique;
		double res_denom = (double)CountBactRes[Timestop - 1] / CountBact[Timestop - 1];
		resstat = min(abs(TargetRes - res_denom) / TargetRes, 1);

		//sumstats for mges
		double chisq = 0.0;
		double totalprev = (double)CountMGE[(Timestop - 1)] / CountBact[(Timestop - 1)];
		double mgeprevstat = 0.0;
		double chistat = -999.9;


		//return big numbers if no MGEs at all
		//#probably need error catcher for each MGE
		if (totalprev == 0) {
			mgeprevstat = 0.0;
			chistat = 1.0; //gets logged later


		}
		else {

			for (k = 0; k < nummges; ++k) {
				double mgeprev = (double)vCountMGE[k][(Timestop - 1)] / CountBact[(Timestop - 1)];
				double mgestat = abs(vMGETargetPrev[k] - mgeprev) / vMGETargetPrev[k];
				//now censor
				mgestat = min(mgestat, 1.0);
				//add to overall stat
				mgeprevstat += log10(mgestat);

				//***this is not secure for multiple MGEs
				for (i = 0; i < numstrains; ++i) {
					double straintot = vCountStrain[i][(Timestop - 1)];
					double mgestrainstat = -999.9;
					if (straintot == 0) {
						mgestrainstat = 0;
					}
					else {
						mgestrainstat = pow(mgeprev * straintot - ((double)vvStrainMGE[i][k][(Timestop - 1)]), 2) / (mgeprev * straintot);

					}
					chisq += mgestrainstat;

				}
			}

			rawchi = chisq;

			//now take chi squared and divide by total number of samples
			chisq /= CountBact[(Timestop - 1)];

			samplechi = chisq;

			//then make this figure relative to expected chi squared
			chistat = abs(MGETargetChi - chisq) / MGETargetChi;
			relvchi = chistat;



			chistat = min(chistat, 1.0);

			censchi = chistat;


		}
		//censor
		chistat = log10(chistat);


		double SumStat = mgeprevstat + chistat;
	nostrains:;
		printf("%.6f", SumStat);


		ofstream ChiStats;
		ChiStats.open(chiname, ios::app);
		ChiStats << "\n" << (double)CountBactRes[(Timestop - 1)] / CountBact[(Timestop - 1)] << "," << rawchi << "," << samplechi << "," << relvchi << "," << censchi << "," << chistat << "," << mgeprevstat << "," << SumStat;
		ChiStats.close();

		//profiling
		auto clock_end = chrono::steady_clock::now();

		ofstream OutputStats;
		OutputStats.open(statname, ios::app); //appends to existing file
		OutputStats << '\n' << Page << "," << ParmsR0 << "," << ParmsIntraPcomp << "," << ParmsInterPcomp << "," << ParmsMGECompetitionCost << "," << ParmsMGEGammaCost << "," <<
			ParmsMGEBetaCost << "," << strainsurvive << "," << strainstat1 << "," << (double)CountHostsInf[(Timestop - 1)] / CountHostsAlive[(Timestop - 1)] << "," << tempstat1 << "," <<
			(double)CountGetInfUnique[0][(Timestop - 1)] / CountHostsInf[(Timestop - 1)] << "," << tempstat2 << "," << (double)CountBactRes[(Timestop - 1)] / CountBact[(Timestop - 1)] << "," <<
			resstat << "," << mgeprevstat << "," << chistat << "," << SumStat << "," << chrono::duration_cast<chrono::seconds>(clock_end - clock_start).count() << "";
		OutputStats.close();

	} //end runs
//} //end pragma omp parallel


	return 0;
}
