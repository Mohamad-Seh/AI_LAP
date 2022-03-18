#pragma once
#pragma warning(disable:4786)		// disable debug warning
#define GA_TARGET		std::string("Hello world!")
#include <string>					// for string class
#include <iostream>					// for cout etc.
using namespace std;				// polluting global namespace, but hey...

//from class 
#define C1	2						
#define C2	2					    
#define W	0.5						


class PSO {

public:
	PSO();
	int calc_fitness_particle(string citizenStr);
	//set & set str
	void set_str(string str);
	string get_str();
	// fitness
	void set_fitness(int fitness);
	int get_fitness();
	//local best
	void set_localBest(string localBest);
	string get_localBest();
	//velocity
	void set_velocity(string velocity);
	string get_velocity();

private:
	string str;						  // the string
	unsigned int fitness;			  // its fitness
	string localBest;                 // its local best
	string velocity;			      // its velocity
};


