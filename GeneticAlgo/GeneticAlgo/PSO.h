#pragma once
#pragma warning(disable:4786)		// disable debug warning
#define GA_TARGET		std::string("Hello world!")

#include <string>					// for string class
#include <iostream>					// for cout etc.
using namespace std;				// polluting global namespace, but hey...


#define C1	2						// for exploration
#define C2	2					    // for exploitation
#define W	0.5						// for Inertia


class PSO {

public:
	PSO();
	int calc_fitness_particle(string citizenStr);

	//setters
	void set_str(string str);
	void set_fitness(int fitness);
	void set_localBest(string localBest);
	void set_velocity(string velocity);


	//getters
	string get_str();
	int get_fitness();
	string get_localBest();
	string get_velocity();




private:
	string str;							// the string
	unsigned int fitness;			  // its fitness
	string localBest;                // its local best
	string velocity;			     // its velocity
};


