#include <iostream>					// for cout etc.
#include <vector>					// for vector class
#include <string>					// for string class
#include <algorithm>				// for sort algorithm
#include <time.h>					// for random seed
#include <math.h>					// for abs()
#include <chrono>					// for elapsed time
#include <ctime>					// for clock ticks

#define GA_POPSIZE		1600		// ga population size
#define GA_MAXITER		16384		// maximum iterations
#define GA_ELITRATE		0.10f		// elitism rate
#define GA_MUTATIONRATE	0.50f		// mutation rate
#define GA_MUTATION		RAND_MAX * GA_MUTATIONRATE
#define GA_TARGET		std::string("Hello world!")
#define heuristic	1				// The Givin Hueristic or Bull's Eye Hueristic
#define SELECTION	3				// for selecting parents method
#define	K	5						// Tournament size 
#define MAX_AGE	10					// Maximum age of a citizen
#define N	100						// for the size of the board
#define MUTATION	2				// for mutation method (exchange mutation or insertion mutation)
using namespace std;

struct nq_struct
{
	int* board = new int[N];	   // The chess size N
	unsigned int fitness;		   // The Fitness;
	unsigned int age;
};