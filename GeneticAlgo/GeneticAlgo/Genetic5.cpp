// Genetic5.cpp : Defines the entry point for the console application.
//


#pragma warning(disable:4786)		// disable debug warning

#include <iostream>					// for cout etc.
#include <vector>					// for vector class
#include <string>					// for string class
#include <algorithm>				// for sort algorithm
#include <time.h>					// for random seed
#include <math.h>					// for abs()
#include <chrono>					// for elapsed time
#include <ctime>					// for clock ticks


#define GA_POPSIZE		2048		// ga population size
#define GA_MAXITER		16384		// maximum iterations
#define GA_ELITRATE		0.10f		// elitism rate
#define GA_MUTATIONRATE	0.25f		// mutation rate
#define GA_MUTATION		RAND_MAX * GA_MUTATIONRATE
#define GA_TARGET		std::string("Hello world!")
#define heuristic	2				// The Givin Hueristic or Bull's Eye Hueristic


using namespace std;				// polluting global namespace, but hey...

struct ga_struct
{
	string str;						// the string
	unsigned int fitness;			// its fitness
	int average = 0;
	int deviation = 0;
};

typedef vector<ga_struct> ga_vector;// for brevity

void init_population(ga_vector &population,
	ga_vector &buffer)
{
	int tsize = GA_TARGET.size();

	for (int i = 0; i < GA_POPSIZE; i++) {
		ga_struct citizen;

		citizen.fitness = 0;
		citizen.str.erase();

		for (int j = 0; j < tsize; j++)
			citizen.str += (rand() % 90) + 32;

		population.push_back(citizen);
	}

	buffer.resize(GA_POPSIZE);
}

void BullsEye_calc_fitness(ga_vector &population)
{
	string target = GA_TARGET;
	int tsize = target.size();
	unsigned int fitness;
	float average = 0;
	float deviation = 0;
	for (int i = 0; i < GA_POPSIZE; i++)
	{
		fitness = tsize * 10;
		for(int j = 0 ;j < tsize; j++)
		{
			if (population[i].str[j] == target[j])
				fitness -= 10;
			else
			{
				for (int k = 0; k < tsize; k++)
				{
					if(population[i].str[j] == target[k])
					{
						fitness--;
						break;
					}
				}
			}
		}
		population[i].fitness = fitness;
		average += fitness;
	}
	average = average / GA_POPSIZE;
	for (int i = 0; i < GA_POPSIZE; i++)
		deviation += pow(population[i].fitness - average, 2);
	deviation = sqrt(deviation / GA_POPSIZE);
	for (int i = 0; i < GA_POPSIZE; i++) // update the population avg. and devi.
	{
		population[i].average = average;
		population[i].deviation = deviation;
	}

}

void calc_fitness(ga_vector &population, int method)
{
	if (method == 0)
	{
		string target = GA_TARGET;
		int tsize = target.size();
		unsigned int fitness;

		float average = 0;
		float deviation = 0;

		for (int i = 0; i < GA_POPSIZE; i++) {
			fitness = 0;
			for (int j = 0; j < tsize; j++) {
				fitness += abs(int(population[i].str[j] - target[j]));
			}

			population[i].fitness = fitness;
			average += fitness; // sum of the fitness for all the GA_POPSIZE
		}
		average = average / GA_POPSIZE;
		for (int i = 0; i < GA_POPSIZE; i++)
			deviation += pow(population[i].fitness - average, 2);
		deviation = sqrt(deviation / GA_POPSIZE);
		for (int i = 0; i < GA_POPSIZE; i++) // update the population avg. and devi.
		{
			population[i].average = average;
			population[i].deviation = deviation;
		}
	}
	else
		BullsEye_calc_fitness(population);
	
}


bool fitness_sort(ga_struct x, ga_struct y)
{
	return (x.fitness < y.fitness);
}

inline void sort_by_fitness(ga_vector &population)
{
	sort(population.begin(), population.end(), fitness_sort);
}

void elitism(ga_vector &population,
	ga_vector &buffer, int esize)
{
	for (int i = 0; i < esize; i++) {
		buffer[i].str = population[i].str;
		buffer[i].fitness = population[i].fitness;
	}
}

void mutate(ga_struct &member)
{
	int tsize = GA_TARGET.size();
	int ipos = rand() % tsize;
	int delta = (rand() % 90) + 32;

	member.str[ipos] = ((member.str[ipos] + delta) % 122);
}

void mate(ga_vector &population, ga_vector &buffer, int PointOp)
{
	int esize = GA_POPSIZE * GA_ELITRATE;
	int tsize = GA_TARGET.size(), spos, i1, i2;
	int spos2;
	string str1, str2, str3;
	elitism(population, buffer, esize);

	// Mate the rest

	for (int i = esize; i < GA_POPSIZE; i++) {
		i1 = rand() % (GA_POPSIZE / 2);
		i2 = rand() % (GA_POPSIZE / 2);
		if(PointOp == 0) // one point
		{
			spos = rand() % tsize;

			buffer[i].str = population[i1].str.substr(0, spos) +
				population[i2].str.substr(spos, tsize - spos);
		}
		if (PointOp == 1) // two point
		{
			int curr_max, curr_min;
			spos = rand() % tsize;
			spos2 = rand() % tsize;
			curr_max = std::max(spos, spos2);
			curr_min = std::min(spos, spos2);
			str1 = population[i1].str.substr(0, curr_min);
			str2 = population[i2].str.substr(curr_min, curr_max - curr_min);
			str3 = population[i1].str.substr(curr_max, tsize - curr_max);
			buffer[i].str = str1 + str2 + str3;
		}
		if (PointOp == -1) // uniform 
		{
			string new_str;
			new_str.erase();
			for (int j = 0; j < tsize; j++)
			{
				spos = rand() % 2;
				if (spos == 0)
					new_str += population[i1].str.substr(j, 1);
				else
					new_str += population[i2].str.substr(j, 1);
			}
			buffer[i].str = new_str;
		}
		if (rand() < GA_MUTATION) mutate(buffer[i]);
	}
}

inline void print_best(ga_vector &gav)
{
	cout << "Best: " << gav[0].str << " (" << gav[0].fitness << ")" << endl;
}

inline void swap(ga_vector *&population,
	ga_vector *&buffer)
{
	ga_vector *temp = population; population = buffer; buffer = temp;
}

int main()
{
	srand(unsigned(time(NULL)));

	ga_vector pop_alpha, pop_beta;
	ga_vector *population, *buffer;

	init_population(pop_alpha, pop_beta);
	population = &pop_alpha;
	buffer = &pop_beta;

	using clock = std::chrono::system_clock;
	using sec = std::chrono::duration<double>;
	const auto before = clock::now();				// for elapsed time
	int numOfGenerations = 0;

	for (int i = 0; i < GA_MAXITER; i++) {
		clock_t start = std::clock();
		calc_fitness(*population, 3);		// calculate fitness
		sort_by_fitness(*population);	// sort them
		print_best(*population);		// print the best one

		if ((*population)[0].fitness == 0) break;

		mate(*population, *buffer, -1);		// mate the population together
		swap(population, buffer);		// swap buffers

		clock_t end = std::clock();
		float toatl_time = (float)(end - start) / CLOCKS_PER_SEC;
		numOfGenerations++;
		std::cout << "Average: " << (*population)[0].average << std::endl;
		std::cout << "Standard Deviation: " << (*population)[0].deviation << std::endl;
		std::cout << "Clock Ticks: " << toatl_time << "s" << std::endl;
	}
	const sec duration = clock::now() - before;
	std::cout << "Time Elapsed: " << duration.count() << std::endl;
	return 0;
}
