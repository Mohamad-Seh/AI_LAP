// Genetic5.cpp : Defines the entry point for the console application.
#pragma warning(disable:4786)		// disable debug warning

#include"Genetic5.h"
using namespace std;				// polluting global namespace, but hey...

typedef vector<ga_struct> ga_vector;// for brevity

//****************************************************************************************************

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

//****************************************************************************************************

void BullsEye(ga_vector &population)
{
	//same struction as normal method 
	string target = GA_TARGET;
	int tsize = target.size();
	unsigned int fitness;
	float average = 0, devi = 0;
	for (int i = 0; i < GA_POPSIZE; i++)
	{
		fitness = tsize * 10;
		for(int k = 0 ;k < tsize; k++)
		{ 
			//if correct assemption
			if (population[i].str[k] == target[k])
				// bonus for correct assemption 
				fitness -= 10;
			else
			{
				for (int k = 0; k < tsize; k++)
				{
					// if correct assemption but wrong place
					if(population[i].str[k] == target[k])
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
	//calculating average and deviation
	average = average / GA_POPSIZE;
	for (int i = 0; i < GA_POPSIZE; i++)
		devi += pow(population[i].fitness - average, 2);
	devi= sqrt(devi / GA_POPSIZE);
	// update the average and deviation
	for (int i = 0; i < GA_POPSIZE; i++) 
	{
		population[i].average = (int)average;
		population[i].devi =(int) devi;
	}

}

//****************************************************************************************************

void calc_fitness(ga_vector &population, int method)
{
	if (method == 0) //normal method
	{
		string target = GA_TARGET;
		int tsize = target.size();
		unsigned int fitness;
		//parametrs
		float average = 0;
		float devi = 0;

		for (int i = 0; i < GA_POPSIZE; i++) {
			fitness = 0;
			for (int j = 0; j < tsize; j++) {
				fitness += abs(int(population[i].str[j] - target[j]));
			}
			//sum fitness
			population[i].fitness = fitness;
			average += fitness; 
		}
		//calculate average 
		average = average / GA_POPSIZE;
		for (int i = 0; i < GA_POPSIZE; i++)
			//calculate deviation
			devi += pow(population[i].fitness - average, 2);
		devi = sqrt(devi / GA_POPSIZE);
		//updating average and deviation
		for (int i = 0; i < GA_POPSIZE; i++) 
		{
			population[i].average = (int)average;
			population[i].devi = (int)devi;
		}
	}
	else // bullseye method
		BullsEye(population);
	
}

//****************************************************************************************************

bool fitness_sort(ga_struct x, ga_struct y)
{
	return (x.fitness < y.fitness);
}

//***************************************************************************************************

inline void sort_by_fitness(ga_vector &population)
{
	sort(population.begin(), population.end(), fitness_sort);
}

//***************************************************************************************************

void elitism(ga_vector &population,
	ga_vector &buffer, int esize)
{
	for (int i = 0; i < esize; i++) {
		buffer[i].str = population[i].str;
		buffer[i].fitness = population[i].fitness;
	}
}

//**************************************************************************************************

void mutate(ga_struct &member)
{
	int tsize = GA_TARGET.size();
	int ipos = rand() % tsize;
	int delta = (rand() % 90) + 32;

	member.str[ipos] = ((member.str[ipos] + delta) % 122);
}

//***************************************************************************************************

void mate(ga_vector &population, ga_vector &buffer, int PointOp)
{
	int esize = (int) GA_POPSIZE * (int)GA_ELITRATE;
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

//****************************************************************************************************

inline void print_best(ga_vector &gav)
{
	cout << "Best: " << gav[0].str << " (" << gav[0].fitness << ")" << endl;
}

//****************************************************************************************************

inline void swap(ga_vector *&population,
	ga_vector *&buffer)
{
	ga_vector *temp = population; population = buffer; buffer = temp;
}
//****************************************************************************************************
// PSO algorithim
typedef vector<ga_struct> ga_vector;
vector<PSO> particle_vector;	        // for particles
string globalBest;						   // for global best

void PSOrun()
{
	//initlizing global
	int tsize = GA_TARGET.size();
	globalBest.erase();

	for (int j = 0; j < tsize; j++)
		globalBest += (rand() % 90) + 32;

	// for initlizing the global , local and velocity randomly

	for (int i = 0; i < GA_POPSIZE; i++) {
		PSO particle;
		particle_vector.push_back(particle);

		if (particle.get_fitness() < particle.calc_fitness_particle(globalBest))
		{
			globalBest = particle.get_str();
		}
	}						

	PSO global_particle;			// for calculating the fitness of the global particle
	int k = 0;
	while (k < GA_MAXITER && global_particle.calc_fitness_particle(globalBest) != 0)
	{
		clock_t start = std::clock();
		for (int i = 0; i < GA_POPSIZE; i++) {
			string myVelocity;
			string myStr;
			myVelocity.erase();
			myStr.erase();
			//Particle Position update
			for (int j = 0; j < tsize; j++) {						
			 //based on equation that we leanded in class
				double r1 = (double)rand() / (RAND_MAX);
				double r2 = (double)rand() / (RAND_MAX);
				double inertia = W * particle_vector[i].get_velocity()[j];
				double cognivtive = C1 * r1 * (particle_vector[i].get_localBest()[j] - particle_vector[i].get_str()[j]);
			    double social = C2 * r2 * (globalBest[j] - particle_vector[i].get_str()[j]);
				double ics = inertia + cognivtive + social;
				myVelocity += ics;
				myStr += particle_vector[i].get_str()[j] + myVelocity[j];

			}
			// updating the velocity and the fitness
			particle_vector[i].set_velocity(myVelocity);					
			particle_vector[i].set_str(myStr);
			particle_vector[i].set_fitness(particle_vector[i].calc_fitness_particle(myStr));
			// updating the local and global best
			if (particle_vector[i].get_fitness()							
				< particle_vector[i].calc_fitness_particle(particle_vector[i].get_localBest()))
			{
				particle_vector[i].set_localBest(particle_vector[i].get_str());


				if (particle_vector[i].calc_fitness_particle(particle_vector[i].get_localBest())
					< particle_vector[i].calc_fitness_particle(globalBest))
				{
					globalBest = particle_vector[i].get_localBest();
				}
			}
		}
		cout << "Best: " << globalBest << " (" << global_particle.calc_fitness_particle(globalBest) << ")" << endl;
		clock_t end = std::clock();
		float toatl_time = (float)(end - start) / CLOCKS_PER_SEC;
		cout << "Clock Ticks: " << toatl_time << "s" << std::endl;
		k++;
	}

}

//*******************************************************************************

int* RWS(ga_vector &population, int *points, int* newFitness)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - esize);
	int *parents = new int[numOfParents];									// the selected parents
	int* sumFitness = new int[GA_POPSIZE];								    // sum of the fitness till index i
	sumFitness[0] = newFitness[0];
	for (int i = 1; i < GA_POPSIZE; i++) {
		sumFitness[i] = sumFitness[i - 1] + newFitness[i];

	}
	for (int i = 0; i < numOfParents; i++) {		// Roulette Wheel Selection 
		int j = 0;
		while (1)
		{
			if (sumFitness[j] >= *(points + i))
				break;
			j++;

			if (j >= GA_POPSIZE)

			{
				j = GA_POPSIZE - 1;
				break;
			}
		}
		if (j >= GA_POPSIZE)
			j = GA_POPSIZE - 1;
		*(parents + i) = j;
	}
	return parents;
}

//*******************************************************************************
//SUS function 

int* SUS(ga_vector &population, long totalFitness, int* newFitness)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - esize);
	int *parents = new int[numOfParents];			    // the selected parents
	int distance = totalFitness / numOfParents;	        // distance between the pointers
	int start = rand() % (distance + 1);			    // random number between 0 and distance
	int *points = new int[numOfParents];		         // list of (sorted)random numbers from 0 to the total fitness

	for (int i = 0; i < numOfParents; i++) {				// points is a (sorted) list of	random numbers														   
		*(points + i) = start + (i * distance);		        // from 0 to total fitness with constant steps. 
	}

	return RWS(population, points, newFitness);

}

//*******************************************************************************

int* Tournament(ga_vector &population, int size)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - esize);
	int *parents = new int[numOfParents];			    // the selected parents
	int* players = new int[size];						// K players in the tournament

	for (int i = 0; i < numOfParents; i++) {
		for (int j = 0; j < size; j++) {
			players[j] = rand() % GA_POPSIZE;
		}
		//playing tournament
		int winner = players[0];
		for (int i = 0; i < size; i++) {
			if (population[players[i]].fitness < population[winner].fitness)
				winner = players[i];
		}
		parents[i] = winner;
	}

	return parents;
}

//*******************************************************************************

void Scaling(ga_vector &population)
{
	for (int i = 0; i < GA_POPSIZE; i++) {
		// based on lecture a*f+b 
		population[i].fitness = static_cast<unsigned int>(0.5 * population[i].fitness + 20);  // linear transformation
	}

}

//*******************************************************************************

// for this we added age in ga struct
#define MAX_AGE	10					// Maximum age of a citizen

int* Aging(ga_vector &population)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int tsize = GA_TARGET.size();
	int numOfParents = 2 * (GA_POPSIZE - esize);
	int *parents = new int[numOfParents];			// the selected parents

	for (int i = esize; i < GA_POPSIZE; i++) {
		if (population[i].age > MAX_AGE) {
			ga_struct citizen;
			citizen.fitness = 0;
			citizen.age = 0;			         	 //reset age
			citizen.str.erase();
			for (int j = 0; j < tsize; j++)
				citizen.str += (rand() % 90) + 32;
			population[i] = citizen;
		}
	}

	for (int i = 0; i < numOfParents; i++) {
		parents[i] = rand() % GA_POPSIZE;
		while (population[parents[i]].age <= 0)
		{
			parents[i] = rand() % GA_POPSIZE;
		}
	}
	return parents;
}

//*******************************************************************************

//*******************************************************************************
//main function 
 
#define PSOflag          0  // if psoflag == 1 pso algorithim will run else GA algorithim will run

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
	if (PSOflag == 1)
	{
		clock_t start = std::clock();
		PSOrun();
		clock_t end = std::clock();
		const sec duration = clock::now() - before;
		cout << "Time Elapsed: " << duration.count() << std::endl;
		return 0;
	}
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
		cout << "Average: " << (*population)[0].average << std::endl;
		cout << "Standard Deviation: " << (*population)[0].devi << std::endl;
		cout << "Clock Ticks: " << toatl_time << "s" << std::endl;
	}
	const sec duration = clock::now() - before;
	cout << "Time Elapsed: " << duration.count() << std::endl;
	return 0;
}
