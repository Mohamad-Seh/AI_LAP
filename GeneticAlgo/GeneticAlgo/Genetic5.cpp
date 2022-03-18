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
	//initlizing paramaters and global
	int tsize = GA_TARGET.size();
	globalBest.erase();
	for (int j = 0; j < tsize; j++)
		globalBest += (rand() % 90) + 32;
	for (int i = 0; i < GA_POPSIZE; i++) {
		PSO particle;
		particle_vector.push_back(particle);

		if (particle.get_fitness() < particle.calc_fitness_particle(globalBest))
		{
			globalBest = particle.get_str();
		}
	}						

	PSO global_particle;		
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
	//parent 
	int *parents = new int[numOfParents];								
	int* sumFitness = new int[GA_POPSIZE];								   
	sumFitness[0] = newFitness[0];
	// summing fitness
	for (int i = 1; i < GA_POPSIZE; i++) {
		sumFitness[i] = sumFitness[i - 1] + newFitness[i];
	}
	//parent select
	for (int i = 0; i < numOfParents; i++) {		
		int j = 0;
		while (sumFitness[j] < *(points + i) && j < GA_POPSIZE) {
			j++;
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
	int *parents = new int[numOfParents];	
	// randnum is list of random integers
	int *randnum = new int[numOfParents];
	for (int i = 0; i < numOfParents; i++) {		    													   
		*(randnum + i) = rand() % (totalFitness / numOfParents + 1) + (i * totalFitness / numOfParents);
	}
	return RWS(population, randnum, newFitness);
}

//*******************************************************************************
#define	Psize	7						// Tournament size 
int* Tournament(ga_vector &population)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - esize);
	//parents
	int *parents = new int[numOfParents];		
	//players of the tournament
	int* players = new int[Psize];					
	//loop
	for (int i = 0; i < numOfParents; i++) {
		for (int j = 0; j < Psize; j++) {
			//random players
			players[j] = rand() % GA_POPSIZE;
		}
		//playing the tournament
		int win = players[0];
		for (int s = 0; s < Psize; i++) {
			if (population[players[i]].fitness < population[win].fitness)
				win = players[i];
		}
		//parent = the winner player
		parents[i] = win;
	}
	return parents;
}


//*******************************************************************************

void Scaling(ga_vector &population)
{
	//based on lecture scale = a*f+b a=0.2 , b= 10
	for (int i = 0; i < GA_POPSIZE; i++) {

		population[i].fitness = static_cast<unsigned int>(0.2 * population[i].fitness + 10);// linear transformation
	}
}

//*******************************************************************************

// for this we added age in ga struct
void updateAge(ga_vector &population,int i)
{
	population[i].age += 1;
}

//*******************************************************************************

#define parentype    4  //parent type 

int* selectParents(ga_vector &population)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - esize);
	int *parents = new int[numOfParents];			
	//sorted numbers
	int *points = new int[numOfParents];			
	//calculating new fitness
	int* newFitness = new int[GA_POPSIZE];								
	for (int i = 0; i < GA_POPSIZE; i++) {

		newFitness[i] = ((-1)*(population[i].fitness) + population[GA_POPSIZE - 1].fitness);
	}
	//total fitness
	long totalFitness = 0;
	for (int i = 0; i < GA_POPSIZE; i++) {
		totalFitness += newFitness[i];
	}
	//parent type:
	if (parentype ==1)                        
	{
		//RWS + Scaling 
		Scaling(population);		          	   
		for (int i = 0; i < numOfParents; i++) {				  															   
			*(points + i) = rand() % (totalFitness + 1);		  
		}
		//sorting
		sort(points, points + numOfParents);		
		//calling RWS
		parents = RWS(population, points, newFitness);			    
	}
	if (parentype == 2)
	{
		// SUS method
		parents = SUS(population, totalFitness, newFitness);     
	}
	if (parentype == 3)
	{
		// Tournament selection
		parents = Tournament(population);           
	}
	return parents;
}

//*******************************************************************************

void mate(ga_vector &population, ga_vector &buffer, int PointOp)
{
	int esize = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int tsize = GA_TARGET.size(), spos, spos2, i1, i2;
	int *parents = selectParents(population);
	elitism(population, buffer, esize);
	// Mate the rest
	for (int i = esize; i < GA_POPSIZE; i++) {
		i1 = parents[2 * (i - esize)];
		i2 = parents[2 * (i - esize) + 1];
		if (PointOp == 0) // Single Point Operator

		{
				spos = rand() % tsize;
				buffer[i].str = population[i1].str.substr(0, spos) +
				population[i2].str.substr(spos, tsize - spos);
		}
		if (PointOp == 1)  // Two Point Operator
		{
				spos = rand() % tsize;
				spos2 = rand() % tsize;
				buffer[i].str = population[i1].str.substr(0, std::min(spos, spos2)) +
				population[i2].str.substr(std::min(spos, spos2), std::max(spos, spos2) - std::min(spos, spos2)) +
				population[i1].str.substr(std::max(spos, spos2), tsize - std::max(spos, spos2));
		}
		if (PointOp == -1)  //// Uniform Point Operator
		{
				string myStr;
				myStr.erase();
				for (int j = 0; j < tsize; j++) {
					spos = rand() % 2;

					if (spos == 0)
					{
						myStr += population[i1].str.substr(j, 1);
					}
					if (spos == 1)
					{
						myStr += population[i2].str.substr(j, 1);
					}
				}
				buffer[i].str = myStr;
		}
		if (rand() < GA_MUTATION) mutate(buffer[i]);
		buffer[i].age = 0;
	}
}

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
		// // calculate fitness 
		calc_fitness(*population,1);		 // population and method 0 for nromal else for Bullseye method
		// sort them
		sort_by_fitness(*population);	
		print_best(*population);		// print the best one

		if ((*population)[0].fitness == 0) break;
		// mate the population together
		mate(*population, *buffer,1);		// population,buffer, PointOp : -1 is uniform , 0 for one point , 1 for two points
		swap(population, buffer);		// swap buffers
		//clock calculation
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
