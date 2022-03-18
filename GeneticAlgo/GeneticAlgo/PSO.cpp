
#include "PSO.h"

PSO::PSO()
{
	int tsize = GA_TARGET.size();
	str.erase();
	velocity.erase();


	for (int j = 0; j < tsize; j++)
	{
		str += (rand() % 90) + 32;
		velocity += (rand() % 90) + 32;
	}

	fitness = calc_fitness_particle(str);
	localBest = str;
}


void PSO::set_str(string str)
{
	this->str = str;
}

void PSO::set_fitness(int fitness)
{
	this->fitness = fitness;
}

void PSO::set_localBest(string localBest)
{
	this->localBest = localBest;
}

void PSO::set_velocity(string velocity)
{
	this->velocity = velocity;
}


string PSO::get_str()
{
	return str;
}

int PSO::get_fitness()
{
	return fitness;
}

string PSO::get_localBest()
{
	return localBest;
}

string PSO::get_velocity()
{
	return velocity;
}


int PSO::calc_fitness_particle(string citizenStr)
{
	string target = GA_TARGET;
	int tsize = target.size();
	unsigned int fitness;
	fitness = 0;
		fitness = tsize * 10;
		for (int j = 0; j < tsize; j++) {
			if (citizenStr[j] == target[j]) {
				fitness -= 10;
			}

			else {
				for (int k = 0; k < tsize; k++) {
					if (citizenStr[j] == target[k]) {
						fitness -= 1;
						break;
					}
				}
			}
		}
		return fitness;
	}

