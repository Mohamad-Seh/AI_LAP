// NQueens.cpp : Defines the entry point for the console application.
//
#include"NQueens.h"
typedef vector<nq_struct> nq_vector;
vector<float> avgs; // averages vector
vector<float> dvts; // deviations vector


void init_game(nq_struct& game)	// used in the minConflict algorithm
{
	for (int i = 0; i < N; i++)
		game.board[i] = i;
	random_shuffle(game.board, game.board + N);
	game.age = 0;
}

void init_population(nq_vector& population, nq_vector& buffer)
{
	for (int i = 0; i < GA_POPSIZE; i++)
	{
		nq_struct citizen;
		citizen.fitness = 0;
		citizen.age = 1;
		for (int i = 0; i < N; i++) // checking that no more 2 queens in the same slot
			citizen.board[i] = i;
		random_shuffle(citizen.board, citizen.board + N);
		population.push_back(citizen);
	}
	buffer.resize(GA_POPSIZE);
}


void printBoard(int* board) {
	for (int i = 0; i < N; i++)
		cout << board[i] << " ";
}

void calc_fitness(nq_vector& population)
{
	unsigned int fitness;
	float average = 0;
	float deviation = 0;
	for (int i = 0; i < GA_POPSIZE; i++)
	{
		fitness = 0;
		for (int j = 0; j < N - 1; j++)
			for (int k = j + 1; k < N; k++)
				if ((k - j) == abs(population[i].board[j] - population[i].board[k]))
					fitness += 1;
		population[i].fitness = fitness;
		average += fitness;
	}
	average = average / GA_POPSIZE;	//the average
	avgs.push_back(average);
	for (int i = 0; i < GA_POPSIZE; i++)
		deviation += pow(population[i].fitness - average, 2);
	deviation = sqrt(deviation / GA_POPSIZE); //the std deviation
	dvts.push_back(deviation);
}

void game_calc_fitness(nq_struct& game)
{
	unsigned int fitness = 0;
	for (int j = 0; j < N - 1; j++)
	{
		for (int k = j + 1; k < N; k++)
		{
			if (game.board[j] == game.board[k])
				fitness++;
			if ((k - j) == abs(game.board[j] - game.board[k]))
				fitness += 1;
		}
	}
	game.fitness = fitness;
}


template <class  S>
bool fitness_sort(S x, S y)
{
	return (x.fitness < y.fitness);
}

template <class  T, class S>
inline void sort_by_fitness(T& population)
{
	sort(population.begin(), population.end(), fitness_sort<S>);
}

void elitism(nq_vector& population, nq_vector& buffer, int esize)
{
	for (int i = 0; i < esize; i++)
	{
		buffer[i].board = population[i].board;
		buffer[i].fitness = population[i].fitness;
		buffer[i].age += 1;
	}
}

int* swapIndexArea(int* arr, int i, int j)
{
	int temp = *(arr + i);
	*(arr + i) = *(arr + j);
	*(arr + j) = temp;
	return arr;
}

void exchangeMutation(nq_struct& member)
{
	int i = rand() % N;
	int j = rand() % N;
	member.board = swapIndexArea(member.board, i, j);
}

void insertionMutation(nq_struct& member)
{
	int index_i = rand() % N;
	int temp = index_i;
	int index_j = rand() % N;
	int* myArr = new int[N];
	index_i = min(index_i, index_j);
	index_j = max(temp, index_j);
	for (int i = 0; i < index_i; i++)
		*(myArr + i) = member.board[i];
	temp = *(myArr + index_i);
	for (int i = index_i; i < index_j; i++)
		*(myArr + i) = member.board[i + 1];
	*(myArr + index_j) = temp;
	for (int i = index_j + 1; i < N; i++)
		*(myArr + i) = member.board[i];
	member.board = myArr;
}

void mutate(nq_struct& member)
{
	if (MUTATION == 1) // swapping
		exchangeMutation(member);
	else if (MUTATION == 2) // insert 
		insertionMutation(member);
}

int* Naive()
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int parentsNum = 2 * (GA_POPSIZE - n_size);
	// creating another child, if  we have already one
		parentsNum = (GA_POPSIZE - n_size);
	int* parents = new int[parentsNum];
	for (int i = 0; i < parentsNum; i++)
		*(parents + i) = rand() % (GA_POPSIZE / 2);
	return parents;
}

int* RWS(nq_vector& population, int* points, int* newFitness)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int parentsNum = 2 * (GA_POPSIZE - n_size);
		// creating another child, if  we have already one
		parentsNum = (GA_POPSIZE - n_size);
	int* parents = new int[parentsNum];	// selected value
	int* fitnessSum = new int[GA_POPSIZE]; // the sum till index i
	fitnessSum[0] = newFitness[0];
	for (int i = 1; i < GA_POPSIZE; i++)
		fitnessSum[i] = fitnessSum[i - 1] + newFitness[i];
	for (int i = 0; i < parentsNum; i++)
	{
		int flag = 1;
		int j = 0;
		while (flag == 1)
		{
			if (fitnessSum[j] >= *(points + i))
				flag = 0;
			j++;
			if (j >= GA_POPSIZE)
			{
				j = GA_POPSIZE - 1;
				flag = 0;
			}
		}
		if (j >= GA_POPSIZE)
			j = GA_POPSIZE - 1;
		*(parents + i) = j;
	}
	return parents;
}

int* SUS(nq_vector& population, long totalFitness, int* newFitness)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int numOfParents = 2 * (GA_POPSIZE - n_size);
		// creating another child, if  we have already one
		numOfParents = (GA_POPSIZE - n_size);
	int* parents = new int[numOfParents]; // the selected parents
	int dist = totalFitness / numOfParents;	// distance between pointers
	int srt = rand() % (dist + 1); // random value between 0 and distance
	int* points = new int[numOfParents]; // list of numbers from 0 till total fit
	for (int i = 0; i < numOfParents; i++) // list of random numbers from 0 till tot. fitness with const. steps														   
		*(points + i) = srt + (i * dist);
	return RWS(population, points, newFitness);
}


void Scaling(nq_vector& population)
{
	for (int i = 0; i < GA_POPSIZE; i++)
		population[i].fitness = static_cast<unsigned int>(0.2 * population[i].fitness + 10); // linear transformation
}


int PlayTournament(nq_vector& population, int* players) {

	int winner = players[0];
	for (int i = 0; i < K; i++)
		if (population[players[i]].fitness < population[winner].fitness)
			winner = players[i];
	return winner;
}

int* Tournament(nq_vector& population)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int parentsNum = 2 * (GA_POPSIZE - n_size);
		// creating another child, if  we have already one
		parentsNum = (GA_POPSIZE - n_size);
	int* parents = new int[parentsNum];	// the currect parents
	int* players = new int[K];	// k players gonna play the tournament 
	for (int i = 0; i < parentsNum; i++)
	{
		for (int j = 0; j < K; j++)
			players[j] = rand() % GA_POPSIZE;
		parents[i] = PlayTournament(population, players);
	}
	return parents;
}


int* Aging(nq_vector& population)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int t_size = N;
	int numOfParents = 2 * (GA_POPSIZE - n_size);
		// creating another child, if  we have already one
		numOfParents = (GA_POPSIZE - n_size);
	int* parents = new int[numOfParents];// the currect parents
	for (int i = n_size; i < GA_POPSIZE; i++)
	{
		if (population[i].age > MAX_AGE)
		{
			nq_struct citizen;
			citizen.fitness = 0;
			citizen.age = 0;
			for (int j = 0; j < t_size; j++)
				citizen.board[j] = j;
			random_shuffle(citizen.board, citizen.board + N);
			population[i] = citizen;
		}
	}
	for (int i = 0; i < numOfParents; i++)
	{
		int flag = 1;
		while (flag == 1)
		{
			parents[i] = rand() % GA_POPSIZE;
			if (population[parents[i]].age > 0)
				flag = 0;
		}
	}
	return parents;
}



int* parentSelector(nq_vector& population)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int parentsNum = 2 * (GA_POPSIZE - n_size);
		// creating another child, if  we have already one
		parentsNum = (GA_POPSIZE - n_size);
	int* parents = new int[parentsNum];	// the selected parents
	int* points = new int[parentsNum]; // list of numbers from 0 to the total fitness
	int* newFit = new int[GA_POPSIZE]; // the original fitness not good enough
	for (int i = 0; i < GA_POPSIZE; i++)
		newFit[i] = ((-1) * (population[i].fitness) + population[GA_POPSIZE - 1].fitness);
	long totalFit = 0;
	for (int i = 0; i < GA_POPSIZE; i++)
		totalFit += newFit[i];
	if (SELECTION == 1)
		parents = Naive(); // naive method
	if (SELECTION == 2) // RWS + Scaling
	{
		Scaling(population); // Scaling

		for (int i = 0; i < parentsNum; i++) // points is a list of	random nums from 0 till fitness															   
			*(points + i) = rand() % (totalFit + 1);
		sort(points, points + parentsNum);
		parents = RWS(population, points, newFit); // RWS use here
	}
	if (SELECTION == 3)
		parents = SUS(population, totalFit, newFit); // SUS method
	if (SELECTION == 4)
		parents = Tournament(population); // start tournament
	if (SELECTION == 5)
		parents = Aging(population); // againg
	return parents;
}

void PMX(nq_vector& population, nq_struct& memb1, nq_struct& memb2, int ind, int i1, int i2)
{

	for (int i = 0; i < N; i++)
	{
		memb1.board[i] = population[i1].board[i];
		memb2.board[i] = population[i2].board[i];
	}
	int temp1 = memb1.board[ind];
	int temp2 = memb2.board[ind];
	for (int i = 0; i < N; i++)	
	{
		if (memb1.board[i] == temp2)
			memb1.board[i] = temp1;
		if (memb2.board[i] == temp1)
			memb2.board[i] = temp2;
	}
	memb1.board[ind] = temp2;
	memb2.board[ind] = temp1;
}

void mate(nq_vector& population, nq_vector& buffer)
{
	int n_size = static_cast<int>(GA_POPSIZE * GA_ELITRATE);
	int t_size = N, spos, i1, i2;
	int* parents = parentSelector(population);
	elitism(population, buffer, n_size);
		for (int i = n_size; i < GA_POPSIZE; i = i + 2)
		{
			i1 = parents[i - n_size];
			i2 = parents[i - n_size + 1];
			spos = rand() % t_size;
			PMX(population, buffer[i], buffer[i + 1], spos, i1, i2);
			if (rand() < GA_MUTATION)
				mutate(buffer[i]);
			if (rand() < GA_MUTATION)
				mutate(buffer[i + 1]);
			buffer[i].age = 0;
			buffer[i + 1].age = 0;
		}
}

inline void print_best(nq_vector& gav)
{
	cout << "Best: ";
	printBoard(gav[0].board);
	cout << " (" << gav[0].fitness << ")" << endl;
}

inline void swap(nq_vector*& population, nq_vector*& buffer)
{
	nq_vector* temp = population; population = buffer; buffer = temp;
}

void minimalConflict(nq_struct& game)
{
	int* randQueens = new int[N];
	int queen = (rand() % N);
	int prevFitness = game.fitness;
	int prevQueen = game.board[queen];
	for (int i = 0; i < N; i++)
		randQueens[i] = i;
	random_shuffle(randQueens, randQueens + N);	//moving queens
	game.board[queen] = randQueens[0];
	game_calc_fitness(game);
	int currentFit = game.fitness;
	int chQueen = randQueens[0];
	for (int i = 1; i < N; i++)
	{
		game.board[queen] = randQueens[i];
		game_calc_fitness(game);
		if (game.fitness < currentFit) // good move
		{
			currentFit = game.fitness;
			chQueen = randQueens[i];
		}
	}
	if (currentFit >= prevFitness) // bad move
	{
		game.board[queen] = prevQueen;
		game_calc_fitness(game);
	}
	game.board[queen] = chQueen;
	game_calc_fitness(game); // updating fitness
}


#define conflict	1				// The given algorithm or minimal conflict

int main()
{
	using clock = std::chrono::system_clock;
	using sec = std::chrono::duration<double>;
	const auto before = clock::now();				// for elapsed time
	int gensNum = 0;
	nq_vector pop_alpha, pop_beta;
	nq_vector* population, *buffer;
	srand(unsigned(time(NULL)));
		if (conflict == 0)
		{

			init_population(pop_alpha, pop_beta);
			population = &pop_alpha;
			buffer = &pop_beta;
			for (int i = 0; i < GA_MAXITER; i++) {

				clock_t begin = std::clock();		// for clock ticks
				calc_fitness(*population);		// calculate fitness using the given heuristic
				sort_by_fitness<vector<nq_struct>, nq_struct>(*population);	// sort them
				print_best(*population);		// print the best one
				if ((*population)[0].fitness == 0)
					break;
				mate(*population, *buffer);		// mate the population together		
				swap(population, buffer);		// swap buffers
				clock_t end = std::clock();
				float time_spent = (float)(end - begin) / CLOCKS_PER_SEC;
				gensNum += 1;
				std::cout << "Average: " << avgs[i] << std::endl;
				std::cout << "Standard Deviation: " << dvts[i] << std::endl;
				std::cout << "Clock Ticks: " << time_spent << "s" << std::endl;
			}
			std::cout << "Number of Generations: " << gensNum << std::endl;
		}
		if (conflict == 1)
		{
			nq_struct game;
			init_game(game);
			for (int i = 0; i < GA_MAXITER; i++)
			{
				if (game.fitness == 0)
					break;
				minimalConflict(game);
				gensNum += 1;
			}
			printBoard(game.board);
			std::cout << "Number of Generations: " << gensNum << std::endl;
		}
	const sec duration = clock::now() - before;
	std::cout << "Time Elapsed: " << duration.count() << "s" << std::endl;
	return 0;
}
