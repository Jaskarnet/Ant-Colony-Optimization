#include <iostream>
#include "AdjacencyMatrix.h"
#include <vector>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <algorithm>
#include <iterator>
#include <windows.h>
#include <iomanip>
#include <iostream>
#include <iomanip>
#include <string>
#include <windows.h>
#include <iomanip>
#include <profileapi.h>
#include <Windows.h>
#include <chrono>
#include <queue>
#include <list>
#include <regex>
#include <thread>
#include <mutex>
#include <condition_variable>


/*
α = 1,
β od 2 do 5,
ρ = 0.5,
m = n (m - liczba mrówek, n - liczba miast)
τ0 = m/Cnn , gdzie Cnn jest szacowaną długością trasy
*/

#define ALPHA 1
#define BETA 2
#define RHO 0.5
#define MAX_THREADS 4
//#define C 10000  // parameter for CAS schema

using namespace std;

string configValues[3];



class Ant {
public:
	vector<unsigned int> visitedCities;
	unsigned int current_city;
	unsigned int next_city;
	int current_cost;

	Ant(unsigned int _current_city) {
		visitedCities.clear();
		visitedCities.push_back(_current_city);
		current_city = _current_city;
		next_city = 999999;
		current_cost = 0;
	}

	// choosing the next city across all the cities with a probability computed for each of ant using the roulette wheel selection method.
	void chooseNextCity(adjacencyMatrix<double>* pheromone_matrix, adjacencyMatrix<unsigned int>* adjMatrix, double alpha, double beta) {
		double probSum = 0;
		vector<double> prob;
		vector<double> probCumulative;
		for (int i = 0; i < pheromone_matrix->getVertexNumber(); i++) {
			//checks if a specific city, represented by the variable i, is not in the list of visited cities.
			if (!(find(visitedCities.begin(), visitedCities.end(), i) != visitedCities.end())) {
				double p = (pow(pheromone_matrix->getA()[current_city][i], alpha) * (double) pow((1/(double)adjMatrix->getA()[current_city][i]), beta));
				//*******\ test
					/*cout << setprecision(2) << pheromone_matrix->getA()[current_city][i] << "^" << alpha << " * " << (1 / (double)adjMatrix->getA()[current_city][i]) << "^" << beta << " = " << (pow(pheromone_matrix->getA()[current_city][i], alpha) * (double)pow((double)adjMatrix->getA()[current_city][i], beta)) << endl;
					cout << "p = " << p << endl;*/
				//*******/
				probSum += p;
				prob.push_back(p);
				probCumulative.push_back(probSum);
			}
			else {
				prob.push_back(0);
				probCumulative.push_back(probSum);
			}
		}
		//a random number is generated between 0 and probSum, and this number is used to determine the next city by checking which cumulative probability is greater than the random number.
		
		double randNum = ((double)rand() / (RAND_MAX)) * probSum;
		//*******\ test
			//cout << "randNum: " << randNum << endl;
			//cout << "probCumulative: " << probCumulative.back() << endl;
		//*******/

		for (int i = 0; i < probCumulative.size(); i++) {
			//*******\ test
				//cout << "probCumulative[" << i << "] = " << probCumulative[i] << " ";
			//*******/
			if (randNum <= probCumulative[i]) {
				
				next_city = i;
				//*******\ test
					//cout << endl << "wybrane miasto: " << next_city << endl << endl;
				//*******/
				break;
			} else if (i == probCumulative.size() - 1) next_city = probCumulative.size() - 1;  //jezeli jakims cudem nie znalazlo (bo double ma zipe) to we ostatni
		}
	}

	void visitCity(int _next_city, adjacencyMatrix<unsigned int>* adjMatrix) {
		current_cost += adjMatrix->getA()[current_city][_next_city];
		current_city = _next_city;
		next_city = _next_city;
		visitedCities.push_back(current_city);
	}
};

bool file_read_line(ifstream& file, string tab[], int size)
{
	string s;
	getline(file, s);

	if (file.fail() || s.empty())
		return(false);

	istringstream in_ss(s);

	for (int i = 0; i < size; i++)
	{
		in_ss >> tab[i];
		if (in_ss.fail())
			return(false);
	}
	return(true);
}

void readConfig(std::string file_name) {
	ifstream file;
	string tab[2];
	file.open(file_name.c_str());

	if (file.is_open())
	{
		for (int i = 0; i < 3; i++)
			if (file_read_line(file, tab, 2))
			{
				// wczytywanie configa: tab[0] - nazwa, tab[1] - wartosc
				if (tab[1] == "same_as_plik_do_otworzenia") configValues[i] = std::regex_replace(configValues[0], std::regex("data/"), "results/");
				else configValues[i] = tab[1];
			}
			else {
				std::cout << "File error - READ INFO" << endl;
				file.close();
			}
	}
	else
		std::cout << "File error - OPEN" << endl;
}

int readSolution(string filename) {
	std::string result = "";
	std::ifstream fin(filename);

	if (fin.is_open()) {
		fin.seekg(0, std::ios_base::end);      //Start at end of file
		char ch = ' ';                        //Init ch not equal to '\n'
		while (ch != '\n') {
			fin.seekg(-2, std::ios_base::cur); //Two steps back, this means we
											  //will NOT check the last character
			if ((int)fin.tellg() <= 0) {        //If passed the start of the file,
				fin.seekg(0);                 //this is the start of the line
				break;
			}
			fin.get(ch);                      //Check the next character
		}

		std::getline(fin, result);
		fin.close();

		//std::cout << "final line length: " << result.size() << std::endl;
		/*std::cout << "final line character codes: ";
		for (size_t i = 0; i < result.size(); i++) {
			std::cout << (int)result[i] << " ";
		}*/
		//std::cout << std::endl;
		//std::cout << endl << endl << "solution_from_file: " << result;
	}
	return stoi(result);
}

void updatePheromoneDecay(adjacencyMatrix<double>* pheromone_matrix, double rho) {
	for (int i = 0; i < pheromone_matrix->getVertexNumber(); i++) {
		for (int j = 0; j < pheromone_matrix->getVertexNumber(); j++) {
			if (i != j) {
				pheromone_matrix->getA()[i][j] = pheromone_matrix->getA()[i][j] * (double) rho;
			}
		}
	}
}

void updatePheromoneTrailCAS(vector<Ant*> ants, adjacencyMatrix<double>* pheromone_matrix, double pheromone_quantity_to_add) {
	for (Ant* ant : ants) {
		double delta = (double) pheromone_quantity_to_add / (double) ant->current_cost;
		for (unsigned int i = 0; i < ant->visitedCities.size() - 1; i++) { //nie dodaje feromonu do sciezki powrotnej MOZLIWE ZE TO ŹLE - juz dodaje
			unsigned int current_city = ant->visitedCities[i];
			unsigned int next_city = ant->visitedCities[i + 1];
			pheromone_matrix->getA()[current_city][next_city] += delta;
		}
	}
}

void updatePheromoneTrailQAS(vector<Ant*> ants, adjacencyMatrix<double>* pheromone_matrix, double pheromone_quantity_to_add, adjacencyMatrix<unsigned int>* adjMatrix) {
	for (Ant* ant : ants) {
		double delta;
		for (unsigned int i = 0; i < ant->visitedCities.size() - 1; i++) { //nie dodaje feromonu do sciezki powrotnej MOZLIWE ZE TO ŹLE - juz dodaje
			
			unsigned int current_city = ant->visitedCities[i];
			unsigned int next_city = ant->visitedCities[i + 1];
			delta = (double)pheromone_quantity_to_add / (double)adjMatrix->getA()[current_city][next_city];
			pheromone_matrix->getA()[current_city][next_city] += delta;
		}
	}
}

void updateBestSolution(vector<Ant*> ants, int& best_cost, vector<unsigned int>& best_path) {
	for (Ant* ant : ants) {
		if (ant->current_cost < best_cost) {
			best_cost = ant->current_cost;
			best_path = ant->visitedCities;
			//*******\ test
				//for (int k : best_path) std::cout << k;
				cout << " (" << best_cost << ") " << endl;
				cout << endl;
			//*******/
		}
	}

}

vector<unsigned int> nearestNeighbor(adjacencyMatrix<unsigned int>* adjMatrix, int startCity) {
	vector<unsigned int> path;
	vector<bool> visited(adjMatrix->getVertexNumber(), false);
	int currentCity = startCity;
	path.push_back(currentCity);
	visited[currentCity] = true;
	for (int i = 0; i < adjMatrix->getVertexNumber() - 1; i++) {
		int nextCity = -1;
		unsigned int shortestDistance = INT_MAX;
		for (int j = 0; j < adjMatrix->getVertexNumber(); j++) {
			if (!visited[j] && adjMatrix->getA()[currentCity][j] < shortestDistance) {
				nextCity = j;
				shortestDistance = adjMatrix->getA()[currentCity][j];
			}
		}
		currentCity = nextCity;
		path.push_back(currentCity);
		visited[currentCity] = true;
	}
	path.push_back(path[0]);
	return path;
}

int computeCycleCost(vector<unsigned int> cycle, adjacencyMatrix<unsigned int>* adjMatrix) {
	int totalCost = 0;
	for (unsigned int i = 0; i < cycle.size() - 1; i++) {
		int currentCity = cycle[i];
		int nextCity = cycle[i + 1];
		totalCost += adjMatrix->getA()[currentCity][nextCity];
	}
	return totalCost;
}

bool stoppingCondition(int bestCandidateSolution, int solution, int error) {
	//cout << std::dec << bestCandidateSolution << " <= " << solution << " + (" << error << " / 100) * " << solution << " (" << (solution + ((float)error / 100) * (float)solution) << ")";
	if (bestCandidateSolution <= (solution + ((float)error / 100) * (float)solution)) return false;
	else return true;
}


void antThread(Ant* ant, int num_cities, adjacencyMatrix<double>* pheromoneMatrix, adjacencyMatrix<unsigned int>* adjMatrix, int alpha, int beta) {
	for (int j = 0; j < num_cities - 1; j++) {
		ant->chooseNextCity(pheromoneMatrix, adjMatrix, alpha, beta);
		ant->visitCity(ant->next_city, adjMatrix);
	}
	ant->visitCity(ant->visitedCities.front(), adjMatrix); //Solution include the edge that completes the cycle.
	//for (int k : ant->visitedCities) std::cout << k << " ";
	//cout << "(" << computeCycleCost(ant->visitedCities, adjMatrix) << ") " << endl;
}

vector<unsigned int> AntColonyOptimization(int num_cities, int num_ants, int alpha, int beta, double pheromone_evaporation_rate, double initial_pheromone, double pheromone_quantity_to_add, int best_in_tsp, adjacencyMatrix<unsigned int>* adjMatrix) {  //heurystyka visibility
	int best_cost = 999999;
	vector<unsigned int> best_path;
	best_path.clear();

	adjacencyMatrix<double>* pheromoneMatrix = new adjacencyMatrix<double>(num_cities);		//initialize pheromone matrix
	for (int i = 0; i < num_cities; i++) {
		for (int j = 0; j < num_cities; j++) {
			if (i != j) pheromoneMatrix->getA()[i][j] = initial_pheromone;
		}
	}
	//*******\ test
		//pheromoneMatrix->showMatrix();
	//*******/
	while (stoppingCondition(best_cost, best_in_tsp, 20)) {		//Repeat for a given number of iterations
		vector<Ant*> ants;							//Create a set of ants
		ants.clear();
		for (int j = 0; j < num_ants; j++) {
			Ant* ant = new Ant(j);					//Assaigning each ant to different city. Note - it only works because num_ants = num_cities
			ants.push_back(ant);
		}

		std::vector<std::thread> threadPool;
		//Let ants construct a solution

		// Creating and using threads in AntColonyOptimization
		for (Ant* ant : ants) {
			if (threadPool.size() < MAX_THREADS) {
				threadPool.emplace_back(antThread, ant, num_cities, pheromoneMatrix, adjMatrix, alpha, beta);
			}
			else {
				// Wait for a thread to finish before adding a new one
				threadPool[0].join();
				threadPool.erase(threadPool.begin());
				threadPool.emplace_back(antThread, ant, num_cities, pheromoneMatrix, adjMatrix, alpha, beta);
			}
		}

		// Join remaining threads
		for (std::thread& thread : threadPool) {
			thread.join();
		}

		// Clear the thread pool for the next iteration
		threadPool.clear();
		//Update the pheromone trail
		updatePheromoneDecay(pheromoneMatrix, RHO);		//pheromone evaporation
		updatePheromoneTrailCAS(ants, pheromoneMatrix, pheromone_quantity_to_add);	//CAS pheromone update

		//*******\ test
			//pheromoneMatrix->showMatrix();
		//*******/

		//Update the best solution found so far
		updateBestSolution(ants, best_cost, best_path);
	}
	return best_path;
}

//void antThread(Ant* ant, int num_cities, adjacencyMatrix<double>* pheromoneMatrix, adjacencyMatrix<unsigned int>* adjMatrix, int alpha, int beta) {
//	for (int j = 0; j < num_cities - 1; j++) {
//		ant->chooseNextCity(pheromoneMatrix, adjMatrix, alpha, beta);
//		ant->visitCity(ant->next_city, adjMatrix);
//	}
//	ant->visitCity(ant->visitedCities.front(), adjMatrix); //Solution include the edge that completes the cycle.
//	//for (int k : ant->visitedCities) std::cout << k << " ";
//	//cout << "(" << computeCycleCost(ant->visitedCities, adjMatrix) << ") " << endl;
//}
//
//vector<unsigned int> AntColonyOptimization(int num_cities, int num_ants, int alpha, int beta, double pheromone_evaporation_rate, double initial_pheromone, double pheromone_quantity_to_add, int best_in_tsp, adjacencyMatrix<unsigned int>* adjMatrix) {  //heurystyka visibility
//	int best_cost = 999999;
//	vector<unsigned int> best_path;
//	best_path.clear();
//
//	adjacencyMatrix<double>* pheromoneMatrix = new adjacencyMatrix<double>(num_cities);		//initialize pheromone matrix
//	for (int i = 0; i < num_cities; i++) {
//		for (int j = 0; j < num_cities; j++) {
//			if (i != j) pheromoneMatrix->getA()[i][j] = initial_pheromone;  
//		}
//	}
//	//*******\ test
//		//pheromoneMatrix->showMatrix();
//	//*******/
//	while (stoppingCondition(best_cost, best_in_tsp, 20)) {		//Repeat for a given number of iterations
//		vector<Ant*> ants;							//Create a set of ants
//		ants.clear();
//		for (int j = 0; j < num_ants; j++) {
//			Ant* ant = new Ant(j);					//Assaigning each ant to different city. Note - it only works because num_ants = num_cities
//			ants.push_back(ant);
//		}
//
//		vector<std::thread> threads;
//		//Let ants construct a solution
//
//		for (Ant* ant : ants) {		
//			threads.emplace_back(antThread, ant, num_cities, pheromoneMatrix, adjMatrix, alpha, beta);
//			//for (int j = 0; j < num_cities - 1; j++) {
//			//	ant->chooseNextCity(pheromoneMatrix, adjMatrix, alpha, beta);
//			//	ant->visitCity(ant->next_city, adjMatrix);
//			//}
//			//ant->visitCity(ant->visitedCities.front(), adjMatrix); //Solution include the edge that completes the cycle.
//
//			//*******\ test pokaz mrowke
//				/*for (int k : ant->visitedCities) std::cout << k << " ";
//				cout << "(" << computeCycleCost(ant->visitedCities, adjMatrix) << ") " << endl;*/
//			//*******/
//		}
//
//		//for (Ant* ant : ants) {
//
//		//	//*******\ test pokaz mrowke
//		//		for (int k : ant->visitedCities) std::cout << k << " ";
//		//		cout << "(" << computeCycleCost(ant->visitedCities, adjMatrix) << ") " << endl;
//		//	//*******/
//		//}
//
//		for (std::thread& thread : threads) {
//			thread.join(); // Czekaj na zakończenie wszystkich wątków
//		}
//		//Update the pheromone trail
//		updatePheromoneDecay(pheromoneMatrix, RHO);		//pheromone evaporation
//		updatePheromoneTrailCAS(ants, pheromoneMatrix, pheromone_quantity_to_add);	//CAS pheromone update
//
//		//*******\ test
//			//pheromoneMatrix->showMatrix();
//		//*******/
//
//		//Update the best solution found so far
//		updateBestSolution(ants, best_cost, best_path);
//	}
//	return best_path;
//}

//--------------------------------------without threads-----------------------------------------------------

//vector<unsigned int> AntColonyOptimization(int num_cities, int num_ants, int alpha, int beta, double pheromone_evaporation_rate, double initial_pheromone, double pheromone_quantity_to_add, int best_in_tsp, adjacencyMatrix<unsigned int>* adjMatrix) {  //heurystyka visibility
//	int best_cost = 999999;
//	vector<unsigned int> best_path;
//	best_path.clear();
//
//	adjacencyMatrix<double>* pheromoneMatrix = new adjacencyMatrix<double>(num_cities);		//initialize pheromone matrix
//	for (int i = 0; i < num_cities; i++) {
//		for (int j = 0; j < num_cities; j++) {
//			if (i != j) pheromoneMatrix->getA()[i][j] = initial_pheromone;
//		}
//	}
//	//*******\ test
//		//pheromoneMatrix->showMatrix();
//	//*******/
//	while (stoppingCondition(best_cost, best_in_tsp, 20)) {		//Repeat for a given number of iterations
//		vector<Ant*> ants;							//Create a set of ants
//		ants.clear();
//		for (int j = 0; j < num_ants; j++) {
//			Ant* ant = new Ant(j);					//Assaigning each ant to different city. Note - it only works because num_ants = num_cities
//			ants.push_back(ant);
//		}
//		//Let ants construct a solution
//		for (Ant* ant : ants) {
//			for (int j = 0; j < num_cities - 1; j++) {
//				ant->chooseNextCity(pheromoneMatrix, adjMatrix, alpha, beta);
//				ant->visitCity(ant->next_city, adjMatrix);
//			}
//			ant->visitCity(ant->visitedCities.front(), adjMatrix); //Solution include the edge that completes the cycle.
//
//			//*******\ test pokaz mrowke
//				/*for (int k : ant->visitedCities) std::cout << k << " ";
//				cout << "(" << computeCycleCost(ant->visitedCities, adjMatrix) << ") " << endl;*/
//				//*******/
//		}
//		//Update the pheromone trail
//		updatePheromoneDecay(pheromoneMatrix, RHO);		//pheromone evaporation
//		updatePheromoneTrailCAS(ants, pheromoneMatrix, pheromone_quantity_to_add);	//CAS pheromone update
//
//		//*******\ test
//			//pheromoneMatrix->showMatrix();
//		//*******/
//
//		//Update the best solution found so far
//		updateBestSolution(ants, best_cost, best_path);
//	}
//	return best_path;
//}

int main() {
	//---------------------------------------------------------------\ potrzebne zainicjowania i wywołania funkcji
	srand(time(NULL));
	readConfig("config.ini");
	string filename(configValues[2]);
	fstream file;
	file.open(filename, std::ofstream::out | std::ofstream::trunc);

	double totalTime = 0;
	vector<short int> vectorOfVisitedVertexes;
	adjacencyMatrix<unsigned int>* adjMatrix = new adjacencyMatrix<unsigned int>();
	adjMatrix = adjMatrix->file_read_graphDirected(configValues[0]);
	if (adjMatrix->getVertexNumber() < 30) adjMatrix->showMatrix();
	//---------------------------------------------------------------/

	//===============================================================\ zmienne do AC
	vector<unsigned int> bestPath;
	int bestWeight;
	int num_cities = adjMatrix->getVertexNumber();
	int num_ants = num_cities;
	int estimatedCost = computeCycleCost(nearestNeighbor(adjMatrix, 0), adjMatrix);  // szacowana długosc trasy uzywajac najblizszego sasiada
	double t0 = (double) num_cities / (double) estimatedCost;						 // poczatkowa ilosc feromonu
	double Q = estimatedCost;
	int bestInTSP = readSolution(configValues[0]);
	std::vector<std::thread> threads;
	//===============================================================/

	
	for (int i = 0; i < stoi(configValues[1]); i++) {
		auto start = std::chrono::high_resolution_clock::now();
		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\ tu kod którego czas chcemy zmierzyć

		bestPath = AntColonyOptimization(num_cities, num_ants, ALPHA, BETA, RHO, t0, 100, bestInTSP, adjMatrix);

		//^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^/
		auto finish = std::chrono::high_resolution_clock::now();
		std::chrono::duration<double> elapsed = finish - start;
		totalTime += elapsed.count();
		for (int i : bestPath) file << i << " ";
		file << "(" << computeCycleCost(bestPath, adjMatrix) << ") ";
		if (file.is_open()) file << "Elapsed time: " << elapsed.count() << " s\n";
		//wyswietlenie wyniku
		bestWeight = computeCycleCost(bestPath, adjMatrix);
		std::cout << endl << "Najlepsza sciezka: ";
		for (int i : bestPath) std::cout << i << " ";
		std::cout << "(" << bestWeight << ") ";
		std::cout << endl;
	}
	file << "sredni czas: " << totalTime / stoi(configValues[1]) << " s";
	
	

	

	

}