// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include "ackleyEEGenotype.h"
#include <ctime>
#include <algorithm>
#include <numeric>
#include <math.h>
#include <time.h>
#include <set>
using namespace std;

//Executa EE funcao de ackley
void ackleyOptimizationEE(int pop, int children, int generations);
void ackleyOptimizationEECheat(int pop, int children, int generations);
//Gera aleatoriamente uma populacao de genes.
vector<ackleyEEGenotype> generatePopulation(int size_population);
//Matematica
double average(vector<double> const& v);
double standard_variation(vector<double> const& v, double average);
double average_distance(vector<ackleyEEGenotype> k);

int _tmain(int argc, _TCHAR* argv[])
{
	int test;
	srand(time(NULL));

	ackleyOptimizationEE(30, 200, 2000);
	ackleyOptimizationEECheat(30, 200, 2000);


	cin >> test;

}

bool compare(const ackleyEEGenotype& geno1, const ackleyEEGenotype& geno2)
{
	ackleyEEGenotype k = geno1;
	ackleyEEGenotype p = geno2;
	vector<double> first = k.getGenome();
	vector<double> second = p.getGenome();

	double difference = 0;
	for (int i = 0; i < 30; i++)
	{
		difference += abs(first[i] - second[i]);
	}

	if (difference >= 2 + 3){
		return false;
	}
	else{
		return true;
	}
}


vector<double> get_fitnesses(vector < ackleyEEGenotype> populacao)
{
	vector<double> fitnesses;
	for (int i = 0; i < populacao.size(); i++){
		fitnesses.push_back(populacao[i].getFitness());
	}
	return fitnesses;
}

vector<double> get_cheat_fitnesses(vector < ackleyEEGenotype> populacao)
{
	vector<double> fitnesses;
	for (int i = 0; i < populacao.size(); i++){
		fitnesses.push_back(populacao[i].distanceToBest());
	}
	return fitnesses;
}

double average_fitnesses(vector<double> v) {
	return 1.0 * accumulate(v.begin(), v.end(), 0.0) / v.size();
}

double standard_variation_fitnesses(vector<double> v, double average)
{
	double temp = 0;
	for (int i = 0; i < v.size(); i++)
	{
		temp += (v[i] - average) * (v[i] - average);
	}
	return sqrt(temp / v.size());
}


//Algoritmo que faz EE//
void ackleyOptimizationEE(int pop, int children, int generationIterations)
{
	cout.precision(13);

	//Impressao de Resultados
	//Fitness médio alcançado nas 30 execuções (média e desvio padrão);
	vector<double> avgFitnesses;
	vector<double> uniques;
	vector<double> maxFitnesses;
	vector<double> distances;

	double solucoes = 0;

	
	//Quantas geracoes//
	double generations = generationIterations;



	//comeco do algoritmo
	for (int exec = 0; exec < 30; exec++)
	{
		cout << exec;
		bool solucao = 0;
		//Chance para mutacao//	
		double mutationChance = 100;
		//Populacao
		double pop_size = 40;
		//cout << "Gen " << exec << endl;
		
		//Cria a populacao//
		vector<ackleyEEGenotype> populacao = generatePopulation(pop);

		double maxFitness = 9999;

		for (int generation = 1; generation < generations; generation++)
		{
			
			//inicializa a proxima geracao//
			vector<ackleyEEGenotype> next_generation;

			//Golden Ratio//
			//Se um quinto dos filhos sao melhores//
			

			if (mutationChance > 40)
			{
				mutationChance -= (30*2)/generations;
			}

			//190 dos filhos novos serao feitos do cruzamento de dois pais aleatorios//
			for (int local_children = 0; local_children < pop_size * 3; local_children++)
			{
				ackleyEEGenotype child = populacao[rand() % 30].binary_local_child(populacao[rand() % 30]);
				child.applyMutation(mutationChance);
				next_generation.push_back(child);
			}
			//10 filhos sao feitos da media de 2 pais aleatorios para cada gene, chamados de filhos globais//
			for (int local_children = 0; local_children < pop_size / (2+3); local_children++)
			{
				ackleyEEGenotype child = populacao[0].global_child(populacao);
				child.applyMutation(mutationChance);
				next_generation.push_back(child);
			}



			double maxFound = populacao[0].fitness;

			if (maxFitness > maxFound){
				maxFitness = maxFound;
			}

			next_generation.insert(next_generation.end(), populacao.begin(), populacao.end());
			sort(next_generation.begin(), next_generation.end(), less_than_key());
			populacao.swap(next_generation);
			populacao.resize(pop_size);
			if (populacao[0].isGlobalMin()){
				solucao = true;
			}
		}


		if (solucao){
			solucoes++;
		}
		vector<double> fitnesses = get_fitnesses(populacao);
	
		//populacao[0].debugGene();
		vector<double> distanc = get_cheat_fitnesses(populacao);
		

		double average = average_fitnesses(fitnesses);
		double stdrddv = standard_variation_fitnesses(fitnesses, average);

		/*cout << "Average: " << average << endl;
		cout << "Max Fitness: " << maxFitness << endl;
		cout << "Best Distance: " << populacao[0].distanceToBest() << endl;
		cout << "Average Distance: " << average_distance(populacao) << endl;*/

		avgFitnesses.push_back(average);
		uniques.push_back((double) populacao.size());
		maxFitnesses.push_back(maxFitness);
		distances.push_back(average_fitnesses(distanc));
	}


	//Debugging, saber se esta melhorando//
	
	double avgFitnessOnAllExec = average_fitnesses(avgFitnesses);
	double stdrdFitn = standard_variation_fitnesses(avgFitnesses, avgFitnessOnAllExec);
	double averageMaxes = average_fitnesses(maxFitnesses);
	double stdrdFitnMaxes = standard_variation_fitnesses(maxFitnesses, averageMaxes);
	double averageDistances = average_fitnesses(distances);
	double stdrdFitnDist = standard_variation_fitnesses(distances, averageDistances);

	cout << "Test Results [EE]:" << endl;
	cout << "Maximum Fitness Found: [AVERAGE: " << averageMaxes <<" STANDARD DEVIATION: " << stdrdFitnMaxes << "]" << endl;
	cout << "Fitness [AVERAGE: " << avgFitnessOnAllExec << " STANDARD DEVIATION: " << stdrdFitn << "]" << endl;
	cout << "Distances [AVERAGE: " << averageDistances << " STANDARD DEVIATION: " << stdrdFitnDist << "]" << endl;
	cout << "Solucoes " << solucoes / 30 << endl;
}

void ackleyOptimizationEECheat(int pop, int children, int generationIterations)
{
	cout.precision(13);

	//Impressao de Resultados
	//Fitness médio alcançado nas 30 execuções (média e desvio padrão);
	vector<double> avgFitnesses;
	vector<double> fitnesses;
	vector<double> minFitnesses;
	vector<double> averageDistances;

	vector<double> averageOldFitnesses;
	vector<double> stndrdOldFitnesses;

	//Quantas geracoes//
	double generations = generationIterations;


	double solucoes = 0;
	


	//comeco do algoritmo
	for (int exec = 0; exec < 30; exec++)
	{
		cout << exec;
		//Populacao
		double pop_size = 40;
		//cout << "Gen " << exec << endl;

		bool solucao = false;

		//Cria a populacao//
		vector<ackleyEEGenotype> populacao = generatePopulation(pop);

		double minFitness = 999;

		for (int generation = 1; generation < generations; generation++)
		{

			//inicializa a proxima geracao//
			vector<ackleyEEGenotype> next_generation;

			//Golden Ratio//
			//Se um quinto dos filhos sao melhores//


			//190 dos filhos novos serao feitos do cruzamento de dois pais aleatorios//
			for (int local_children = 0; local_children < pop_size * 3; local_children++)
			{
				ackleyEEGenotype father1 = populacao[rand() % 30];
				ackleyEEGenotype father2 = populacao[rand() % 30];
				ackleyEEGenotype child = father1.binary_local_child(father2);
				child.applyMutationCheat();

				if (child.getFitness() < father1.getFitness() || child.getFitness() < child.getFitness())
				{
					next_generation.push_back(child);
				}
			}
			//10 filhos sao feitos da media de 2 pais aleatorios para cada gene, chamados de filhos globais//
			for (int local_children = 0; local_children < pop_size / (2 + 3); local_children++)
			{
				ackleyEEGenotype child = populacao[0].global_child(populacao);
				child.applyMutationCheat();
				if (child.getFitness() < populacao[29].getFitness())
				{
					next_generation.push_back(child);
				}
			}


			double maxFound = populacao[0].getFitness();

			if (minFitness > maxFound){
				minFitness = maxFound;
			}

			next_generation.insert(next_generation.end(), populacao.begin(), populacao.end());
			sort(next_generation.begin(), next_generation.end(), less_than_key());
			next_generation.resize(pop_size);


			populacao.swap(next_generation);
			if (populacao[0].isGlobalMin())
			{
				solucao = true;
			}
		}

		vector<double> fitnesses = get_cheat_fitnesses(populacao);
		vector<double> fitnessesOld = get_fitnesses(populacao);
		//populacao[0].debugGene();

		if (solucao)
		{
			solucoes++;
		}

		double average = average_fitnesses(fitnesses);
		double oldFitness = average_fitnesses(fitnessesOld);
		double stdrddv = standard_variation_fitnesses(fitnesses, oldFitness);

		/*cout << "Average: " << oldFitness << endl;
		cout << "Max Fitness: " << minFitness << endl;
		cout << "Best Distance: " << populacao[0].distanceToBest() << endl;
		cout << "Average Distance: " << average_distance(populacao) << endl;*/


		avgFitnesses.push_back(average);
		minFitnesses.push_back(minFitness);
		averageOldFitnesses.push_back(oldFitness);
	}


	double avgFitnessOnAllExec = average_fitnesses(avgFitnesses);
	double stdrdFitn = standard_variation_fitnesses(avgFitnesses, avgFitnessOnAllExec);
	double averageMaxes = average_fitnesses(minFitnesses);
	double stdrdFitnMaxes = standard_variation_fitnesses(minFitnesses, averageMaxes);
	double avg_old_average = average_fitnesses(averageOldFitnesses);
	double avg_old_stdrd = standard_variation_fitnesses(averageOldFitnesses,avg_old_average);


	cout << "Test Results [EE]:" << endl;
	cout << "Minimum Fitness Found: [AVERAGE: " << averageMaxes << " STANDARD DEVIATION: " << stdrdFitnMaxes << "]" << endl;
	cout << "Distances [AVERAGE: " << avgFitnessOnAllExec << " STANDARD DEVIATION: " << stdrdFitn << "]" << endl;
	cout << "Fitnesses [AVERAGE:  " << avg_old_average << " STANDARD: DEVIATION: " << avg_old_stdrd << "]" << endl;
	cout << "Solucoes: " << solucoes / generations << endl;
}

double average_distance(vector<ackleyEEGenotype> k)
{
	double distancesum = 0;
	for (int i = 0; i < k.size(); i++)
	{
		distancesum += k[i].distanceToBest();
	}
	return distancesum /= k.size();
}

vector<ackleyEEGenotype> generatePopulation(int size_population){
	vector<ackleyEEGenotype> population;

	for (int i = 0; i < size_population; i++)
	{
		ackleyEEGenotype k = ackleyEEGenotype();
		k.getAckley();
		population.push_back(k);
	
	}

	return population;
}