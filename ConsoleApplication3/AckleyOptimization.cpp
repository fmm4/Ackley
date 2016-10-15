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
using namespace std;

//Executa EE funcao de ackley
void ackleyOptimizationEE(int pop, int children);
//Gera aleatoriamente uma populacao de genes.
vector<ackleyEEGenotype> generatePopulation(int size_population);

int _tmain(int argc, _TCHAR* argv[])
{
	int test;
	srand(time(NULL));

	ackleyOptimizationEE(30, 200);


	cin >> test;

}


//Algoritmo que faz EE//
void ackleyOptimizationEE(int pop, int children)
{
	//Cria a populacao//
	vector<ackleyEEGenotype> populacao = generatePopulation(pop);

	//Reordena para que sempre o menor esteja mais abaixo//
	sort(populacao.begin(), populacao.end(), less_than_key());

	//Chance para mutacao//	
	double mutationChance = 80;

	//comeco do algoritmo
	for (int generation = 1; generation < 2000 + 3000; generation++)
	{
		//inicializa a proxima geracao//
		vector<ackleyEEGenotype> next_generation;

		//190 dos filhos novos serao feitos do cruzamento de dois pais aleatorios//
		for (int local_children = 0; local_children < 190; local_children++)
		{
			ackleyEEGenotype child = populacao[rand() % 30].binary_local_child(populacao[rand() % 30]);
			child.applyMutation(mutationChance);
			next_generation.push_back(child);
		}
		//10 filhos sao feitos da media de 2 pais aleatorios para cada gene, chamados de filhos globais//
		for (int local_children = 0; local_children < 10; local_children++)
		{
			ackleyEEGenotype child = populacao[0].global_child(populacao);
			child.applyMutation(mutationChance);
			next_generation.push_back(child);
		}

		//Forma a nova populacao. Essa nova populacao e formado somente dos filhos, uma vez que os slides das aulas indicam que e ruim usar os pais devido a elitismo//
		sort(next_generation.begin(), next_generation.end(), less_than_key());

		next_generation.resize(30);

		populacao.swap(next_generation);

		cout << "Maximum Fitness: " << populacao[0].getFitness() << endl;
	}


	//Debugging, saber se esta melhorando//
	cout << "Best Fitness: " << populacao[0].getFitness() << endl;
	
	vector<double> gene = populacao[0].getGenome();
	cout << "Mean: " << accumulate(gene.begin(), gene.end(), 0.0 / gene.size());
	cout << endl << "[";
	for (int i = 0; i < gene.size(); i++)
	{
		cout << gene[i] << ",";
	}
	cout << "]";

}


vector<ackleyEEGenotype> generatePopulation(int size_population){
	vector<ackleyEEGenotype> population;

	for (int i = 0; i < size_population; i++)
	{
		population.push_back(ackleyEEGenotype());
	}

	return population;
}