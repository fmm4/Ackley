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
//Matematica
double average(vector<double> const& v);
double standard_variation(vector<double> const& v, double average);

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
	cout.precision(4);
	//Cria a populacao//
	vector<ackleyEEGenotype> populacao = generatePopulation(pop);

	//Reordena para que sempre o menor esteja mais abaixo//
	sort(populacao.begin(), populacao.end(), less_than_key());

	//Chance para mutacao//	
	double mutationChance = 100;
	//Diminui o maximo que a mutacao pode chegar, para no final ficar mais focado em explotacao//
	double mutationCeiling = 100;
	//Quantas geracoes//
	double generations = 2000;

	//comeco do algoritmo
	for (int generation = 1; generation < generations; generation++)
	{
		//inicializa a proxima geracao//
		vector<ackleyEEGenotype> next_generation;

		int pop_size = 20 + 30;

		mutationChance--;

		if (mutationChance < 10){
			mutationChance = mutationCeiling;
			mutationCeiling-= 90/generations;
		}

		//190 dos filhos novos serao feitos do cruzamento de dois pais aleatorios//
		for (int local_children = 0; local_children < pop_size*(4); local_children++)
		{
			ackleyEEGenotype child = populacao[rand() % 30].binary_local_child(populacao[rand() % 30]);
			child.applyMutation(mutationChance);
			next_generation.push_back(child);
		}
		//10 filhos sao feitos da media de 2 pais aleatorios para cada gene, chamados de filhos globais//
		for (int local_children = 0; local_children < pop_size; local_children++)
		{
			ackleyEEGenotype child = populacao[0].global_child(populacao);
			child.applyMutation(mutationChance);
			next_generation.push_back(child);
		}

		//Forma a nova populacao. Essa nova populacao e formado somente dos filhos, uma vez que os slides das aulas indicam que e ruim usar os pais devido a elitismo//
		sort(next_generation.begin(), next_generation.end(), less_than_key());


		next_generation.resize(pop_size);


		populacao.swap(next_generation);
		
		vector<double> fitnesses;
		for (int i = 0; i < populacao.size(); i++){
			fitnesses.push_back(populacao[i].getFitness());
		}

		double averageFit = average(fitnesses);
		double standard_deviationFit = standard_variation(fitnesses,averageFit);
		cout << "Maximum Fitness: " << populacao[0].getFitness() << endl;
		cout << "Average Fitness: " << averageFit << endl;
		cout << "Standard Deviation: " << standard_deviationFit << endl;
	}


	//Debugging, saber se esta melhorando//
	
	cout << "Best Fitness: " << populacao[0].getFitness() << endl;

	vector<double> gene = populacao[0].getGenome();
	cout << "Mean: " << accumulate(gene.begin(), gene.end(), 0.0 / gene.size()) << endl;

	for (int i = 0; i < 30; i++){
		populacao[i].debugGene();
	}
}

double average(vector<double> const& v) {
	return 1.0 * accumulate(v.begin(), v.end(),0.0) / v.size();
}

double standard_variation(vector<double> const& v, double average)
{
	return 1.0 * (accumulate(v.begin(), v.end(),0.0-average*v.size()))/ v.size();
}


vector<ackleyEEGenotype> generatePopulation(int size_population){
	vector<ackleyEEGenotype> population;

	for (int i = 0; i < size_population; i++)
	{
		population.push_back(ackleyEEGenotype());
	}

	return population;
}