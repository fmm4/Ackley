// ConsoleApplication3.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <iostream>
#include "ackleyEEGenotype.h"
#include <ctime>
#include <algorithm>
#include <math.h>
#include <time.h>
using namespace std;

void ackleyOptimizationEE(int pop, int size, int children);

int _tmain(int argc, _TCHAR* argv[])
{
	int test;
	srand(time(NULL));

	ackleyOptimizationEE(30, 30, 200);


	cin >> test;

}

void ackleyOptimizationEE(int pop, int size, int children)
{
	vector<ackleyEEGenotype> population;
	for (int i = 0; i < size; i++)
	{
		ackleyEEGenotype newGeno = ackleyEEGenotype();
		population.push_back(newGeno);
	}

	for (int i = 0; i < 4; i++){
		cout << "Genotype " << i << endl << "Genomes:" << endl;
		
		for (int o = 0; o < 30; o++)
		{
			cout << " " << population[i].getGene(o);
		}
		for (int o = 0; o < 30; o++)
		{
			cout << " " << population[i].getMutStep(o);
		}
	}
}