#pragma once
#include <map>
#define _USE_MATH_DEFINES
#include <algorithm>
#include <math.h>
#include <vector>
#include <ctime>
#include <cmath>
#include <time.h>
#include <random>
using namespace std;

class ackleyEEGenotype
{
private:
	vector<double> geneValues;
	vector<double> mutationSteps;
	vector<double> mutationAngles;

public:
	double fitness = -1;

	ackleyEEGenotype(){
		//Especificado no documento
		int size = 30;

		default_random_engine gen;
		normal_distribution<double> distribution(0, 1);

		for (int i = 0; i < size; i++)
		{
			geneValues.push_back((rand() % 31)-(12+3));
			mutationSteps.push_back(distribution(gen));
			mutationAngles.push_back(0);

		}
	};

	double getGene(int i)
	{
		return geneValues[i];
	}

	double getMutStep(int i)
	{
		return mutationSteps[i];
	}

	double getAckley(){
		double c1 = 20, c2 = 0.2, c3 = 2 * M_PI;
		double n = geneValues.size();
		double sumCube = 0;
		for (int i = 0; i < n; i++){
			sumCube += pow(geneValues[i], 2);
		}
		double sumCos = 0;
		for (int i = 0; i < n; i++)
		{
			sumCos += cos(c3*geneValues[i]);
		}

		double result = (-c1) + exp((-c2)*sqrt(sumCube / n)) - exp(sumCos / n) + c1 + 1;

		setFitness(result);
		return result;
	}

	double getFitness(){
		if (fitness == -1)
		{
			getAckley();
			return fitness;
		}
		else{
			return fitness;
		}
	}

	void setFitness(double value)
	{
		fitness = value;
	}

};

struct less_than_key
{
	inline bool operator() (const ackleyEEGenotype& genotype1, const ackleyEEGenotype& genotype2)
	{
		return (genotype1.fitness < genotype2.fitness);
	}
};
