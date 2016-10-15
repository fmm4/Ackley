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
	//Representacao//
	//geneValues = <X1,X2,...,X30>, valores do gene que vao ser utilizados pela funcao de ackley.
	//mutationSteps = <Q1,Q2,...,Q30>, valores de desvio padrao que vao ser utilizados para modificar os genes.
	//mutationAngle = nao usado, iria indicar o angulo se nossas variaveis fossem relacionadas. Nao sei se e o caso.
	//globalLR = LR global (aula 10, pg9)
	//localLR = LR local (aula 10, pg9)
	//c = deprecated, nao tem utilidade mas pode vir a ser util mais tarde
	//gen = gerador de valores aleatorios

	vector<double> geneValues;
	vector<double> mutationSteps;
	vector<double> mutationAngles;
	double globalLR;
	double localLR;
	double c = 0.8;
	default_random_engine gen;
	
	//deprecated//
	//Achava que era usando isso, ignorar. pode ser util mais tarde//
	void raise_mutationStep(int mutationStepPosition)
	{
		mutationSteps[mutationStepPosition] = mutationSteps[mutationStepPosition] / c;
	}
	void lower_mutationStep(int mutationStepPosition)
	{
		mutationSteps[mutationStepPosition] = mutationSteps[mutationStepPosition] * c;
	}
	//

	//Mutacao//
	//Varia o valor da mutacao baseados nos slides do professor. Mutaciona primeiro o passo de mutacao para depois modificar o gene//
	//Verificar aula 10 pg 9 para saber como funciona//
	void mutateMutationSteps(int position)
	{
		normal_distribution<double> distribution(0, 1);
		double value = mutationSteps[position] * exp(globalLR*distribution(gen) + localLR*distribution(gen));
		if (value >= 13 + 2){
			value = 13 + 2;
		}
		else if (value <= -13 - 2)
		{
			value = -13 - 2;
		}
		mutationSteps[position] = value;
	}

	//mutaciona os valores do gene especificado, verificar slides da aula 10, pg 9
	void mutateGenes(int i)
	{
		normal_distribution<double> distribution(0, 1);
		geneValues[i] = geneValues[i] + mutationSteps[i] * distribution(gen);
	}

public:
	double fitness = -1;


	//construtor de genoma//
	//automaticamente gera aleatoriamente todos os valores//
	ackleyEEGenotype(){
		//Especificado no documento, tamanho 30//
		int size = 30;

		normal_distribution<double> distribution(0, 1);

		for (int i = 0; i < size; i++)
		{
			geneValues.push_back((rand() % 31)-(12+3));
			mutationSteps.push_back(distribution(gen));
			mutationAngles.push_back(0);
		}

		//ajusta a taxa de aprendizado, ver slides da aula 10, pg 9
		startLearningRatio(30);
	};

	//retorna os gene da posicao no individuo//
	double getGene(int i)
	{
		return geneValues[i];
	}

	//retorna todos os genes como um array de doubles//
	vector<double> getGenome()
	{
		return geneValues;
	}

	//retorna o passo de mutacao do usuario//
	double getMutStep(int i)
	{
		return mutationSteps[i];
	}

	//funcao de ackley especificada pelos slides do professor.//
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

	//nome descritivo//
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

	//Fitness ainda nao especificado, no momento usa o proprio valor do ackley para isso//
	//talvez especificar um algoritmo especifico?//
	void setFitness(double value)
	{
		fitness = value;
	}

	//Filho de dois pais//
	//Media dos valores de cada um//
	ackleyEEGenotype binary_local_child(ackleyEEGenotype parent)
	{
		ackleyEEGenotype midpoint = ackleyEEGenotype();
		for (int i = 0; i < geneValues.size(); i++)
		{
			midpoint.geneValues[i] = (parent.geneValues[i] + geneValues[i]) / 2;
			midpoint.mutationSteps[i] = (parent.mutationSteps[i] + mutationSteps[i]) / 2;

		}

		return midpoint;
	}


	//Prole Global//
	//Tirado dos slides: Escolhe dois pais aleatoriamente para criar um filho onde cada gene e formado da media do gene de dois pais aleatorios.//
	ackleyEEGenotype global_child(vector<ackleyEEGenotype> population)
	{
		ackleyEEGenotype global_child = ackleyEEGenotype();

		for (int i = 0; i < geneValues.size(); i++)
		{
			global_child.geneValues[i] = (population[rand() % 30].geneValues[i] + population[rand() % 30].geneValues[i]) / 2;
			global_child.mutationSteps[i] = (population[rand() % 30].mutationSteps[i] + population[rand() % 30].mutationSteps[i]) / 2;

		}

		return global_child;
	}


	//Mutacao//
	//Mutaciona cada gene e seu respectivo passo de mutacao (simbolo desvio padrao), checar slide 9 do professor.//
	void applyMutation(int chance){
		for (int i = 0; i < mutationSteps.size(); i++)
		{
			if ((rand() % 100) < chance )
			{
				mutateMutationSteps(i);
				mutateGenes(i);
			}
		}
		getAckley();
	}

	//Learning Ratio//
	//Especificado nos slides da aula 11, nao sei se esta correto, mas foi dessa forma que entendi. Utilizado para mutacao unicas para cada gene.//
	void startLearningRatio(double n)
	{
		globalLR = 1/ (pow(2*n,1/2));
		localLR = 1 / (pow((pow(2 * n, 1 / 2)), 1 / 2));
	}

};

//Estrutura que permite a utilizacao de sort com vetores de objetos especificados pelo usuario, no caso, ackleyEEGenotype
struct less_than_key
{
	inline bool operator() (const ackleyEEGenotype& genotype1, const ackleyEEGenotype& genotype2)
	{
		return (genotype1.fitness < genotype2.fitness);
	}
};
