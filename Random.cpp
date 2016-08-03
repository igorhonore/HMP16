#define _USE_MATH_DEFINES

#include <cmath>
#include <iostream>
#include <vector>
#include "Random.h"

using namespace::std;
	

double Random::Unif(double a, double b)
{
	return (b-a)*rand() /(double)RAND_MAX + a;
}

double Random::Normal(double mu, double sigma)
{
	double U1, U2;
	double result;
	U1 = Unif();
	U2 = Unif();
	result = std::sqrt(-2 * std::log(U1))*std::cos(2 * M_PI*U2);
	return mu+sigma*result;
}

std::vector<double> Random::Normalstd_d(int d)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normal());
	}
	return result;
}

std::vector<double> Random::Normaliid_d(int d, double mu, double sigma)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normalstd_d(d)[i]*sigma+mu);
	}
	return result;
}

