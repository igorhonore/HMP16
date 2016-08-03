#define _USE_MATH_DEFINES


#include <cmath>
#include <iostream>
#include <vector>
#include "Bernoulli.h"

using namespace::std;
double Bernoulli::get_value()
{
	int i;
	double U;
	double temp;
		temp=Unif(0,1);
		U=(temp<=0.5)-(temp>0.5);
	return U;
	
}

std::vector<double> Bernoulli::get_value (int d)
{
	std::vector<double> Ud;
	double U;
	int i;
	for(i=1;i<d+1;i++)
	{
		U=Unif(0,1);
		Ud.push_back((U<=0.5)-(U>0.5));
	}
	return Ud;
}
	
double Bernoulli::Unif(double a, double b)
{
	return (b-a)*rand() /(double)RAND_MAX + a;
}

double Bernoulli::Normal(double mu, double sigma)
{
	double U1, U2;
	double result;
	U1 = Unif();
	U2 = Unif();
	result = std::sqrt(-2 * std::log(U1))*std::cos(2 * M_PI*U2);
	return result;
}

std::vector<double> Bernoulli::Normalstd_d(int d)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normal());
	}
	return result;
}

std::vector<double> Bernoulli::Normaliid_d(int d, double mu, double sigma)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normalstd_d(d)[i]*sigma+mu);
	}
	return result;
}
