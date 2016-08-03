#define _USE_MATH_DEFINES

//#include "Algorithm.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "Gaussian.h"

using namespace::std;
double Gaussian::get_value(double mu, double sigma)
{
	double U;
	U=Normal(mu,sigma);
	return U;
	
}
double Gaussian::get_value()
{
	double U;
	U=Normal(0.,1.);
	return U;
	
}
std::vector<double> Gaussian::get_value (int d)
{
	std::vector<double> Ud;
	Ud=Normalstd_d(d);
	return Ud;
}
	
double Gaussian::Unif(double a, double b)
{
	return (b-a)*rand() /(double)RAND_MAX + a;
}

double Gaussian::Normal(double mu, double sigma)
{
	double U1, U2;
	double result;
	U1 = Unif();
	U2 = Unif();
	result = std::sqrt(-2 * std::log(U1))*std::cos(2 * M_PI*U2);
	return result;
}

std::vector<double> Gaussian::Normalstd_d(int d)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normal());
	}
	return result;
}

std::vector<double> Gaussian::Normaliid_d(int d, double mu, double sigma)
{
	std::vector<double> result;
	for (int i = 0; i < d; i++)
	{
		result.push_back(Normalstd_d(d)[i]*sigma+mu);
	}
	return result;
}

