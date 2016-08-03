#ifndef GAUSSIAN_H
#define GAUSSIAN_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include <cstdlib>
#include "Random.h"

class Gaussian : public Random
{
protected :

public:
	Gaussian(){};
	double get_value();// return Normal(mu,sigma)
	double get_value(double mu, double sigma);// return Normalstd_d(d);
	std::vector<double> get_value (int d);//return Normal(0.,1.)
	double Unif(double a = 0.0, double b = 1.0);// return a uniform random variable defined on [a,b]
	std::vector<double> Normalstd_d(int d = 1);// return a Standard Gaussian vector 
	double Normal(double mu = 0, double sigma = 1);// return a Gaussian with mean=mu, and variance = sigma
	std::vector<double> Normaliid_d(int d, double mu=0, double sigma=1);// return a Gaussian vector mean=mu*I_d, and variance = sigma*I_d
	 
};

#endif
