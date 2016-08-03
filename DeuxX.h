#ifndef DEUXX_H
#define DEUXX_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include "Fonction.h"



	class DeuxX : public Fonction
	{
	protected:

	public:
				DeuxX() {};
		// value of f(x)
		 double get_value( double x);
		 double get_value( std::vector<double>  x);
		 std::vector<double> get_value( std::vector<double>  x, int d);
		 std::vector<std::vector<double> > get_value( std::vector<double>  x, int d, int r);
		// value of \nabla f(x)
		 std::vector<double> grad(std::vector<double> x);
 		double grad(double x);
		 // value of \D2 f(x)
		 std::vector<std::vector<double> > laplacian(std::vector<double> x);
		 double laplacian(double x);
	};

#endif

