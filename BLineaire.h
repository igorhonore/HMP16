#ifndef BLINEAIRE_H
#define BLINEAIRE_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include "Fonction.h"



	class BLineaire : public Fonction
	{
	protected:

	public:
		BLineaire() {};
		// value of f(x)
		double get_value( double x);
		double get_value( std::vector<double>  x);
		 double* get_value (double*  x, int d);		
		 std::vector<double> get_value( std::vector<double>&  x, int d);
		std::vector<double> get_value( std::vector<double>  x, int d);
		 std::vector<std::vector<double> > get_value( std::vector<double>  x, int d, int r);
		virtual double ** get_value( double * x, int r,int d);
		virtual double get_value( double *  x);
		// value of \nabla f(x)
		virtual double *grad(double * x);
		 double grad(double x);
		std::vector<double> grad(std::vector<double> x);
		// value of D2f(x) 
		 std::vector<std::vector<double> > laplacian(std::vector<double> x);		
		 double laplacian(double x);
		virtual double ** laplacian(double * x);
		// value of D3f(x)
		double D3(double x);
		virtual double *** D3(double* x);
		// value of D4f(x)
		virtual double D4(double x);
	};

#endif

