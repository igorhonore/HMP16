#ifndef FONCTION_H
#define FONCTION_H

#include <cmath>
#include <iostream>
#include <vector>
#include <iterator>
#include <stdlib.h>


	class Fonction
	{
	protected:

	public:
		// value of f(x)
		virtual double get_value( double x)=0; 
		virtual double get_value( std::vector<double>  x)=0;
		virtual std::vector<double> get_value( std::vector<double>  x, int d)=0;
		virtual std::vector<std::vector<double> > get_value( std::vector<double>  x, int r,int d)=0;
		virtual double get_value( double *  x)=0;
		virtual double * get_value( double *  x, int d)=0;
		virtual double ** get_value( double * x, int r,int d)=0;
		// value of \nabla f(x)
		virtual double grad(double x)=0;
		virtual std::vector<double> grad(std::vector<double> x)=0;
		virtual double *grad(double * x)=0;
		// value of D2f(x) 
		virtual double laplacian(double x)=0;
		virtual std::vector<std::vector<double> > laplacian(std::vector<double>x)=0;
		virtual double ** laplacian(double * x)=0;
		// value of D3f(x)
		virtual double D3(double x)=0;
		virtual double *** D3(double* x)=0;
		// value of D4f(x)
		virtual double D4(double x)=0;




	};


#endif
