#include "BLineaire.h"
#include <cmath>
#include <iostream>
#include <vector>
#include "Fonction.h"
#include <stdlib.h>

using namespace::std;

double BLineaire::D4(double x)
{
	return 0;
}
	
double * BLineaire::get_value( double*  x, int d)
{	
	int i;
	//d=sizeof(x);
	double* temp;
	temp= (double*)malloc(d*sizeof(double));
	for(i=0;i<d;i++)
	{
		//for(j=1;j<d+1;j++)
			//temp[i][j]=0;
		temp[i]=-0.5*x[i];
	}
	return temp;
//	free(temp);
}

	double BLineaire::get_value(double x)
{
	return -0.5*x;
}

	double BLineaire::grad(double x)
	{
		return -0.5;
	}

	double BLineaire::laplacian(double x)
	{
		return 0;
	}

	double BLineaire::get_value(std::vector<double> x)
	 {
		 double result;
		  double temp=0;
		  		 int i;
		int n;
			 n=sizeof(x);
		for(i=1;i<n+1;i++)
			temp+=fabs(x[i]);
			result=-0.5*temp;
		 return result;
	 }
		 std::vector<double> BLineaire::grad(std::vector<double> x)
	{
		 std::vector<double> result=x;
		 		  //double temp=0;
		//for(i=1;i<n+1;i++)
			//temp+=fabs(x[i]);
			//result=-0.5*temp;
		 return result;
	 }
		 std::vector<std::vector<double> > BLineaire::laplacian(std::vector<double> x)
	{
		 std::vector< std::vector<double> > result;
		 return result;
	 }
		 std::vector<std::vector<double> > BLineaire::get_value( std::vector<double>  x, int d, int r)
		 {
			 std::vector<std::vector<double> > result;
			 return result;
		 }
		 
	std::vector<double > BLineaire::get_value( std::vector<double> & x, int d)
		 {
			 std::vector<double> result;
			 double temp=0;
			 int i;
			 int n;
			 //n=sizeof(x);
			 result=x;
		for(i=1;i<d+1;i++)
			temp+=fabs(x[i]);
		for(i=1;i<d+1;i++)
			result[i]=-0.5*temp;//*(2+cos(x[i]));
			 return result;
		 }
		 	std::vector<double > BLineaire::get_value( std::vector<double>  x, int d)
		 {
			 std::vector<double> result;
			 double temp=0;
			 int i;
			 int n;
			 n=sizeof(x);
			 result=x;
		for(i=1;i<d+1;i++)
			temp+=fabs(x[i]);
		for(i=1;i<d+1;i++)
			result[i]=-0.5*temp;//*(2+cos(x[i]));
			 return result;
		 }
		 //std::vector<std::vector<double> > XBeta::get_value( std::vector<double>  x, int d, int r)
	 //{
		 //std::vector<std::vector<double> > result;
		 //return result;
	 //}
	 //
//std::vector<double > BLineaire::get_value( std::vector<double>  x, int d)
	 //{
		 //std::vector<double> result;
		 //return result;
	 //}
	 
	double ** BLineaire::get_value( double * x, int r,int d)
	{
		int i;
		
			double ** result;
			return result;
		for (i = 0; i < d; i++)
		{
			free(result[i]);
		} 
		free(result);
	}
		double *BLineaire::grad(double * x)
		{
				double * result;
				result= (double*)malloc(sizeof(x)*sizeof(double));
				int i;
				for(i=0;i<sizeof(x)/sizeof(double);i++)
					result[i]=-0.5*x[i];
				return result;
		//for (i = 0; i < sizeof(x); i++)
		//{
			//free(result[i]);
		//} 
		free(result);
		}

		double ** BLineaire::laplacian(double * x)
		{
			int i;
			double ** result;
			return result;
				for (i = 0; i < sizeof(x)/sizeof(double); i++)
	{
		free(result[i]);
	} 
	free(result);
		}
		double BLineaire::get_value( double *  x)
		{
			double result;
			return result;
		}
		double BLineaire::D3(double x)
		{
			return 0;
		}
		double *** BLineaire::D3(double * x)
		{
			int i,j;
			double *** result;
			return result;
	for (i = 0; i < sizeof(x)/sizeof(double); i++)
	{
		for (j = 0; j < sizeof(x)/sizeof(double); j++)
		free(result[i][j]);
		free(result[i]);
	} 
	free(result);
		}
