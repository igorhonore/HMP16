#include "Sinus.h"
#include <cmath>
#include <iostream>
#include <vector>
//#include "Fonction.h"
#include <stdlib.h>

using namespace::std;

	double ** Sinus::get_value( double*  x, int d, int r)
{	
	int i,j, m;
	m=min(d,r);
	d=sizeof(x);
	double** temp;
	temp= (double**)malloc(d*sizeof(double*));
	for (i = 1; i < m+1; i++)
	{
		temp[i] = (double*) calloc(d , sizeof(double));
	} 
	for(i=1;i<d+1;i++)
	{
		temp[i][i]=sin(x[i]);
	}
	return temp;
}

	double Sinus::get_value(double x)
{
	return pow(0.01*sin(2.*x)/2.+x*sin(x),2);
}

	double Sinus::grad(double x)
	{
		return cos(x);
	}

	double Sinus::laplacian(double x)
	{
		return -sin(x);
	}

	 double Sinus::get_value(std::vector<double> x)
	 {
		 double result, temp=0;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			temp+=x[i];
		result=sin(temp);
		 return result;
	 }
		 std::vector<double> Sinus::grad(std::vector<double> x)
	{
		 std::vector<double> result;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			result.push_back(cos(x[i]));
		 return result;
	 }
		 std::vector<std::vector<double> > Sinus::laplacian(std::vector<double> x)
	{
		 std::vector< std::vector<double> > result;
		 int n,i,j;
		 n=sizeof(x);
			result.resize(n);
		 for(i=1;i<n+1;i++)
		 	result[i].resize(n);
		 for(i=1;i<n+1;i++)
		 {
			 for(j=1;j<n+1;j++)
			 {
				result[i][j]=0;
			}
			result[i][i]=-sin(x[i]);
		 }
		 return result;
	 }
		 std::vector<std::vector<double> > Sinus::get_value( std::vector<double>  x, int d, int r)
		 {
		std::vector<std::vector<double> > result;
		 int i,j;
		 result.resize(d);
		 for(i=1;i<d+1;i++)
		 	result[i].resize(r);
		 for(i=1;i<d+1;i++)
		 {
			 for(j=1;j<r+1;j++)
				{
					result[i][j]=1;
				}
			result[i][i]=sin(x[i]);
		}
		 return result;
		 }
		 
	std::vector<double > Sinus::get_value( std::vector<double>  x, int d)
		 {
			 std::vector<double> result;
			 int i;
			 for(i=1;i<d+1;i++)
				result.push_back(sin(x[i]));
			 return result;
		 }
		 
	double Sinus::D3(double x)
	{
		return -cos(x);
	}

		double Sinus::D4(double x)
	{
		return sin(x);
	}
	
double Sinus::get_value(double* x)
	 {
		 double result, temp=0;
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			temp+=x[i];
		result=sin(temp);
		 return result;
	 }
	 
double * Sinus::get_value(double* x, int d)
	 {
		 double temp=0;
		 double * result;
		 result=(double*)calloc(d,sizeof(double));
		 int i;
		 for(i=1;i<sizeof(x)+1;i++)
			result[i]=sin(x[i]);
		 return result;
	 }
double * Sinus::grad(double*)
	 {
		 double* result;
		 return result;
	 }
double** Sinus::laplacian(double*)
	 {
		 double** result;
		 return result;
	 }
double*** Sinus::D3(double*)
	 {
		 double*** result;
		 return result;
	 }

