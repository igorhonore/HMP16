#include "Constant.h"
#include <cmath>
#include <iostream>
#include <vector>
//#include "Fonction.h"
#include <stdlib.h>

using namespace::std;

double ** Constant::get_value( double*  x, int d, int r)
{	
	int i,j;
	double** temp;
	temp= (double**)malloc(d*sizeof(double*));
	for (i = 0; i < d; i++)
	{
		temp[i] = (double*) calloc(d , sizeof(double));
	} 
	for(i=0;i<d;i++)
	{
		for(j=0;j<d;j++)
			temp[i][j]=0;
		temp[i][i]=1;
	}
	return temp;
}

double Constant::get_value(double x)
{
	return 1;
}

double Constant::grad(double x)
{
	return 0;
}

double Constant::laplacian(double x)
{
	return 0;
}

double Constant::get_value(std::vector<double> x)
 {
	 return 1;
 }
std::vector<double> Constant::grad(std::vector<double> x)
{
	int n;
	n=sizeof(x);
	 std::vector<double> result(n);
	 int i;
	 for(i=0;i<n;i++)
		{
			result[i]=0;
		}
	 return result;
 }
std::vector<std::vector<double> > Constant::laplacian(std::vector<double> x)
{
	int i,j,n;
	 n=sizeof(x);
	 std::vector< std::vector<double> > result(n);
	 for(i=0;i<n;i++)
		for(j=0;j<n;j++)
			result[j].push_back(1);
	 return result;
 }
std::vector<std::vector<double> > Constant::get_value( std::vector<double>  x, int d, int r)
 {
	 std::vector<std::vector<double> > result(d);
	 int i,j;
	 for(i=0;i<d;i++)
	 {
		result[i].resize(r);
	}
	 for(i=0;i<d;i++)
	 {
		 for(j=0;j<r;j++)
		 {
			result[i][j]=0;
			}
		result[i][i]=1;
	}
	 return result;
 }
 
std::vector<std::vector<double> > Constant::get_value( std::vector<double>  &x, int d, int r)
{
	 std::vector<std::vector<double> > result;
	 int i,j;
	 result.resize(d);
	 for(i=0;i<d;i++)
		result[i].resize(r);
	 for(i=0;i<d;i++)
	 {
		 for(j=0;j<r;j++)
			result[i][j]=0;
		result[i][i]=1;
	}
	 return result;
}
	 
std::vector<double > Constant::get_value( std::vector<double>  x, int d)
 {
	 std::vector<double> result;
	 int i;
	 int n;
	 result.resize(d);
	 for(i=0;i<d;i++)
	 {
		result[i]=1;
	}
	 return result;
 }
	 
double Constant::get_value(double* x)
{
	return 1;
}
double * Constant::get_value(double* x, int d)
{
	double * result;
	result=(double *)calloc(d,sizeof(double));
	return result;
	free(result);
}

double* Constant::grad(double* x)
{
	double * result;
	int i,n;
	n=sizeof(x)/sizeof(double);
	result=(double *)calloc(n,sizeof(double));
	return result;
	free(result);
}
double** Constant::laplacian(double* x)
{
	double ** result;
	int i,n;
	n=sizeof(x)/sizeof(double);
	result=(double **)malloc(n*sizeof(double*));
	for(i=0;i<n;i++)
	{
		result[i]=(double *)calloc(n,sizeof(double));
	}
	return result;
			for (i = 0; i < n; i++)
	{
		free(result[i]);
	} 
	free(result);
}
double Constant::D3(double x)
{
	return 0;
}
double*** Constant::D3(double* x)
{
	double*** result;
	int i,j,k,n;
	n=sizeof(x)/sizeof(double);
	result=(double ***)malloc(n*sizeof(double**));
	for(i=0;i<n;i++)
	{
		result[i]=(double **)malloc(n*sizeof(double*));
		for(j=0;j<n;j++)
			result[i][j]=(double*)calloc(n,sizeof(double));
		
	}
	return result;
		for(i=0;i<n;i++)
	{
		for(j=0;j<n;j++)
			free(result[i][j]);
		free(result[i]);
	}
	free(result);
}
double Constant::D4(double x)
{
	return 0;
}

